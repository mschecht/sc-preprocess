"""Cell Ranger workflow rules for all modalities (GEX, ATAC, ARC)."""

# Standard library imports
import os
import sys
import shutil

import pandas as pd

from pathlib import Path
from tempfile import gettempdir
from collections import defaultdict
from sc_preprocess import utils as u
from sc_preprocess.config_validator import parse_output_directories, parse_cellranger_config


# ============================================================================
# CELL RANGER GEX
# ============================================================================

# Get centralized output directories
OUTPUT_DIRS = parse_output_directories(config)

if config.get("cellranger_gex"):
    gex_cfg = parse_cellranger_config(config, "cellranger_gex", OUTPUT_DIRS, has_chemistry=True)
    
    GEX_REFERENCE = gex_cfg["reference"]
    GEX_LIBRARIES = gex_cfg["libraries"]
    GEX_CHEMISTRY = gex_cfg["chemistry"]
    GEX_CREATE_BAM = gex_cfg["create_bam"]
    GEX_NORMALIZE = gex_cfg["normalize"]
    GEX_LOGS_DIR = gex_cfg["logs_dir"]
    GEX_COUNT_DIR = gex_cfg["count_dir"]
    GEX_AGGR_DIR = gex_cfg["aggr_dir"]
    
    # Parse libraries file
    gex_df = u.sanity_check_libraries_list_tsv(
        GEX_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    # Build sample mappings - convert batch to string for consistency
    gex_samples = gex_df["capture"].unique().tolist()
    gex_batch_to_samples = {str(k): v for k, v in gex_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger GEX: Found {len(gex_samples)} sample(s) across {len(gex_batch_to_samples)} batch(es)")

    _gex_cm = config["cellranger_gex"].get("cluster-mode") or {}
    GEX_CLUSTER_JOBMODE     = _gex_cm.get("jobmode") if _gex_cm.get("enabled") else None
    GEX_CLUSTER_MEMPERCORE  = _gex_cm.get("mempercore") if _gex_cm.get("enabled") else None
    GEX_CLUSTER_MAXJOBS     = _gex_cm.get("maxjobs") if _gex_cm.get("enabled") else None
    GEX_CLUSTER_JOBINTERVAL = _gex_cm.get("jobinterval") if _gex_cm.get("enabled") else None


    rule cellranger_gex_count:
        """Run Cell Ranger count for gene expression data."""
        input:
            reference = GEX_REFERENCE,
            fastqs = lambda wc: gex_df[gex_df["capture"] == wc.capture]["fastqs"].iloc[0].split(",")
        output:
            h5 = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
            summary = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            bam = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "possorted_genome_bam.bam"),
            barcodes = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
            done = touch(os.path.join(GEX_LOGS_DIR, "{batch}_{capture}_gex_count.done"))
        params:
            outdir = GEX_COUNT_DIR,
            chemistry = GEX_CHEMISTRY,
            create_bam = GEX_CREATE_BAM,
            sample_name = lambda wc: gex_df[gex_df["capture"] == wc.capture]["sample"].iloc[0],
            jobmode     = GEX_CLUSTER_JOBMODE,
            mempercore  = GEX_CLUSTER_MEMPERCORE,
            maxjobs     = GEX_CLUSTER_MAXJOBS,
            jobinterval = GEX_CLUSTER_JOBINTERVAL
        threads: config["cellranger_gex"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_gex"].get("mem_gb", 64) * 1024,
            runtime = config["cellranger_gex"].get("runtime_minutes", 720),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(GEX_LOGS_DIR, "{batch}_{capture}_gex_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"
            create_bam_str = str(params.create_bam).lower() # Convert Python boolean to lowercase string for Cell Ranger

            # Remove stale _lock file so Cell Ranger can resume a partial run
            lock_file = os.path.join(output_id, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            if params.jobmode and params.jobmode != "local":
                cluster_flags = f"--jobmode={params.jobmode}"
                if params.mempercore:
                    cluster_flags += f" \\\n                    --mempercore={params.mempercore}"
                if params.maxjobs:
                    cluster_flags += f" \\\n                    --maxjobs={params.maxjobs}"
                if params.jobinterval:
                    cluster_flags += f" \\\n                    --jobinterval={params.jobinterval}"
            else:
                cluster_flags = f"--localcores={threads} \\\n                    --localmem={resources.mem_mb // 1024}"

            shell(
                f"""
                cellranger count \\
                    --id={output_id} \\
                    --transcriptome={input.reference} \\
                    --fastqs={",".join(input.fastqs)} \\
                    --sample={params.sample_name} \\
                    --create-bam={create_bam_str} \\
                    --chemistry={params.chemistry} \\
                    {cluster_flags} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_gex_aggr:
        """Aggregate multiple GEX samples."""
        input:
            samples = lambda wc: expand(
                os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
                batch=wc.batch,
                capture=gex_batch_to_samples[wc.batch]
            )
        output:
            done = touch(os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.done")),
            csv = os.path.join(GEX_AGGR_DIR, "{batch}_aggregation.csv")
        params:
            outdir = GEX_AGGR_DIR,
            normalize = GEX_NORMALIZE,
        threads: config["cellranger_gex"].get("aggr", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_gex"].get("aggr", {}).get("mem_gb", 64) * 1024,
            runtime = config["cellranger_gex"].get("aggr", {}).get("runtime_minutes", 240),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.log")
        run:
            lock_file = os.path.join(wildcards.batch, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            # Create aggregation CSV
            captures = gex_batch_to_samples[wildcards.batch]
            aggr_data = []
            for capture in captures:
                molecule_h5 = os.path.join(
                    GEX_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "molecule_info.h5"
                )
                aggr_data.append({"sample_id": f"{wildcards.batch}_{capture}", "molecule_h5": molecule_h5})
            
            pd.DataFrame(aggr_data).to_csv(output.csv, index=False)

            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger aggr \\
                        --id={wildcards.batch} \\
                        --csv={output.csv} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger aggr step.")


# ============================================================================
# CELL RANGER ATAC
# ============================================================================

if config.get("cellranger_atac"):
    atac_cfg = parse_cellranger_config(config, "cellranger_atac", OUTPUT_DIRS, has_chemistry=True)
    
    ATAC_REFERENCE = atac_cfg["reference"]
    ATAC_LIBRARIES = atac_cfg["libraries"]
    ATAC_CHEMISTRY = atac_cfg["chemistry"]
    ATAC_NORMALIZE = atac_cfg["normalize"]
    ATAC_LOGS_DIR = atac_cfg["logs_dir"]
    ATAC_COUNT_DIR = atac_cfg["count_dir"]
    ATAC_AGGR_DIR = atac_cfg["aggr_dir"]
    
    # Parse libraries file
    atac_df = u.sanity_check_libraries_list_tsv(
        ATAC_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    atac_samples = atac_df["capture"].unique().tolist()
    atac_batch_to_samples = {str(k): v for k, v in atac_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger ATAC: Found {len(atac_samples)} sample(s) across {len(atac_batch_to_samples)} batch(es)")

    _atac_cm = config["cellranger_atac"].get("cluster-mode") or {}
    ATAC_CLUSTER_JOBMODE     = _atac_cm.get("jobmode") if _atac_cm.get("enabled") else None
    ATAC_CLUSTER_MEMPERCORE  = _atac_cm.get("mempercore") if _atac_cm.get("enabled") else None
    ATAC_CLUSTER_MAXJOBS     = _atac_cm.get("maxjobs") if _atac_cm.get("enabled") else None
    ATAC_CLUSTER_JOBINTERVAL = _atac_cm.get("jobinterval") if _atac_cm.get("enabled") else None


    rule cellranger_atac_count:
        """Run Cell Ranger ATAC count."""
        input:
            reference = ATAC_REFERENCE,
            fastqs = lambda wc: atac_df[atac_df["capture"] == wc.capture]["fastqs"].iloc[0].split(",")
        output:
            h5 = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix.h5"),
            fragments = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "fragments.tsv.gz"),
            bam = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "possorted_bam.bam"),
            barcodes = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix", "barcodes.tsv"),
            summary = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            done = touch(os.path.join(ATAC_LOGS_DIR, "{batch}_{capture}_atac_count.done"))
        params:
            outdir = ATAC_COUNT_DIR,
            sample_name = lambda wc: atac_df[atac_df["capture"] == wc.capture]["sample"].iloc[0],
            jobmode     = ATAC_CLUSTER_JOBMODE,
            mempercore  = ATAC_CLUSTER_MEMPERCORE,
            maxjobs     = ATAC_CLUSTER_MAXJOBS,
            jobinterval = ATAC_CLUSTER_JOBINTERVAL
        threads: config["cellranger_atac"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_atac"].get("mem_gb", 64) * 1024,
            runtime = config["cellranger_atac"].get("runtime_minutes", 720),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ATAC_LOGS_DIR, "{batch}_{capture}_atac_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"

            # Remove stale _lock file so Cell Ranger ATAC can resume a partial run
            lock_file = os.path.join(output_id, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            if params.jobmode and params.jobmode != "local":
                cluster_flags = f"--jobmode={params.jobmode}"
                if params.mempercore:
                    cluster_flags += f" \\\n                    --mempercore={params.mempercore}"
                if params.maxjobs:
                    cluster_flags += f" \\\n                    --maxjobs={params.maxjobs}"
                if params.jobinterval:
                    cluster_flags += f" \\\n                    --jobinterval={params.jobinterval}"
            else:
                cluster_flags = f"--localcores={threads} \\\n                    --localmem={resources.mem_mb // 1024}"

            shell(
                f"""
                cellranger-atac count \\
                    --id={output_id} \\
                    --reference={input.reference} \\
                    --fastqs={",".join(input.fastqs)} \\
                    --sample={params.sample_name} \\
                    {cluster_flags} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_atac_aggr:
        """Aggregate multiple ATAC samples."""
        input:
            reference = ATAC_REFERENCE,
            samples = lambda wc: expand(
                os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix.h5"),
                batch=wc.batch,
                capture=atac_batch_to_samples[wc.batch]
            )
        output:
            done = touch(os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.done")),
            csv = os.path.join(ATAC_AGGR_DIR, "{batch}_aggregation.csv")
        params:
            outdir = ATAC_AGGR_DIR,
            normalize = ATAC_NORMALIZE,
        threads: config["cellranger_atac"].get("aggr", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_atac"].get("aggr", {}).get("mem_gb", 64) * 1024,
            runtime = config["cellranger_atac"].get("aggr", {}).get("runtime_minutes", 240),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.log")
        run:
            lock_file = os.path.join(wildcards.batch, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            # Create aggregation CSV
            captures = atac_batch_to_samples[wildcards.batch]
            aggr_data = []
            for capture in captures:
                fragments = os.path.abspath(os.path.join(
                    ATAC_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "fragments.tsv.gz"
                ))
                singlecell = os.path.abspath(os.path.join(
                    ATAC_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "singlecell.csv"
                ))
                aggr_data.append({"library_id": f"{wildcards.batch}_{capture}", "fragments": fragments, "cells": singlecell})
            
            pd.DataFrame(aggr_data).to_csv(output.csv, index=False)

            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger-atac aggr \\
                        --id={wildcards.batch} \\
                        --csv={output.csv} \\
                        --reference={input.reference} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger-atac aggr step.")


# ============================================================================
# CELL RANGER ARC
# ============================================================================

if config.get("cellranger_arc"):
    arc_cfg = parse_cellranger_config(config, "cellranger_arc", OUTPUT_DIRS, has_chemistry=False)
    
    ARC_REFERENCE = arc_cfg["reference"]
    ARC_LIBRARIES = arc_cfg["libraries"]
    ARC_NORMALIZE = arc_cfg["normalize"]
    ARC_LOGS_DIR = arc_cfg["logs_dir"]
    ARC_COUNT_DIR = arc_cfg["count_dir"]
    ARC_AGGR_DIR = arc_cfg["aggr_dir"]
    
    # Parse libraries file (ARC uses CSV column)
    arc_df = u.sanity_check_libraries_list_tsv(
        ARC_LIBRARIES,
        expected_columns={"batch", "capture", "CSV"},
        path_column="CSV",
        file_extension=".csv"
    )
    
    arc_captures = arc_df["capture"].unique().tolist()
    arc_batch_to_captures = {str(k): v for k, v in arc_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger ARC: Found {len(arc_captures)} capture(s) across {len(arc_batch_to_captures)} batch(es)")

    _arc_cm = config["cellranger_arc"].get("cluster-mode") or {}
    ARC_CLUSTER_JOBMODE     = _arc_cm.get("jobmode") if _arc_cm.get("enabled") else None
    ARC_CLUSTER_MEMPERCORE  = _arc_cm.get("mempercore") if _arc_cm.get("enabled") else None
    ARC_CLUSTER_MAXJOBS     = _arc_cm.get("maxjobs") if _arc_cm.get("enabled") else None
    ARC_CLUSTER_JOBINTERVAL = _arc_cm.get("jobinterval") if _arc_cm.get("enabled") else None


    rule cellranger_arc_count:
        """Run Cell Ranger ARC count for multiome data."""
        input:
            reference = ARC_REFERENCE,
            libraries_csv = lambda wc: arc_df[arc_df["capture"] == wc.capture]["CSV"].iloc[0]
        output:
            h5 = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
            fragments = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "atac_fragments.tsv.gz"),
            gex_bam = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "gex_possorted_bam.bam"),
            barcodes = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
            summary = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            done = touch(os.path.join(ARC_LOGS_DIR, "{batch}_{capture}_arc_count.done"))
        params:
            outdir       = ARC_COUNT_DIR,
            jobmode      = ARC_CLUSTER_JOBMODE,
            mempercore   = ARC_CLUSTER_MEMPERCORE,
            maxjobs      = ARC_CLUSTER_MAXJOBS,
            jobinterval  = ARC_CLUSTER_JOBINTERVAL
        threads: config["cellranger_arc"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_arc"].get("mem_gb", 64) * 1024,
            runtime = config["cellranger_arc"].get("runtime_minutes", 720),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ARC_LOGS_DIR, "{batch}_{capture}_arc_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"

            # Remove stale _lock file so Cell Ranger ARC can resume a partial run
            lock_file = os.path.join(output_id, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            if params.jobmode and params.jobmode != "local":
                cluster_flags = f"--jobmode={params.jobmode}"
                if params.mempercore:
                    cluster_flags += f" \\\n                    --mempercore={params.mempercore}"
                if params.maxjobs:
                    cluster_flags += f" \\\n                    --maxjobs={params.maxjobs}"
                if params.jobinterval:
                    cluster_flags += f" \\\n                    --jobinterval={params.jobinterval}"
            else:
                cluster_flags = f"--localcores={threads} \\\n                    --localmem={resources.mem_mb // 1024}"

            shell(
                f"""
                cellranger-arc count \\
                    --id={output_id} \\
                    --reference={input.reference} \\
                    --libraries={input.libraries_csv} \\
                    {cluster_flags} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_arc_aggr:
        """Aggregate multiple ARC captures."""
        input:
            captures = lambda wc: expand(
                os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
                batch=wc.batch,
                capture=arc_batch_to_captures[wc.batch]
            )
        output:
            done = touch(os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.done")),
            csv = os.path.join(ARC_AGGR_DIR, "{batch}_aggregation.csv")
        params:
            outdir = ARC_AGGR_DIR,
            normalize = ARC_NORMALIZE,
            reference = ARC_REFERENCE
        threads: config["cellranger_arc"].get("aggr", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_arc"].get("aggr", {}).get("mem_gb", 64) * 1024,
            runtime = config["cellranger_arc"].get("aggr", {}).get("runtime_minutes", 240),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.log")
        run:
            lock_file = os.path.join(wildcards.batch, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            # Create aggregation CSV for cellranger-arc aggr
            captures = arc_batch_to_captures[wildcards.batch]
            aggr_data = []
            for capture in captures:
                library_id = f"{wildcards.batch}_{capture}"
                atac_fragments = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "atac_fragments.tsv.gz"
                ))
                per_barcode_metrics = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "per_barcode_metrics.csv"
                ))
                gex_molecule_info = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "gex_molecule_info.h5"
                ))
                aggr_data.append({
                    "library_id": library_id,
                    "atac_fragments": atac_fragments,
                    "per_barcode_metrics": per_barcode_metrics,
                    "gex_molecule_info": gex_molecule_info
                })
            
            pd.DataFrame(aggr_data).to_csv(output.csv, index=False)

            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger-arc aggr \\
                        --id={wildcards.batch} \\
                        --csv={output.csv} \\
                        --reference={params.reference} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger-arc aggr step.")


# ============================================================================
# CELL RANGER MULTI
# ============================================================================

if config.get("cellranger_multi"):
    multi_cfg = parse_cellranger_config(config, "cellranger_multi", OUTPUT_DIRS, has_chemistry=False)

    MULTI_LIBRARIES = multi_cfg["libraries"]
    MULTI_LOGS_DIR = multi_cfg["logs_dir"]
    MULTI_COUNT_DIR = multi_cfg["count_dir"]
    # Parse libraries file (batch, capture, CSV — multi config CSV is not tabular, skip content validation)
    multi_df = u.sanity_check_libraries_list_tsv(
        MULTI_LIBRARIES,
        expected_columns={"batch", "capture", "CSV"},
        path_column="CSV",
        file_extension=".csv",
        validate_csv_contents=False
    )

    multi_captures = multi_df["capture"].unique().tolist()
    multi_batch_to_captures = {str(k): v for k, v in multi_df.groupby("batch")["capture"].apply(list).to_dict().items()}

    u.custom_logger.info(f"Cell Ranger multi: Found {len(multi_captures)} capture(s) across {len(multi_batch_to_captures)} batch(es)")

    _multi_cm = config["cellranger_multi"].get("cluster-mode") or {}
    MULTI_CLUSTER_JOBMODE     = _multi_cm.get("jobmode") if _multi_cm.get("enabled") else None
    MULTI_CLUSTER_MEMPERCORE  = _multi_cm.get("mempercore") if _multi_cm.get("enabled") else None
    MULTI_CLUSTER_MAXJOBS     = _multi_cm.get("maxjobs") if _multi_cm.get("enabled") else None
    MULTI_CLUSTER_JOBINTERVAL = _multi_cm.get("jobinterval") if _multi_cm.get("enabled") else None


    rule cellranger_multi_count:
        """Run Cell Ranger multi for a single capture (5' immune profiling or Flex).

        The multi config CSV (user-supplied per capture) declares all library types,
        references, and probe sets. Supports GEX + VDJ (TCR/BCR) + Feature Barcoding
        and Flex (Fixed RNA Profiling) in the same rule.
        """
        input:
            libraries_csv = lambda wc: multi_df[multi_df["capture"] == wc.capture]["CSV"].iloc[0]
        output:
            done = touch(os.path.join(MULTI_LOGS_DIR, "{batch}_{capture}_multi_count.done"))
        params:
            outdir      = MULTI_COUNT_DIR,
            jobmode     = MULTI_CLUSTER_JOBMODE,
            mempercore  = MULTI_CLUSTER_MEMPERCORE,
            maxjobs     = MULTI_CLUSTER_MAXJOBS,
            jobinterval = MULTI_CLUSTER_JOBINTERVAL
        threads: config["cellranger_multi"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_multi"].get("mem_gb", 64) * 1024,
            runtime = config["cellranger_multi"].get("runtime_minutes", 720),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(MULTI_LOGS_DIR, "{batch}_{capture}_multi_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"

            # Remove stale _lock file so Cell Ranger multi can resume a partial run
            lock_file = os.path.join(output_id, "_lock")
            if os.path.exists(lock_file):
                os.remove(lock_file)

            if params.jobmode and params.jobmode != "local":
                cluster_flags = f"--jobmode={params.jobmode}"
                if params.mempercore:
                    cluster_flags += f" \\\n                    --mempercore={params.mempercore}"
                if params.maxjobs:
                    cluster_flags += f" \\\n                    --maxjobs={params.maxjobs}"
                if params.jobinterval:
                    cluster_flags += f" \\\n                    --jobinterval={params.jobinterval}"
            else:
                cluster_flags = f"--localcores={threads} \\\n                    --localmem={resources.mem_mb // 1024}"

            shell(
                f"""
                cellranger multi \\
                    --id={output_id} \\
                    --csv={input.libraries_csv} \\
                    {cluster_flags} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)

