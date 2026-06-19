"""Build target files for rule all based on pipeline configuration."""

import os
import pandas as pd

from pathlib import Path
from sc_preprocess.config_validator import parse_output_directories

def build_all_targets(config, enabled_steps):
    """
    Build list of final output files based on enabled steps.
    
    Args:
        config: Snakemake config dictionary
        enabled_steps: List of enabled pipeline step names
        
    Returns:
        list: Paths to final output files
    """
    # Defensive check - ensure enabled_steps is a list
    if not enabled_steps:
        return []
    
    if not isinstance(enabled_steps, list):
        enabled_steps = list(enabled_steps)
    
    targets = []
    
    # Cell Ranger outputs
    if "cellranger_gex" in enabled_steps:
        targets.extend(get_cellranger_gex_outputs(config))
    if "cellranger_atac" in enabled_steps:
        targets.extend(get_cellranger_atac_outputs(config))
    if "cellranger_arc" in enabled_steps:
        targets.extend(get_cellranger_arc_outputs(config))
    if "cellranger_multi" in enabled_steps:
        targets.extend(get_cellranger_multi_outputs(config))

    _cr_steps = ["cellranger_gex", "cellranger_atac", "cellranger_arc", "cellranger_multi"]

    # Phase 1: Object creation outputs (after cellranger)
    if any(step in enabled_steps for step in _cr_steps):
        targets.extend(get_object_creation_outputs(config))

    # Phase 5: Batch aggregation outputs (merge per-capture objects)
    if any(step in enabled_steps for step in _cr_steps):
        targets.extend(get_batch_aggregation_outputs(config))

    # Demux outputs
    if "demultiplexing" in enabled_steps:
        targets.extend(get_demux_outputs(config))

    # Doublet detection outputs
    if "doublet_detection" in enabled_steps:
        targets.extend(get_doublet_outputs(config))

    # Final enriched objects (merge all metadata into batch objects)
    if any(step in enabled_steps for step in _cr_steps):
        targets.extend(get_enriched_object_outputs(config))

    return targets


def parse_libraries_file(libraries_path):
    """
    Parse libraries TSV file to get sample information.
    
    Args:
        libraries_path: Path to libraries TSV file
        
    Returns:
        list: Sample/capture names
    """
    df = pd.read_csv(libraries_path, sep="\t")
    
    # For GEX/ATAC, samples are in 'capture' or 'sample' column
    if 'capture' in df.columns:
        return df['capture'].unique().tolist()
    elif 'sample' in df.columns:
        return df['sample'].unique().tolist()
    else:
        raise ValueError(f"Libraries file must have 'capture' or 'sample' column: {libraries_path}")


def get_cellranger_gex_outputs(config):
    """
    Get GEX output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to GEX done files
    """
    if not config.get("cellranger_gex"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    gex_config = config["cellranger_gex"]
    
    # Parse libraries to get batches
    df = pd.read_csv(gex_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_gex_aggr.done"))
    
    return outputs


def get_cellranger_atac_outputs(config):
    """
    Get ATAC output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to ATAC done files
    """
    if not config.get("cellranger_atac"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    atac_config = config["cellranger_atac"]
    
    # Parse libraries to get batches
    df = pd.read_csv(atac_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_atac_aggr.done"))
    
    return outputs

def get_cellranger_arc_outputs(config):
    """
    Get ARC output file paths.

    Args:
        config: Snakemake config dictionary

    Returns:
        list: Paths to ARC done files
    """
    if not config.get("cellranger_arc"):
        return []

    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    arc_config = config["cellranger_arc"]

    # Parse libraries to get batches
    df = pd.read_csv(arc_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()

    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_arc_aggr.done"))

    return outputs


def get_cellranger_multi_outputs(config):
    """Get multi output file paths (per-batch aggr done files)."""
    if not config.get("cellranger_multi"):
        return []

    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    multi_config = config["cellranger_multi"]

    df = pd.read_csv(multi_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()

    return [os.path.join(logs_dir, f"{batch}_multi_aggr.done") for batch in batches]


def get_object_creation_outputs(config):
    """
    Get AnnData/MuData object creation output file paths.
    Phase 1: Object creation after cellranger count (per-capture).

    Args:
        config: Snakemake config dictionary

    Returns:
        list: Paths to object creation done files
    """
    outputs = []
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]

    # GEX AnnData creation (per-capture)
    if config.get("cellranger_gex"):
        df = pd.read_csv(config["cellranger_gex"]["libraries"], sep="\t")
        for _, row in df.iterrows():
            batch = row['batch']
            capture = row['capture']
            outputs.append(os.path.join(logs_dir, f"{batch}_{capture}_gex_anndata.done"))

    # ATAC AnnData creation (per-capture)
    if config.get("cellranger_atac"):
        df = pd.read_csv(config["cellranger_atac"]["libraries"], sep="\t")
        for _, row in df.iterrows():
            batch = row['batch']
            capture = row['capture']
            outputs.append(os.path.join(logs_dir, f"{batch}_{capture}_atac_anndata.done"))

    # ARC MuData creation (per-capture)
    if config.get("cellranger_arc"):
        df = pd.read_csv(config["cellranger_arc"]["libraries"], sep="\t")
        for _, row in df.iterrows():
            batch = row['batch']
            capture = row['capture']
            outputs.append(os.path.join(logs_dir, f"{batch}_{capture}_arc_mudata.done"))

    # Multi MuData creation (per-capture)
    if config.get("cellranger_multi"):
        df = pd.read_csv(config["cellranger_multi"]["libraries"], sep="\t")
        for _, row in df.iterrows():
            batch = row['batch']
            capture = row['capture']
            outputs.append(os.path.join(logs_dir, f"{batch}_{capture}_multi_mudata.done"))

    return outputs


def get_batch_aggregation_outputs(config):
    """
    Get batch aggregation output file paths.
    Phase 5: Aggregate per-capture objects into batch-level objects.

    Args:
        config: Snakemake config dictionary

    Returns:
        list: Paths to batch aggregation done files
    """
    outputs = []
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]

    # GEX batch aggregation
    if config.get("cellranger_gex"):
        df = pd.read_csv(config["cellranger_gex"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_gex_batch_aggregation.done"))

    # ATAC batch aggregation
    if config.get("cellranger_atac"):
        df = pd.read_csv(config["cellranger_atac"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_atac_batch_aggregation.done"))

    # ARC batch aggregation
    if config.get("cellranger_arc"):
        df = pd.read_csv(config["cellranger_arc"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_arc_batch_aggregation.done"))

    # Multi batch aggregation
    if config.get("cellranger_multi"):
        df = pd.read_csv(config["cellranger_multi"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_multi_batch_aggregation.done"))

    return outputs


def get_enriched_object_outputs(config):
    """
    Get enriched object output file paths.
    Final step: Merge all analysis metadata into batch objects.

    Args:
        config: Snakemake config dictionary

    Returns:
        list: Paths to enriched object done files
    """
    outputs = []
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]

    # GEX enriched objects
    if config.get("cellranger_gex"):
        df = pd.read_csv(config["cellranger_gex"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_gex_enrichment.done"))

    # ATAC enriched objects
    if config.get("cellranger_atac"):
        df = pd.read_csv(config["cellranger_atac"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_atac_enrichment.done"))

    # ARC enriched objects
    if config.get("cellranger_arc"):
        df = pd.read_csv(config["cellranger_arc"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_arc_enrichment.done"))

    # Multi enriched objects
    if config.get("cellranger_multi"):
        df = pd.read_csv(config["cellranger_multi"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        for batch in batches:
            outputs.append(os.path.join(logs_dir, f"{batch}_multi_enrichment.done"))

    return outputs


def get_demux_outputs(config):
    """
    Get demultiplexing output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to demux output files
    """
    if not config.get("demultiplexing"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    demux_config = config["demultiplexing"]
    method = demux_config["method"]
    
    outputs = []
    
    # Get libraries from whichever modality is enabled
    libraries_path = None
    if config.get("cellranger_gex"):
        libraries_path = config["cellranger_gex"]["libraries"]
    elif config.get("cellranger_atac"):
        libraries_path = config["cellranger_atac"]["libraries"]
    elif config.get("cellranger_arc"):
        libraries_path = config["cellranger_arc"]["libraries"]
    elif config.get("cellranger_multi"):
        libraries_path = config["cellranger_multi"]["libraries"]

    if libraries_path:
        df = pd.read_csv(libraries_path, sep="\t")
        batches = df['batch'].unique().tolist()
        captures = df['capture'].unique().tolist()

        if method == "vireo":
            for batch in batches:
                for capture in captures:
                    outputs.append(os.path.join(logs_dir, f"vireo_output_{batch}_{capture}.done"))

        if method == "demuxalot":
            for batch in batches:
                for capture in captures:
                    outputs.append(os.path.join(logs_dir, f"demuxalot_output_{batch}_{capture}.done"))

    return outputs


def get_doublet_outputs(config):
    """
    Get doublet detection output file paths (per-capture sidecar TSVs).

    Args:
        config: Snakemake config dictionary

    Returns:
        list: Paths to doublet detection .done files
    """
    if not config.get("doublet_detection"):
        return []

    output_dirs = parse_output_directories(config)
    doublet_config = config["doublet_detection"]
    method = doublet_config["method"]
    logs_dir = output_dirs["logs_dir"]

    outputs = []

    # Get per-capture outputs from whichever modality is enabled
    libraries_path = None
    if config.get("cellranger_gex"):
        libraries_path = config["cellranger_gex"]["libraries"]
    elif config.get("cellranger_atac"):
        libraries_path = config["cellranger_atac"]["libraries"]
    elif config.get("cellranger_arc"):
        libraries_path = config["cellranger_arc"]["libraries"]
    elif config.get("cellranger_multi"):
        libraries_path = config["cellranger_multi"]["libraries"]

    if libraries_path:
        df = pd.read_csv(libraries_path, sep="\t")
        for _, row in df.iterrows():
            batch = row['batch']
            capture = row['capture']
            # Output pattern: {batch}_{capture}_{method}.done
            outputs.append(os.path.join(logs_dir, f"{batch}_{capture}_{method}.done"))

    return outputs


