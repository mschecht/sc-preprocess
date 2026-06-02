"""Create AnnData object from ATAC Cell Ranger count output using SnapATAC2."""

import os
import sys
import gzip
import subprocess

import pandas as pd
import scanpy as sc
import snapatac2 as snap

# Access snakemake object
fragments_path = snakemake.input.fragments
peak_matrix_path = snakemake.input.peak_matrix
batch_id = snakemake.params.batch
capture_id = snakemake.params.capture
output_h5ad = snakemake.output.h5ad
output_snap_h5ad = snakemake.output.snap_h5ad
log_file = snakemake.log[0]
tmp_dir = snakemake.resources.tmpdir
n_threads = snakemake.threads

# Redirect output to log
sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout


def get_chrom_sizes_from_reference(fragments_path):
    """Extract chromosome sizes from reference path in fragments file header."""
    print("Extracting reference path from fragments file...")

    # Read fragments file header to get reference path
    reference_path = None
    with gzip.open(fragments_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# reference_path='):
                reference_path = line.strip().split('=', 1)[1]
                break

    if not reference_path:
        raise ValueError("Could not find reference_path in fragments file header")

    print(f"Reference path: {reference_path}")

    # Read chromosome sizes from reference .fai file
    fai_file = os.path.join(reference_path, "fasta", "genome.fa.fai")

    if not os.path.exists(fai_file):
        raise FileNotFoundError(f"Reference .fai file not found: {fai_file}")

    print(f"Reading chromosome sizes from: {fai_file}")

    chrom_sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom_name = fields[0]
            chrom_length = int(fields[1])
            chrom_sizes[chrom_name] = chrom_length

    print(f"✓ Found {len(chrom_sizes)} chromosomes: {list(chrom_sizes.keys())}")

    return chrom_sizes


try:
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)

    # ---- Step 1: Read peak matrix (Cell Ranger filtered cells x peaks) ----
    print(f"Reading Cell Ranger peak matrix from: {peak_matrix_path}")
    adata = sc.read_10x_h5(peak_matrix_path, gex_only=False)
    # Cell Ranger ATAC h5 stores peaks in var; rename for clarity
    adata.var_names_make_unique()
    print(f"✓ Loaded peak matrix: {adata.n_obs} cells × {adata.n_vars} peaks")

    # Keep filtered barcodes for SnapATAC2 import
    filtered_barcodes = set(adata.obs.index)

    # ---- Step 2: Compute QC metrics via SnapATAC2 import_fragments ----
    print(f"\nImporting ATAC fragments from: {fragments_path}")

    # Sort fragments by barcode for efficient import
    sorted_fragments = fragments_path.replace('.tsv.gz', '.sorted_by_barcode.tsv.gz')

    if not os.path.exists(sorted_fragments) or os.path.getsize(sorted_fragments) == 0:
        print(f"Sorting fragments by barcode for efficient import...")
        mem_gb = snakemake.resources.get("mem_mb", 32768) // 1024
        sort_mem = max(1, mem_gb // 2)
        sort_cmd = f"""
        zcat {fragments_path} | \
        grep -v '^#' | \
        sort -k4,4 -k1,1 -k2,2n -T {tmp_dir} --parallel={n_threads} -S {sort_mem}G | \
        bgzip --threads {n_threads} > {sorted_fragments}
        """
        subprocess.run(sort_cmd, shell=True, check=True, executable='/bin/bash')
        print(f"✓ Fragments sorted and saved to: {sorted_fragments}")
    else:
        print(f"Using existing sorted fragments: {sorted_fragments}")

    chrom_sizes = get_chrom_sizes_from_reference(fragments_path)

    # Import only Cell Ranger-filtered barcodes to get QC metrics
    snap_adata = snap.pp.import_fragments(
        sorted_fragments,
        chrom_sizes=chrom_sizes,
        file=output_snap_h5ad,
        sorted_by_barcode=True,
        whitelist=filtered_barcodes,
        min_num_fragments=1,
    )
    print(f"✓ SnapATAC2 dataset saved to: {output_snap_h5ad}")
    print(f"✓ SnapATAC2 QC computed for {snap_adata.n_obs} cells")

    # ---- Step 3: Transfer SnapATAC2 QC metrics to peak-matrix AnnData ----
    # When import_fragments writes to file=, .obs is PyDataFrameElem (backed),
    # not a pandas DataFrame — access known columns directly by name.
    qc_cols = ['n_fragment', 'frac_dup', 'frac_mito']
    print(f"Transferring QC metrics: {qc_cols}")

    snap_obs = pd.DataFrame(
        {col: list(snap_adata.obs[col]) for col in qc_cols},
        index=list(snap_adata.obs_names),
    ).reindex(adata.obs.index)
    for col in qc_cols:
        adata.obs[col] = snap_obs[col].values

    # ---- Step 4: Add traceability metadata ----
    adata.obs['batch_id'] = str(batch_id)
    adata.obs['capture_id'] = str(capture_id)
    adata.obs['cell_id'] = (
        adata.obs['batch_id'].astype(str) + '_' +
        adata.obs['capture_id'].astype(str) + '_' +
        adata.obs.index.astype(str)
    )

    if not adata.obs['cell_id'].is_unique:
        raise ValueError("cell_id is not unique!")
    adata.obs['barcode'] = adata.obs_names
    adata.obs_names = adata.obs['cell_id']

    print(f"\nATAC QC metrics: {qc_cols}")
    print(f"Cell metadata columns: {adata.obs.columns.tolist()}")
    print(f"Batch: {batch_id}, Capture: {capture_id}")
    print(f"First few cells:\n{adata.obs.head()}")

    # ---- Step 5: Compute QC metrics ----
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)

    # ---- Step 6: Write to disk ----
    print(f"\nWriting ATAC AnnData to: {output_h5ad}")
    adata.write_h5ad(output_h5ad)
    print(f"✓ ATAC AnnData creation complete: {adata.n_obs} cells × {adata.n_vars} peaks")

except Exception as e:
    print(f"ERROR creating ATAC AnnData: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
