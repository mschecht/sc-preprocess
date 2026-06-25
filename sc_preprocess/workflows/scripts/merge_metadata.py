"""Merge analysis metadata into batch-level AnnData/MuData objects."""

import os
import sys
import glob
import pandas as pd
import scanpy as sc

# Access snakemake object
batch_id = snakemake.params.batch
modality = snakemake.params.modality
input_object = snakemake.input.batch_object
output_object = snakemake.output.enriched_object
obs_table_path = snakemake.output.obs_table
log_file = snakemake.log[0]

# Directories for metadata files
demux_dir = snakemake.params.get("demux_dir")
doublet_dir = snakemake.params.get("doublet_dir")

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout


def _join_by_cell_id(obs_df, metadata_df):
    """Left-join metadata onto obs using cell_id index, preserving obs index."""
    return obs_df.join(metadata_df, how='left')


def _parse_capture_and_method(basename, batch_id, suffix=".tsv.gz"):
    """
    Extract capture_id and method from a flat per-capture result filename.
    Expected format: {batch}_{capture}_{method}.tsv.gz
    e.g. '1_L001_scrublet.tsv.gz' -> ('L001', 'scrublet')
    Assumes capture IDs contain no underscores (L001, L002, etc.).
    """
    stem = basename[len(f"{batch_id}_"):].replace(suffix, "")  # "L001_scrublet"
    parts = stem.split("_")
    capture_id = parts[0]
    method = "_".join(parts[1:])
    return capture_id, method


def _build_cell_id_index(df, batch_id, capture_id):
    """Replace df index (barcodes) with cell_id = {batch}_{capture}_{barcode}."""
    df.index = df.index.map(lambda b: f"{batch_id}_{capture_id}_{b}")
    df.index.name = "cell_id"
    return df


def collect_step_metadata(step_dir, batch_id, col_prefix):
    """
    Collect all per-capture result files from a step directory.

    Convention: {step_dir}/{batch}_{capture}_{method}.tsv.gz
      - col 0 = plain barcodes (index)
      - remaining cols = method-specific results

    Returns a DataFrame indexed by cell_id, or None if no files found.
    """
    if not step_dir or not os.path.exists(step_dir):
        return None

    files = glob.glob(os.path.join(step_dir, f"{batch_id}_*.tsv.gz"))
    if not files:
        return None

    all_dfs = []
    for f in files:
        print(f"  - {os.path.basename(f)}")
        capture_id, method = _parse_capture_and_method(os.path.basename(f), batch_id)
        df = pd.read_csv(f, sep="\t", index_col=0)
        df = _build_cell_id_index(df, batch_id, capture_id)
        df.columns = [f"{col_prefix}_{method}_{col}" for col in df.columns]
        all_dfs.append(df)

    return pd.concat(all_dfs, axis=0)


try:
    print(f"{'='*70}")
    print(f"Enriching batch {batch_id} {modality.upper()} object with metadata")
    print(f"{'='*70}\n")

    # Read batch object
    print(f"Reading batch object: {input_object}")
    if modality in ["arc", "multi"]:
        import muon as mu
        adata = mu.read(input_object)
        is_mudata = True
    else:
        adata = sc.read_h5ad(input_object)
        is_mudata = False

    print(f"✓ Loaded {adata.n_obs} cells\n")

    metadata_added = []

    # ========================================================================
    # Merge per-step metadata
    # Convention for all steps: {step_dir}/{batch}_{capture}_{method}.tsv.gz
    #   col 0 = plain barcodes, joined to obs via cell_id
    # ========================================================================
    steps = [
        (demux_dir,   "demux",   "demultiplexing"),
        (doublet_dir, "doublet", "doublet_detection"),
    ]

    for step_dir, col_prefix, label in steps:
        print(f"Searching for {label} metadata...")
        metadata = collect_step_metadata(step_dir, batch_id, col_prefix)

        if metadata is not None:
            print(f"  Found {len(metadata.columns)} column(s) across {len(metadata)} cells")
            if is_mudata:
                for mod_name, mod_data in adata.mod.items():
                    mod_data.obs = _join_by_cell_id(mod_data.obs, metadata)
                    matched = mod_data.obs[metadata.columns].notna().any(axis=1).sum()
                    print(f"✓ Added {label} metadata to {mod_name} ({matched} cells matched)")
            else:
                adata.obs = _join_by_cell_id(adata.obs, metadata)
                matched = adata.obs[metadata.columns].notna().any(axis=1).sum()
                print(f"✓ Added {label} metadata ({matched} cells matched)")
            metadata_added.append(label)
        else:
            print(f"  No {label} files found")

        print()

    # ========================================================================
    # Summary
    # ========================================================================
    print(f"{'='*70}")
    print(f"Metadata enrichment summary:")
    print(f"{'='*70}")
    print(f"Batch: {batch_id}")
    print(f"Modality: {modality.upper()}")
    print(f"Cells: {adata.n_obs}")

    if is_mudata:
        print(f"\nMetadata columns per modality:")
        for mod_name, mod_data in adata.mod.items():
            print(f"  {mod_name}: {mod_data.obs.columns.tolist()}")
    else:
        print(f"\nMetadata columns: {adata.obs.columns.tolist()}")

    print(f"\nMetadata added: {', '.join(metadata_added) if metadata_added else 'None (no analysis completed yet)'}")
    print(f"{'='*70}\n")

    # Write enriched object
    print(f"Writing enriched object to: {output_object}")
    if is_mudata:
        adata.write(output_object)
    else:
        adata.write_h5ad(output_object)

    # Export obs table
    print(f"Writing obs table to: {obs_table_path}")
    if is_mudata:
        obs_frames = []
        for mod_name, mod_data in adata.mod.items():
            df = mod_data.obs.copy()
            df.insert(0, "modality", mod_name)
            obs_frames.append(df)
        obs_df = pd.concat(obs_frames, axis=0)
    else:
        obs_df = adata.obs.copy()
    obs_df.to_csv(obs_table_path, sep="\t", index=False, compression="gzip")

    print("✓ Metadata enrichment complete!\n")

    # QC summary per batch
    qc_summary_path = obs_table_path.replace('.tsv.gz', '_summary.tsv.gz')
    print(f'Writing QC summary to: {qc_summary_path}')
    summary_frames = []

    if is_mudata:
        for mod_name, mod_data in adata.mod.items():
            if 'total_counts' not in mod_data.obs.columns or 'n_genes_by_counts' not in mod_data.obs.columns:
                continue
            is_gex = mod_name in ('rna', 'gex')
            agg_dict = {
                'n_cells': ('total_counts', 'count'),
                'median_umi' if is_gex else 'median_fragments': ('total_counts', 'median'),
                'median_genes' if is_gex else 'median_peaks': ('n_genes_by_counts', 'median'),
            }
            df = mod_data.obs.groupby('batch_id').agg(**agg_dict).reset_index()
            df.insert(0, 'modality', mod_name)
            summary_frames.append(df)
        summary_df = pd.concat(summary_frames, axis=0)
    elif modality == 'gex':
        summary_df = adata.obs.groupby('batch_id').agg(
            n_cells=('total_counts', 'count'),
            median_umi=('total_counts', 'median'),
            median_genes=('n_genes_by_counts', 'median'),
        ).reset_index()
    else:
        summary_df = adata.obs.groupby('batch_id').agg(
            n_cells=('n_fragment', 'count'),
            median_fragments=('n_fragment', 'median'),
            median_peaks=('n_genes_by_counts', 'median'),
        ).reset_index()

    summary_df.to_csv(qc_summary_path, sep='\t', index=False, compression='gzip')
    print(summary_df)

except Exception as e:
    print(f"ERROR in metadata enrichment: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
