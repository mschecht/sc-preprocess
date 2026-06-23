"""Create MuData from Cell Ranger multi output (5' immune profiling or Flex).

Output structure detection:
- per_sample_outs/ present → modern layout (both 5' multi and Flex)
    - single sample dir   → 5' multi singleplex (cell_id = {batch}_{capture}_{bc})
    - multiple sample dirs → Flex multiplex (cell_id = {batch}_{capture}_{sample_id}_{bc})
- fallback: top-level filtered_feature_bc_matrix.h5 → legacy 5' multi

Per-sample matrix is located with find_sample_matrix(), which tries:
    per_sample/{sample}/count/sample_filtered_feature_bc_matrix.h5  (5' multi, newer CR)
    per_sample/{sample}/sample_filtered_feature_bc_matrix.h5        (Flex / older CR)

VDJ-T and VDJ-B are combined into a single `airr` modality via ir.pp.merge_airr()
if either (or both) are present. ir.pp.index_chains() is called on the MuData after
the airr modality is inserted.

Output is always one .h5mu per capture.
"""

import os
import sys
import warnings

import muon as mu
import scirpy as ir
import scanpy as sc
import anndata as ad

# Suppress known, handled third-party warnings that clutter the log
warnings.filterwarnings("ignore", message="Variable names are not unique", category=UserWarning)
warnings.filterwarnings("ignore", message="var_names.*not unique", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")
warnings.filterwarnings("ignore", message="Support for Awkward Arrays", category=Warning)
mu.set_options(pull_on_update=False)  # adopt mudata 0.4 behaviour proactively


def find_sample_matrix(sample_dir):
    """Locate the filtered matrix inside a per-sample output directory.

    Tries the 5' multi layout (count/ subdir) first, then the Flex flat layout.
    """
    candidates = [
        os.path.join(sample_dir, "count", "sample_filtered_feature_bc_matrix.h5"),
        os.path.join(sample_dir, "sample_filtered_feature_bc_matrix.h5"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


def _find_vdj_csv(vdj_dir):
    """Return path to best available VDJ annotation CSV in a Cell Ranger VDJ output directory."""
    for fname in ("filtered_contig_annotations.csv", "all_contig_annotations.csv"):
        path = os.path.join(vdj_dir, fname)
        if os.path.exists(path):
            return path
    return None


def _load_gex_prot_non_multiplexed(h5_path):
    """Load GEX + Antibody Capture from a single top-level h5 file (legacy 5' multi)."""
    mdata = mu.read_10x_h5(h5_path, extended=False)
    if "rna" in mdata.mod:
        mdata.mod["gex"] = mdata.mod.pop("rna")
    for mod_name in mdata.mod:
        mdata[mod_name].var_names_make_unique()
    return mdata


def _load_gex_prot_per_sample(per_sample_outs_dir):
    """Load GEX + Antibody Capture from per_sample_outs/.

    Returns (mdata, is_multiplexed). is_multiplexed is True when multiple sample
    directories are present (Flex multiplex); sample_id is added to obs only then.
    """
    sample_dirs = sorted([d for d in os.scandir(per_sample_outs_dir) if d.is_dir()])
    if not sample_dirs:
        raise FileNotFoundError(f"No sample directories found in {per_sample_outs_dir}")

    is_multiplexed = len(sample_dirs) > 1

    per_sample_mdatas = []
    for sample_entry in sample_dirs:
        sample_id = sample_entry.name
        h5_path = find_sample_matrix(sample_entry.path)
        if h5_path is None:
            print(
                f"  Warning: no filtered matrix under {sample_entry.path}, skipping {sample_id}",
                file=sys.stderr,
            )
            continue
        sample_mdata = mu.read_10x_h5(h5_path, extended=False)
        if "rna" in sample_mdata.mod:
            sample_mdata.mod["gex"] = sample_mdata.mod.pop("rna")
        for mod_name in sample_mdata.mod:
            sample_mdata[mod_name].var_names_make_unique()
        # Tag with sample_id only for multiplexed Flex so it surfaces in cell_id
        if is_multiplexed:
            for mod_name in sample_mdata.mod:
                sample_mdata[mod_name].obs["sample_id"] = sample_id
        per_sample_mdatas.append(sample_mdata)

    if not per_sample_mdatas:
        raise FileNotFoundError(f"No usable sample matrices found under {per_sample_outs_dir}")

    if len(per_sample_mdatas) == 1:
        mdata = per_sample_mdatas[0]
        mdata.update()
        return mdata, is_multiplexed

    # Concatenate per modality, then reconstruct MuData
    mod_names = list(per_sample_mdatas[0].mod.keys())
    merged_mods = {
        mod_name: ad.concat(
            [m.mod[mod_name] for m in per_sample_mdatas],
            axis=0,
            join="outer",
            merge="same",
            index_unique=None,
        )
        for mod_name in mod_names
        if all(mod_name in m.mod for m in per_sample_mdatas)
    }
    mdata = mu.MuData(merged_mods)
    mdata.update()
    return mdata, is_multiplexed


batch = snakemake.params.batch
capture = snakemake.params.capture
outs_dir = snakemake.params.outs_dir
vdj_t_dir = snakemake.params.vdj_t_dir
vdj_b_dir = snakemake.params.vdj_b_dir
output_h5mu = snakemake.output.h5mu
log_path = snakemake.log[0]

with open(log_path, "w") as log_fh:
    sys.stdout = log_fh
    sys.stderr = log_fh

    print(f"Creating MuData for {batch}_{capture}")

    # Detect output structure: per_sample_outs takes priority (modern CR multi);
    # fall back to top-level h5 for legacy non-demuxed 5' multi.
    per_sample_dir = os.path.join(outs_dir, "per_sample_outs")
    count_h5 = os.path.join(outs_dir, "filtered_feature_bc_matrix.h5")

    if os.path.exists(per_sample_dir):
        print(f"  Detected per-sample output: {per_sample_dir}")
        mdata, is_multiplexed = _load_gex_prot_per_sample(per_sample_dir)
        layout = "multiplexed Flex" if is_multiplexed else "5' multi (singleplex)"
        print(f"  Layout: {layout}")
    elif os.path.exists(count_h5):
        print(f"  Detected legacy top-level output: {count_h5}")
        mdata = _load_gex_prot_non_multiplexed(count_h5)
        is_multiplexed = False
    else:
        raise FileNotFoundError(
            f"Cannot find multi output at {per_sample_dir} or {count_h5}. "
            "Ensure cellranger multi completed successfully."
        )

    modalities_loaded = list(mdata.mod.keys())
    print(f"  Loaded modalities: {modalities_loaded}")

    # Compute GEX QC metrics
    adata_gex = mdata["gex"]
    qc_vars = []
    if adata_gex.var_names.str.startswith(("MT-", "MT.")).any():
        adata_gex.var["mt"] = adata_gex.var_names.str.startswith(("MT-", "MT."))
        qc_vars.append("mt")
    if adata_gex.var_names.str.startswith(("RPS", "RPL")).any():
        adata_gex.var["ribo"] = adata_gex.var_names.str.startswith(("RPS", "RPL"))
        qc_vars.append("ribo")
    sc.pp.calculate_qc_metrics(adata_gex, qc_vars=qc_vars, percent_top=None, inplace=True)

    # Load VDJ-T (TCR) and VDJ-B (BCR); combine into single airr modality
    loaded_vdj = {}
    for chain_type, vdj_dir in (("t", vdj_t_dir), ("b", vdj_b_dir)):
        if os.path.exists(vdj_dir):
            vdj_csv = _find_vdj_csv(vdj_dir)
            if vdj_csv:
                print(f"  Loading VDJ-{chain_type.upper()} from {vdj_csv}")
                loaded_vdj[chain_type] = ir.io.read_10x_vdj(vdj_csv)
            else:
                print(f"  No VDJ-{chain_type.upper()} annotation CSV in {vdj_dir}, skipping")
        else:
            print(f"  No VDJ-{chain_type.upper()} directory at {vdj_dir}, skipping")

    if len(loaded_vdj) == 2:
        airr_adata = ir.pp.merge_airr(loaded_vdj["t"], loaded_vdj["b"])
    elif len(loaded_vdj) == 1:
        airr_adata = next(iter(loaded_vdj.values()))
    else:
        airr_adata = None

    if airr_adata is not None:
        mdata.mod["airr"] = airr_adata
        print(f"  Added airr modality ({airr_adata.n_obs} cells with receptor data)")
    else:
        print("  No VDJ data found — airr modality omitted")

    # Add traceability metadata to every modality
    for mod_name, adata in mdata.mod.items():
        adata.obs.index.name = None  # scirpy sets obs index name to "cell_id"; clear it
        adata.obs["batch_id"] = str(batch)
        adata.obs["capture_id"] = str(capture)
        # For multiplexed Flex: include sample_id in cell_id to guarantee uniqueness
        if "sample_id" in adata.obs.columns:
            adata.obs["cell_id"] = [
                f"{batch}_{capture}_{row.sample_id}_{barcode}"
                for barcode, row in adata.obs.iterrows()
            ]
        else:
            adata.obs["cell_id"] = [
                f"{batch}_{capture}_{barcode}"
                for barcode in adata.obs_names
            ]

    mdata.update()

    # Expose traceability at the top-level MuData obs (driven by gex)
    mdata.obs["batch_id"] = mdata["gex"].obs["batch_id"]
    mdata.obs["capture_id"] = mdata["gex"].obs["capture_id"]
    mdata.obs["cell_id"] = mdata["gex"].obs["cell_id"]

    # Index chains and run chain QC after mdata.update() so scirpy columns
    # are written to mdata.obs last and not wiped by update().
    if "airr" in mdata.mod:
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)

    assert mdata["gex"].obs["cell_id"].is_unique, "cell_id must be unique across all cells!"

    print(f"  Total cells (gex): {mdata['gex'].n_obs}")
    print(f"  Final modalities: {list(mdata.mod.keys())}")

    os.makedirs(os.path.dirname(output_h5mu), exist_ok=True)
    mdata.write_h5mu(output_h5mu)
    print(f"  Written to {output_h5mu}")
