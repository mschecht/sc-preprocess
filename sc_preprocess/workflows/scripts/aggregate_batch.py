"""Aggregate per-capture objects into batch-level objects."""

import sys
import os

import muon as mu
import scanpy as sc
import anndata as ad

# Access snakemake object
batch_id = snakemake.params.batch
modality = snakemake.params.modality
output_file = snakemake.output.h5ad if modality in ["gex", "atac"] else snakemake.output.h5mu
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Aggregating batch {batch_id} for {modality.upper()} modality")

    if modality in ["arc", "multi"]:
        input_files = snakemake.input.h5mus
        print(f"Found {len(input_files)} per-capture MuData files to aggregate:")
    else:
        input_files = snakemake.input.h5ads
        print(f"Found {len(input_files)} per-capture AnnData files to aggregate:")

    for f in input_files:
        print(f"  - {f}")

    # Read all objects
    print("\nReading per-capture objects...")
    objects = []
    for i, file_path in enumerate(input_files):
        print(f"  [{i+1}/{len(input_files)}] Reading: {os.path.basename(file_path)}")

        if modality in ["arc", "multi"]:
            obj = mu.read(file_path)
        else:
            obj = sc.read_h5ad(file_path)

        objects.append(obj)
        print(f"      Loaded {obj.n_obs} cells")

    # Concatenate objects
    print(f"\nConcatenating {len(objects)} objects...")

    if modality in ["arc", "multi"]:
        mod_names = list(dict.fromkeys(k for obj in objects for k in obj.mod.keys()))
        merged_mods = {}
        for mod_name in mod_names:
            mod_list = [obj.mod[mod_name] for obj in objects if mod_name in obj.mod]
            merged_mods[mod_name] = ad.concat(mod_list, axis=0, join='outer', merge='same')
        batch_obj = mu.MuData(merged_mods)
        batch_obj.update()
        # Promote traceability metadata to top-level obs (muon prefixes these with modality name by default)
        first_mod = list(batch_obj.mod.values())[0]
        batch_obj.obs['batch_id'] = first_mod.obs['batch_id']
        batch_obj.obs['capture_id'] = first_mod.obs['capture_id']
        batch_obj.obs['cell_id'] = first_mod.obs['cell_id']
    else:
        batch_obj = ad.concat(objects, axis=0, join='outer', merge='same')

    print(f"✓ Concatenated to {batch_obj.n_obs} cells total")

    # Verify cell_id uniqueness
    if modality in ["arc", "multi"]:
        # For MuData, check each modality and top-level obs
        for mod_name, mod_data in batch_obj.mod.items():
            if not mod_data.obs['cell_id'].is_unique:
                raise ValueError(f"cell_id is not unique in {mod_name} modality after aggregation!")
        if not batch_obj.obs['cell_id'].is_unique:
            raise ValueError("cell_id is not unique at MuData level after aggregation!")
        print(f"✓ cell_id is unique in all modalities")
    else:
        # For AnnData
        if not batch_obj.obs['cell_id'].is_unique:
            raise ValueError("cell_id is not unique after aggregation!")
        print(f"✓ cell_id is unique")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Batch {batch_id} {modality.upper()} aggregation summary:")
    print(f"{'='*60}")
    print(f"Total cells: {batch_obj.n_obs}")

    if modality in ["arc", "multi"]:
        print(f"Modalities: {list(batch_obj.mod.keys())}")
        for mod_name, mod_data in batch_obj.mod.items():
            print(f"  {mod_name}: {mod_data.n_obs} cells x {mod_data.n_vars} features")
    else:
        print(f"Features: {batch_obj.n_vars}")

    if modality in ["arc", "multi"]:
        # For MuData, metadata lives in modality-level obs
        first_mod = list(batch_obj.mod.values())[0]
        print(f"Metadata columns (per-modality): {first_mod.obs.columns.tolist()}")
        print(f"\nCapture distribution:")
        print(first_mod.obs['capture_id'].value_counts())
    else:
        print(f"Metadata columns: {batch_obj.obs.columns.tolist()}")
        print(f"\nCapture distribution:")
        print(batch_obj.obs['capture_id'].value_counts())

    print(f"\n{'='*60}")

    # Write batch object
    print(f"\nWriting batch-level object to: {output_file}")

    if modality in ["arc", "multi"]:
        batch_obj.write(output_file)
    else:
        batch_obj.write_h5ad(output_file)

    print(f"✓ Batch aggregation complete!")

except Exception as e:
    print(f"ERROR in batch aggregation: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
