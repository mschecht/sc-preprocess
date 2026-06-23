"""Config generator for single-cell preprocessing pipeline."""

_GEX_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32  # default memory; individual steps override this with their own mem_gb
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_gex:
  enabled: true
  reference: /path/to/cellranger/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, sample, fastqs
  chemistry: auto  # options: auto, SC3Pv1, SC3Pv2, SC3Pv3, fiveprime, SC5P-PE, SC5P-R2, ARC-v1
  normalize: none  # options: none, mapped, depth
  create-bam: false
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  cluster-mode:          # cellranger cluster-mode: see https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode
    enabled: false       # set true to submit Cell Ranger jobs via a cluster scheduler
    jobmode: slurm       # slurm | lsf | sge | /path/to/custom.template
    mempercore: null     # SLURM users: leave null — memory is requested directly via --mem=__MRO_MEM_GB__G
    maxjobs: 64          # max concurrent cluster subjobs
    jobinterval: null    # delay between submissions in ms; increase if cluster rate-limits
  # cellranger_gex_aggr
  aggr:
    threads: 16
    mem_gb: 64
    runtime_minutes: 240
  # create_gex_anndata
  anndata:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120
  # aggregate_gex_batch
  batch_aggregation:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: demuxalot  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_ATAC_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32  # default memory; individual steps override this with their own mem_gb
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_atac:
  enabled: true
  reference: /path/to/cellranger-atac/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, sample, fastqs
  chemistry: auto  # options: auto, ARC-v1
  normalize: none  # options: none, depth
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  cluster-mode:          # cellranger-atac cluster-mode: see https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/advanced/cluster-mode
    enabled: false       # set true to submit Cell Ranger ATAC jobs via a cluster scheduler
    jobmode: slurm       # slurm | lsf | sge | /path/to/custom.template
    mempercore: null     # SLURM users: leave null — memory is requested directly via --mem=__MRO_MEM_GB__G
    maxjobs: 64          # max concurrent cluster subjobs
    jobinterval: null    # delay between submissions in ms; increase if cluster rate-limits
  # cellranger_atac_aggr
  aggr:
    threads: 16
    mem_gb: 64
    runtime_minutes: 240
  # create_atac_anndata (includes fragment sort — can take 10–30 min on real data)
  anndata:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120
  # aggregate_atac_batch
  batch_aggregation:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: vireo  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_ARC_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32  # default memory; individual steps override this with their own mem_gb
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_arc:
  enabled: true
  reference: /path/to/cellranger-arc/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, CSV (per-capture CSV path)
  normalize: none  # options: none, depth
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  cluster-mode:          # cellranger-arc cluster-mode: see https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/advanced/cluster-mode
    enabled: false       # set true to submit Cell Ranger ARC jobs via a cluster scheduler
    jobmode: slurm       # slurm | lsf | sge | /path/to/custom.template
    mempercore: null     # SLURM users: leave null — memory is requested directly via --mem=__MRO_MEM_GB__G
    maxjobs: 64          # max concurrent cluster subjobs
    jobinterval: null    # delay between submissions in ms; increase if cluster rate-limits
  # cellranger_arc_aggr
  aggr:
    threads: 16
    mem_gb: 64
    runtime_minutes: 240
  # create_arc_mudata
  anndata:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120
  # aggregate_arc_batch
  batch_aggregation:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: vireo  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_MULTI_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32  # default memory; individual steps override this with their own mem_gb
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_multi:
  enabled: true
  libraries: /path/to/libraries_list.tsv  # REQUIRED: TSV with columns: batch, capture, CSV
  # No reference field — declare reference(s) inside each capture's multi config CSV.
  # The multi config CSV also specifies: probe-set (Flex), [vdj] reference, [feature] reference.
  # See: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi
  # Flex: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  cluster-mode:          # cellranger cluster-mode: see https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode
    enabled: false       # set true to submit Cell Ranger jobs via a cluster scheduler
    jobmode: slurm       # slurm | lsf | sge | /path/to/custom.template
    mempercore: null     # SLURM users: leave null — memory is requested directly via --mem=__MRO_MEM_GB__G
    maxjobs: 64          # max concurrent cluster subjobs
    jobinterval: null    # delay between submissions in ms; increase if cluster rate-limits
  # create_multi_mudata (VDJ + protein loading can be memory-intensive)
  anndata:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120
  # aggregate_multi_batch
  batch_aggregation:
    threads: 16
    mem_gb: 32
    runtime_minutes: 120

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: vireo  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_TEMPLATES = {
    "gex": _GEX_TEMPLATE,
    "atac": _ATAC_TEMPLATE,
    "arc": _ARC_TEMPLATE,
    "multi": _MULTI_TEMPLATE,
}


def generate_config_yaml(modality: str) -> str:
    """Return a fully-commented YAML config string for the given modality."""
    return _TEMPLATES[modality]
