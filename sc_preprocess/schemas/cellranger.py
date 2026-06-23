"""Cell Ranger configuration schemas for all modalities."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field, field_validator
from .base import BaseStepConfig, DirectoryConfig, ToolMeta


class AnndataConfig(BaseModel):
    """Resources for the create_{modality}_anndata rule."""
    threads: int = Field(default=16, ge=1, description="CPU threads")
    mem_gb: int = Field(default=32, gt=0, description="Memory (GB)")
    runtime_minutes: int = Field(default=120, gt=0, description="Wall-clock limit (minutes)")

    class Config:
        extra = "forbid"


class AggrConfig(BaseModel):
    """Resources for the cellranger_{modality}_aggr rule."""
    threads: int = Field(default=16, ge=1, description="CPU threads")
    mem_gb: int = Field(default=64, gt=0, description="Memory (GB)")
    runtime_minutes: int = Field(default=240, gt=0, description="Wall-clock limit (minutes)")

    class Config:
        extra = "forbid"


class BatchAggregationConfig(BaseModel):
    """Resources for the aggregate_{modality}_batch rule."""
    threads: int = Field(default=16, ge=1, description="CPU threads")
    mem_gb: int = Field(default=32, gt=0, description="Memory (GB)")
    runtime_minutes: int = Field(default=120, gt=0, description="Wall-clock limit (minutes)")

    class Config:
        extra = "forbid"


class EnrichmentConfig(BaseModel):
    """Resources for the enrich_{modality}_metadata rule."""
    threads: int = Field(default=4, ge=1, description="CPU threads")
    mem_gb: int = Field(default=16, gt=0, description="Memory (GB)")
    runtime_minutes: int = Field(default=60, gt=0, description="Wall-clock limit (minutes)")

    class Config:
        extra = "forbid"


class ClusterModeConfig(BaseModel):
    """Cell Ranger cluster mode — submits compute jobs via a cluster scheduler."""

    enabled: bool = Field(
        default=False,
        description="Enable cluster mode"
    )
    jobmode: str = Field(
        description=(
            "Job manager to use. Valid options: local (default), sge, lsf, slurm, "
            "or path to a custom .template file. All three Cell Ranger tools (GEX, ATAC, ARC) "
            "support 'slurm' natively. The custom path option allows a site-specific template "
            "with partition/account/time directives."
        )
    )
    mempercore: Optional[int] = Field(
        default=None,
        description=(
            "Reserve enough threads per job to ensure enough memory will be available, "
            "assuming each core on your cluster has at least this much memory (GB). "
            "Only applies to cluster jobmodes. "
            "SLURM users: leave null — the slurm.template requests memory directly via "
            "--mem=__MRO_MEM_GB__G, so mempercore has no effect on SLURM. "
            "Set this only for SGE/LSF clusters that cannot allocate memory independently."
        )
    )
    maxjobs: Optional[int] = Field(
        default=10,
        description="Maximum concurrent cluster jobs."
    )
    jobinterval: Optional[int] = Field(
        default=None,
        description=(
            "Delay between submitting jobs to cluster, in ms (10x default: 100ms). "
            "Increase on clusters with strict job submission rate limits."
        )
    )

    class Config:
        extra = "forbid"


class CellRangerGEXConfig(BaseStepConfig):
    """Cell Ranger GEX (Gene Expression) configuration."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger",
        url="https://www.10xgenomics.com/support/software/cell-ranger/latest",
        shell_version_cmd="cellranger --version",
    )

    reference: str = Field(
        description="Path to Cell Ranger GEX reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, sample, fastqs"
    )
    chemistry: Literal["auto", "ARC-v1", "threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2"] = Field(
        default="auto",
        description="Assay chemistry configuration"
    )
    normalize: Literal["none", "mapped", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    create_bam: bool = Field(
        default=False,
        description="Whether to create a BAM file",
        alias="create-bam"
    )
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger count"
    )
    runtime_minutes: int = Field(
        default=720,
        gt=0,
        description="Maximum runtime in minutes for the SLURM job"
    )
    
    cluster_mode: Optional[ClusterModeConfig] = Field(
        default=None,
        alias="cluster-mode",
        description="Cell Ranger cluster mode — see https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode"
    )

    anndata: AnndataConfig = Field(default_factory=AnndataConfig)
    aggr: AggrConfig = Field(default_factory=AggrConfig)
    batch_aggregation: BatchAggregationConfig = Field(default_factory=BatchAggregationConfig)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERGEX_COUNT_DIR="01_CELLRANGERGEX_COUNT",
            CELLRANGERGEX_AGGR_DIR="02_CELLRANGERGEX_AGGR"
        )
    )

    class DirectoryConfig(DirectoryConfig):
        """GEX-specific directories."""
        CELLRANGERGEX_COUNT_DIR: str = Field(default="01_CELLRANGERGEX_COUNT")
        CELLRANGERGEX_AGGR_DIR: str = Field(default="02_CELLRANGERGEX_AGGR")


class CellRangerATACConfig(BaseStepConfig):
    """Cell Ranger ATAC configuration."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger-atac",
        url="https://www.10xgenomics.com/support/software/cell-ranger-atac/latest",
        shell_version_cmd="cellranger-atac --version",
    )

    reference: str = Field(
        description="Path to Cell Ranger ATAC reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, sample, fastqs"
    )
    chemistry: Literal["auto", "ARC-v1"] = Field(
        default="auto",
        description="Assay chemistry configuration"
    )
    normalize: Literal["none", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger-atac count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger-atac count"
    )
    runtime_minutes: int = Field(
        default=720,
        gt=0,
        description="Maximum runtime in minutes for the SLURM job"
    )
    
    cluster_mode: Optional[ClusterModeConfig] = Field(
        default=None,
        alias="cluster-mode",
        description="Cell Ranger cluster mode — see https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/advanced/cluster-mode"
    )

    anndata: AnndataConfig = Field(default_factory=AnndataConfig)
    aggr: AggrConfig = Field(default_factory=AggrConfig)
    batch_aggregation: BatchAggregationConfig = Field(default_factory=BatchAggregationConfig)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERATAC_COUNT_DIR="01_CELLRANGERATAC_COUNT",
            CELLRANGERATAC_AGGR_DIR="02_CELLRANGERATAC_AGGR"
        )
    )

    class DirectoryConfig(DirectoryConfig):
        """ATAC-specific directories."""
        CELLRANGERATAC_COUNT_DIR: str = Field(default="01_CELLRANGERATAC_COUNT")
        CELLRANGERATAC_AGGR_DIR: str = Field(default="02_CELLRANGERATAC_AGGR")


class CellRangerARCConfig(BaseStepConfig):
    """Cell Ranger ARC (multiome: ATAC + RNA) configuration."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger-arc",
        url="https://www.10xgenomics.com/support/software/cell-ranger-arc/latest",
        shell_version_cmd="cellranger-arc --version",
    )

    reference: str = Field(
        description="Path to Cell Ranger ARC reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, CSV"
    )
    normalize: Literal["none", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    anndata: AnndataConfig = Field(default_factory=AnndataConfig)
    aggr: AggrConfig = Field(default_factory=AggrConfig)
    batch_aggregation: BatchAggregationConfig = Field(default_factory=BatchAggregationConfig)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERARC_COUNT_DIR="01_CELLRANGERARC_COUNT",
            CELLRANGERARC_AGGR_DIR="02_CELLRANGERARC_AGGR"
        )
    )
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger-arc count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger-arc count"
    )
    runtime_minutes: int = Field(
        default=720,
        gt=0,
        description="Maximum runtime in minutes for the SLURM job"
    )
    
    cluster_mode: Optional[ClusterModeConfig] = Field(
        default=None,
        alias="cluster-mode",
        description="Cell Ranger cluster mode — see https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/advanced/cluster-mode"
    )

    class DirectoryConfig(DirectoryConfig):
        """ARC-specific directories."""
        CELLRANGERARC_COUNT_DIR: str = Field(default="01_CELLRANGERARC_COUNT")
        CELLRANGERARC_AGGR_DIR: str = Field(default="02_CELLRANGERARC_AGGR")


class CellRangerMultiConfig(BaseStepConfig):
    """Cell Ranger multi configuration (5' immune profiling: GEX + VDJ + Feature Barcoding, and Flex).

    References and probe sets are declared inside each capture's multi config CSV,
    not as CLI arguments. This means all library types (GEX, VDJ-T, VDJ-B, Antibody Capture,
    CRISPR, Flex probe set) are configured per-capture in the CSV.
    """

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger",
        url="https://www.10xgenomics.com/support/software/cell-ranger/latest",
        shell_version_cmd="cellranger --version",
    )

    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, CSV (per-capture multi config CSV path)"
    )
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger multi"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger multi"
    )
    runtime_minutes: int = Field(
        default=720,
        gt=0,
        description="Maximum runtime in minutes for the SLURM job"
    )
    cluster_mode: Optional[ClusterModeConfig] = Field(
        default=None,
        alias="cluster-mode",
        description="Cell Ranger cluster mode — see https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode"
    )
    anndata: AnndataConfig = Field(default_factory=AnndataConfig)
    batch_aggregation: BatchAggregationConfig = Field(default_factory=BatchAggregationConfig)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)

    class DirectoryConfig(DirectoryConfig):
        """Multi-specific directories."""
        CELLRANGERMULTI_COUNT_DIR: str = Field(default="01_CELLRANGERMULTI_COUNT")
