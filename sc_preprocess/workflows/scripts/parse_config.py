"""Configuration parsing utilities for Snakemake workflows."""

import sys
from pathlib import Path


def get_enabled_steps(config):
    """
    Extract enabled pipeline steps from config.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Names of enabled pipeline steps
    """
    steps = []
    
    # Check each step
    for step in ["cellranger_gex", "cellranger_atac", "cellranger_arc", "cellranger_multi",
                 "demultiplexing", "doublet_detection", "celltype_annotation"]:
        if config.get(step) and config[step].get("enabled", True):
            steps.append(step)
    
    return steps


def get_step_method(config, step):
    """
    Get the method selected for a pipeline step.
    
    Args:
        config: Snakemake config dictionary
        step: Name of the pipeline step
        
    Returns:
        str: Method name or None
    """
    step_config = config.get(step, {})
    return step_config.get("method")


def get_step_config(config, step):
    """
    Get configuration for a specific pipeline step.
    
    Args:
        config: Snakemake config dictionary
        step: Name of the pipeline step
        
    Returns:
        dict: Step configuration
    """
    return config.get(step, {})


def get_samples(config):
    """
    Get sample information from config.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: Sample metadata dictionary
    """
    return config.get("samples", {})


def get_sample_ids(config):
    """
    Get list of sample IDs.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Sample IDs
    """
    return list(config.get("samples", {}).keys())


def get_output_dir(config):
    """
    Get base output directory.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        str: Output directory path
    """
    return config.get("output_dir", "output")


def get_resources(config):
    """
    Get resource configuration.
    - [cellranger count args](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments#count)
    - `--localcores` and `--localmem`: Specify the number of cores and memory to use for local execution.
    - How resources flow from config to cluster:

        1. **YAML config**:
           ```yaml
           resources:
             mem_gb: 77    # Memory Cell Ranger will use
           ```

        2. **Snakemake rule** (`cellranger.smk`):
           ```python
           resources:
             mem_gb = RESOURCES.get("mem_gb", 64)  # Available to cluster command
           threads: 8                                # Available to cluster command
           ```

        3. **Cluster submission** (use Snakemake wildcards):
           ```bash
           --cluster 'sbatch --mem={resources.mem_gb}G --cpus-per-task={threads}'
           ```
           - `{resources.mem_gb}` → expands to `77` (from YAML)
           - `{threads}` → expands to `8` (from rule)
           - Cluster allocates what Cell Ranger will actually use

    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: Resource configuration
    """
    return config.get("resources", {})


def get_hpc_config(config):
    """
    Get HPC configuration.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: HPC configuration
    """
    return config.get("hpc", {})
