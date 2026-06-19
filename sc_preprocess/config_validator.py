"""Configuration validation and parameter documentation utilities."""

import os
import sys
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, Type
from pydantic import ValidationError

from sc_preprocess.schemas.config import PipelineConfig
from sc_preprocess.schemas.cellranger import (
    CellRangerGEXConfig, CellRangerATACConfig, CellRangerARCConfig, CellRangerMultiConfig
)
from sc_preprocess.schemas.demultiplexing import (
    DemuxalotConfig, VireoConfig
)
from sc_preprocess.schemas.doublet_detection import ScrubletConfig


# Pipeline directory structure - update this when adding new analysis steps
PIPELINE_DIRECTORIES = [
    ("logs", "00_LOGS"),
    ("cellrangergex_count", "01_CELLRANGERGEX_COUNT"),
    ("cellrangergex_aggr", "02_CELLRANGERGEX_AGGR"),
    ("cellrangeratac_count", "01_CELLRANGERATAC_COUNT"),
    ("cellrangeratac_aggr", "02_CELLRANGERATAC_AGGR"),
    ("cellrangerarc_count", "01_CELLRANGERARC_COUNT"),
    ("cellrangerarc_aggr", "02_CELLRANGERARC_AGGR"),
    ("cellrangermulti_count", "01_CELLRANGERMULTI_COUNT"),
    ("cellrangermulti_aggr", "02_CELLRANGERMULTI_AGGR"),
    ("anndata", "03_ANNDATA"),  # Phase 1: Per-capture AnnData/MuData objects
    ("batch_objects", "04_BATCH_OBJECTS"),  # Phase 5: Batch-level aggregated objects
    ("demultiplexing", "05_DEMULTIPLEXING"),  # Renumbered
    ("doublet_detection", "06_DOUBLET_DETECTION"),
    ("final", "07_FINAL"),  # Future: Final merged objects with all metadata
]


def parse_output_directories(config: Dict[str, Any]) -> Dict[str, str]:
    """
    Parse output directory structure for all pipeline steps.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: Directory paths for all steps (keys follow pattern: {step_name}_dir)
    """
    output_dir = config.get("output_dir", "output")
    
    return {
        f"{step_name}_dir": os.path.join(output_dir, dir_path)
        for step_name, dir_path in PIPELINE_DIRECTORIES
    }


def parse_cellranger_config(config: Dict[str, Any], modality_key: str, output_dirs: Dict[str, str], has_chemistry: bool = True) -> Dict[str, Any]:
    """
    Parse Cell Ranger configuration for a specific modality.
    
    Args:
        config: Snakemake config dictionary
        modality_key: Key for the modality (e.g., "cellranger_gex")
        output_dirs: Pre-computed output directories dict
        has_chemistry: Whether this modality supports chemistry parameter
        
    Returns:
        dict: Parsed configuration with standardized keys
    """
    mod_config = config[modality_key]
    
    # Determine directory keys based on modality
    if "gex" in modality_key:
        count_key = "cellrangergex_count_dir"
        aggr_key = "cellrangergex_aggr_dir"
    elif "atac" in modality_key:
        count_key = "cellrangeratac_count_dir"
        aggr_key = "cellrangeratac_aggr_dir"
    elif "arc" in modality_key:
        count_key = "cellrangerarc_count_dir"
        aggr_key = "cellrangerarc_aggr_dir"
    elif "multi" in modality_key:
        count_key = "cellrangermulti_count_dir"
        aggr_key = "cellrangermulti_aggr_dir"

    is_multi = "multi" in modality_key

    parsed = {
        "libraries": mod_config["libraries"],
        "create_bam": mod_config.get("create-bam", False),
        "logs_dir": output_dirs["logs_dir"],
        "count_dir": output_dirs[count_key],
        "aggr_dir": output_dirs[aggr_key],
    }

    if not is_multi:
        parsed["reference"] = mod_config["reference"]
        parsed["normalize"] = mod_config.get("normalize", "none")

    if has_chemistry:
        parsed["chemistry"] = mod_config.get("chemistry", "auto")
    
    return parsed


class ConfigValidator:
    """Configuration validation and help utilities."""
    
    # Method-specific parameter schemas
    METHOD_SCHEMAS = {
        "cellranger": {
            "gex": CellRangerGEXConfig,
            "atac": CellRangerATACConfig,
            "arc": CellRangerARCConfig,
            "multi": CellRangerMultiConfig,
        },
        "demultiplexing": {
            "demuxalot": DemuxalotConfig,
            "vireo": VireoConfig,
        },
        "doublet_detection": {
            "scrublet": ScrubletConfig,
        },
    }
    
    @classmethod
    def validate_config_file(cls, config_path: str) -> PipelineConfig:
        """
        Validate a configuration file.

        Note: The 'run' command automatically validates configs before execution.
        This method is useful for:
        - Quick validation during config development
        - CI/CD pipelines that need standalone validation
        - Checking configs without triggering Snakemake setup

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Validated PipelineConfig object

        Raises:
            ValidationError: If configuration is invalid
            SystemExit: If any configured paths do not exist on disk
        """
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)

        try:
            config = PipelineConfig(**config_dict)
        except ValidationError as e:
            print(f"\n❌ Configuration validation failed for: {config_path}\n")
            cls._print_validation_errors(e)
            raise

        # Check that all user-provided paths exist on disk
        path_errors = cls._validate_paths(config)
        if path_errors:
            print(f"\n❌ Path validation failed for: {config_path}\n")
            for err in path_errors:
                print(f"  • {err}")
            print()
            raise SystemExit(1)

        return config
    
    @classmethod
    def _print_validation_errors(cls, error: ValidationError):
        """Print detailed validation errors with helpful messages."""
        for err in error.errors():
            location = " -> ".join(str(loc) for loc in err['loc'])
            msg = err['msg']
            err_type = err['type']
            
            print(f"\n  Location: {location}")
            print(f"  Error: {msg}")
            
            # Provide helpful context for common errors
            if 'required' in err_type.lower():
                print(f"  💡 This field is required. Please provide a value.")
            elif 'literal' in err_type.lower():
                # Extract allowed values from error context
                if 'expected' in err:
                    print(f"  💡 Allowed values: {err['expected']}")
            elif 'extra' in err_type.lower():
                print(f"  💡 This parameter is not recognized. Check for typos or see available parameters.")

    @classmethod
    def _validate_paths(cls, config: PipelineConfig) -> list[str]:
        """
        Check that all user-provided file/directory paths in the config exist on disk.

        Only checks paths for enabled pipeline steps.

        Args:
            config: Validated PipelineConfig object

        Returns:
            List of error messages for paths that do not exist (empty if all valid)
        """
        errors = []

        def _check(path_value: str, location: str):
            if not Path(path_value).exists():
                errors.append(f"{location}: path does not exist: {path_value}")

        # Cell Ranger GEX
        if config.cellranger_gex and config.cellranger_gex.enabled:
            _check(config.cellranger_gex.reference, "cellranger_gex.reference")
            _check(config.cellranger_gex.libraries, "cellranger_gex.libraries")

        # Cell Ranger ATAC
        if config.cellranger_atac and config.cellranger_atac.enabled:
            _check(config.cellranger_atac.reference, "cellranger_atac.reference")
            _check(config.cellranger_atac.libraries, "cellranger_atac.libraries")

        # Cell Ranger ARC
        if config.cellranger_arc and config.cellranger_arc.enabled:
            _check(config.cellranger_arc.reference, "cellranger_arc.reference")
            _check(config.cellranger_arc.libraries, "cellranger_arc.libraries")

        # Cell Ranger multi
        if config.cellranger_multi and config.cellranger_multi.enabled:
            _check(config.cellranger_multi.libraries, "cellranger_multi.libraries")

        # Demultiplexing
        if config.demultiplexing and config.demultiplexing.enabled:
            if config.demultiplexing.demuxalot:
                _check(config.demultiplexing.demuxalot.vcf, "demultiplexing.demuxalot.vcf")
                _check(config.demultiplexing.demuxalot.genome_names, "demultiplexing.demuxalot.genome_names")
            if config.demultiplexing.vireo:
                _check(config.demultiplexing.vireo.cellsnp.vcf, "demultiplexing.vireo.cellsnp.vcf")

        # Sample-level paths
        for sample_id, sample in config.samples.items():
            if sample.gex_fastqs:
                _check(sample.gex_fastqs, f"samples.{sample_id}.gex_fastqs")
            if sample.atac_fastqs:
                _check(sample.atac_fastqs, f"samples.{sample_id}.atac_fastqs")
            if sample.arc_csv:
                _check(sample.arc_csv, f"samples.{sample_id}.arc_csv")

        return errors

    @classmethod
    def show_method_params(cls, step: str | None = None, method: str | None = None):
        """
        Show available parameters for a specific method, or list available steps/methods.

        Args:
            step: Pipeline step (e.g., 'cellranger', 'demultiplexing', 'doublet_detection', 'celltype_annotation')
            method: Method name (e.g., 'gex', 'atac', 'arc', 'scrublet', 'demuxalot', etc.)

        Behavior:
            - No args: List all available steps
            - --step only: List methods for that step
            - --step and --method: Show parameters for that method
        """
        # Case 1: No step provided - list all available steps
        if step is None:
            print("\nAvailable steps:")
            for step_name in cls.METHOD_SCHEMAS.keys():
                print(f"  - {step_name}")
            print("\nUsage:")
            print("  sc-preprocess show-params --step <STEP>")
            print("  sc-preprocess show-params --step <STEP> --method <METHOD>")
            return

        # Validate step
        if step not in cls.METHOD_SCHEMAS:
            print(f"Unknown step: {step}")
            print(f"Available steps: {', '.join(cls.METHOD_SCHEMAS.keys())}")
            return

        # Case 2: Step provided but no method - list methods for that step
        if method is None:
            print(f"\nAvailable methods for '{step}':")
            for method_name in cls.METHOD_SCHEMAS[step].keys():
                print(f"  - {method_name}")
            print("\nUsage:")
            print(f"  sc-preprocess show-params --step {step} --method <METHOD>")
            return

        # Case 3: Both step and method provided - show parameters
        if method not in cls.METHOD_SCHEMAS[step]:
            print(f"Unknown method '{method}' for step '{step}'")
            print(f"Available methods: {', '.join(cls.METHOD_SCHEMAS[step].keys())}")
            return

        schema_class = cls.METHOD_SCHEMAS[step][method]
        cls._print_schema_params(schema_class, f"{step}.{method}")
    
    @staticmethod
    def _get_tool_version(tool_meta) -> str:
        """Get the installed version of a tool from its ToolMeta.

        Tries importlib.metadata first (Python packages), then falls back
        to a shell command if shell_version_cmd is set.

        Returns:
            Version string or 'not installed' if lookup fails.
        """
        import importlib.metadata
        import subprocess

        # Try Python package lookup first
        try:
            return importlib.metadata.version(tool_meta.package)
        except importlib.metadata.PackageNotFoundError:
            pass

        # Fall back to shell command for non-Python tools (e.g., Cell Ranger)
        if tool_meta.shell_version_cmd:
            from sc_preprocess.utils.version_check import CellRangerVersionChecker
            command = tool_meta.shell_version_cmd.split()[0]
            version = CellRangerVersionChecker.get_version(command)
            if version:
                return version

        return "not installed"

    @classmethod
    def _print_schema_params(cls, schema_class: Type, title: str):
        """Print parameters from a pydantic schema."""
        print(f"\n{'='*60}")
        print(f"Parameters for: {title}")

        # Display tool metadata if available
        tool_meta = getattr(schema_class, 'tool_meta', None)
        if tool_meta is not None:
            version = cls._get_tool_version(tool_meta)
            version_str = f"v{version}" if version != "not installed" else version
            print(f"Tool: {tool_meta.package} {version_str}")
            print(f"URL: {tool_meta.url}")

        print(f"{'='*60}\n")
        
        schema = schema_class.model_json_schema()
        properties = schema.get('properties', {})
        required = schema.get('required', [])
        
        for param_name, param_info in properties.items():
            is_required = param_name in required
            param_type = param_info.get('type', 'unknown')
            default = param_info.get('default', 'N/A')
            description = param_info.get('description', '')
            
            # Handle enum/literal types
            if 'enum' in param_info:
                param_type = f"Literal[{', '.join(repr(v) for v in param_info['enum'])}]"
            elif 'anyOf' in param_info:
                types = [t.get('type', 'unknown') for t in param_info['anyOf']]
                param_type = ' | '.join(types)
            
            # Format parameter info
            required_str = "REQUIRED" if is_required else f"optional (default: {default})"
            
            print(f"  {param_name}")
            print(f"    Type: {param_type}")
            print(f"    Status: {required_str}")
            if description:
                print(f"    Description: {description}")
            
            # Show constraints if any
            if 'minimum' in param_info or 'maximum' in param_info:
                constraints = []
                if 'minimum' in param_info:
                    constraints.append(f">= {param_info['minimum']}")
                if 'maximum' in param_info:
                    constraints.append(f"<= {param_info['maximum']}")
                print(f"    Constraints: {', '.join(constraints)}")
            
            print()
    
    @classmethod
    def list_all_methods(cls):
        """List all available methods for each step."""
        print("\n" + "="*60)
        print("Available Methods by Pipeline Step")
        print("="*60 + "\n")
        
        for step, methods in cls.METHOD_SCHEMAS.items():
            print(f"\n{step.replace('_', ' ').title()}:")
            for method in methods.keys():
                print(f"  • {method}")


def main():
    """CLI for config validation and help."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Configuration validation and parameter documentation"
    )
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate a configuration file')
    validate_parser.add_argument('config', help='Path to configuration file')
    
    # Show params command
    params_parser = subparsers.add_parser('show-params', help='Show parameters for a method')
    params_parser.add_argument('--step', required=True, 
                               help='Pipeline step (cellranger, demultiplexing, doublet_detection, celltype_annotation)')
    params_parser.add_argument('--method', required=True, help='Method name')
    
    # List methods command
    subparsers.add_parser('list-methods', help='List all available methods')
    
    # Generate example command
    example_parser = subparsers.add_parser('generate-example', help='Generate example configuration')
    example_parser.add_argument('--output', default='example_pipeline_config.yaml',
                                help='Output path for example config')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        try:
            config = ConfigValidator.validate_config_file(args.config)
            print(f"\n✓ Configuration is valid!")
            print(f"\nEnabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            sys.exit(1)
    
    elif args.command == 'show-params':
        ConfigValidator.show_method_params(args.step, args.method)
    
    elif args.command == 'list-methods':
        ConfigValidator.list_all_methods()
    
    elif args.command == 'generate-example':
        ConfigValidator.generate_example_config(args.output)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
