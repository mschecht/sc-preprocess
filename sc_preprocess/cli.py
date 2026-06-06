#!/usr/bin/env python

import sys
import os
import shlex
import shutil
import socket
import argparse
import subprocess

from pathlib import Path

from sc_preprocess.utils.utils import validate_cores
from sc_preprocess.config_generator import generate_config_yaml
from sc_preprocess.config_validator import ConfigValidator
from sc_preprocess.utils.custom_logger import custom_logger
from sc_preprocess.utils.version_check import CellRangerVersionChecker
from sc_preprocess.utils.test_data_generator import TestDataGenerator

__version__ = "0.1.0-dev"
__description__ = "Single-cell preprocessing pipeline for 10X Genomics data"

def main():
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate config for a GEX pipeline (print to stdout)
  sc-preprocess init-config --modality gex

  # Write config directly to a file
  sc-preprocess init-config --modality gex --output pipeline_config.yaml

  # Quick config validation (optional - 'run' auto-validates)
  sc-preprocess validate-config --config-file pipeline_config.yaml

  # Check installed Cell Ranger versions
  sc-preprocess check-versions
  sc-preprocess check-versions --workflow GEX

  # Run the pipeline (dry-run, auto-validates config first)
  sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dry-run

  # Run the pipeline (real run with 8 cores, auto-validates config first)
  sc-preprocess run --config-file pipeline_config.yaml --cores 8

  # Generate DAG visualization
  sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag.png

  # Show available parameters for a Cell Ranger modality
  sc-preprocess show-params --step cellranger --method gex

  # Show available parameters for a processing method
  sc-preprocess show-params --step doublet_detection --method scrublet

  # List all available methods
  sc-preprocess list-methods

  # Generate example test data files e.g. for ATAC workflow
  sc-preprocess generate-test-data ATAC --output-dir 00_TEST_DATA

  # Build and serve documentation locally
  sc-preprocess render-docs --port 8000
        """
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Subcommands
    subparsers = parser.add_subparsers(dest='subcommand', help='Available commands')
    
    # Run workflow subcommand
    run_parser = subparsers.add_parser(
        'run',
        help='Run the Snakemake workflow (auto-validates config)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Note: Config is automatically validated before running. No need to run validate-config first.

Additional Snakemake arguments can be passed via --snakemake-args.
IMPORTANT: --snakemake-args must be the LAST argument in your command.

Common Snakemake arguments:
  --forceall, -F             Force re-run of all rules
  --unlock                   Unlock working directory
  --profile DIR              Use a Snakemake profile for HPC/cloud execution
  --jobs N                   Use at most N CPU cluster/cloud jobs in parallel
  --printshellcmds, -p       Print shell commands that are executed

Examples:
  # Dry run to check what will be executed
  sc-preprocess run --config-file config.yaml --cores 1 --dry-run
  
  # Run with 8 cores
  sc-preprocess run --config-file config.yaml --cores 8
  
  # Use all available cores
  sc-preprocess run --config-file config.yaml --cores all

  # Unlock working directory
  sc-preprocess run --config-file config.yaml --cores 1 --snakemake-args --unlock

  # Pass additional arguments to Snakemake (e.g., limit parallel jobs)
  sc-preprocess run --config-file config.yaml --cores all --snakemake-args --jobs 10

  # Run with HPC execution via a Snakemake profile (note: --snakemake-args is LAST)
  sc-preprocess run --config-file config.yaml --cores all --snakemake-args --profile HPC_profiles
        """
    )
    run_parser.add_argument(
        '--config-file',
        required=True,
        help='Path to pipeline configuration file'
    )
    run_parser.add_argument(
        '--cores',
        required=True,
        type=validate_cores,
        help='Number of CPU cores for Snakemake to use (REQUIRED BY SNAKEMAKE). '
             'Specify a positive integer (e.g., 8) or "all" to use all available cores. '
             'This controls local parallelization. For cluster execution, also use --snakemake-args with --jobs flag.'
    )
    run_parser.add_argument(
        '--dry-run',
        '-n',
        action='store_true',
        help='Perform a dry run without executing actions'
    )
    run_parser.add_argument(
        '--dag',
        action='store_true',
        help='Generate workflow DAG (directed acyclic graph) visualization. Pipe to dot: | dot -Tpng > dag.png'
    )
    run_parser.add_argument(
        '--snakemake-args',
        nargs=argparse.REMAINDER,
        default=[],
        help='Additional arguments to pass to Snakemake. '
             'MUST be the last argument in the command. '
             'Example: --snakemake-args --forceall --jobs 10 --cluster "sbatch -J {rule} --ntasks=1 --cpus-per-task=12 --mem=40G"'
    )
    
    # Init config subcommand
    init_parser = subparsers.add_parser(
        'init-config',
        help='Generate a configuration file with all parameters for a given modality'
    )
    init_parser.add_argument(
        '--modality',
        required=True,
        choices=['gex', 'atac', 'arc'],
        help='Pipeline modality: gex (Gene Expression), atac (Chromatin Accessibility), arc (Multiome)'
    )
    init_parser.add_argument(
        '--output',
        default=None,
        help='Write config to this file (default: print to stdout)'
    )
    
    # Validate config subcommand
    validate_parser = subparsers.add_parser(
        'validate-config',
        help='Validate a configuration file (quick check without running workflow)'
    )
    validate_parser.add_argument(
        '--config-file',
        required=True,
        help='Path to configuration file to validate'
    )
    
    # Show params subcommand
    params_parser = subparsers.add_parser(
        'show-params',
        help='Show available parameters for a method (or list steps/methods if not specified)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all available steps
  sc-preprocess show-params

  # List methods for a specific step
  sc-preprocess show-params --step demultiplexing

  # Show parameters for a specific method
  sc-preprocess show-params --step demultiplexing --method vireo
        """
    )
    params_parser.add_argument(
        '--step',
        help='Pipeline step (cellranger, demultiplexing, doublet_detection, celltype_annotation)'
    )
    params_parser.add_argument(
        '--method',
        help='Method name (requires --step)'
    )
    
    # List methods subcommand
    subparsers.add_parser(
        'list-methods',
        help='List all available methods for each pipeline step'
    )
    
    # Check versions subcommand
    versions_parser = subparsers.add_parser(
        'check-versions',
        help='Check installed Cell Ranger versions'
    )
    versions_parser.add_argument(
        '--workflow',
        choices=['GEX', 'ATAC', 'ARC'],
        help='Check versions for specific workflow (default: check all)'
    )
    
    # Generate test data subcommand
    testdata_parser = subparsers.add_parser(
        'generate-test-data',
        help='Generate example test data files and configuration for pipeline testing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        This command creates a complete set of test files for the specified workflow, including:
          - libraries_list.tsv: TSV containing the metadata and paths for your Cell Ranger libraries
          - Reference genome path file
          - Test configuration YAML file with example settings

        These files provide a sandbox environment to test and develop your workflow.

        NOTE: The test data is derived from the Cell Ranger `testrun` command, which leverages test data 
        that comes with `cellranger`. For more information, see:
        https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in#testrun

        Examples:
          # Generate GEX workflow test files
          sc-preprocess generate-test-data GEX --output-dir 00_TEST_DATA_GEX

          # Generate ATAC workflow test files in custom directory
          sc-preprocess generate-test-data ATAC --output-dir 00_TEST_DATA_ATAC

          # Generate ARC workflow test files
          sc-preprocess generate-test-data ARC --output-dir 00_TEST_DATA_ARC
        """
    )
    testdata_parser.add_argument(
        'workflow',
        choices=['GEX', 'ATAC', 'ARC'],
        help='Workflow type: GEX (Gene Expression), ATAC (Chromatin Accessibility), or ARC (Multiome)'
    )
    testdata_parser.add_argument(
        '--output-dir',
        default='00_TEST_DATA',
        help='Directory where test data files will be generated (default: 00_TEST_DATA). '
             'Creates the directory if it does not exist.'
    )

    # Render docs subcommand
    render_docs_parser = subparsers.add_parser(
        'render-docs',
        help='Build and serve the documentation locally with live reload',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build docs and serve on default port 8000
  sc-preprocess render-docs

  # Use a custom port
  sc-preprocess render-docs --port 9000
        """
    )
    render_docs_parser.add_argument(
        '--port',
        type=int,
        default=8000,
        help='Port to serve documentation on (default: 8000)'
    )

    # Parse arguments
    args = parser.parse_args()

    # Handle subcommands
    if args.subcommand == 'run':
        # Get the path to the main.smk file
        package_dir = Path(__file__).parent
        snakefile = os.path.join(package_dir, "workflows", "main.smk")
        
        if not os.path.exists(snakefile):
            custom_logger.error(f"Snakefile not found: {snakefile}")
            sys.exit(1)
        
        # Validate config file exists
        if not os.path.exists(args.config_file):
            custom_logger.error(f"Config file not found: {args.config_file}")
            sys.exit(1)
        
        # Validate config before running
        try:
            config = ConfigValidator.validate_config_file(args.config_file)
            custom_logger.info(f"Config validated. Enabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            custom_logger.error("Config validation failed. Fix errors before running.")
            sys.exit(1)
        
        # Validate that reserved flags are not in snakemake-args
        if args.snakemake_args:
            snakemake_args_str = ' '.join(args.snakemake_args)
            if '--cores' in snakemake_args_str or '-c' in args.snakemake_args:
                custom_logger.error("Do not pass --cores or -c to --snakemake-args. Use the --cores flag instead.")
                sys.exit(1)
            if '--dry-run' in snakemake_args_str or '-n' in args.snakemake_args:
                custom_logger.error("Do not pass --dry-run or -n to --snakemake-args. Use the --dry-run flag instead.")
                sys.exit(1)
            if '--dag' in args.snakemake_args:
                custom_logger.error("Do not pass --dag to --snakemake-args. Use the --dag flag instead.")
                sys.exit(1)
        
        # Build snakemake command
        snakemake_cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", args.config_file,
            "--cores", args.cores,
            "--use-conda"
        ]
        
        # Add dry-run flag if specified
        if args.dry_run:
            snakemake_cmd.append("--dry-run")
        
        # Add dag flag if specified
        if args.dag:
            snakemake_cmd.append("--dag")
        
        # Add additional parameters if provided (already a list, no need to split)
        if args.snakemake_args:
            snakemake_cmd.extend(args.snakemake_args)

        custom_logger.info(f"Running Snakemake with command: {' '.join(snakemake_cmd)}")

        # Ensure conda is in PATH for snakemake's --use-conda integration.
        # When a conda env is active but the base conda dir is not in PATH,
        # snakemake fails to run `conda info`. Fix by adding the base conda bin.
        run_env = os.environ.copy()
        if not shutil.which("conda"):
            # Standard conda structure: CONDA_BASE/envs/ENV_NAME/bin/python -> CONDA_BASE/bin
            conda_base_bin = Path(sys.executable).parent.parent.parent.parent / "bin"
            if (conda_base_bin / "conda").exists():
                run_env["PATH"] = f"{conda_base_bin}:{run_env['PATH']}"

        # Run snakemake
        try:
            result = subprocess.run(snakemake_cmd, env=run_env)
            sys.exit(result.returncode)
        except FileNotFoundError:
            custom_logger.error("Snakemake not found. Please install snakemake: pip install snakemake")
            sys.exit(1)
        except KeyboardInterrupt:
            custom_logger.warning("Pipeline execution cancelled.")
            sys.exit(130)
    
    elif args.subcommand == 'init-config':
        yaml_str = generate_config_yaml(args.modality)
        if args.output is None:
            print(yaml_str, end="")
        else:
            output_path = Path(args.output)
            if output_path.exists():
                custom_logger.error(f"Output file already exists: {args.output}. Remove it first or choose a different name.")
                sys.exit(1)
            output_path.write_text(yaml_str)
            custom_logger.info(f"Config written to {args.output}")
    
    elif args.subcommand == 'validate-config':
        try:
            config = ConfigValidator.validate_config_file(args.config_file)
            custom_logger.info("Configuration is valid!")
            custom_logger.info(f"Enabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            sys.exit(1)
    
    elif args.subcommand == 'show-params':
        ConfigValidator.show_method_params(args.step, args.method)
    
    elif args.subcommand == 'list-methods':
        ConfigValidator.list_all_methods()
    
    elif args.subcommand == 'generate-test-data':
        workflow_type = args.workflow  # 'GEX', 'ATAC', or 'ARC'
        output_dir = Path(args.output_dir)
        
        custom_logger.info(f"Generating {workflow_type} test data in {output_dir}")
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate test data files (libraries TSV, reference file)
        test_data_paths = TestDataGenerator.generate_test_data(workflow_type, output_dir)
        
        if not test_data_paths:
            custom_logger.error("Failed to generate test data")
            sys.exit(1)
        
        # Generate test configuration file with actual paths
        config_path = os.path.join(output_dir, f"test_config_{workflow_type.lower()}.yaml")
        TestDataGenerator.generate_test_config(
            workflow=workflow_type,
            test_data_dir=str(output_dir),
            reference_path=test_data_paths['reference_path'],
            output_path=str(config_path)
        )
        
        custom_logger.info(f"\n{'='*60}")
        custom_logger.info(f"Test data generation complete!")
        custom_logger.info(f"{'='*60}")
        custom_logger.info(f"Libraries file: {test_data_paths['libraries_file']}")
        custom_logger.info(f"Reference path: {test_data_paths['reference_path']}")
        custom_logger.info(f"Config file: {config_path}")
        custom_logger.info(f"HPC profile: {output_dir}/HPC_profiles/config.yaml")
        custom_logger.info(f"\nTo run the {workflow_type} test pipeline:")
        custom_logger.info(f"  # Local execution:")
        custom_logger.info(f"  sc-preprocess run --config-file {config_path} --cores 8")
        custom_logger.info(f"\n  # Cluster execution (edit HPC_profiles/config.yaml first):")
        custom_logger.info(f"  sc-preprocess run --config-file {config_path} --cores 1 --snakemake-args --profile {output_dir}/HPC_profiles")
    
    elif args.subcommand == 'check-versions':
        checker = CellRangerVersionChecker()
        
        if args.workflow:
            # Check specific workflow
            is_valid, messages = checker.validate_for_workflow(args.workflow)
            
            for msg in messages:
                if "not found" in msg.lower() or "below minimum" in msg.lower():
                    custom_logger.error(msg)
                else:
                    custom_logger.info(msg)
            
            if is_valid:
                custom_logger.info(f"✓ All required tools for {args.workflow} workflow are installed and meet minimum version requirements")
                sys.exit(0)
            else:
                custom_logger.error(f"✗ Version check failed for {args.workflow} workflow")
                sys.exit(1)
        else:
            # Check all tools
            checker.print_all_versions()
            versions = checker.get_all_versions()
            
            if all(v is not None for v in versions.values()):
                custom_logger.info("\n✓ All Cell Ranger tools are installed")
                sys.exit(0)
            else:
                custom_logger.warning("\n⚠ Some Cell Ranger tools are not installed")
                sys.exit(0)  # Don't exit with error if just checking all
    
    elif args.subcommand == 'render-docs':
        docs_dir = Path(__file__).parent.parent / "docs"

        if not docs_dir.exists():
            custom_logger.error(f"Docs directory not found: {docs_dir}")
            sys.exit(1)

        if not shutil.which("sphinx-autobuild"):
            custom_logger.error("sphinx-autobuild not found. Install it with: pip install sphinx-autobuild")
            sys.exit(1)

        default_port = 8000
        port = args.port
        user_specified = port != default_port

        if not user_specified:
            max_attempts = 10
            for attempt in range(max_attempts):
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                    try:
                        s.bind(("127.0.0.1", port))
                        break
                    except OSError:
                        custom_logger.warning(f"Port {port} is in use, trying {port + 1}...")
                        port += 1
            else:
                custom_logger.error(
                    f"Could not find a free port in range {default_port}–{default_port + max_attempts - 1}. "
                    f"Use --port to specify one explicitly."
                )
                sys.exit(1)

        custom_logger.info(f"Serving docs at http://localhost:{port} (Ctrl+C to stop)")
        try:
            subprocess.run(
                ["sphinx-autobuild", "source", "build/html", "--port", str(port)],
                cwd=docs_dir
            )
        except KeyboardInterrupt:
            custom_logger.info("Docs server stopped.")

    else:
        parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    main()