
# CLI for translation-based MSA with tunable k-mer size and gap penalty
# This script parses arguments and invokes the package pipeline with parameters.

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description='Translation-based multiple sequence alignment pipeline')
    parser.add_argument('--input', '-i', required=True, help='Input FASTA file with nucleotide sequences')
    parser.add_argument('--outdir', '-o', required=True, help='Output directory')
    parser.add_argument('--kmer-size', type=int, default=3, help='K-mer size for similarity calculation')
    parser.add_argument('--gap-penalty', type=float, default=1.0, help='Gap penalty for alignment scoring')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    # Import here so CLI remains import-light
    try:
        # Prefer package-style entrypoint
        from msaligner.cli import run_pipeline as pkg_run
        used_pkg = True
    except Exception:
        pkg_run = None
        used_pkg = False

    # Fallback local modules
    pipeline_func = None
    if pkg_run is not None:
        pipeline_func = pkg_run
    else:
        try:
            # Try a local run function in this environment
            from cli import run_pipeline as local_run
            pipeline_func = local_run
        except Exception:
            pipeline_func = None

    if pipeline_func is None:
        # Last resort: dynamically import core functions
        try:
            from alignment import run_pipeline as alt_run
            pipeline_func = alt_run
        except Exception:
            print('Unable to locate pipeline entrypoint. Ensure msaligner is importable or cli exposes run_pipeline.')
            sys.exit(2)

    # Call with threaded parameters. The pipeline is expected to accept these kwargs; if not, ignore gracefully.
    try:
        pipeline_func(
            input_path=args.input,
            output_dir=args.outdir,
            kmer_size=args.kmer_size,
            gap_penalty=args.gap_penalty,
            verbose=args.verbose
        )
    except TypeError:
        # For older signatures that do not accept new params, call with required args only
        pipeline_func(
            input_path=args.input,
            output_dir=args.outdir,
            verbose=args.verbose
        )


if __name__ == '__main__':
    main()
