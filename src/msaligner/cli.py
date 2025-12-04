#!/usr/bin/env python3
"""
Command-line interface for the Translation-Based Multiple Sequence Aligner.
"""

import argparse
import sys
import os
from pathlib import Path
import logging

from msaligner.orf_detection import process_sequences
from msaligner.kmer_similarity import sort_sequences_by_similarity
from msaligner.alignment import perform_progressive_alignment
from msaligner.back_translation import create_codon_alignment
from msaligner.statistics import generate_statistics
from msaligner.visualization import create_visualizations

def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Translation-Based Multiple Sequence Aligner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  msaligner --input sequences.fasta --outdir results
  msaligner --input sequences.fasta --outdir results --kmer-size 6 --verbose
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Input FASTA file containing nucleotide sequences'
    )
    
    parser.add_argument(
        '--outdir', '-o',
        type=str,
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '--kmer-size', '-k',
        type=int,
        default=6,
        help='k-mer size for similarity calculation (default: 6)'
    )
    
    parser.add_argument(
        '--gap-penalty', '-g',
        type=int,
        default=-8,
        help='Gap penalty for alignment (default: -8)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def main():
    """Main entry point for the CLI."""
    args = parse_arguments()
    setup_logging(args.verbose)
    
    logging.info("Starting Translation-Based Multiple Sequence Alignment")
    
    # Validate input file
    if not os.path.isfile(args.input):
        logging.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: ORF detection and translation
        logging.info("Step 1: Detecting ORFs and translating sequences")
        sequences_data = process_sequences(args.input)
        
        # Step 2: k-mer similarity and sorting
        logging.info("Step 2: Calculating k-mer similarity and sorting sequences")
        sorted_sequences = sort_sequences_by_similarity(
            sequences_data, args.kmer_size
        )
        
        # Step 3: Progressive alignment
        logging.info("Step 3: Performing progressive alignment")
        aa_alignment = perform_progressive_alignment(
            sorted_sequences, args.gap_penalty
        )
        
        # Step 4: Back-translation
        logging.info("Step 4: Back-translating to nucleotide alignment")
        nt_alignment = create_codon_alignment(aa_alignment, sequences_data)
        
        # Step 5: Generate statistics
        logging.info("Step 5: Generating codon position statistics")
        stats_df = generate_statistics(nt_alignment)
        
        # Step 6: Create visualizations
        logging.info("Step 6: Creating visualizations")
        create_visualizations(stats_df, outdir)
        
        # Save outputs
        logging.info("Saving output files")
        
        # Save amino acid alignment
        aa_output = outdir / "aligned_amino_acids.fasta"
        with open(aa_output, 'w') as f:
            for seq_id, alignment in aa_alignment.items():
                f.write(f">{seq_id}\n{alignment}\n")
        
        # Save nucleotide alignment
        nt_output = outdir / "aligned_nucleotides.fasta"
        with open(nt_output, 'w') as f:
            for seq_id, alignment in nt_alignment.items():
                f.write(f">{seq_id}\n{alignment}\n")
        
        # Save statistics
        stats_output = outdir / "codon_statistics.csv"
        stats_df.to_csv(stats_output, index=False)
        
        logging.info(f"Pipeline completed successfully!")
        logging.info(f"Results saved to: {outdir}")
        logging.info(f"  - Amino acid alignment: {aa_output}")
        logging.info(f"  - Nucleotide alignment: {nt_output}")
        logging.info(f"  - Statistics: {stats_output}")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()