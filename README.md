# BINF-6112 Final Project

Ian Jones <br>
ijones20@charlotte.edu <br>
801317591 <br>


##  Translation Based Multiple Sequence Aligner

This bioinformatics pipleline uses python to perform translation based multiple aequence alignments through open reading frame (ORF) detection, k-mer similarity scoring, and alignment using the Needleman-Wunsch algorithm.

### Project Overview
This package includes the following functionality:

  1. Detects open reading frames (ORFs) in nucleotide sequences
  2. Translates sequences into amino acids for alignment
  3. Uses k-mer similarity to sort sequences for progressive alignment
  4. Performs pairwise alignment with the Needleman-Wunsch algorithm
  5. Back-translates amino acid alignments into codon-aware nucleotide alignments
  6. Generates codon position statistics and visualizations

### Installation
```bash
git clone https://github.com/ianjones02/BINF_6112_Final
cd BINF_6112_Final
pip install -e .
```
### Usage

```bash
python -m msaligner.cli --input data/example_orfs_cds.fasta --outdir Results/
python -m msaligner.cli --input data/example_orfs_cds.fasta --outdir Results/ --kmer-size 5 --gap-penalty -2.5
python -m msaligner.cli --input data/example_orfs_cds.fasta --outdir outputs/run_verbose --kmer-size 4 --gap-penalty -1.5 --verbose
```

### Inputs/Outputs
- Input: nucleotide FASTA with CDS.
- Outputs: aligned_protein.fasta, aligned_codon.fasta, codon_stats.csv, variability_plot.png.

  ### AI Disclosure
  This project was completed in part using Julius AI (GPT-5) including preliminary scaffolding and code base, as well as generated FASTA files.

