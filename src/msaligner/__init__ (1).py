"""
Translation-Based Multiple Sequence Aligner
"""
__version__ = "1.0.1"
__author__ = "Ian Jones"

from .alignment import needleman_wunsch, progressive_alignment, get_blosum_score
from .kmer_similarity import generate_kmers, kmer_frequency, cosine_similarity, sort_sequences_by_similarity
from .orf_detection import detect_orf, translate_sequence, read_fasta, has_ambiguous_bases, process_sequences
from .back_translation import create_codon_alignment, back_translate_single_sequence, validate_codon_alignment
from .statistics import calculate_codon_statistics, calculate_overall_statistics
