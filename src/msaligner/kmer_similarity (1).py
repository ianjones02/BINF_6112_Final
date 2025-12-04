"""
k-mer similarity calculation and sequence sorting module.
"""

import logging
from typing import Dict, List, Tuple
import numpy as np
from collections import defaultdict

def generate_kmers(sequence: str, k: int) -> List[str]:
    """
    Generate all k-mers from a sequence.
    
    Args:
        sequence: Input sequence
        k: k-mer size
        
    Returns:
        List of k-mers
    """
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def kmer_frequency(sequence: str, k: int) -> Dict[str, int]:
    """
    Calculate k-mer frequency distribution.
    
    Args:
        sequence: Input sequence
        k: k-mer size
        
    Returns:
        Dictionary mapping k-mers to frequencies
    """
    kmers = generate_kmers(sequence, k)
    freq_dict = defaultdict(int)
    for kmer in kmers:
        freq_dict[kmer] += 1
    return freq_dict

def cosine_similarity(vec1: Dict[str, int], vec2: Dict[str, int]) -> float:
    """
    Calculate cosine similarity between two k-mer frequency vectors.
    
    Args:
        vec1: First k-mer frequency vector
        vec2: Second k-mer frequency vector
        
    Returns:
        Cosine similarity score
    """
    all_kmers = set(vec1.keys()) | set(vec2.keys())
    
    dot_product = 0
    norm1 = 0
    norm2 = 0
    
    for kmer in all_kmers:
        v1 = vec1.get(kmer, 0)
        v2 = vec2.get(kmer, 0)
        dot_product += v1 * v2
        norm1 += v1 ** 2
        norm2 += v2 ** 2
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    return dot_product / (np.sqrt(norm1) * np.sqrt(norm2))

def compute_kmer_similarity(sequences: Dict[str, dict], k: int = 6) -> np.ndarray:
    """
    Compute pairwise k-mer similarity matrix.
    
    Args:
        sequences: Dictionary of processed sequences
        k: k-mer size
        
    Returns:
        Similarity matrix
    """
    seq_ids = list(sequences.keys())
    n_seqs = len(seq_ids)
    similarity_matrix = np.zeros((n_seqs, n_seqs))
    
    # Precompute k-mer frequencies
    kmer_freqs = {}
    for seq_id in seq_ids:
        aa_seq = sequences[seq_id]['amino_acid']
        kmer_freqs[seq_id] = kmer_frequency(aa_seq, k)
    
    # Compute pairwise similarities
    for i in range(n_seqs):
        for j in range(i, n_seqs):
            if i == j:
                similarity_matrix[i, j] = 1.0
            else:
                seq1_id = seq_ids[i]
                seq2_id = seq_ids[j]
                sim = cosine_similarity(
                    kmer_freqs[seq1_id], 
                    kmer_freqs[seq2_id]
                )
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
    
    return similarity_matrix, seq_ids

def sort_sequences_by_similarity(sequences: Dict[str, dict], k: int=6) -> Dict[str, dict]:
    """
    Sort sequences for progressive alignment based on k-mer similarity.
    
    Args:
        sequences: Dictionary of processed sequences
        k: k-mer size
        
    Returns:
        List of sequence IDs in optimal order for progressive alignment
    """
    if len(sequences) <= 1:
        return sequences
    
    similarity_matrix, seq_ids = compute_kmer_similarity(sequences, k)
    
    # Use hierarchical clustering-like approach
    sorted_ids = [seq_ids[0]]  # Start with first sequence
    remaining_ids = seq_ids[1:]
    
    while remaining_ids:
        # Find sequence most similar to current set
        best_similarity = -1
        best_seq = None
        
        for candidate in remaining_ids:
            candidate_idx = seq_ids.index(candidate)
            avg_similarity = np.mean([
                similarity_matrix[candidate_idx, seq_ids.index(sorted_id)]
                for sorted_id in sorted_ids
            ])
            
            if avg_similarity > best_similarity:
                best_similarity = avg_similarity
                best_seq = candidate
        
        sorted_ids.append(best_seq)
        remaining_ids.remove(best_seq)
    
    logging.info(f"Sorted {len(sorted_ids)} sequences by k-mer similarity")
    return {seq_id: sequences[seq_id] for seq_id in sorted_ids}