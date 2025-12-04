"""
Codon position statistics calculation module.
"""

import logging
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
from collections import Counter

def calculate_codon_statistics(nt_alignment: Dict[str, str]) -> pd.DataFrame:
    """
    Calculate codon position statistics from nucleotide alignment.
    
    Args:
        nt_alignment: Aligned nucleotide sequences
        
    Returns:
        DataFrame with codon position statistics
    """
    if not nt_alignment:
        return pd.DataFrame()
    
    # Get alignment length (should be divisible by 3)
    first_seq = next(iter(nt_alignment.values()))
    alignment_length = len(first_seq)
    num_codons = alignment_length // 3
    
    # Initialize statistics arrays
    positions = []
    variability_scores = []
    gap_fractions = []
    conserved_positions = []
    
    for codon_pos in range(num_codons):
        # Extract nucleotides for this codon position across all sequences
        pos_nucleotides = []
        gap_count = 0
        total_sequences = 0
        
        for seq_id, aligned_nt in nt_alignment.items():
            codon_start = codon_pos * 3
            codon_end = codon_start + 3
            
            if codon_end > len(aligned_nt):
                continue
                
            codon = aligned_nt[codon_start:codon_end]
            
            # Check for gaps
            if '-' in codon:
                gap_count += 1
            else:
                pos_nucleotides.extend(list(codon))
            
            total_sequences += 1
        
        # Calculate statistics
        if total_sequences > 0:
            gap_fraction = gap_count / total_sequences
            
            # Variability score (1 - conservation)
            if pos_nucleotides:
                nucleotide_counts = Counter(pos_nucleotides)
                most_common_count = nucleotide_counts.most_common(1)[0][1]
                conservation = most_common_count / len(pos_nucleotides)
                variability = 1 - conservation
            else:
                variability = 0.0
            
            # Determine if position is conserved
            is_conserved = variability < 0.1 and gap_fraction < 0.5
            
            positions.append(codon_pos + 1)  # 1-based indexing
            variability_scores.append(variability)
            gap_fractions.append(gap_fraction)
            conserved_positions.append(is_conserved)
    
    # Create DataFrame
    stats_df = pd.DataFrame({
        'codon_position': positions,
        'variability_score': variability_scores,
        'gap_fraction': gap_fractions,
        'is_conserved': conserved_positions
    })
    
    logging.info(f"Calculated statistics for {len(positions)} codon positions")
    return stats_df

def generate_statistics(nt_alignment: Dict[str, str]) -> pd.DataFrame:
    """
    Main function to generate codon position statistics.
    
    Args:
        nt_alignment: Nucleotide alignment dictionary
        
    Returns:
        DataFrame with comprehensive statistics
    """
    return calculate_codon_statistics(nt_alignment)

def calculate_overall_statistics(stats_df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate overall alignment statistics.
    
    Args:
        stats_df: Codon position statistics DataFrame
        
    Returns:
        Dictionary of overall statistics
    """
    if stats_df.empty:
        return {}
    
    overall_stats = {
        'total_codon_positions': len(stats_df),
        'conserved_positions': stats_df['is_conserved'].sum(),
        'average_variability': stats_df['variability_score'].mean(),
        'average_gap_fraction': stats_df['gap_fraction'].mean(),
        'max_variability': stats_df['variability_score'].max(),
        'min_variability': stats_df['variability_score'].min()
    }
    
    return overall_stats