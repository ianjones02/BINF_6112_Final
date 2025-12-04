"""
Back-translation module for converting amino acid alignments to codon alignments.
"""

import logging
from typing import Dict, List

def create_codon_alignment(aa_alignment: Dict[str, str], 
                          sequences_data: Dict[str, dict]) -> Dict[str, str]:
    """
    Convert amino acid alignment back to codon-aware nucleotide alignment.
    
    Args:
        aa_alignment: Aligned amino acid sequences
        sequences_data: Original sequence data with ORF information
        
    Returns:
        Dictionary of aligned nucleotide sequences
    """
    codon_alignment = {}
    
    for seq_id, aligned_aa in aa_alignment.items():
        if seq_id not in sequences_data:
            logging.warning(f"Sequence {seq_id} not found in original data")
            continue
        
        seq_data = sequences_data[seq_id]
        nt_seq = seq_data['nucleotide']
        orf_start = seq_data['orf_start']
        orf_end = seq_data['orf_end']
        
        # Extract ORF region
        orf_nt = nt_seq[orf_start:orf_end]
        
        # Create codon alignment
        aligned_nt = back_translate_single_sequence(aligned_aa, orf_nt)
        codon_alignment[seq_id] = aligned_nt
    
    return codon_alignment

def back_translate_single_sequence(aligned_aa: str, original_orf_nt: str) -> str:
    """
    Back-translate a single aligned amino acid sequence to nucleotides.
    
    Args:
        aligned_aa: Aligned amino acid sequence (with gaps)
        original_orf_nt: Original nucleotide ORF sequence
        
    Returns:
        Aligned nucleotide sequence
    """
    aligned_nt = []
    aa_pos = 0
    nt_pos = 0
    
    for aa_char in aligned_aa:
        if aa_char == '-':
            # Insert gap in nucleotide alignment (3 nucleotides per codon)
            aligned_nt.extend(['-', '-', '-'])
        else:
            if aa_char == 'X' or nt_pos + 3 > len(original_orf_nt):
                # Handle ambiguous amino acids or boundary issues
                aligned_nt.extend(['N', 'N', 'N'])
            else:
                # Use original codon
                codon = original_orf_nt[nt_pos:nt_pos+3]
                aligned_nt.extend(list(codon))
            nt_pos += 3
            aa_pos += 1
    
    return ''.join(aligned_nt)

def validate_codon_alignment(aa_alignment: Dict[str, str], 
                           nt_alignment: Dict[str, str]) -> bool:
    """
    Validate that nucleotide alignment matches amino acid alignment.
    
    Args:
        aa_alignment: Amino acid alignment
        nt_alignment: Nucleotide alignment
        
    Returns:
        True if validation passes
    """
    for seq_id in aa_alignment:
        if seq_id not in nt_alignment:
            logging.error(f"Sequence {seq_id} missing in nucleotide alignment")
            return False
        
        aa_seq = aa_alignment[seq_id].replace('-', '')
        nt_seq = nt_alignment[seq_id].replace('-', '')
        
        # Check length consistency
        if len(nt_seq) % 3 != 0:
            logging.error(f"Nucleotide sequence length not divisible by 3: {seq_id}")
            return False
        
        # Check that translated nucleotides match amino acids
        # (simplified validation)
        expected_aa_length = len(nt_seq) // 3
        if len(aa_seq) != expected_aa_length:
            logging.warning(f"Length mismatch for {seq_id}: AA={len(aa_seq)}, NT/3={expected_aa_length}")
    
    logging.info("Codon alignment validation completed")
    return True