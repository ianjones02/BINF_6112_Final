"""
ORF detection and sequence translation module.
"""

import logging
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Standard genetic code
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

# Ambiguous base handling
AMBIGUOUS_BASES = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Read sequences from a FASTA file.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq).upper()
        logging.info(f"Read {len(sequences)} sequences from {file_path}")
    except Exception as e:
        logging.error(f"Error reading FASTA file {file_path}: {str(e)}")
        raise
    
    return sequences

def has_ambiguous_bases(sequence: str) -> bool:
    """
    Check if sequence contains ambiguous bases.
    
    Args:
        sequence: Nucleotide sequence
        
    Returns:
        True if sequence contains ambiguous bases
    """
    return any(base in AMBIGUOUS_BASES for base in sequence)

def resolve_ambiguous_codon(codon: str) -> str:
    """
    Resolve ambiguous codon to possible amino acids.
    
    Args:
        codon: 3-nucleotide codon possibly containing ambiguous bases
        
    Returns:
        Amino acid or 'X' if ambiguous
    """
    if not has_ambiguous_bases(codon):
        return GENETIC_CODE.get(codon, 'X')
    
    # For ambiguous codons, return 'X' (unknown amino acid)
    logging.warning(f"Ambiguous codon detected: {codon}")
    return 'X'

def detect_orf(sequence: str, min_length: int = 30) -> Optional[Tuple[int, int]]:
    """
    Detect the longest open reading frame in a nucleotide sequence.
    
    Args:
        sequence: Nucleotide sequence
        min_length: Minimum ORF length in nucleotides
        
    Returns:
        Tuple of (start_position, end_position) or None if no ORF found
    """
    sequence = sequence.upper()
    best_orf = None
    max_length = 0
    
    # Check all reading frames
    for frame in range(3):
        # Look for start codon
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':  # Start codon
                # Look for stop codon
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in ['TAA', 'TAG', 'TGA']:
                        orf_length = j + 3 - i
                        if orf_length >= min_length and orf_length > max_length:
                            max_length = orf_length
                            best_orf = (i, j + 3)
                        break
    
    return best_orf

def translate_sequence(sequence: str, orf_coords: Tuple[int, int]) -> str:
    """
    Translate nucleotide sequence to amino acid sequence.
    
    Args:
        sequence: Nucleotide sequence
        orf_coords: ORF coordinates (start, end)
        
    Returns:
        Amino acid sequence
    """
    start, end = orf_coords
    orf_sequence = sequence[start:end]
    
    amino_acids = []
    for i in range(0, len(orf_sequence) - 2, 3):
        codon = orf_sequence[i:i+3]
        aa = resolve_ambiguous_codon(codon)
        amino_acids.append(aa)
    
    return ''.join(amino_acids)

def process_sequences(fasta_file: str) -> Dict[str, dict]:
    """
    Process all sequences in a FASTA file: detect ORFs and translate.
    
    Args:
        fasta_file: Path to input FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequence data
    """
    nucleotide_sequences = read_fasta(fasta_file)
    processed_sequences = {}
    
    for seq_id, nt_seq in nucleotide_sequences.items():
        # Check for ambiguous bases
        if has_ambiguous_bases(nt_seq):
            logging.warning(f"Sequence {seq_id} contains ambiguous bases")
        
        # Detect ORF
        orf_coords = detect_orf(nt_seq)
        
        if orf_coords is None:
            logging.warning(f"No ORF found in sequence {seq_id}")
            continue
        
        # Translate to amino acids
        aa_seq = translate_sequence(nt_seq, orf_coords)
        
        processed_sequences[seq_id] = {
            'nucleotide': nt_seq,
            'amino_acid': aa_seq,
            'orf_start': orf_coords[0],
            'orf_end': orf_coords[1],
            'has_ambiguous': has_ambiguous_bases(nt_seq)
        }
    
    logging.info(f"Successfully processed {len(processed_sequences)} sequences")
    return processed_sequences