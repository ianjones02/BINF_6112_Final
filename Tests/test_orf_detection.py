"""
Unit tests for ORF detection module.
"""

import pytest
import tempfile
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from msaligner.orf_detection import (
    detect_orf, translate_sequence, read_fasta, 
    has_ambiguous_bases, process_sequences
)

def test_detect_orf_basic():
    """Test basic ORF detection."""
    # Simple sequence with ORF
    sequence = "ATG" + "GCC" * 10 + "TAA"  # Start + 10 codons + stop
    orf_coords = detect_orf(sequence)
    assert orf_coords is not None
    assert orf_coords[0] == 0
    assert orf_coords[1] == len(sequence)

def test_detect_orf_no_start():
    """Test ORF detection with no start codon."""
    sequence = "GCC" * 20  # No start codon
    orf_coords = detect_orf(sequence)
    assert orf_coords is None

def test_translate_sequence():
    """Test sequence translation."""
    sequence = "ATGGCTTGA"  # M A *
    orf_coords = (0, 9)
    aa_seq = translate_sequence(sequence, orf_coords)
    assert aa_seq == "MA*"

def test_has_ambiguous_bases():
    """Test ambiguous base detection."""
    assert has_ambiguous_bases("ATCG") == False
    assert has_ambiguous_bases("ATCN") == True
    assert has_ambiguous_bases("ATCR") == True

def test_read_fasta():
    """Test FASTA file reading."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">seq1\nATCG\n>seq2\nGCTA\n")
        temp_file = f.name
    
    try:
        sequences = read_fasta(temp_file)
        assert len(sequences) == 2
        assert "seq1" in sequences
        assert sequences["seq1"] == "ATCG"
    finally:
        os.unlink(temp_file)