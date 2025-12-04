"""
Unit tests for alignment module.
"""

import pytest
from msaligner.alignment import (
    needleman_wunsch, progressive_alignment, get_blosum_score
)

def test_get_blosum_score():
    """Test BLOSUM score retrieval."""
    # Same amino acid should have positive score
    assert get_blosum_score("A", "A") > 0
    # Different amino acids should have lower scores
    assert get_blosum_score("A", "A") > get_blosum_score("A", "W")

def test_needleman_wunsch_identical():
    """Test alignment of identical sequences."""
    seq1 = "MAG"
    seq2 = "MAG"
    aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
    assert aligned1 == aligned2 == "MAG"
    assert score > 0

def test_needleman_wunsch_gaps():
    """Test alignment with gaps."""
    seq1 = "MAG"
    seq2 = "MG"
    aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
    assert len(aligned1) == len(aligned2)
    assert "-" in aligned1 or "-" in aligned2

def test_progressive_alignment_single():
    """Test progressive alignment with single sequence."""
    sequences = {"seq1": {"amino_acid": "MAG"}}
    sorted_ids = ["seq1"]
    aligned = progressive_alignment(sequences, sorted_ids)
    assert aligned["seq1"] == "MAG"