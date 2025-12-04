"""
Unit tests for back-translation module.
"""

import pytest
from msaligner.back_translation import (
    create_codon_alignment, back_translate_single_sequence, validate_codon_alignment
)

def test_back_translate_single_sequence():
    """Test single sequence back-translation."""
    aligned_aa = "MA-G"
    original_orf_nt = "ATGGCTGGT"  # MAG
    aligned_nt = back_translate_single_sequence(aligned_aa, original_orf_nt)
    # Should have gaps for the gap in AA alignment
    assert "---" in aligned_nt
    assert len(aligned_nt) == 12  # 4 codons * 3 nucleotides

def test_create_codon_alignment():
    """Test codon alignment creation."""
    aa_alignment = {"seq1": "MAG", "seq2": "M-G"}
    sequences_data = {
        "seq1": {
            "nucleotide": "ATGGCTGGTTAA",
            "orf_start": 0,
            "orf_end": 9,
            "has_ambiguous": False
        },
        "seq2": {
            "nucleotide": "ATGGGTTAA",
            "orf_start": 0,
            "orf_end": 6,
            "has_ambiguous": False
        }
    }
    codon_alignment = create_codon_alignment(aa_alignment, sequences_data)
    assert "seq1" in codon_alignment
    assert "seq2" in codon_alignment
    assert len(codon_alignment["seq1"]) == len(codon_alignment["seq2"])