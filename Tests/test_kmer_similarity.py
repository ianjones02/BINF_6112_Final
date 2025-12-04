"""
Unit tests for k-mer similarity module.
"""

import pytest
from msaligner.kmer_similarity import (
    generate_kmers, kmer_frequency, cosine_similarity, sort_sequences_by_similarity
)

def test_generate_kmers():
    """Test k-mer generation."""
    sequence = "ABCDEF"
    kmers = generate_kmers(sequence, 3)
    assert kmers == ["ABC", "BCD", "CDE", "DEF"]

def test_kmer_frequency():
    """Test k-mer frequency calculation."""
    sequence = "ABCABC"
    freq = kmer_frequency(sequence, 3)
    assert freq["ABC"] == 2

def test_cosine_similarity():
    """Test cosine similarity calculation."""
    vec1 = {"ABC": 2, "DEF": 1}
    vec2 = {"ABC": 1, "DEF": 2}
    similarity = cosine_similarity(vec1, vec2)
    assert 0 <= similarity <= 1

def test_sort_sequences_by_similarity():
    """Test sequence sorting by similarity."""
    sequences = {
        "seq1": {"amino_acid": "MAAAAA"},
        "seq2": {"amino_acid": "MAAAAA"},  # Similar to seq1
        "seq3": {"amino_acid": "MGGGGG"},  # Different
    }
    sorted_ids = sort_sequences_by_similarity(sequences, k=3)
    assert len(sorted_ids) == 3
    assert "seq1" in sorted_ids
    assert "seq2" in sorted_ids
    assert "seq3" in sorted_ids