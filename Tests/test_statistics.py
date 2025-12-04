"""
Unit tests for statistics module.
"""

import pytest
import pandas as pd
from msaligner.statistics import calculate_codon_statistics, calculate_overall_statistics

def test_calculate_codon_statistics():
    """Test codon statistics calculation."""
    nt_alignment = {
        "seq1": "ATGGCT---TAA",  # 4 codons with gap
        "seq2": "ATGGCTTAGTAA",  # 4 complete codons
    }
    stats_df = calculate_codon_statistics(nt_alignment)
    assert not stats_df.empty
    assert "codon_position" in stats_df.columns
    assert "variability_score" in stats_df.columns
    assert "gap_fraction" in stats_df.columns

def test_calculate_overall_statistics():
    """Test overall statistics calculation."""
    stats_df = pd.DataFrame({
        'codon_position': [1, 2, 3],
        'variability_score': [0.1, 0.2, 0.3],
        'gap_fraction': [0.0, 0.1, 0.0],
        'is_conserved': [True, False, False]
    })
    overall_stats = calculate_overall_statistics(stats_df)
    assert "total_codon_positions" in overall_stats
    assert "conserved_positions" in overall_stats
    assert overall_stats["conserved_positions"] == 1