"""
Visualization module for codon statistics and alignment results.
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def plot_codon_variability(stats_df: pd.DataFrame, output_path: Path):
    """
    Plot codon position variability.
    
    Args:
        stats_df: DataFrame with codon statistics
        output_path: Path to save the plot
    """
    if stats_df.empty:
        logging.warning("No data available for visualization")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot variability scores
    positions = stats_df['codon_position']
    variability = stats_df['variability_score']
    gap_fraction = stats_df['gap_fraction']
    
    ax1.plot(positions, variability, 'b-', linewidth=2, label='Variability')
    ax1.fill_between(positions, variability, alpha=0.3)
    ax1.set_xlabel('Codon Position')
    ax1.set_ylabel('Variability Score')
    ax1.set_title('Codon Position Variability')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot gap fractions
    ax2.bar(positions, gap_fraction, alpha=0.7, color='red', label='Gap Fraction')
    ax2.set_xlabel('Codon Position')
    ax2.set_ylabel('Gap Fraction')
    ax2.set_title('Gap Distribution Across Codon Positions')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_path / 'codon_variability_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Saved codon variability plot to {output_path / 'codon_variability_plot.png'}")

def plot_conservation_heatmap(stats_df: pd.DataFrame, output_path: Path):
    """
    Plot conservation pattern across codon positions.
    
    Args:
        stats_df: DataFrame with codon statistics
        output_path: Path to save the plot
    """
    if stats_df.empty:
        return
    
    fig, ax = plt.subplots(figsize=(15, 4))
    
    # Create binary conservation data
    conserved = stats_df['is_conserved'].astype(int).values
    positions = len(conserved)
    
    # Create heatmap-like plot
    im = ax.imshow([conserved], cmap='RdYlGn', aspect='auto', 
                   extent=[1, positions, 0, 1])
    
    ax.set_xlabel('Codon Position')
    ax.set_yticks([])
    ax.set_title('Conservation Pattern (Green = Conserved, Red = Variable)')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(['Variable', 'Conserved'])
    
    plt.tight_layout()
    plt.savefig(output_path / 'conservation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Saved conservation heatmap to {output_path / 'conservation_heatmap.png'}")

def create_visualizations(stats_df: pd.DataFrame, output_dir: Path):
    """
    Create all visualizations for the project.
    
    Args:
        stats_df: DataFrame with codon statistics
        output_dir: Output directory for saving plots
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plot_codon_variability(stats_df, output_dir)
    plot_conservation_heatmap(stats_df, output_dir)
    
    logging.info("All visualizations created successfully")