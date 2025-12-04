
import matplotlib.pyplot as plt

def plot_codon_position_variability(stats_dict, save_path=None):
    positions = sorted(stats_dict.keys())
    variability = [stats_dict[p].get('variable', 0) for p in positions]
    plt.figure(figsize=(8,4))
    plt.plot(positions, variability, marker='o')
    plt.xlabel('Codon position')
    plt.ylabel('Variability')
    plt.title('Variability by codon position')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    plt.show()
