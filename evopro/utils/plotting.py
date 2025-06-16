import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from string import ascii_uppercase, ascii_lowercase

alphabet_list = list(ascii_uppercase+ascii_lowercase)

def create_line_plot(plot_name, vals):
    flat_vals = [item for sublist in vals for item in sublist]
    plt.rcParams.update({'font.size': 25})
    plt.figure(figsize=(15, 12), dpi=1200)
    for l, i in zip(vals, range(len(vals))):
        plt.plot(l, label="Score "+str(i))
    plt.xlabel("Iteration number")
    plt.ylabel("Weighted score value")
    plt.ylim([int(min(flat_vals) - abs(min(flat_vals)*0.1)), int(max(flat_vals) + abs(max(flat_vals)*0.1))])
    plt.legend()
    plt.savefig(plot_name)
    
def save_csv(csv_name, vals):
    iter_list = list(range(1, len(vals[0])+1))
    d = {'iter':iter_list}
    for l, i in zip(vals, range(len(vals))):
        d['score '+str(i)] = l
    df = pd.DataFrame(data = d)
    df.to_csv(csv_name)

def plot_scores(data, plot_name, stat_type='average'):
    """
    Plot score trajectories using either mean, max, or median statistics per iteration.
    
    Args:
        data: List of iterations containing score sets
        stat_type: String indicating statistic to use ('average', 'top', or 'median')
    """
    if stat_type not in ['average', 'top', 'median']:
        raise ValueError("stat_type must be either 'average', 'top', or 'median'")
    
    # Function to compute the requested statistic
    def compute_stat(scores):
        if stat_type == 'average':
            return np.mean(scores)
        elif stat_type == 'median':
            return np.median(scores)
        elif stat_type == 'top':
            return np.max(scores)
    
    # Extract scores and create iteration numbers
    iterations = range(len(data))
    
    # Get the number of components
    num_components = len(data[0][0][0][1:])  # Exclude the total score
    
    # Initialize lists to store statistical scores
    total_scores_stat = []
    component_scores_stat = [[] for _ in range(num_components)]
    weighted_scores_stat = [[] for _ in range(num_components)]
    
    # Calculate statistics for each iteration
    for iteration in data:
        # Calculate statistic for total scores
        all_scores = [score_set[0] for score_set in iteration]
        sorted_scores = sorted(all_scores, key=lambda x: x[0])
        
        # Take only the top half of the sorted scores
        top_half_count = max(1, len(sorted_scores) // 2)  # Ensure at least one score is used
        top_half_scores = sorted_scores[:top_half_count]
        
        total_scores = [score[0] for score in top_half_scores]
        total_scores_stat.append(compute_stat(total_scores))
        
        # Calculate statistic for component scores
        for i in range(num_components):
            component_scores = [score[i+1][0] for score in top_half_scores]
            component_scores_stat[i].append(compute_stat(component_scores))
            
            weighted_scores = [score[i+1][0] * score[i+1][1] for score in top_half_scores]
            weighted_scores_stat[i].append(compute_stat(weighted_scores))
    
    plot_name_weighted = plot_name.replace('.png', '_weighted.png')
    create_line_plot(plot_name_weighted, [total_scores_stat, *weighted_scores_stat])
    # Create the plot
    plt.rcParams.update({'font.size': 25})
    plt.figure(figsize=(15, 12), dpi=1200)   
       
    # Plot total score statistic
    if stat_type == 'average':
        stat_label = 'Average'
    elif stat_type == 'median':
        stat_label = 'Median'
    else:
        stat_label = 'Top'
    
    plt.plot(iterations, total_scores_stat, 'b-', 
            label=f'{stat_label} Total Score', linewidth=2)
    
    # Colors for components
    colors = ['r', 'g', 'm', 'c', 'y', 'k', 'b']
    
    # Plot component score statistics
    for i in range(num_components):
        color_idx = i % len(colors)  # Cycle through colors if more components than colors
        plt.plot(iterations, component_scores_stat[i], f'{colors[color_idx]}--', 
                label=f'Component {i+1} (weight={data[0][0][0][i+1][1]})', 
                alpha=0.7)
    
    # Customize the plot
    plt.xlabel('Iteration')
    plt.ylabel(f'{stat_label} Score')
    plt.title(f'{stat_label} Score Trajectories Across Iterations')
    plt.legend(loc='center right')

    
    # Adjust layout to prevent legend cutoff
    plt.savefig(plot_name)
    plt.close()
    
    return plot_name

def plot_ticks(Ls):
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color='black')
        plt.plot([L, L], [0, Ln], color='black')
    ticks = np.cumsum([0]+Ls)
    ticks = (ticks[1:] + ticks[:-1])/2
    plt.yticks(ticks, alphabet_list[:len(ticks)])


def plot_pae(pae, Ls=None, dpi=300, fig=True, ptm=None, iptm=None):
    if fig: plt.figure(figsize=(15,10), dpi=dpi)
    pae_title = 'Predicted Aligned Error.'
    if ptm:
        pae_title += f' ptm = {ptm:.3f}.'
    if iptm:
        pae_title += f' iptm = {iptm:.3f}.'
    plt.title(pae_title)
    Ln = pae.shape[0]
    plt.imshow(pae, cmap='bwr', vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
    if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
    plt.colorbar()
    plt.ylabel('Positions')
    plt.xlabel('Positions')
    return plt


def plot_plddt(plddt, Ls=None, dpi=100, fig=True):
    if fig: plt.figure(figsize=(15,10), dpi=300)
    plt.title(f'Predicted lDDT per position. Mean = {np.mean(plddt):.3f}')
    plt.plot(plddt)
    if Ls is not None:
        L_prev = 0
        for L_i in Ls[:-1]:
            L = L_prev + L_i
            L_prev += L_i
            plt.plot([L, L], [0, 100], color='black')
    plt.ylim(0, 100)
    plt.ylabel('Predicted lDDT')
    plt.xlabel('Positions')
    return plt