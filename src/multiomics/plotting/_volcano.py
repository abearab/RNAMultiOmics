import numpy as np
import matplotlib.pyplot as plt
from ._utils import run_adjust_text


def plot_volcano(df, ax=None, 
                 up_threshold=1, down_threshold=-1, 
                 delta_label='log2FoldChange', 
                 pv_label='pvalue', 
                 pv_threshold=0.01, 
                 label_list=None, 
                 label_column='gene_name', 
                 xlabel='auto', ylabel='auto',
                 xlims="auto", ylims="auto",
                 dot_size=10, 
                 label_fontsize=8,
                 figsize=(10, 8)):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = df.copy()

    if xlabel == 'auto': xlabel = delta_label
    if ylabel == 'auto': ylabel=f'-log10({pv_label})'

    # Calculate -log10(pv)
    df[f'-log10({pv_label})'] = -np.log10(df[pv_label])
    
    # Define upregulated and downregulated genes
    df['label'] = 'non-significant'
    df.loc[(df[delta_label] > up_threshold) & (df[pv_label] < pv_threshold), 'label'] = 'up'
    df.loc[(df[delta_label] < down_threshold) & (df[pv_label] < pv_threshold), 'label'] = 'down'
    
    # Plot non-significant points
    non_sig = df[df['label'] == 'non-significant']
    ax.scatter(non_sig[delta_label], non_sig[f'-log10({pv_label})'], color='grey', s=dot_size, alpha=0.5, label='Non-significant')
    
    # Plot upregulated points
    up = df[df['label'] == 'up']
    ax.scatter(up[delta_label], up[f'-log10({pv_label})'], color='red', s=dot_size, label='Upregulated')
    
    # Plot downregulated points
    down = df[df['label'] == 'down']
    ax.scatter(down[delta_label], down[f'-log10({pv_label})'], color='blue', s=dot_size, label='Downregulated')
    
    # Set labels and limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Determine x-axis and y-axis limits if set to "auto"
    if xlims == "auto":
        xlims = (df[delta_label].min() - 1, df[delta_label].max() + 1)
    if ylims == "auto":
        ylims = (0, df[f'-log10({pv_label})'].max() + 1)

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
    # Add legend
    ax.legend()


    if label_list == "up":
        for i, row in df.query('label == "up"').iterrows():
            run_adjust_text(row[delta_label], row[f'-log10({pv_label})'], row[label_column], ax, font_size=label_fontsize)
    elif label_list == "down":
        for i, row in df.query('label == "down"').iterrows():
            run_adjust_text(row[delta_label], row[f'-log10({pv_label})'], row[label_column], ax, font_size=label_fontsize)
    elif label_list == "both":
        for i, row in df.query('label in ["down","up"]').iterrows():
            run_adjust_text(row[delta_label], row[f'-log10({pv_label})'], row[label_column], ax, font_size=label_fontsize)
    elif type(label_list) is list:
        for i, row in df.query(f'{label_column} in {label_list}').iterrows():
            run_adjust_text(row[delta_label], row[f'-log10({pv_label})'], row[label_column], ax, font_size=label_fontsize)

    return ax