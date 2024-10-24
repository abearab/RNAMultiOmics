import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text


def run_adjust_text(x, y, label, ax=None, use_arrow=False, font_weight='bold', font_size=8):
    texts = [
        plt.text(
            x, y, 
            label,
            fontdict={'weight': font_weight, 'size': font_size},
            ha='center', va='center'
        )
    ]
    
    if use_arrow:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax = ax)
    else:
        adjust_text(texts, ax = ax,force_pull=(0.1,0.1))

def plot_volcano(df, ax=None, 
                 up_threshold=1, down_threshold=-1, 
                 delta_label='log2FoldChange', 
                 pval_label='pvalue', 
                 xlabel='log2(Fold Change)',
                 ylabel='-log10(pvalue)',
                 pvalue_threshold=0.01, 
                 label_list=None, 
                 label_column='gene_name', label_fontsize=8,
                 dot_size=10, xlims="auto", ylims="auto",figsize=(10, 8)):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = df.copy()
    # Calculate -log10(pvalue)
    df['-log10(pvalue)'] = -np.log10(df[pval_label])
    
    # Define upregulated and downregulated genes
    df['label'] = 'non-significant'
    df.loc[(df[delta_label] > up_threshold) & (df[pval_label] < pvalue_threshold), 'label'] = 'up'
    df.loc[(df[delta_label] < down_threshold) & (df[pval_label] < pvalue_threshold), 'label'] = 'down'
    
    # Plot non-significant points
    non_sig = df[df['label'] == 'non-significant']
    ax.scatter(non_sig[delta_label], non_sig['-log10(pvalue)'], color='grey', s=dot_size, alpha=0.5, label='Non-significant')
    
    # Plot upregulated points
    up = df[df['label'] == 'up']
    ax.scatter(up[delta_label], up['-log10(pvalue)'], color='red', s=dot_size, label='Upregulated')
    
    # Plot downregulated points
    down = df[df['label'] == 'down']
    ax.scatter(down[delta_label], down['-log10(pvalue)'], color='blue', s=dot_size, label='Downregulated')
    
    # Set labels and limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Determine x-axis and y-axis limits if set to "auto"
    if xlims == "auto":
        xlims = (df[delta_label].min() - 1, df[delta_label].max() + 1)
    if ylims == "auto":
        ylims = (0, df['-log10(pvalue)'].max() + 1)

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
    # Add legend
    ax.legend()


    if label_list == "up":
        for i, row in df.query('label == "up"').iterrows():
            run_adjust_text(row[delta_label], row['-log10(pvalue)'], row[label_column], ax, font_size=label_fontsize)
    elif label_list == "down":
        for i, row in df.query('label == "down"').iterrows():
            run_adjust_text(row[delta_label], row['-log10(pvalue)'], row[label_column], ax, font_size=label_fontsize)
    elif label_list == "both":
        for i, row in df.query('label in ["down","up"]').iterrows():
            run_adjust_text(row[delta_label], row['-log10(pvalue)'], row[label_column], ax, font_size=label_fontsize)
    elif type(label_list) is list:
        for i, row in df.query(f'{label_column} in {label_list}').iterrows():
            run_adjust_text(row[delta_label], row['-log10(pvalue)'], row[label_column], ax, font_size=label_fontsize)

    return ax