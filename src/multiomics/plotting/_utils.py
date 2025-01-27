import pandas as pd
from matplotlib import pyplot as plt
from adjustText import adjust_text


def run_adjust_text(x, y, labels, ax=None, use_arrow=True, font_weight='bold', font_size=8):
    texts = [
        plt.text(
            x[i], y[i], 
            labels[i],
            fontdict={'weight': font_weight, 'size': font_size},
            ha='center', va='center'
        ) for i in range(len(x))
    ]
    
    if use_arrow:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax = ax)
    else:
        adjust_text(texts, ax = ax)