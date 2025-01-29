import pandas as pd
from matplotlib import pyplot as plt
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
