'''
Common functions
'''

def plot_corr_matrix(pcov):
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Compute correlation matrix
    perr = np.sqrt(np.diag(pcov))
    corr_matrix = pcov / np.outer(perr, perr)
    
    # Plot
    fig, ax = plt.subplots()
    im = ax.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax, label='Correlation coefficient')
    
    # Label axes
    param_names = ["t0","f_base","delta","t_b","t_w"]
    ax.set_xticks(np.arange(len(param_names)))
    ax.set_yticks(np.arange(len(param_names)))
    ax.set_xticklabels(param_names)
    ax.set_yticklabels(param_names)
    plt.title("Correlation Matrix of Fit Parameters")
    plt.show()
