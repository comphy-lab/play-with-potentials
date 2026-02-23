#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


def u_attractive(r, a_cc, d):
    term1 = 2.0 * d**2 / (r**2 - 4.0 * d**2)
    term2 = 2.0 * d**2 / r**2
    term3 = np.log(1.0 - (2.0 * d / r) ** 2)
    return -a_cc / 6.0 * (term1 + term2 + term3)


def u_repulsive(r, a_cc, d, sigma):
    prefactor = a_cc * sigma**6 / (37800.0 * r)
    t1 = (r**2 - 14.0 * r * d + 54.0 * d**2) / (r - 2.0 * d) ** 7
    t2 = (r**2 + 14.0 * r * d + 54.0 * d**2) / (r + 2.0 * d) ** 7
    t3 = (2.0 * r**2 - 60.0 * d**2) / r**7
    return prefactor * (t1 + t2 - t3)


def lennard_jones(r, epsilon, sigma):
    """Lennard-Jones potential: E = 4ε[(σ/r)^12 - (σ/r)^6]"""
    return 4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)


def main():
    sigma = 1.0
    epsilon = 1.0
    
    rho = 1.0
    # Range for d values: from sigma/10 to sigma (avoiding very large d)
    # Large d values would require r > 2d, pushing the valid range beyond our plot window
    d_values = np.logspace(np.log10(sigma / 10.0), np.log10(sigma), 8)
    
    # Use viridis colormap for colloidal potentials
    cmap = plt.get_cmap('viridis')
    colors = [cmap(i / (len(d_values) - 1)) for i in range(len(d_values))]

    # Global r range - needs to accommodate the largest d value
    r_min_global = 0.1 * sigma
    r_max = 10.0 * sigma
    r_global = np.linspace(r_min_global, r_max, 2000)

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot colloidal potentials for different d values
    for i, d in enumerate(d_values):
        a_cc = 4.0 * np.pi**2 * epsilon * (rho * sigma**3) ** 2
        
        # Only compute where r > 2d + small buffer to avoid singularity
        r_cutoff = 2.0 * d + 1.0e-2
        valid_mask = r_global > r_cutoff
        r = r_global[valid_mask]
        
        # Compute colloidal potential only for valid r values
        u = u_attractive(r, a_cc, d) + u_repulsive(r, a_cc, d, sigma)
        
        ax.plot(r, u, color=colors[i], linewidth=3, 
                label=r'$d = {:.2f}\sigma$'.format(d / sigma))
    
    # Compute Lennard-Jones potential for reference
    r_lj_min = 0.8 * sigma
    r_lj = np.linspace(r_lj_min, r_max, 2000)
    e_lj = lennard_jones(r_lj, epsilon, sigma)
    ax.plot(r_lj, e_lj, color='k', linewidth=3, linestyle='--', label=r'LJ')

    ax.set_xlabel(r'$r/\sigma$', fontsize=40, labelpad=15)
    ax.set_ylabel(r'$U(r)$', fontsize=40, labelpad=15)

    ax.tick_params(which='both', direction='out', width=3, labelsize=30, pad=10)
    ax.tick_params(which='major', length=12)
    ax.tick_params(which='minor', length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    ax.set_xlim((0.0, r_max))
    ax.set_ylim((-1.0, 1.0))
    
    # Add legend with best location
    ax.legend(fontsize=24, frameon=True, fancybox=False, edgecolor='k', 
              framealpha=1.0, loc='best', ncol=2)

    plt.tight_layout()
    plt.savefig('u_of_r.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
