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
    
    # Create figure with 1 row, 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
    
    # Use viridis colormap for both plots
    cmap = plt.get_cmap('viridis')
    
    # Global r range
    r_min_global = 0.1 * sigma
    r_max = 10.0 * sigma
    r_global = np.linspace(r_min_global, r_max, 2000)
    
    # ========== Panel 1: Varying d ==========
    rho = 1.0
    d_values = np.logspace(np.log10(sigma / 10.0), np.log10(sigma), 8)
    colors_d = [cmap(i / (len(d_values) - 1)) for i in range(len(d_values))]
    
    for i, d in enumerate(d_values):
        a_cc = 4.0 * np.pi**2 * epsilon * (rho * sigma**3) ** 2
        
        # Only compute where r > 2d + small buffer to avoid singularity
        r_cutoff = 2.0 * d + 1.0e-2
        valid_mask = r_global > r_cutoff
        r = r_global[valid_mask]
        
        # Compute colloidal potential only for valid r values
        u = u_attractive(r, a_cc, d) + u_repulsive(r, a_cc, d, sigma)
        
        ax1.plot(r, u, color=colors_d[i], linewidth=3, 
                 label=r'$d = {:.2f}\sigma$'.format(d / sigma))
    
    # Compute Lennard-Jones potential for reference
    r_lj_min = 0.8 * sigma
    r_lj = np.linspace(r_lj_min, r_max, 2000)
    e_lj = lennard_jones(r_lj, epsilon, sigma)
    ax1.plot(r_lj, e_lj, color='k', linewidth=3, linestyle='--', label=r'LJ')
    
    ax1.set_xlabel(r'$r/\sigma$', fontsize=40, labelpad=15)
    ax1.set_ylabel(r'$U(r)$', fontsize=40, labelpad=15)
    ax1.tick_params(which='both', direction='out', width=3, labelsize=30, pad=10)
    ax1.tick_params(which='major', length=12)
    ax1.tick_params(which='minor', length=6)
    for spine in ax1.spines.values():
        spine.set_linewidth(3)
    ax1.minorticks_on()
    ax1.set_xlim((0.0, r_max))
    ax1.set_ylim((-1.0, 1.0))
    ax1.legend(fontsize=20, frameon=True, fancybox=False, edgecolor='k', 
               framealpha=1.0, loc='best', ncol=2)
    ax1.set_title(r'Varying $d$ ($\rho = 1.0$)', fontsize=40, pad=20)
    
    # ========== Panel 2: Varying rho ==========
    d = 0.25 * sigma
    rho_values = np.logspace(np.log10(10.0), np.log10(20.0), 8)
    colors_rho = [cmap(i / (len(rho_values) - 1)) for i in range(len(rho_values))]
    
    r_cutoff = 2.0 * d + 1.0e-2
    valid_mask = r_global > r_cutoff
    r = r_global[valid_mask]
    
    for i, rho in enumerate(rho_values):
        a_cc = 4.0 * np.pi**2 * epsilon * (rho * sigma**3) ** 2
        
        # Compute colloidal potential
        u = u_attractive(r, a_cc, d) + u_repulsive(r, a_cc, d, sigma)
        
        ax2.plot(r, u, color=colors_rho[i], linewidth=3, 
                 label=r'$\rho = {:.2f}$'.format(rho))
    
    # Add Lennard-Jones reference
    ax2.plot(r_lj, e_lj, color='k', linewidth=3, linestyle='--', label=r'LJ')
    
    ax2.set_xlabel(r'$r/\sigma$', fontsize=40, labelpad=15)
    ax2.set_ylabel(r'$U(r)$', fontsize=40, labelpad=15)
    ax2.tick_params(which='both', direction='out', width=3, labelsize=30, pad=10)
    ax2.tick_params(which='major', length=12)
    ax2.tick_params(which='minor', length=6)
    for spine in ax2.spines.values():
        spine.set_linewidth(3)
    ax2.minorticks_on()
    ax2.set_xlim((0.0, r_max))
    ax2.set_ylim((-1.0, 1.0))
    ax2.legend(fontsize=20, frameon=True, fancybox=False, edgecolor='k', 
               framealpha=1.0, loc='best', ncol=2)
    ax2.set_title(r'Varying $\rho$ ($d = 0.25\sigma$)', fontsize=40, pad=20)

    plt.tight_layout()
    plt.savefig('u_of_r.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
