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
    d = sigma

    a_cc = 4.0 * np.pi**2 * epsilon * (rho * sigma**3) ** 2

    r_min = 2.0 * d + 1.0e-3
    r_max = 10.0 * sigma
    r = np.linspace(r_min, r_max, 2000)

    # Compute colloidal potential
    u = u_attractive(r, a_cc, d) + u_repulsive(r, a_cc, d, sigma)
    
    # Compute Lennard-Jones potential
    # Use a slightly different r range to avoid singularity at r=0
    r_lj_min = 0.8 * sigma
    r_lj = np.linspace(r_lj_min, r_max, 2000)
    e_lj = lennard_jones(r_lj, epsilon, sigma)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(r, u, color='k', linewidth=3, label=r'$U(r)$ (Colloidal)')
    ax.plot(r_lj, e_lj, color='b', linewidth=3, linestyle='--', label=r'$E(r)$ (LJ)')

    ax.set_xlabel(r'$r$', fontsize=40, labelpad=15)
    ax.set_ylabel(r'$U(r)$', fontsize=40, labelpad=15)

    ax.tick_params(which='both', direction='out', width=3, labelsize=30, pad=10)
    ax.tick_params(which='major', length=12)
    ax.tick_params(which='minor', length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    ax.set_xlim((0.0, r_max))
    ax.set_ylim((-1.0, 1.0))
    
    # Add legend
    ax.legend(fontsize=30, frameon=True, fancybox=False, edgecolor='k', framealpha=1.0, loc='best')

    plt.tight_layout()
    plt.savefig('u_of_r.pdf', bbox_inches='tight', dpi=300)
    plt.savefig('u_of_r.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
