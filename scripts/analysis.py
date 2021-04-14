import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os.path
from os import path
import decimal

TC_EXACT = 0.89294

params = {
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 10,
    'legend.fontsize': 10,
    'axes.labelsize': 10,
    'xtick.labelsize':10,
    'ytick.labelsize':10,
    'lines.linewidth':1,
    "patch.edgecolor": "black"
}

plt.rcParams.update(params)
plt.style.use('seaborn-deep')

def observables_json(L, beta, seed=1234):
    # template: L=4_beta=0.9_dtheta=0.01_seed=1234_observables
    if int(beta) == beta:
        beta == int(beta)
        beta = int(beta)

    dtheta = dtheta_from_beta(beta)
    path = '../examples/data/L={}_beta={}_dtheta={}_seed={}_observables.json'.format(L, beta, dtheta, seed)
    return path

def dtheta_from_beta(beta):
    t = (beta**(-1) - TC_EXACT) / TC_EXACT
    '''
    if -1 <= t <= -0.5:
        dtheta = 0.1
    elif -0.5 < t <= 0.5:
        dtheta = 0.2
    elif 0.5 < t <= 1.0:
        dtheta = 0.5
    else:
        dtheta = 0.8
    '''
    if -1 <= t <= -0.5:
        dtheta = 0.05
    elif -0.5 < t <= 0.5:
        dtheta = 0.1
    elif 0.5 < t <= 1.0:
        dtheta = 0.2
    else:
        dtheta = 0.5

    return dtheta

def colors_from_Ls(Ls, c="Blues"):
    cmap = plt.get_cmap(c)
    colors = [cmap(i) for i in np.linspace(0.3,1,len(Ls))]
    return colors
    
def plots_nocollapse(c="Blues"):
    #betas = np.arange(0.9, 1.355, 0.005)\
    betas = np.arange(0.8, 3.0, 0.1)
    betas *= 100
    betas = np.round(betas)/100
    T = betas**(-1)
    #t = (T - Tc_exact) / TC_EXACT

    #Ls = [10, 12, 16, 20, 24, 30, 40, 50]
    Ls = [4,10]
    colors = colors_from_Ls(Ls, c=c)

    rhos = np.zeros((len(Ls), len(betas)))
    energies = np.zeros((len(Ls), len(betas)))
    specific_heats = np.zeros((len(Ls), len(betas)))

    fig, ax = plt.subplots(3,1, sharex=True)
    for (j,L) in enumerate(Ls):
        for (i,beta) in enumerate(betas):
            path = observables_json(L, beta)
            with open(path) as f:
                contents = json.load(f)
            
            rhos[j,i] = contents["spin_stiffness"]["mean"]
            energies[j,i] = contents["energy"]["mean"]
            energy_sqr = contents["sqr_energy"]["mean"]
            specific_heats[j,i] = (energy_sqr - energies[j,i]**2)*(beta**2.)*(L**4)
        
        ax[0].plot(T, energies[j,:], ls='--', label=r'$L = {}$'.format(L), color=colors[j], marker='o')
        ax[1].plot(T, specific_heats[j,:], ls='--', color=colors[j], marker='o')
        ax[2].plot(T, rhos[j,:], ls='--', color=colors[j], marker='o')

    ax[0].set_ylabel(r'$\frac{\langle E \rangle}{N}$', rotation=0)
    ax[0].legend()

    ax[1].set_ylabel(r'$c_v$', rotation=0)

    #ax[2].plot(T, 2*T/np.pi, color='black')
    ax[2].vlines(x=TC_EXACT, ymax = 0, ymin = 0.5, color='red')
    ax[2].set_ylabel(r'$\rho_s$', rotation=0)
    ax[2].set_xlabel(r'$T$')
    
    fig.tight_layout()
    fig.savefig("../examples/figures/observables_nocollapse.pdf", dpi=500, bbox_inches='tight')

plots_nocollapse()