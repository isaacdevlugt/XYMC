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
    
def plots(c="Blues"):
    betas_coarse = np.array([0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 1.4, 1.45, 1.5, 2.0, 3.0])
    betas1 = np.arange(0.9, 1.105, 0.005)
    betas2 = np.arange(1.105, 1.2, 0.005)
    betas3 = np.arange(1.201, 1.351, 0.005)
    betas = np.sort(np.concatenate([betas_coarse, betas1, betas2, betas3]))
    betas *= 1000
    betas = np.round(betas)/1000
    T = betas**(-1)
    #t = (T - Tc_exact) / TC_EXACT

    #Ls = [10, 12, 16, 20, 24, 30, 40, 50]
    Ls = [12, 16, 20, 24]
    colors = colors_from_Ls(Ls, c=c)

    rhos = np.zeros((len(Ls), len(betas)))
    energies = np.zeros((len(Ls), len(betas)))
    specific_heats = np.zeros((len(Ls), len(betas)))
    T_KT_estimates = np.zeros(len(Ls))
    inv_log_sqr_Ls = (1./ np.log(Ls))**2.
    rhos_linear_intercept_idx = np.zeros(len(Ls), dtype=int)

    fig, ax1 = plt.subplots(2,1)
    fig.subplots_adjust(hspace=0.3)
    left, bottom, width, height = [0.6, 0.65, 0.22, 0.22]
    ax2 = fig.add_axes([left, bottom, width, height])
    for (j,L) in enumerate(Ls):
        for (i,beta) in enumerate(betas):
            path = observables_json(L, beta)
            with open(path) as f:
                contents = json.load(f)
            
            rhos[j,i] = contents["spin_stiffness"]["mean"]

        reduced_rhos = abs(rhos[j,:] - 2*T/np.pi)
        T_KT_estimates[j] = T[np.where(reduced_rhos == np.amin(reduced_rhos))]
        rhos_linear_intercept_idx[j] = int(np.where(reduced_rhos == np.amin(reduced_rhos))[0])
        print("For L = ",L,": T_KT = ", T_KT_estimates[j])

        ax1[0].plot(T, rhos[j,:], label=r'$L = {}$'.format(L), color=colors[j])
        ax1[1].scatter(inv_log_sqr_Ls[j], T_KT_estimates[j], color=colors[j])
        ax2.plot(T, rhos[j,:], label=r'$L = {}$'.format(L), color=colors[j])
        ax2.scatter(T_KT_estimates[j], rhos[j, rhos_linear_intercept_idx[j]], marker='x', color='red')

    ax1[0].axvline(x=TC_EXACT, ymin=0, ymax=1, ls='--', label=r'$T_{\mathrm{KT}, \mathrm{exact}}$', color='grey')
    ax1[0].plot(T, 2*T/np.pi, color='orange', label=r'$\rho_s = \frac{2}{\pi}T$', ls='dotted')
    
    ax1[0].set_ylabel(r'$\rho_s$', rotation=0)
    ax1[0].set_xlabel(r'$T$')
    ax1[0].set_xlim(0.7, 1.5)
    ax1[0].set_ylim(0, 1.0)
    ax1[0].yaxis.set_label_coords(-0.1,0.5)

    ax1[1].set_ylabel(r'$T_{\mathrm{KT}}(L)$', rotation=0)
    ax1[1].set_xlabel(r'$(\ln(L))^{-2}$')
    ax1[1].yaxis.set_label_coords(-0.12,0.5)

    p, V = np.polyfit(inv_log_sqr_Ls, T_KT_estimates, 1, cov=True)
    print("T_KT numerical = ", p[1]," +\- ",np.sqrt(V[1][1]))

    ax1[0].axvline(x=p[1], ymin=-1, ymax=1, label=r'$T_{\mathrm{KT}, \mathrm{MC}}$', color='purple')
    ax1[0].axvline(x=p[1]+np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)
    ax1[0].axvline(x=p[1]-np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)

    ax1[1].plot(inv_log_sqr_Ls, p[0]*inv_log_sqr_Ls + p[1], color = 'black')
    ax2.plot(T, 2*T/np.pi, color='orange', label=r'$\rho_s = \frac{2}{\pi}T$', ls='dotted')
    ax2.set_xlim(0.91, 0.95)
    ax2.set_ylim(0.575, 0.625)
    ax1[0].legend(loc=(1.01, 0.1))
    
    fig.savefig("../examples/figures/rhos.pdf", dpi=500, bbox_inches='tight')

    fig, ax = plt.subplots(2,1,sharex=True)
    for (j,L) in enumerate(Ls):
        for (i,beta) in enumerate(betas):
            path = observables_json(L, beta)
            with open(path) as f:
                contents = json.load(f)            
            energies[j,i] = contents["energy"]["mean"]
            energy_sqr = contents["sqr_energy"]["mean"]
            specific_heats[j,i] = (energy_sqr - energies[j,i]**2)*(beta**2.)*(L**4)
        
        ax[0].plot(T, energies[j,:], label=r'$L = {}$'.format(L), color=colors[j])
        ax[1].plot(T, specific_heats[j,:], label=r'$L = {}$'.format(L), color=colors[j])

    ax[0].axvline(x=TC_EXACT, ymin=-1, ymax=1, ls='--', label=r'$T_{\mathrm{KT}, \mathrm{exact}}$', color='grey')
    ax[0].set_ylabel(r'$\frac{\langle E \rangle}{N}$', rotation=0)
    ax[0].yaxis.set_label_coords(-0.1,0.5)
    ax[0].set_xlim(0.6, 1.6)
    ax[0].axvline(x=p[1], ymin=-1, ymax=1, label=r'$T_{\mathrm{KT}, \mathrm{MC}}$', color='purple')
    ax[0].axvline(x=p[1]+np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)
    ax[0].axvline(x=p[1]-np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)
    ax[0].legend()

    ax[1].axvline(x=TC_EXACT, ymin=-1, ymax=1, ls='--', label=r'$T_{\mathrm{KT}, \mathrm{exact}}$', color='grey')
    ax[1].axvline(x=p[1], ymin=-1, ymax=1, label=r'$T_{\mathrm{KT}, \mathrm{MC}}$', color='purple')
    ax[1].axvline(x=p[1]+np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)
    ax[1].axvline(x=p[1]-np.sqrt(V[1][1]), ymin=-1, ymax=1, color='purple', alpha=0.5)
    ax[1].set_ylabel(r'$c_v$', rotation=0)
    ax[1].yaxis.set_label_coords(-0.1,0.5)
    ax[1].set_xlim(0.6, 1.6)
    ax[1].set_xlabel(r'$T$')

    plt.savefig("../examples/figures/energy_cv.pdf", dpi=500, bbox_inches='tight')

plots()