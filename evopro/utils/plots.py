import os
from argparse import ArgumentParser
import numpy as np
from string import ascii_uppercase, ascii_lowercase
import matplotlib.pyplot as plt

from utils import utils


alphabet_list = list(ascii_uppercase+ascii_lowercase)


def plot_ticks(Ls):
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color='black')
        plt.plot([L, L], [0, Ln], color='black')
    ticks = np.cumsum([0]+Ls)
    ticks = (ticks[1:] + ticks[:-1])/2
    plt.yticks(ticks, alphabet_list[:len(ticks)])


def plot_pae(pae, Ls=None, dpi=300, fig=True, ptm=None, iptm=None):
    if fig: plt.figure(figsize=(15,10), dpi=dpi)
    pae_title = 'Predicted Aligned Error.'
    if ptm:
        pae_title += f' ptm = {ptm:.3f}.'
    if iptm:
        pae_title += f' iptm = {iptm:.3f}.'
    plt.title(pae_title)
    Ln = pae.shape[0]
    plt.imshow(pae, cmap='bwr', vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
    if Ls is not None and len(Ls) > 1: plot_ticks(Ls)
    plt.colorbar()
    plt.ylabel('Positions')
    plt.xlabel('Positions')
    return plt


def plot_plddt(plddt, Ls=None, dpi=100, fig=True):
    if fig: plt.figure(figsize=(15,10), dpi=300)
    plt.title(f'Predicted lDDT per position. Mean = {np.mean(plddt):.3f}')
    plt.plot(plddt)
    if Ls is not None:
        L_prev = 0
        for L_i in Ls[:-1]:
            L = L_prev + L_i
            L_prev += L_i
            plt.plot([L, L], [0, 100], color='black')
    plt.ylim(0, 100)
    plt.ylabel('Predicted lDDT')
    plt.xlabel('Positions')
    return plt


def get_chain_lengths(results_dict):
    protein = results_dict['unrelaxed_protein']
    chain_index = protein.chain_index
    residue_index = protein.residue_index

    Ls = []
    chain_idx = 0
    for res_idx in range(len(chain_index)):
        if chain_index[res_idx] == chain_idx:
            if len(Ls) < chain_idx + 1:
                Ls.append(0)
        else:
            Ls.append(0)
            chain_idx += 1

        Ls[chain_idx] += 1

    # Double check that we didn't run a monomer model on a multimer.
    if len(Ls) == 1:
        # We ran a monomer.
        if residue_index[-1] == len(residue_index):
            pass
        # We ran a multimer.
        else:
            # Create new container.
            Ls = []

            # Determine where chain breaks occur.
            chain_breaks = []
            adder = 0
            x = iter(range(len(residue_index)))
            for res_idx in range(len(residue_index)):
                differ = abs(next(x) + 1 + adder - residue_index[res_idx])
                if differ != 0:
                    adder += differ
                    chain_breaks.append(res_idx)

            # Determine chain lengths.
            remainder = len(residue_index)
            for breaks in chain_breaks:
                if remainder >= 0:
                    Ls.append(breaks)
                    remainder -= breaks
            if remainder >= 0:
                Ls.append(remainder)

    return Ls
