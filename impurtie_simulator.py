'''
date: 2022-09-15
Author: Michael Bedard

intro: This code was created to simulate single quantum dots of diffrent
types of atoms. it is rather basic
'''


import set
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst
import warnings
import scipy.linalg as linalg


# build simulation (levels, caps) returns simulation set
def build_simulation(Cd, Cs, Cg, levels, degens, Gd=1.0, Gs=1.0):
    '''
    Cd, Cs, Cg : are (in order) drain, source and gate coupling capacitance with the dot (in aF)
    levels : are the energy levels of the dot (in mV)
    degens : are the degenerencies of each levels
    Gd, Gs : unknown at the time

    returns new_set created'''
    new_set = set.SET()

    new_set.add_quantum_dot('dot', list(levels), degens)

    # Add components to the dot to form the structure
    new_set.add_lead('source')
    new_set.add_lead('drain')
    new_set.add_gate('gate')
    new_set.add_link('dl', 'dot', 'drain', Cd * 1e-18, Gd)
    new_set.add_link('dl', 'dot', 'source', Cs * 1e-18, Gs)
    new_set.add_link('dg', 'dot', 'gate', Cg * 1e-18)
    return new_set


def build_simulation2(Cd, Cs, Cg, max_nb_e, offset=0, Gd=1.0, Gs=1.0):
    '''
    builds a simulation with a metalic dot
    Cd, Cs, Cg : are (in order) drain, source and gate coupling capacitance with the dot (in aF)
    max_nb_e : maximum number of electrons on the dot
    offset : an offset energy on the dot
    Gd, Gs : unknown at the time

    returns new_set created'''
    new_set = set.SET()

    new_set.add_metallic_dot('dot', max_nb_e, 0, offset)

    # Add components to the dot to form the structure
    new_set.add_lead('source')
    new_set.add_lead('drain')
    new_set.add_gate('gate')
    new_set.add_link('dl', 'dot', 'drain', Cd * 1e-18, Gd)
    new_set.add_link('dl', 'dot', 'source', Cs * 1e-18, Gs)
    new_set.add_link('dg', 'dot', 'gate', Cg * 1e-18)
    return new_set


# run simulation (voltages (Vs=0), simulation set, T) returns currents
def simulate_current(myset, Vg, Vd, T):
    """
    myset: a set we want to simulate
    Vg: an array of gate voltages (in mV)
    Vd: an array of drain voltages (in mV)
        NOTE: here, we assume Vsource = 0V
    T: temperature in K

    returns a 2D matrix of currents (first indices iterates
        over Vd and second indicies iterates over Vg)
    """
    # initialising simulation
    myset.set_temperature(T)
    myset.pre_processing()

    # running the simulation
    I = []
    for vd in Vd:
        temp = []
        for vg in Vg:
            myset.tunnel_rate([0, vd, vg])
            myset.solver()
            temp.append(myset.current('drain', 'dot'))
        I.append(temp)
    I = np.array(I)
    return I


# plt data (voltages, current, title)
def plt_current(Vg, Vd, I, title='no name'):
    """
    Vg: array of Gate voltages, I use the first element as the min Vg value
        and the last as the max value. every other values are unused
    Vd: array of drain voltages, I use the first element as the min Vd value
        and the last as the max value. every other values are unused
    I: Vds current matrix (in ????????????????????)
    title: graph title"""
    plt.title(title)
    plt.imshow(I, extent=[Vg[0],Vg[-1],Vd[0],Vd[-1]], aspect='auto')
    cbar = plt.colorbar(label='current in ??????')
    #cbar.set_label('current in ??????', rotation=270)
    plt.xlabel(r'$V_g$ in mV')
    plt.ylabel(r'$V_{ds}$ in mV')
    plt.axhline(0, color='k', alpha=0.3)
    #plt.axvline(0, color='k', alpha=0.3)


def main():
    levels = np.array([0.0, 10.0, 15.0, 17.0, 18.0])

    # fig, axs = plt.subplots(3, 1)
    set1 = build_simulation(0.86, 0.87, 3.52, levels, [1, 1, 1, 1, 1])
    # set1 = build_simulation2(0.86, 0.87, 3.52, 10, -10)
    Vg = np.linspace(-50, 350, 100)
    Vd = np.linspace(-70, 70, 100)
    warnings.filterwarnings(action='ignore', category=linalg.LinAlgWarning)  # manually supressing linalg warnings
    I1 = simulate_current(set1, Vg, Vd, 40)
    # plt.sca(axs[0])
    plt_current(Vg, Vd, I1, 'Simulated stability diagrams')

    plt.tight_layout()
    plt.show()
    return


if __name__ == '__main__':
    main()




