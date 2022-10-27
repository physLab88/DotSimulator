import set
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst
from scipy.stats import beta
import yaml  # to save python objects
import time  # to generate filenames
import warnings  # to supress linalg warnings
import scipy.linalg as linalg
import copy


# ====================== DECLARING CONSTANTS ======================
FILEPATH = "Data/single_dot/train/"


# ========================== VISUAL TOOLS =========================
def pltBeta(a, b, loc=0.0, scale=1.0):
    ''' this function allows you to visualise a distribution function
    (not by using a statistical aproach)'''
    x = np.linspace(loc + 0.001*scale, loc + 0.999*scale, 100)
    plt.title(r"Beta distribution $\alpha$=%s, $\beta$=%s, loc=%s, scale=%s" %
              ('{:.2f}'.format(a), '{:.2f}'.format(b), '{:.2f}'.format(loc), '{:.2f}'.format(scale)))
    plt.plot(x, beta.pdf(x, a, b, loc=loc, scale=scale),
            'r-', lw=3, alpha=0.6, label='beta pdf')
    plt.ylim(bottom=0)
    plt.show()


def plt_current(I, Vg, Vd, title='no name', multi_plt=False):
    """
    Vg: array of Gate voltages, I use the first element as the min Vg value
        and the last as the max value. every other values are unused
    Vd: array of drain voltages, I use the first element as the min Vd value
        and the last as the max value. every other values are unused
    I: Vds current matrix (in ????????????????????)
    title: graph title"""
    plt.title(title)
    plt.imshow(np.abs(I), extent=[Vg[0],Vg[-1],Vd[0],Vd[-1]], aspect='auto', cmap='hot')
    if not multi_plt:
        cbar = plt.colorbar(label='current in ??????')
        plt.xlabel(r'$V_g$ in mV')
        plt.ylabel(r'$V_{ds}$ in mV')
    plt.axhline(0, color='k', alpha=0.3)
    # plt.axvline(0, color='k', alpha=0.3)


def plt_conduct(G, Vg, Vd, title='no name', multi_plt=False):
    """
    Vg: array of Gate voltages, I use the first element as the min Vg value
        and the last as the max value. every other values are unused
    Vd: array of drain voltages, I use the first element as the min Vd value
        and the last as the max value. every other values are unused
    G: Gds conductance matrix (in ????????????????????)
    title: graph title"""
    plt.imshow(G, extent=[Vg[0],Vg[-1],Vd[0],Vd[-1]], aspect='auto', cmap='RdPu')
    if not multi_plt:
        cbar = plt.colorbar(label='conductance in ??????')
        plt.xlabel(r'$V_g$ in mV')
        plt.ylabel(r'$V_{ds}$ in mV')
    plt.axhline(0, color='k', alpha=0.3)
    plt.title(title)
    # plt.axvline(0, color='k', alpha=0.3)


def plt_file(fileIndex, multi_plt=False):
    f = open(FILEPATH + '_data_indexer.yaml', 'r')
    info = yaml.load(f, Loader=yaml.FullLoader)[fileIndex]
    data = np.load(FILEPATH + info['f'] + ".npy")
    Vg = info['Vg_range']
    Vds = info['Vds_range']
    mesure = info['mesure']

    title = ("Ec:" + '{:.2f}'.format(info['Ec']) + "meV  T:" + '{:.2f}'.format(info['T']) + r"K    $\alpha$:" +
             '{:.2f}'.format(info['Cg']/(info['Cg'] + info['Cs'] + info['Cd'])) + "    Cs:" + '{:.2f}'.format(info['Cs']) +
              "aF  Cd:" + '{:.2f}'.format(info['Cd']) + "aF")
    plt.annotate('LV: %s' % ['{:.1f}'.format(level) for level in info['levels']], (15, 20), xycoords='axes pixels', color='b')
    plt.annotate('degen: %s' % info['degens'], (15, 5), xycoords='axes pixels', color='b')
    if mesure == 'I':
        plt_current(data, Vg, Vds, title, multi_plt=multi_plt)
    elif mesure == 'G':
        plt_conduct(data, Vg, Vds, title, multi_plt=multi_plt)
    if not multi_plt:
        plt.show()
    return


# ========================= SUB FUNCTIONS =========================
def build_simulation(Cd, Cs, Cg, levels, degens, Gd=1.0, Gs=1.0):
    '''
    Cd, Cs, Cg : are (in order) drain, source and gate coupling capacitance with the dot (in aF)
    levels : are the energy levels of the dot (in meV)
    degens : are the degenerencies of each levels
    Gd, Gs : unknown at the time

    returns new_set created'''
    new_set = set.SET()

    new_set.add_quantum_dot('dot', list(levels), list(degens))

    # Add components to the dot to form the structure
    new_set.add_lead('source')
    new_set.add_lead('drain')
    new_set.add_gate('gate')
    new_set.add_link('dl', 'dot', 'drain', Cd * 1e-18, Gd)
    new_set.add_link('dl', 'dot', 'source', Cs * 1e-18, Gs)
    new_set.add_link('dg', 'dot', 'gate', Cg * 1e-18)
    return new_set


def simulate(myset, Vg, Vd, T, mesure='I'):
    """
    myset: a set we want to simulate
    Vg: an array of gate voltages (in mV)
    Vd: an array of drain voltages (in mV)
        NOTE: here, we assume Vsource = 0V
    T: temperature in K
    mesure: the quantity to mesure ('I' for current and 'G' for conductance)

    returns a 2D matrix of currents (first indices iterates
        over Vd and second indicies iterates over Vg)
    """
    # initialising simulation
    myset.set_temperature(T)
    myset.pre_processing()

    # running the simulation
    data = []
    for vd in Vd:
        temp = []
        for vg in Vg:
            myset.tunnel_rate([0, vd, vg])
            myset.solver()
            temp.append(myset.current('drain', 'dot'))
        data.append(temp)
    data = np.array(data)
    if mesure == 'G':
        data = np.gradient(data, (Vd[1] - Vd[0])*1E-3, axis=0)
    return data


def randCapGenerator(Ec_dist, g_ratio, snd_ratio):
    """ This function outputs random capacitances out of the
    distributions provided:
    Ec_dist: probability distribution of the charging energy (in meV)
    g_ratio: the Cg/C ratio probability distribution
    snd_ratio: the Cs/(Cs+Cd) ratio. if we assume the source
      and drain are equivalent, it would make sens to make
      this distribution symetrical around 0.5"""
    Ec = Ec_dist()  # meV
    Ctot = cst.e/(Ec*1E-3) * 1E18  # aF
    Cg = Ctot*g_ratio()
    Cs = snd_ratio()*(Ctot-Cg)
    Cd = Ctot - Cg - Cs
    return Cd, Cs, Cg, Ec


def randLevelGenerator():
    """ This function outputs random energy levels out of the
    distributions provided and random level degenerancies:"""
    # TODO add in level compacting at high energies
    rand_level_numb = lambda: np.random.choice(5, 1, p=[0.21, 0.35, 0.27, 0.12, 0.05]) + 1  # arbitrairy probabilities (what seems good)
    n_levels = int(rand_level_numb())
    rand_level_spacing = lambda: beta.rvs(2.3, 2.0, loc=1.0, scale=35)/(n_levels - 1)  # meV
    if n_levels == 1:
        rand_level_degens = lambda: np.random.choice(5, 1,
                                                     p=[0.03, 0.25, 0.37, 0.25, 0.1])  # arbitrairy probabilities
    elif n_levels == 2:
        rand_level_degens = lambda: np.random.choice(6, 1,
                                                     p=[0.1, 0.2, 0.32, 0.2, 0.12, 0.06])  # arbitrairy probabilities
    elif n_levels == 3:
        rand_level_degens = lambda: np.random.choice(6, 1,
                                                     p=[0.3, 0.22, 0.18, 0.14, 0.1, 0.06])  # arbitrairy probabilities
    elif n_levels == 4:
        rand_level_degens = lambda: np.random.choice(6, 1,
                                                     p=[0.3, 0.22, 0.18, 0.14, 0.1, 0.06])  # arbitrairy probabilities
    else:
        rand_level_degens = lambda: np.random.choice(5, 1,
                                                     p=[0.4, 0.3, 0.17, 0.1, 0.03])  # arbitrairy probabilities

    levels = [0.0]
    degens = [1]
    for i in range(n_levels - 1):
        levels.append(float(rand_level_spacing() + levels[-1]))
        degens.append(1)

    # add random degens
    add_degens = int(rand_level_degens())
    for i in range(add_degens):
        degens[int(np.random.choice(n_levels))] += 1
    print(levels)
    print(degens)
    return levels, degens


def generation_loop(n, T_dist, Ec_dist, g_ratio, snd_ratio, Vg_range, nVg, Vds_range, nVds, mesure='I'):
    global dopants
    try:
        f = open(FILEPATH + '_data_indexer.yaml', 'r')
        data = yaml.load(f, Loader=yaml.FullLoader)
    except IOError:
        f = open(FILEPATH + '_data_indexer.yaml', 'w')  # creating a file if it does not exist
        data = []
    warnings.filterwarnings(action='ignore', category=linalg.LinAlgWarning)  # manually surpressing linalg warnings

    for i in range(n):
        print("GENERATING SAMPLE # %s" % i)
        Cd, Cs, Cg, Ec = randCapGenerator(Ec_dist, g_ratio, snd_ratio)
        Cd, Cs, Cg, Ec = float(Cd), float(Cs), float(Cg), float(Ec)
        print("Capacities: %s" % [Cd, Cs, Cg])
        T = float(T_dist())
        ID = str(time.time()).replace('.', 's')

        # simulation
        levels, degens = randLevelGenerator()
        set1 = build_simulation(Cd, Cs, Cg, levels, degens)
        Vg = np.linspace(Vg_range[0], Vg_range[1], nVg)
        Vd = np.linspace(Vds_range[0], Vds_range[1], nVds)
        diagram = simulate(set1, Vg, Vd, T, mesure=mesure)

        # saving data
        print("saving...")
        temp = {'f': ID,
                'Ec': Ec,
                'Cg': Cg,
                'Cd': Cd,
                'Cs': Cs,
                'T': T,
                'levels': levels,
                'degens': degens,
                'Vds_range': Vds_range,
                'Vg_range': Vg_range,
                'nVds': nVds,
                'nVg': nVg,
                'mesure': mesure,
                }
        data.append(temp)
        f = open(FILEPATH + '_data_indexer.yaml', 'w')
        np.save(FILEPATH + ID + '.npy', diagram)
        yaml.dump(data, f)
    return


def generateFunction(n, mesure='I'):
    global dopants

    g_ratio = lambda: beta.rvs(1.2, 1.6, loc=0.10, scale=0.70)  # aF
    snd_ratio = lambda: beta.rvs(2, 2, loc=0.15, scale=0.7)  # aF
    Ec_dist = lambda: beta.rvs(2, 1.7, loc=12, scale=55)  # meV
    T_dist = lambda: beta.rvs(1.8, 2.1, loc=1.5, scale=20)  # K

    generation_loop(n, T_dist, Ec_dist, g_ratio, snd_ratio, [-10, 290], 100, [-70, 70], 100, mesure=mesure)


# =========================== MAIN ===========================
def main():
    pltBeta(2, 2, loc=0.15, scale=0.7)
    num = 5
    generateFunction(num, mesure='I')
    for i in range(num):
        plt_file(-(i+1))


if __name__ == '__main__':
    main()




