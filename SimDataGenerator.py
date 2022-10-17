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
FILEPATH = "Data/sim_data/"

dopants = []
dopants.append({
    'type': 'P',
    'levels': [-45.59, -32.58, -33.89],  # mV
    'degen': [1, 2, 3],
    'capGenerator': lambda: None,  # here, put cap generation function as lambda
    })
dopants.append({
    'type': 'As',
    'levels': [-53.76, -31.26, -32.67],  # mV
    'degen': [1, 2, 3],
    'capGenerator': lambda: None,  # here, put cap generation function as lambda
    })
dopants.append({
    'type': 'B',
    'levels': [-15.33, -11.18, -7.36],  # mV
    'degen': [1, 2, 3],
    'capGenerator': lambda: None,  # here, put cap generation function as lambda
    })


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
    plt.xlabel(r'$V_g$ in mV')
    plt.ylabel(r'$V_{ds}$ in mV')
    plt.axhline(0, color='k', alpha=0.3)
    # plt.axvline(0, color='k', alpha=0.3)


def plt_file(fileIndex):
    f = open(FILEPATH + '_data_indexer.yaml', 'r')
    info = yaml.load(f, Loader=yaml.FullLoader)[fileIndex]
    data = np.load(FILEPATH + info['f'] + ".npy")
    Vg = info['Vg_range']
    Vds = info['Vds_range']

    plt.title("type:" + info['type'] + "  T:" + '{:.2f}'.format(info['T']) + "K  Index:" + str(fileIndex) +
              "  Cg:" + '{:.2f}'.format(info['Cg']) + "aF  Cs:" + '{:.2f}'.format(info['Cs']) +
              "aF  Cd:" + '{:.2f}'.format(info['Cd']) + "aF")
    plt.imshow(np.abs(data), extent=[Vg[0], Vg[-1], Vds[0], Vds[-1]], aspect='auto')
    cbar = plt.colorbar(label='current in ??????')
    plt.xlabel(r'$V_g$ in mV')
    plt.ylabel(r'$V_{ds}$ in mV')
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

    new_set.add_quantum_dot('dot', list(levels), degens)

    # Add components to the dot to form the structure
    new_set.add_lead('source')
    new_set.add_lead('drain')
    new_set.add_gate('gate')
    new_set.add_link('dl', 'dot', 'drain', Cd * 1e-18, Gd)
    new_set.add_link('dl', 'dot', 'source', Cs * 1e-18, Gs)
    new_set.add_link('dg', 'dot', 'gate', Cg * 1e-18)
    return new_set


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


def randCapGenerator(C_dist, g_ratio, snd_ratio):
    """ This function outputs random capacitances out of the
    distributions provided:
    C_dist: probability distribution of the total capacitance (in aF)
    g_ratio: the Cg/C ratio probability distribution
    snd_ratio: the Cs/(Cs+Cd) ratio. if we assume the source
      and drain are equivalent, it would make sens to make
      this distribution symetrical around 0.5"""
    Ctot = C_dist()
    Cg = Ctot*g_ratio()
    Cs = snd_ratio()*(Ctot-Cg)
    Cd = Ctot - Cg - Cs
    return Cd, Cs, Cg


def generation_loop(n, dop_dist, T_dist, Vg_range, nVg, Vds_range, nVds):
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
        dop_label = dop_dist()
        dop = dopants[dop_label]
        print("Dopant type: %s" % dop['type'])
        Cd, Cs, Cg = dop['capGenerator']()
        Cd, Cs, Cg = float(Cd), float(Cs), float(Cg)
        T = float(T_dist())
        ID = dop['type'] + '_' + str(time.time()).replace('.', 's')

        # simulation
        set1 = build_simulation(Cd, Cs, Cg, dop['levels'], dop['degen'])
        Vg = np.linspace(Vg_range[0], Vg_range[1], nVg)
        Vd = np.linspace(Vds_range[0], Vds_range[1], nVds)
        I = simulate_current(set1, Vg, Vd, T)

        # saving data
        print("saving...")
        temp = {'f': ID,
                'label': dop_label,  # number
                'type': dop['type'],  # string
                'Cg': Cg,
                'Cd': Cd,
                'Cs': Cs,
                'T': T,
                'Vds_range': Vds_range,
                'Vg_range': Vg_range,
                'nVds': nVds,
                'nVg': nVg,
                'mesure': 'I',
                }
        data.append(temp)
        f = open(FILEPATH + '_data_indexer.yaml', 'w')
        np.save(FILEPATH + ID + '.npy', I)
        yaml.dump(data, f)
    return


def generateFunction(n):
    global dopants

    g_ratio = lambda: beta.rvs(1.2, 1.2, loc=0.40, scale=0.40)
    snd_ratio = lambda: beta.rvs(2, 2, loc=0.15, scale=0.7)
    C_dist = lambda: beta.rvs(1.3, 1.3, loc=3.5, scale=6)
    # ----------->>> P impurity
    dopants[0]['capGenerator'] = lambda: randCapGenerator(C_dist=C_dist,
                                                          g_ratio=g_ratio,
                                                          snd_ratio=snd_ratio)
    # ----------->>> As impurity
    dopants[1]['capGenerator'] = lambda: randCapGenerator(C_dist=C_dist,
                                                          g_ratio=g_ratio,
                                                          snd_ratio=snd_ratio)
    # ----------->>> B impurity
    dopants[2]['capGenerator'] = lambda: randCapGenerator(C_dist=C_dist,
                                                          g_ratio=g_ratio,
                                                          snd_ratio=snd_ratio)

    dop_dist = lambda: np.random.randint(0, 2)
    T_dist = lambda: beta.rvs(4, 4, loc=0, scale=6)
    generation_loop(n, dop_dist, T_dist, [-150, 200], 100, [-70, 70], 100)


# =========================== MAIN ===========================
def main():
    # pltBeta(4, 4, loc=0, scale=6)
    num = 50
    generateFunction(num)
    for i in range(num):
        plt_file(-(i+1))


if __name__ == '__main__':
    main()




