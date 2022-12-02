import set
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst
from scipy.stats import beta
import yaml  # to save python objects
import time  # to generate filenames
import warnings  # to supress linalg warnings
import scipy.linalg as linalg
from math import ceil, floor
import copy
# TODO adjust the step size of the data!!

# ====================== DECLARING CONSTANTS ======================
FILEPATH = "Data/sim3_0/train/"
EXP_FILEPATH = "Data/exp_w_labels/"


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
    plt.imshow(np.log10(np.abs(I)), extent=[Vg[0],Vg[-1],Vd[0],Vd[-1]], aspect='auto', cmap='hot')
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
             '{:.2f}'.format(info['ag']) + "    s_ratio:" + '{:.2f}'.format(info['s_ratio']))
    plt.annotate('LV: %s' % ['{:.1f}'.format(level) for level in info['levels']], (15, 20), xycoords='axes pixels', color='b')
    plt.annotate('degen: %s' % info['degens'], (15, 5), xycoords='axes pixels', color='b')
    if mesure == 'I':
        plt_current(data, Vg, Vds, title, multi_plt=multi_plt)
    elif mesure == 'G':
        plt_conduct(data, Vg, Vds, title, multi_plt=multi_plt)
    box = info['box']
    Vg = np.linspace(Vg[0], Vg[1], info["nVg"])
    Vds = np.linspace(Vds[0], Vds[1], info["nVds"])
    plt.gca().add_patch(plt.Rectangle([Vg[box[0][0]], Vds[box[1][1]]], Vg[box[1][0]] - Vg[box[0][0]], Vds[box[0][1]] - Vds[box[1][1]], fc='none', ec="b"))
    if not multi_plt:
        plt.show()
    return


def plt_exp_file(fileIndex, multi_plt=False):
    f = open(EXP_FILEPATH + '_data_indexer.yaml', 'r')
    info = yaml.load(f, Loader=yaml.FullLoader)[fileIndex]
    data = np.load(EXP_FILEPATH + info['f'] + ".npy")
    Vg = info['Vg_range']
    Vds = info['Vds_range']
    mesure = info['mesure']

    title = ("Ec:" + '{:.2f}'.format(info['Ec']))
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
    # TODO back gate
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
                                                     p=[0.17, 0.21, 0.32, 0.21, 0.09])  # arbitrairy probabilities
    elif n_levels == 2:
        rand_level_degens = lambda: np.random.choice(5, 1,
                                                     p=[0.23, 0.36, 0.22, 0.13, 0.06])  # arbitrairy probabilities
    elif n_levels == 3:
        rand_level_degens = lambda: np.random.choice(5, 1,
                                                     p=[0.33, 0.23, 0.19, 0.15, 0.1])  # arbitrairy probabilities
    elif n_levels == 4:
        rand_level_degens = lambda: np.random.choice(5, 1,
                                                     p=[0.35, 0.24, 0.2, 0.15, 0.06])  # arbitrairy probabilities
    else:
        rand_level_degens = lambda: np.random.choice(4, 1,
                                                     p=[0.48, 0.36, 0.12, 0.04])  # arbitrairy probabilities

    # here, it is important to understand that the first level must always have at least a degen of 2
    # because if not, the charging energy of the first diamond will not be Ec, but will be shifted by an
    # amount equal to the level spacing. also, if we only have 1 level, their won't be anny diamond if the
    # degen = 1. thus, this is why we initiate degens = [2] instead of [1]. this is justifiable because in
    # reality, all levels are at least degenerated 2 times because of spin if there is no magnetic field
    # (which is assumed here)
    levels = [0.0]
    degens = [2]
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


def generation_loop(n, T_dist, Ec_dist, g_ratio, snd_ratio, mesure='I'):
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
        Ctot = Cg + Cs + Cd  # aF
        ag = Cg/Ctot  # calculating alpha (the lever arm)
        s_ratio = Cs/(Ctot - Cg)
        T = float(T_dist())
        ID = str(time.time()).replace('.', 's')
        print("Propreties: \t Ec %smeV\t ag %s\t s_ratio %s" % (Ec, ag, s_ratio))

        # simulation
        levels, degens = randLevelGenerator()
        set1 = build_simulation(Cd, Cs, Cg, levels, degens)

        # ----------------->>> calculating Vg and Vds range and resolution
        # diamond dimentions
        height = Ec  # mV
        width = Ec / ag  # mV

        # setting the diamond resolution
        Vg_res = 1/beta.rvs(2.3, 1.5, 5.0, 60.0)
        if Vg_res*height <= 0.5:  # if step size < 0.5mV, clamp it
            Vg_res = 0.5/height
        Vds_res = 1 / beta.rvs(2.3, 1.5, 5.0, 80.0)
        if Vds_res*height <= 0.5:  # if step size < 0.5mV, clamp it
            Vds_res = 0.5/height

        Vg_step = Vg_res*height
        Vds_step = Vds_res*height
        temp = Ec*2
        if temp < 50.0:
            temp = 40.0
        if temp > 80:
            temp = 80
        # TODO add in a minimum number of steps
        temp = ceil(temp/Vds_step)*Vds_step
        Vds_range = [float(-temp), float(temp)]

        temp = width*(4 + 1/2)
        if temp > 200.0:  # maximum volt width
            temp = 200.0
        if temp < 60:  # minimum volt width
            temp = 60
        if temp/Vg_step > 200:  # maximum number of steps
            temp = 200.0*Vg_step
            if temp < 2*width:
                temp = 2.0*width
        # TODO add in a minimum number of steps
        temp = ceil(temp / Vds_step) * Vds_step
        Vg_range = [0, float(temp)]

        nVg = int((Vg_range[1] - Vg_range[0]) // Vg_step) + 1
        nVds = int((Vds_range[1] - Vg_range[0]) // Vds_step) + 1
        Vg = np.linspace(Vg_range[0], Vg_range[1], nVg)
        Vd = np.linspace(Vds_range[0], Vds_range[1], nVds)

        # ------------------->>> simulating
        diagram = simulate(set1, Vg, Vd, T, mesure=mesure)

        # calculating smalest box that frames the first diamond
        # this box is a 2*2 array [UperLeft corner, LowerRight corner] in indicies coordinates
        start = width/2  # mV
        left = int(np.argmin(np.abs(Vg - start)))
        right = int(np.argmin(np.abs(Vg - (start + width))))
        up = int(np.argmin(np.abs(Vd - height)))
        down = int(np.argmin(np.abs(Vd - (-height))))
        box = [[right, up], [left, down]]

        # saving data
        print("saving...")
        temp = {'f': ID,
                'Ec': Ec,
                'ag': ag,
                's_ratio': s_ratio,
                'Ctot': Ctot,
                'box': box,
                'T': T,
                'levels': levels,
                'degens': degens,
                'Vds_range': Vds_range,
                'Vg_range': Vg_range,
                'nVds': nVds,
                'nVg': nVg,
                'mesure': mesure,
                }
        print(temp)
        data.append(temp)
        f = open(FILEPATH + '_data_indexer.yaml', 'w')
        np.save(FILEPATH + ID + '.npy', diagram)
        yaml.dump(data, f)
    return


def generateFunction(n, mesure='I'):
    global dopants

    g_ratio = lambda: beta.rvs(1.2, 1.2, loc=0.40, scale=0.40)  # (1.2, 1.6, loc=0.10, scale=0.70)  # aF
    snd_ratio = lambda: beta.rvs(2, 2, loc=0.15, scale=0.7)  # aF
    Ec_dist = lambda: beta.rvs(1.2, 1.2, loc=15, scale=45)  # (2, 1.7, loc=12, scale=55)  # meV
    T_dist = lambda: beta.rvs(1.8, 2.1, loc=15, scale=50) #loc=1.5, scale=20)  # K

    generation_loop(n, T_dist, Ec_dist, g_ratio, snd_ratio, mesure=mesure)


# =========================== MAIN ===========================
def main():
    # pltBeta(2.3, 1.5, 5.0, 60.0)
    num = 100000
    generateFunction(num, mesure='I')
    # for i in range(num):
    #     plt_file(-(i+1))
    #     # plt_exp_file(-(i+1))


if __name__ == '__main__':
    main()




