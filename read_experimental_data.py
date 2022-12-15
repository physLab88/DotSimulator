import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import yaml  # to save python objects
import scipy.constants as cst

directory = 'Data/exp_data_2'


def to_grid(I, Vg, Vd):
    ''' This function rearanges the data in a matrix so it is easier to work with'''
    Vg = np.unique(Vg)
    nb_Vg = len(Vg)
    Vd = np.unique(Vd)
    nb_Vd = len(Vd)
    I = I.reshape((nb_Vg, nb_Vd))
    # inverse every odd lines

    i = np.indices(I.shape)
    # print(i[0, :, :] % 2 == 0)
    I[i[0, :, :] % 2 == 0] = np.flip(I, 1)[i[0, :, :] % 2 == 0]
    I = I.T
    #I = np.flip(I, 1)
    return I, Vg, Vd


def load_exp(filepath, claude=False):
    data = np.loadtxt(filepath)
    if claude:
        Vg = data[:, 0]
        Vd = data[:, 2]
        I = data[:, 4]
    else:
        Vg = data[:, 0]
        Vd = data[:, 1]
        I = data[:, 2]

    I, Vg, Vd = to_grid(I, Vg, Vd)  # re-arange the data so it is easier to work with
    return I, Vg, Vd


def plot_exp_data(filepath, title=None, mesure='I', log=True, claude=False, clamp=10**-14):
    if title is None:
        title = filepath
    plt.title(title)

    I, Vg, Vd = load_exp(filepath, claude=claude)
    Vg *= 1E3
    Vd *= 1E3
    print(I.shape)
    if mesure == 'I':
        plt.title(title)
        I = np.abs(I)
        I[I<clamp] = clamp
        if log:
            plt.imshow(np.log10(I), extent=[Vg[0], Vg[-1], Vd[0], Vd[-1]], aspect='auto', cmap='hot')
            cbar = plt.colorbar(label='log(|current|)')
        else:
            plt.imshow(I, extent=[Vg[0], Vg[-1], Vd[0], Vd[-1]], aspect='auto', cmap='hot')
            cbar = plt.colorbar(label='current in ??????')
        plt.xlabel(r'$V_g$ in mV')
        plt.ylabel(r'$V_{ds}$ in mV')
        plt.axhline(0, color='k', alpha=0.3)
        #plt.axvline(0, color='k', alpha=0.3)
    elif mesure == 'G':
        G = np.gradient(I, (Vd[1] - Vd[0]) * 1E-3, axis=0)
        plt.title(title)
        if log:
            plt.imshow(np.log10(np.abs(G)), extent=[Vg[0], Vg[-1], Vd[0], Vd[-1]], aspect='auto', cmap='RdPu')
            cbar = plt.colorbar(label='log(|conductance|)')
        else:
            plt.imshow(G, extent=[Vg[0], Vg[-1], Vd[0], Vd[-1]], aspect='auto', cmap='RdPu')
            cbar = plt.colorbar(label='conductance in ??????')
        plt.xlabel(r'$V_g$ in mV')
        plt.ylabel(r'$V_{ds}$ in mV')
        plt.axhline(0, color='k', alpha=0.3)
        #plt.axvline(0, color='k', alpha=0.3)


def plt_I_vs_G(log=True):

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=[12, 5])
            plt.sca(axs[0])
            plot_exp_data(f, f, 'I', log)
            plt.sca(axs[1])
            plot_exp_data(f, 'experimental data', 'G', log)
            plt.tight_layout()
            print(f)
            plt.show()
    return


def plt_all(mesure='I', log=True):
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            plot_exp_data(f, 'experimental data', mesure, log)
            plt.tight_layout()
            plt.show()


def plt_I_analytics(filepath, n_curves, log=False):
    I, Vgs, Vds = load_exp(filepath)
    I = np.abs(I)
    if log:
        I = np.log(I)

    n_vd = Vds.shape[0]
    vd_index = np.linspace(0, n_vd-1, n_curves, dtype=int)
    for i in vd_index:
        vd = Vds[i]
        plt.plot(Vgs*1E3, I[i], label="Vds = {:0.2f} mV".format(vd*1E3))
    plt.xlabel("Vg in mV")
    plt.ylabel("I in A")
    plt.legend()
    # plt.tight_layout()
    return


def plt_I_analytics2(filepath, n_curves, log=False):
    cmap = cm.get_cmap("brg")
    I, Vgs, Vds = load_exp(filepath)
    # I = np.abs(I)
    if log:
        I = np.log(I)
    Imax = I.max()

    n_vg = Vgs.shape[0]
    vg_index = np.linspace(0, n_vg-1, n_curves, dtype=int)
    for k in range(len(vg_index)):
        i = vg_index[k]
        vg = Vgs[i]
        plt.plot(Vds*1E3, I[:,i],'-', label="Vg = {:0.2f} mV".format(vg*1E3), lw=1, color=cmap(float(k)/float(n_curves-1)))
    if log == False:
        #plt.ylim([0, Imax/12])
        pass
    plt.xlabel("Vds in mV")
    plt.ylabel("I in A")
    # plt.legend()
    # plt.tight_layout()
    return


def compare_analytics(filepath, log=False):
    fig, axs = plt.subplots(1, 3, figsize=[14, 5])
    plt.sca(axs[0])
    plot_exp_data(filepath, filepath, 'I', True)
    plt.sca(axs[1])
    plt_I_analytics(filepath, 6, log)
    plt.sca(axs[2])
    plt_I_analytics2(filepath, 18, log)

    #plt.tight_layout()
    plt.show()


conv_directory = "Data/exp_claude/"
files_to_convert = ["ST28_Q05_Center_PCF19A_DUT11_4K_diamants_20220720-160541multi.txt",
                    "ST28_Q05_Center_PCF19A_DUT11_4K_diamants_20220720-165134multi.txt",
                    #"ST28_Q05_Center_PCF19A_DUT11_4K_diamants_20220720-174550multi.txt",
                    "ST28_Q05_Center_PCF20A_DUT2_4K_diamants_20220622-153810.txt_20220622-153810multi.txt",
                    "ST28_Q05_North_PCF19A_DUT8_4K_diamants_20221027-130719multi.txt",
                    "ST28_Q05_North_PCF20A_DUT8_4K_diamants_20221026-160229multi_abs.txt",
                    ]
values_of_files = [{'Ec':cst.e/(8.2228*1E-18)*1E3,},
                   {'Ec':cst.e/(8.2236*1E-18)*1E3,},
                   {'Ec':cst.e/(5.2601*1E-18)*1E3,},
                   {'Ec':cst.e/(2.257*1E-18)*1E3,},
                   {'Ec':cst.e/(6.4581*1E-18)*1E3,},]
files_to_convert2 = ["20220615-095619_Map_STR_VDSR.txt",
                     "20220615-134200_Map_STR_VDSR.txt",
                    ]
values_of_files2 = [{'Ec': 1.6926,},  # mV
                   {'Ec':3.11,},]


def convert_exp_data(files, values, claude=False):
    try:
        f = open(conv_directory + '_data_indexer.yaml', 'r')
        data = yaml.load(f, Loader=yaml.FullLoader)
    except IOError:
        f = open(conv_directory + '_data_indexer.yaml', 'w')  # creating a file if it does not exist
        data = []

    for file, info in zip(files, values):
        diagram, Vg, Vds = load_exp(directory + '/' + file, claude=claude)
        Vg *= 1E3
        Vds *= 1E3
        file = file[:file.find('.')]
        info['f'] = file
        info['Vds_range'] = [float(Vds[0]), float(Vds[-1])]
        info['Vg_range'] = [float(Vg[0]), float(Vg[-1])]
        info['nVds'] = len(Vds)
        info['nVg'] = len(Vg)
        info['mesure'] = 'I'
        data.append(info)
        np.save(conv_directory + file + '.npy', diagram)
    f = open(conv_directory + '_data_indexer.yaml', 'w')
    yaml.dump(data, f)
    return


# ========================== MAIN =============================
def main():
    # plot_exp_data("Data/exp_data/20220615-095619_Map_STR_VDSR.txt", claude=True, clamp=10**-11)
    # plt.show()
    # plot_exp_data("Data/exp_data/20220615-134200_Map_STR_VDSR.txt", claude=True, clamp=10**-11)
    # plt.show()
    # convert_exp_data(files_to_convert2, values_of_files2, claude=True)
    # filepaths = [
    #     "Data/exp_data_2\ST28_Q05_Center_PCF19A_DUT11_4K_diamants_20220720-160541multi.txt",
    #     #"Data/exp_data_2\ST28_Q05_Est_PCF26A_DUT2_2K_diamants_20221013-162619multi.txt",
    #     #"Data/exp_data_2\ST28_Q05_Est_PCF26A_DUT2_2K_diamants_20221013-165319multi_abs.txt",
    #     "Data/exp_data_2\ST28_Q05_Est_PCF26A_DUT2_4K_diamants_20221013-130233multi.txt",
    #     #"Data/exp_data_2\ST28_Q05_Est_PCF26A_DUT4_4K_diamants_20221013-141659multi.txt",
    #     #"Data/exp_data_2\ST28_Q05_Est_PCF26A_DUT5_4K_diamants_20221013-144919multi.txt",
    #     "Data/exp_data_2\ST28_Q05_North_PCF20A_DUT2_4K_diamants_20221026-134517multi.txt",
    #     "Data/exp_data_2\ST28_Q05_North_PCF20A_DUT2_4K_diamants_20221026-134517multi_abs.txt",
    #     "Data/exp_data_2\ST28_Q05_North_PCF20A_DUT3_4K_diamants_20221026-140959multi.txt",
    #     "Data/exp_data_2\ST28_Q05_North_PCF20A_DUT7_4K_diamants_20221026-153023multi.txt",
    #     "Data/exp_data_2\ST28_Q05_West_PCF19A_DUT5_4K_diamants_20220914-165449multi.txt",
    #     "Data/exp_data_2\ST28_Q05_West_PCF19A_DUT6_4K_diamants_20220914-172504multi.txt",
    #     "Data/exp_data_2\ST28_Q05_West_PCF19A_DUT7_4K_diamants_20220914-174924multi.txt",
    #     "Data/exp_data_2\ST28_Q05_West_PCF20A_DUT2_4K_diamants_20220914-144200multi.txt",
    # ]
    # plt_I_vs_G(True)
    plt_all('I', True)
    # for i in range(len(filepaths)):
    #     filepath = filepaths[i]
    #     print(i)
    #     compare_analytics(filepath, False)
    return


if __name__ == '__main__':
    main()




