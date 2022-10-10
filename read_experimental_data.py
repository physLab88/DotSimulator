import numpy as np
import matplotlib.pyplot as plt
import os


def to_grid(I, Vg, Vd):
    ''' This function rearanges the data in a matrix so it is easier to work with'''
    Vg = np.unique(Vg)
    nb_Vg = len(Vg)
    Vd = np.unique(Vd)
    nb_Vd = len(Vd)
    I = I.reshape((nb_Vd, nb_Vg))
    # inverse every odd lines
    i = np.indices(I.shape)
    print(i[0, :, :] % 2 == 0)
    I[i[0, :, :] % 2 == 0] = np.flip(I, 1)[i[0, :, :] % 2 == 0]
    I = I.T
    #I = np.flip(I, 1)
    return I, Vg, Vd


def plot_exp_data(filepath, title=None):
    data = np.loadtxt(filepath)
    Vg = data[:,0]
    Vd = data[:,1]
    I = data[:,2]
    if title is None:
        title = filepath
    plt.title(title)

    I, Vg, Vd = to_grid(I, Vg, Vd)  # re-arange the data so it is easier to work with
    print(I.shape)

    plt.title(title)
    plt.imshow(I, extent=[Vg[0], Vg[-1], Vd[0], Vd[-1]], aspect='auto')
    cbar = plt.colorbar(label='current in ??????')
    plt.xlabel(r'$V_g$ in mV')
    plt.ylabel(r'$V_{ds}$ in mV')
    #plt.axhline(0, color='k', alpha=0.3)
    #plt.axvline(0, color='k', alpha=0.3)


def plt_all():
    directory = 'Data/exp_data'
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            plot_exp_data(f, 'experimental data')
            plt.tight_layout()
            plt.show()

# ========================== MAIN =============================
def main():
    # fig, axs = plt.subplots(1,2)
    # plt.sca(axs[0])
    # plot_exp_data('Data/exp_data/ST28_Q05_Center_PCF19A_DUT11_4K_diamants_20220720-174550multi.txt', 'test')
    # plt.sca(axs[1])
    # plot_exp_data('Data/exp_data/ST28_Q05_Center_PCF20A_DUT2_4K_diamants__20220622-165225multi.txt', 'test')
    # plt.tight_layout()
    # plt.show()
    plt_all()
    return


if __name__ == '__main__':
    main()




