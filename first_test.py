import set
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst
myset = set.SET()

# ============================ CONSTANTES =================================
# Cg, Cbg, Cd, Cs = np.array([100.0, 100.0, 100.0, 100.0])*cst.e  # F (aF = 1E-18 F)
Cg, Cd, Cs = (3.52, 0.86, 0.87)  # aF (aF = 1E-18 F)
Gd, Gs = (1.0, 1.0)
# Ag, Abg, Ad, As = (0.64, 0.16, 0.09, 0.13)


# ============================ SET CREATORS ===============================
def set1():
    print('starting test_set1')
    myset.add_lead('source')
    myset.add_lead('drain')
    myset.add_gate('gate')
    # myset.add_gate('back_gate')
    myset.add_metallic_dot('dot1', 10)
    #myset.add_quantum_dot('dot1', [0, 2, 1, 0.5], [1, 1, 1, 1])

    # adding links
    temp = 1E-1
    myset.add_link('dl', 'dot1', 'source', Cs, temp)
    myset.add_link('dl', 'dot1', 'drain', Cd, temp)
    myset.add_link('dg', 'dot1', 'gate', Cg)

    myset.set_temperature(2.0)

    # --------------- >>> running steps
    size = 100
    Vs = np.linspace(-50, 50, size)  # -35  # mV
    Vd = 0.0  # mV
    Vg = np.linspace(-800, -500, size)  # mV
    myset.pre_processing()
    I=[]
    for i in range(size):
        temp = []
        for j in range(size):
            myset.tunnel_rate([Vs[i], Vd, Vg[j]])
            myset.solver()
            temp.append(myset.current('source', 'dot1'))
        I.append(temp)
    I = np.array(I)
    print('Current is %s' % I)
    return I


def set2():
    print('starting test_set2')
    myset.add_lead('source')
    myset.add_lead('drain')
    myset.add_gate('gate')
    myset.add_metallic_dot('dot1', 10, 0, 0)
    #myset.add_quantum_dot('dot2', [0, 2, 1], [1, 1, 1])

    myset.add_link('dg', 'dot1', 'gate', 16E-18)
    myset.add_link('dl', 'dot1', 'source', 50E-18, 0.1)
    myset.add_link('dl', 'dot1', 'drain', 50E-18, 0.1)

    myset.set_temperature(1.0)

    myset.pre_processing()

    # --------------- >>> running steps
    size = 100
    Vs = np.linspace(-50, 50, size)  # -35  # mV
    Vd = 0.0  # mV
    Vg = np.linspace(-800, 800, 500)  # mV
    I = []
    # for i in range(size):
    #     temp = []
    #     for j in range(500):
    #         myset.tunnel_rate([Vs[i], Vd, Vg[j]])
    #         myset.solver()
    #         temp.append(myset.current('source', 'drain'))
    #     I.append(temp)
    # I = np.array(I)
    print('Current is %s' % I)
    return I


def set3(Vg_span, nb_Vg, Vds_span, nb_Vds, T=1.0):
    myset.add_quantum_dot('dot', [0, 10], [2, 2])  #, 10.0, 15.0, 17.0, 18.0], [1, 1, 1, 1, 1])

    # Add components to the dot to form the structure
    myset.add_lead('source')
    myset.add_lead('drain')
    myset.add_gate('gate')
    myset.add_link('dl', 'dot', 'drain', Cd * 1e-18, Gd * 2 / 77.461)
    myset.add_link('dl', 'dot', 'source', Cs * 1e-18, Gs * 2 / 77.461)
    myset.add_link('dg', 'dot', 'gate', Cg * 1e-18)
    myset.set_temperature(T)

    myset.pre_processing()

    # Calculate current for each volatge couple
    Vg = np.linspace(Vg_span[0], Vg_span[1], nb_Vg)
    Vd = np.linspace(Vds_span[0], Vds_span[1], nb_Vds)
    data_matrix = []
    temp = []
    for (i_vg, vg) in enumerate(Vg):
        I = []
        P = []
        V_dot = []
        for vd in Vd:
            myset.tunnel_rate([0, vd, vg])
            myset.solver()
            I.append(myset.current('drain', 'dot'))
            # I.append(myset.current('source','dot'))
            P.append(myset.proba('dot'))
            # V_dot.append(myset.voltage('dot'))
        # convert lists to scipy arrays
        I = np.array(I)
        temp.append(I)
        P = np.array(P)
        V_dot = np.array(V_dot)
    temp = np.array(temp)
    print(temp)
    return temp


def plt_mat(mat):
    plt.imshow(mat)
    plt.show()


# ====================================== MAIN ==================================
def main():
    global Cg, Cd, Cs, Gd, Gs
    #temp = set2()
    #plt_mat(temp)
    #temp = set1()
    #plt_mat(temp)
    Ec = 25  # meV
    Ctot = cst.e/(Ec*1E-3) * 1E18  # aF
    alpha = 0.9
    slew = 0.8
    Cg = Ctot*alpha
    Cs = Ctot*(1-alpha)*slew
    Cd = Ctot - Cg - Cs
    Gd, Gs = (1.0, 1.0)
    temp = set3([0.0, 100.0], 100, [-50.0, 50.0], 100, 1)
    plt_mat(temp.T)



if __name__ == '__main__':
    main()









'''
self.Gs = self.Gs_sim.value()
self.Gd = self.Gd_sim.value()
self.T = self.T_sim.value()
Vg_start = self.spin_Vg_start.value()
Vg_end = self.spin_Vg_end.value()
Vds_start = self.spin_Vds_start.value()
Vds_end = self.spin_Vds_end.value()
nb_Vg = self.spin_Vg_Pts.value()
nb_Vds = self.spin_Vds_Pts.value()
mode = self.box_mode.currentText()

# Define the set
myset = set.SET()
# Quantum dots or mettallic dot
if self.check_MD.isChecked():
    """
    If check_ag is checked, divide offset energy by lever arm
    Value divided by the lever arm is the offset observed on the absciss of stability diagramm 
    """
    if self.check_ag.isChecked():
        ag = self.Cg / (self.Cg + self.Cd + self.Cs)
        energy_offset = self.spin_E_offset.value() * ag
    else:
        energy_offset = self.spin_E_offset.value()
    myset.add_metallic_dot('dot', self.spin_nb_e.value(), 0, energy_offset)

elif self.check_QD.isChecked():
    # If levels and degeneracy are not put, a message is displayed in the command
    if len(self.line_levels.text()) != 0 and len(self.line_degeneracy.text()) != 0:
        global levels
        global degeneracy
        l = self.line_levels.text()[0::2]
        d = self.line_degeneracy.text()[0::2]
        levels = []
        degeneracy = []
        for i in l: levels.append(float(i))
        for i in d: degeneracy.append(int(i))
    try:
        myset.add_quantum_dot('dot', levels, degeneracy)
    except:
        print 'Select levels and degeneracy'

# Add components to the dot to form the structure
myset.add_lead('source')
myset.add_lead('drain')
myset.add_gate('gate')
myset.add_link('dl', 'dot', 'drain', self.Cd * 1e-18, self.Gd * 2 / 77.461)
myset.add_link('dl', 'dot', 'source', self.Cs * 1e-18, self.Gs * 2 / 77.461)
myset.add_link('dg', 'dot', 'gate', self.Cg * 1e-18)
myset.set_temperature(self.T)

myset.pre_processing()

# Calculate current for each volatge couple
Vg = np.linspace(Vg_start, Vg_end, nb_Vg)
Vd = np.linspace(Vds_start, Vds_end, nb_Vds)
data_matrix = []
for (i_vg, vg) in enumerate(Vg):
    I = []
    P = []
    V_dot = []
    for vd in Vd:
        myset.tunnel_rate([0, vd, vg])
        myset.solver()
        I.append(myset.current('drain', 'dot'))
        # I.append(myset.current('source','dot'))
        P.append(myset.proba('dot'))
        # V_dot.append(myset.voltage('dot'))
    # convert lists to scipy arrays
    I = np.array(I)
    P = np.array(P)
    V_dot = np.array(V_dot)
    # compute the diffential conductance
    if mode == 'Current':
        Y = Vd
        F = I
    elif mode == 'Diff conductance':
        F, Y = self.derive(I, Vd)
        F *= 1e3
    data_matrix.append(F)
data_matrix = np.array(data_matrix)
self.data = np.transpose(data_matrix)
self.plotSim(Vg.min(), Vg.max(), Y.min(), Y.max(), self.data, mode=mode)

self.runned = True'''