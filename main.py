# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import math
import cmath
import random
import numpy as np
import sys


def p_random(arr1, arr2):
    """
    生成指定概率的随机数
    :param arr1: Length does not match.
    :param arr2: Total rate is not 1.
    :return:
    """
    assert len(arr1) == len(arr2), "Length does not match."
    assert sum(arr2) == 1, "Total rate is not 1."

    sup_list = [len(str(i).split(".")[-1]) for i in arr2]
    top = 10 ** max(sup_list)
    new_rate = [int(i * top) for i in arr2]
    rate_arr = []
    for i in range(1, len(new_rate) + 1):
        rate_arr.append(sum(new_rate[:i]))
    rand = random.randint(1, top)
    data = None
    for i in range(len(rate_arr)):
        if rand <= rate_arr[i]:
            data = arr1[i]
            break
    return data


class SimRIS:

    def __init__(self, Environment, Scenario, Frequency, ArrayType, N, Nt, Nr, Nc, Tx_xyz, Rx_xyz, RIS_xyz):
        """

        :param Environment:
        :param Scenario:
        :param Frequency:
        :param ArrayType:
        :param N:
        :param Nt:
        :param Nr:
        :param Tx_xyz:
        :param Rx_xyz:
        :param RIS_xyz:
        """

        self.Environment = Environment
        self.Scenario = Scenario
        self.Frequency = Frequency
        self.ArrayType = ArrayType
        self.N = N
        self.Nt = Nt
        self.Nr = Nr
        self.Nc = Nc
        self.Tx_xyz = Tx_xyz
        self.Rx_xyz = Rx_xyz
        self.RIS_xyz = RIS_xyz

        self.H = []
        self.G = []

        self.G_c_num=0
        self.G_s_list=[]
        self.G_c=[]
        self.G_s=[]
        self.G_i = []

        self.waveLength = (3 * (10 ** 8)) / (self.Frequency * (10 ** 9))  # lambda
        self.k = 2 * math.pi / self.waveLength  # Wave Number
        self.dis = self.waveLength / 2  # RIS element spacing (can be modified to observe its effect on h and g)

        if math.sqrt(self.N) % 1:
            sys.exit("N should be an integer power of 2")

        if Environment == 1:
            '''
            INDOORS
            '''
            # Environment 1 (InH - Office) - NLOS
            self.room_height = 3.5  # indoor room height
            self.n_NLOS = 3.19  # Path Loss Exponent (Indoor Office NLOS)
            self.sigma_NLOS = 8.29  # Shadow Fading Term (dB) (Indoor Office NLOS)
            self.b_NLOS = 0.06  # Path Loss Parameter (Indoor Office NLOS)
            self.f0 = 24.2  # Path Loss Parameter (GHz) (Indoor Office NLOS)
            # Environment 1 (InH - Office) - LOS
            self.n_LOS = 1.73  # Path Loss Exponent (Indoor Office LOS)
            self.sigma_LOS = 3.02  # Shadow Fading Term (dB) (Indoor Office LOS)
            self.b_LOS = 0.0  # Path Loss Parameter (Indoor Office NLOS)
        elif Environment == 2:
            '''
            OUTDOORS
            '''
            # Environment 1 (InH - Office) - NLOS
            self.n_NLOS = 3.19  # Path Loss Exponent (Indoor Office NLOS)
            self.sigma_NLOS = 8.2  # Shadow Fading Term (dB) (Indoor Office NLOS)
            self.b_NLOS = 0.0  # Path Loss Parameter (Indoor Office NLOS)
            self.f0 = 24.2  # Path Loss Parameter (GHz) (Indoor Office NLOS)
            # Environment 1 (InH - Office) - LOS
            self.n_LOS = 1.98  # Path Loss Exponent (Indoor Office LOS)
            self.sigma_LOS = 3.1  # Shadow Fading Term (dB) (Indoor Office LOS)
            self.b_LOS = 0.0  # Path Loss Parameter (Indoor Office NLOS)
        else:
            sys.exit("Environment must 1 (INDOORS) or 2 (OUTDOORS)")

        # Parameter of Number of Clusters
        if Frequency == 28:
            self.lambda_p = 1.8
        elif Frequency == 73:
            self.lambda_p = 1.9

        # Element Radiation Pattern Parameters
        self.q = 0.285
        self.Gain = math.pi

    def generateH(self):
        x_Tx = self.Tx_xyz[0]
        y_Tx = self.Tx_xyz[1]
        z_Tx = self.Tx_xyz[2]

        x_Rx = self.Rx_xyz[0]
        y_Rx = self.Rx_xyz[1]
        z_Rx = self.Rx_xyz[2]

        x_RIS = self.RIS_xyz[0]
        y_RIS = self.RIS_xyz[1]
        z_RIS = self.RIS_xyz[2]

        # Calculate Tx-RIS LOS distance  and Generate LOS Component for h
        d_T_RIS = np.linalg.norm(self.Tx_xyz - self.RIS_xyz)
        # LOS Probability is Relatively Low for Indoors if d_T_RIS> 20

        for repeat in range(self.Nc):
            '''
            STEP 1
            '''

            if self.Environment == 1:
                if z_RIS < z_Tx:
                    # InH LOS Probability
                    if d_T_RIS <= 1.2:
                        p_LOS = 1.0
                    elif 1.2 < d_T_RIS < 6.5:
                        p_LOS = math.exp(-(d_T_RIS - 1.2) / 4.7)
                    else:
                        p_LOS = 0.32 * math.exp(-(d_T_RIS - 6.5) / 32.6)
                    I_LOS = p_random([1, 0], [p_LOS, 1 - p_LOS])
                elif z_RIS >= z_Tx:
                    # for an RIS mounted at a high place (100% LOS)
                    I_LOS = 1
            elif self.Environment == 2:
                # UMi LOS Probability
                p_LOS = min([20 / d_T_RIS, 1]) * (1 - math.exp(-d_T_RIS / 39)) + math.exp(-d_T_RIS / 39)
                I_LOS = p_random([1, 0], [p_LOS, 1 - p_LOS])

            theta_T_RIS_LOS = 0.0
            phi_T_RIS_LOS = 0.0
            phi_Tx_LOS = 0.0
            theta_Tx_LOS = 0.0
            if I_LOS == 1:
                # Calculate Tx Departure and RIS arrival angles to calculate array
                # response vectors
                if self.Scenario == 1:  # side-wall RIS
                    I_phi = np.sign(x_RIS - x_Tx)
                    phi_T_RIS_LOS = I_phi * math.atan(abs(x_RIS - x_Tx) / abs(y_RIS - y_Tx))

                    I_theta = np.sign(z_Tx - z_RIS)
                    theta_T_RIS_LOS = I_theta * math.asin(abs(z_RIS - z_Tx) / d_T_RIS)

                    # Tx departure angles for LOS component
                    I_phi_Tx = np.sign(y_Tx - y_RIS)
                    phi_Tx_LOS = I_phi_Tx * math.atan(abs(y_Tx - y_RIS) / abs(x_Tx - x_RIS))

                    I_theta_Tx = np.sign(z_Tx - z_RIS)
                    theta_Tx_LOS = I_theta_Tx * math.asin(abs(z_RIS - z_Tx) / d_T_RIS)
                elif self.Scenario == 2:  # opposite-wall RIS
                    I_phi = np.sign(y_Tx - y_RIS)  # These are different from Scenario 1
                    phi_T_RIS_LOS = I_phi * math.atan(abs(y_RIS - y_Tx) / abs(x_RIS - x_Tx))

                    I_theta = np.sign(z_Tx - z_RIS)  # Same as Scenario 1
                    theta_T_RIS_LOS = I_theta * math.asin(abs(z_RIS - z_Tx) / d_T_RIS)

                    I_phi_Tx = np.sign(y_Tx - y_RIS)
                    phi_Tx_LOS = I_phi_Tx * math.atan(abs(y_RIS - y_Tx) / abs(x_RIS - x_Tx))

                    I_theta_Tx = np.sign(z_RIS - z_Tx)
                    theta_Tx_LOS = I_theta_Tx * math.asin(abs(z_RIS - z_Tx) / d_T_RIS)

                #  Array Response Calculation (LOS) (Be Careful with sin/cos)
                array_RIS_LOS = np.zeros([1, self.N], dtype=complex)
                counter = 0
                for x in range(0, round(math.sqrt(self.N))):
                    for y in range(0, round(math.sqrt(self.N))):
                        array_RIS_LOS[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                x * math.sin(theta_T_RIS_LOS) + y * math.sin(phi_T_RIS_LOS) * math.cos(
                            theta_T_RIS_LOS)))
                        counter = counter + 1

                array_Tx_LOS = np.zeros([1, self.Nt], dtype=complex)
                if self.ArrayType == 1:
                    counter = 0
                    for x in range(0, self.Nt):
                        array_Tx_LOS[0, counter] = cmath.exp(
                            1j * self.k * self.dis * (x * math.sin(phi_Tx_LOS) * math.cos(theta_Tx_LOS)))
                        counter = counter + 1
                elif self.ArrayType == 2:
                    counter = 0
                    for x in range(0, round(math.sqrt(self.Nt))):
                        for y in range(0, round(math.sqrt(self.Nt))):
                            array_Tx_LOS[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                    x * math.sin(phi_Tx_LOS) * math.cos(theta_Tx_LOS) + y * math.sin(theta_Tx_LOS)))
                            counter = counter + 1

                # Link Attentuation (LOS)
                # Note: This is different than FSPL (Shadowing/Waveguiding effect included with n < 2)
                L_dB_LOS = -20 * math.log10(4 * math.pi / self.waveLength) - 10 * self.n_LOS * (
                        1 + self.b_LOS * ((self.Frequency - self.f0) / self.f0)) * math.log10(d_T_RIS) - math.fabs(
                    np.random.randn()) * self.sigma_LOS
                L_LOS = np.power(10.0, L_dB_LOS / 10.0)

                # LOS COMPONENT GENERATED (with random phase) (practical PL and shadowing)
                # Element radiation pattern considered

                h_LOS = math.sqrt(L_LOS) * np.transpose(array_RIS_LOS) * array_Tx_LOS * cmath.exp(
                    1j * random.random() * 2 * math.pi) * math.sqrt(
                    self.Gain * (math.cos(theta_T_RIS_LOS)) ** (2 * self.q))
            else:
                h_LOS = 0

            '''
            STEP 2
            '''
            # Generate Clusters/Sub-rays, Azimuth/Elevation Departure Angles and Cluster Distances
            for generate in range(100):
                # Number of Clusters - lambda_p was defined earlier: 1.8/1.9
                C = max([1, np.random.poisson(lam=self.lambda_p)])  # Poisson distributed

                # Number of Sub-rays per Cluster
                S = np.random.randint(30, size=C)

                while True:
                    if sum(S) != 0:
                        break
                    else:
                        S = np.random.randint(30, size=C)

                # Azimuth/Elevation Departure Angles
                phi_Tx = []
                theta_Tx = []
                phi_av = np.zeros(C)
                theta_av = np.zeros(C)

                isCoordinatesDone = False
                while True:
                    if isCoordinatesDone:
                        break
                    for counter in range(C):
                        phi_av[counter] = np.random.rand() * 180 - 90  # mean azimuth
                        theta_av[counter] = np.random.rand() * 90 - 45  # mean elevation
                        #  (needs update & most of the clusters are above ceiling)

                        # for outdoors this might not be a big concern
                        # cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
                        new_phi_Tx = np.log(
                            np.random.rand(1, S[counter]) / np.random.rand(1, S[counter])).flatten() * math.sqrt(
                            25 / 2) + phi_av[counter]
                        phi_Tx.extend(new_phi_Tx)
                        new_theta_Tx = np.log(
                            np.random.rand(1, S[counter]) / np.random.rand(1, S[counter])).flatten() * math.sqrt(
                            25 / 2) + theta_av[counter]
                        theta_Tx.extend(new_theta_Tx)

                    # Cluster Distances
                    #  Cluster distances uniform [1,d_T_RIS]
                    #  can be modified later
                    a_c = 1 + np.random.rand(C) * (d_T_RIS - 1)

                    # Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
                    Coordinates = np.zeros(shape=[C, 3])  # for Clusters
                    Coordinates2 = np.zeros(shape=[sum(S), 3])  # for Scatterers

                    # Correction on Cluster Locations for Indoors
                    if self.Environment == 1:
                        #  Room dimensions (Indoor Hotspot)
                        # x-y dimensions recommended by 5G Channel Model, height is assumed as 3.5 m
                        dim = [75, 50, self.room_height]
                        for counter in range(C):
                            loop = 1
                            while True:
                                Coordinates[counter, :] = [
                                    x_Tx + a_c[counter] * math.cos(theta_av[counter]) * math.cos(phi_av[counter]),
                                    y_Tx - a_c[counter] * math.cos(theta_av[counter]) * math.sin(phi_av[counter]),
                                    z_Tx + a_c[counter] * math.sin(theta_av[counter])]
                                if 0 < Coordinates[counter, 2] < dim[2] and 0 < Coordinates[counter, 1] < dim[1] and 0 < \
                                        Coordinates[counter, 0] < dim[0]:
                                    isCoordinatesDone = True
                                    break
                                else:
                                    if loop > 10:
                                        isCoordinatesDone = False
                                        break
                                    a_c[counter] = 0.8 * a_c[counter]
                                    # Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                                    # in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
                                    loop += 1
                    elif self.Environment == 2:
                        # Outdoors
                        for counter in range(C):
                            loop = 1
                            while True:
                                Coordinates[counter, :] = [
                                    x_Tx + a_c[counter] * math.cos(theta_av[counter]) * math.cos(phi_av[counter]),
                                    y_Tx - a_c[counter] * math.cos(theta_av[counter]) * math.sin(phi_av[counter]),
                                    z_Tx + a_c[counter] * math.sin(theta_av[counter])]
                                if Coordinates[counter, 2] < 0:
                                    isCoordinatesDone = True
                                    break
                                else:
                                    if loop > 10:
                                        isCoordinatesDone = False
                                        break
                                    else:
                                        a_c[counter] = 0.8 * a_c[counter]
                                        # Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                                        # in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
                                        loop += 1
                # Plot Scatterers
                a_c_rep = []
                for counter in range(C):
                    a_c_rep.extend(np.full((1, S[counter]), a_c[counter]).flatten())
                for counter2 in range(sum(S)):
                    Coordinates2[counter2, :] = [
                        x_Tx + a_c_rep[counter2] * math.cos(theta_Tx[counter2]) * math.cos(phi_Tx[counter2]),
                        y_Tx - a_c_rep[counter2] * math.cos(theta_Tx[counter2]) * math.sin(phi_Tx[counter2]),
                        z_Tx + a_c_rep[counter2] * math.sin(theta_Tx[counter2])]
                if self.Environment == 1:
                    ignore = []
                    for counter2 in range(sum(S)):
                        if Coordinates2[counter2, 2] > dim[2] or Coordinates2[counter2, 2] < 0 or Coordinates2[
                            counter2, 1] > dim[1] or Coordinates2[counter2, 1] < 0 or Coordinates2[counter2, 0] > dim[
                            0] or \
                                Coordinates2[counter2, 0] < 0:
                            ignore.extend([int(counter2)])
                elif self.Environment == 2:
                    ignore = []
                    for counter2 in range(sum(S)):
                        if Coordinates2[counter2, 2] < 0:  # only underground scatterers
                            ignore.extend([int(counter2)])
                # updated indices
                indices = np.setdiff1d(np.arange(1, sum(S), 1), ignore)
                M_new = len(indices)

                # if you want to revert back Version 1.1 from Version 1.2, simply set
                # indices=1:sum(S);
                # M_new=sum(S)

                # Necessary Loop to have at least one scatterer
                if M_new > 0:
                    break

            '''
            STEP 3
            '''
            # Calculate Arrival Angles for the RIS and the Link Distances
            phi_cs_RIS = np.zeros(sum(S))
            theta_cs_RIS = np.zeros(sum(S))
            phi_Tx_cs = np.zeros(sum(S))
            theta_Tx_cs = np.zeros(sum(S))
            b_cs = np.zeros(sum(S))
            d_cs = np.zeros(sum(S))
            for counter in indices:
                b_cs[counter] = np.linalg.norm(
                    self.RIS_xyz - Coordinates2[counter, :])  # Distance between Scatterer and RIS
                d_cs[counter] = a_c_rep[counter] + b_cs[counter]  # Total distance Tx-Scatterer-RIS
                if self.Environment == 1:
                    # side-wall RIS

                    I_phi = np.sign(x_RIS - Coordinates2[counter, 0] - y_RIS)
                    phi_cs_RIS[counter] = I_phi * math.atan(
                        abs(x_RIS - Coordinates2[counter, 0]) / abs(y_RIS - Coordinates2[counter, 1]))

                    I_theta = np.sign(Coordinates2[counter, 2] - z_RIS)
                    theta_cs_RIS[counter] = I_theta * math.asin(abs(z_RIS - Coordinates2[counter, 2]) / b_cs[counter])

                    I_phi_Tx_cs = np.sign(y_Tx - Coordinates2[counter, 1])
                    phi_Tx_cs[counter] = I_phi_Tx_cs * math.atan(
                        abs(Coordinates2[counter, 1] - y_Tx) / abs(Coordinates2[counter, 0] - x_Tx))

                    I_theta_Tx_cs = np.sign(Coordinates2[counter, 2] - z_Tx)
                    theta_Tx_cs[counter] = I_theta_Tx_cs * math.asin(
                        abs(Coordinates2[counter, 2] - z_Tx) / a_c_rep[counter])

                if self.Environment == 2:
                    # opposite-wall RIS

                    I_phi = np.sign(Coordinates2[counter, 1] - y_RIS)
                    phi_cs_RIS[counter] = I_phi * math.atan(
                        abs(y_RIS - Coordinates2[counter, 1]) / abs(x_RIS - Coordinates2[counter, 0]))

                    I_theta = np.sign(Coordinates2[counter, 2] - z_RIS)
                    theta_cs_RIS[counter] = I_theta * math.asin(abs(z_RIS - Coordinates2[counter, 2]) / b_cs[counter])

                    I_phi_Tx_cs = np.sign(y_Tx - Coordinates2[counter, 1])
                    phi_Tx_cs[counter] = I_phi_Tx_cs * math.atan(
                        abs(Coordinates2[counter, 1] - y_Tx) / abs(Coordinates2[counter, 0] - x_Tx))

                    I_theta_Tx_cs = np.sign(Coordinates2[counter, 2] - z_Tx)
                    theta_Tx_cs[counter] = I_theta_Tx_cs * math.asin(
                        abs(Coordinates2[counter, 2] - z_Tx) / a_c_rep[counter])

            '''
            STEP 4
            '''
            # Array Response Calculation
            array_cs_RIS = np.zeros([sum(S), self.N], dtype=complex)
            for counter in indices:
                counter2 = 0
                for x in range(round(math.sqrt(self.N))):
                    for y in range(round(math.sqrt(self.N))):
                        array_cs_RIS[counter, counter2] = cmath.exp(1j * self.k * self.dis * (
                                x * math.sin(theta_cs_RIS[counter]) + y * math.sin(phi_cs_RIS[counter]) * math.cos(
                            theta_cs_RIS[counter])))
                        counter2 += 1
            array_Tx_cs = np.zeros([sum(S), self.Nt], dtype=complex)
            if self.ArrayType == 1:
                for counter in indices:
                    counter2 = 0
                    for x in range(self.Nt):
                        array_Tx_cs[counter, counter2] = cmath.exp(
                            1j * self.k * self.dis * (
                                    x * math.sin(phi_Tx_cs[counter]) * math.cos(theta_Tx_cs[counter])))
                        counter2 += 1

            elif self.ArrayType == 2:
                for counter in indices:
                    counter2 = 0
                    for x in range(round(math.sqrt(self.Nt))):
                        for y in range(round(math.sqrt(self.Nt))):
                            array_Tx_cs[counter, counter2] = cmath.exp(1j * self.k * self.dis * (
                                    x * math.sin(phi_Tx_cs[counter]) * math.cos(theta_Tx_cs[counter]) + y * math.sin(
                                theta_Tx_cs[counter])))
                            counter2 += 1
            '''
            STEP 5
            '''
            # Common for Enviroments 1 and 2
            # Calculate Link Attenuation and Generate Tx-RIS Channel (h) using 5G Channel Model
            # Ge Updated
            h_NLOS = np.zeros([self.N, self.Nt])
            beta = np.zeros([1, sum(S)], dtype=complex)  # to be reused for shared clusters (Environment 1)
            shadow = beta  # to be reused for shared clusters (Environment 1)
            for counter in indices:
                X_sigma = np.random.randn() * self.sigma_NLOS
                Lcs_dB = -20 * math.log10(4 * math.pi / self.waveLength) - 10 * self.n_NLOS * (
                        1 + self.b_NLOS * ((self.Frequency - self.f0) / self.f0)) * math.log10(d_cs[counter]) - X_sigma
                Lcs = 10 ** (Lcs_dB / 10)
                beta[0, counter] = (np.random.randn() + 1j * np.random.randn()) / math.sqrt(
                    2)  # common complex gain for shared clusters
                shadow[0, counter] = X_sigma  # commun shadow factor for shared clusters

                # consider all scatters
                h_NLOS = h_NLOS + beta[0, counter] * math.sqrt(
                    self.Gain * (math.cos(theta_cs_RIS[counter])) ** (2 * self.q)) * math.sqrt(Lcs) * array_cs_RIS[
                                                                                                      counter,
                                                                                                      :].reshape(
                    len(array_cs_RIS[counter, :]), 1) * array_Tx_cs[counter, :]
            h_NLOS = h_NLOS * math.sqrt(1 / M_new)
            h = h_NLOS + h_LOS  # include the LOS component (if any) h_LOS=0 when there is no LOS
        self.H.append(h)
        return self.H

    def generateG(self):
        x_Tx = self.Tx_xyz[0]
        y_Tx = self.Tx_xyz[1]
        z_Tx = self.Tx_xyz[2]

        x_Rx = self.Rx_xyz[0]
        y_Rx = self.Rx_xyz[1]
        z_Rx = self.Rx_xyz[2]

        x_RIS = self.RIS_xyz[0]
        y_RIS = self.RIS_xyz[1]
        z_RIS = self.RIS_xyz[2]

        # Calculate Tx-RIS LOS distance  and Generate LOS Component for h
        d_T_RIS = np.linalg.norm(self.Tx_xyz - self.RIS_xyz)
        # LOS Probability is Relatively Low for Indoors if d_T_RIS> 20

        G_new = []

        for repeat in range(self.Nc):
            '''
            STEP 6-7
            '''
            # Generation of g (RIS-Rx Channel)
            if self.Environment == 1:
                # GENERATE A LOS CHANNEL
                # Calculate Departure Angles Considering RIS and Rx Coordinates
                d_RIS_R = np.linalg.norm(self.RIS_xyz - self.Rx_xyz)
                # Elevation Departure Angle
                I_theta = np.sign(z_Rx - z_RIS)
                theta_Rx_RIS = I_theta * math.asin(abs(z_Rx - z_RIS) / d_RIS_R)  # AoD of RIS

                # AoA angles of Rx for g_LOS channel in an Indoor
                phi_av_Rx = np.random.rand() * 180 - 90  # mean azimuth
                theta_av_Rx = np.random.rand() * 180 - 90  # mean elevation

                phi_Rx = math.log(np.random.rand() / np.random.rand()) * math.sqrt(25 / 2) + phi_av_Rx
                theta_Rx = math.log(np.random.rand() / np.random.rand()) * math.sqrt(25 / 2) + theta_av_Rx

                # Azimuth Departure Angle
                if self.Scenario == 1:
                    I_phi = np.sign(x_RIS - x_Rx)
                    phi_Rx_RIS = I_phi * math.atan(abs(x_Rx - x_RIS) / abs(y_Rx - y_RIS))

                elif self.Scenario == 2:
                    I_phi = np.sign(y_Rx - y_RIS)
                    phi_Rx_RIS = I_phi * math.atan(abs(y_Rx - y_RIS) / abs(x_Rx - x_RIS))

                # Recalculate Array Response for Two Angles (in Rx direction)
                array_2 = np.zeros([1, self.N], dtype=complex)
                counter = 0
                for x in range(round(math.sqrt(self.N))):
                    for y in range(round(math.sqrt(self.N))):
                        array_2[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                x * math.sin(theta_Rx_RIS) + y * math.sin(phi_Rx_RIS) * math.cos(theta_Rx_RIS)))
                        counter += 1

                array_Rx = np.zeros([1, self.Nr], dtype=complex)
                counter = 0
                if self.ArrayType == 1:
                    for x in range(self.Nr):
                        array_Rx[0, counter] = cmath.exp(
                            1j * self.k * self.dis * (x * math.sin(phi_Rx) * math.cos(theta_Rx)))
                        counter += 1
                elif self.ArrayType == 2:
                    for x in range(round(math.sqrt(self.Nr))):
                        for y in range(round(math.sqrt(self.Nr))):
                            array_Rx[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                    x * math.sin(phi_Rx) * math.cos(theta_Rx) + y * math.sin(theta_Rx)))
                            counter += 1

                # LOS Link Attenuation
                L_dB_LOS_2 = -20 * math.log10(4 * math.pi / self.waveLength) - 10 * self.n_LOS * (
                        1 + self.b_LOS * ((self.Frequency - self.f0) / self.f0)) * math.log10(
                    d_RIS_R) - np.random.randn() * self.sigma_LOS
                L_LOS_2 = 10 ** (L_dB_LOS_2 / 10)
                # Generate g (Pure LOS)
                g = math.sqrt(self.Gain * (math.cos(theta_Rx_RIS)) ** (2 * self.q)) * math.sqrt(L_LOS_2) * np.transpose(
                    array_2) * array_Rx * cmath.exp(1j * np.random.rand() * 2 * math.pi)
            elif self.Environment == 2:
                d_RIS_R = np.linalg.norm(self.RIS_xyz - self.Rx_xyz)
                # UMi LOS Probability
                p_LOS_2 = min([20 / d_RIS_R, 1]) * (1 - math.exp(-d_RIS_R / 39)) + math.exp(-d_RIS_R / 39)
                I_LOS_2 = p_random([1, 0], [p_LOS_2, 1 - p_LOS_2])

                I_theta = np.sign(z_Rx - z_RIS)
                theta_RIS_Rx_LOS = I_theta * math.asin(abs(z_Rx - z_RIS) / d_RIS_R)  # AoD of RIS

                # AoA angles of Rx for g_LOS channel in an Outdoor
                phi_av_Rx_LOS = np.random.rand() * 180 - 90  # mean azimuth
                theta_av_Rx_LOS = np.random.rand() * 180 - 90  # mean elevation

                phi_Rx_LOS = math.log(np.random.rand() / np.random.rand()) * math.sqrt(25 / 2) + phi_av_Rx_LOS
                theta_Rx_LOS = math.log(np.random.rand() / np.random.rand()) * math.sqrt(25 / 2) + theta_av_Rx_LOS

                if I_LOS_2 == 1:
                    if self.Scenario == 1:
                        I_phi = np.sign(x_RIS - x_Rx)
                        phi_RIS_Rx_LOS = I_phi * math.atan(abs(x_RIS - x_Rx) / abs(y_RIS - y_Rx))

                    elif self.Scenario == 2:
                        I_phi = np.sign(y_Rx - y_RIS)
                        phi_RIS_Rx_LOS = I_phi * math.atan(abs(y_Rx - y_RIS) / abs(x_Rx - x_RIS))

                    # Array Response Calculation(LOS)
                    array_RIS_Rx_LOS = np.zeros([1, self.N], dtype=complex)
                    array_Rx_LOS = np.zeros([1, self.Nr], dtype=complex)

                    counter = 0
                    for x in range(round(math.sqrt(self.N))):
                        for y in range(round(math.sqrt(self.N))):
                            array_RIS_Rx_LOS[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                    x * math.sin(theta_RIS_Rx_LOS) + y * math.sin(phi_RIS_Rx_LOS) * math.cos(
                                theta_RIS_Rx_LOS)))
                            counter += 1

                    counter = 0
                    if self.ArrayType == 1:
                        for x in range(self.Nr):
                            array_Rx_LOS[0, counter] = cmath.exp(
                                1j * self.k * self.dis * (x * math.sin(phi_Rx_LOS) * math.cos(theta_Rx_LOS)))
                            counter += 1
                    elif self.ArrayType == 2:
                        for x in range(round(math.sqrt(self.Nr))):
                            for y in range(round(math.sqrt(self.Nr))):
                                array_Rx_LOS[0, counter] = cmath.exp(1j * self.k * self.dis * (
                                        x * math.sin(phi_Rx_LOS) * math.cos(theta_Rx_LOS) + y * math.sin(theta_Rx_LOS)))
                                counter += 1

                    #  Link Attentuation (LOS)
                    L_dB_LOS_2 = -20 * math.log10(4 * math.pi / self.waveLength) - 10 * self.n_LOS * (
                            1 + self.b_LOS * ((self.Frequency - self.f0) / self.f0)) * math.log10(
                        d_RIS_R) - np.random.randn() * self.sigma_LOS
                    L_LOS_2 = 10 ** (L_dB_LOS_2 / 10)
                    # Generate g (Pure LOS)
                    g_LOS = math.sqrt(self.Gain * (math.cos(theta_RIS_Rx_LOS)) ** (2 * self.q)) * math.sqrt(
                        L_LOS_2) * np.transpose(array_RIS_Rx_LOS) * array_Rx_LOS * cmath.exp(
                        1j * np.random.rand() * 2 * math.pi)
                else:
                    g_LOS = 0



                #  Generate New Clusters/Scatters (Step 2) - Ensure that All Scattters are above ground
                # Generate Clusters/Sub-rays, Azimuth/Elevation Departure Angles and Cluster Distances
                # Number of Clusters - lambda_p was defined earlier: 1.8/1.9
                if self.G_c_num==0:
                    C_2 = max([1, np.random.poisson(lam=self.lambda_p)])  # Poisson distributed
                    self.G_c_num=C_2
                else:
                    C_2=self.G_c_num

                if len(self.G_s_list) == 0:
                    # Number of Sub-rays per Cluster
                    S_2 = np.random.randint(30, size=C_2)

                    while True:
                        if sum(S_2) != 0:
                            break
                        else:
                            S_2 = np.random.randint(30, size=C_2)
                    self.G_s_list=S_2
                else:
                    S_2=self.G_s_list

                for generate in range(100):
                    # Azimuth/Elevation Departure Angles
                    phi_Tx_2 = []
                    theta_Tx_2 = []
                    phi_av_2 = np.zeros(C_2)
                    theta_av_2 = np.zeros(C_2)

                    isCoordinatesDone = False
                    while True:
                        if isCoordinatesDone:
                            break
                        for counter in range(C_2):
                            phi_av_2[counter] = np.random.rand() * 90 - 45  # mean azimuth
                            theta_av_2[counter] = np.random.rand() * 90 - 45  # mean elevation
                            #  (needs update & most of the clusters are above ceiling)

                            # for outdoors this might not be a big concern
                            # cluster angles: First S_2(1) belongs to Cluster 1, Next S_2(2) belongs to Cluster 2....
                            new_phi_Tx = np.log(
                                np.random.rand(1, S_2[counter]) / np.random.rand(1,
                                                                                 S_2[
                                                                                     counter])).flatten() * math.sqrt(
                                25 / 2) + phi_av_2[counter]
                            phi_Tx_2.extend(new_phi_Tx)
                            new_theta_Tx = np.log(
                                np.random.rand(1, S_2[counter]) / np.random.rand(1,
                                                                                 S_2[
                                                                                     counter])).flatten() * math.sqrt(
                                25 / 2) + theta_av_2[counter]
                            theta_Tx_2.extend(new_theta_Tx)

                        # Cluster Distances
                        #  Cluster distances uniform [1,d_T_RIS]
                        #  can be modified later
                        a_c = 1 + np.random.rand(C_2) * (d_T_RIS - 1)

                        # Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
                        if len(self.G_c)==0:
                            Coordinates = np.zeros(shape=[C_2, 3])  # for Clusters
                            # Correction on Cluster Locations for Indoors
                            for counter in range(C_2):
                                loop = 1
                                while True:
                                    if self.Scenario == 1:
                                        Coordinates[counter, :] = [
                                            x_RIS - a_c[counter] * math.cos(theta_av_2[counter]) * math.sin(
                                                phi_av_2[counter]),
                                            y_RIS - a_c[counter] * math.cos(theta_av_2[counter]) * math.sin(
                                                phi_av_2[counter]),
                                            z_RIS + a_c[counter] * math.sin(theta_av_2[counter])]
                                    elif self.Scenario == 2:
                                        Coordinates[counter, :] = [
                                            x_RIS - a_c[counter] * math.cos(theta_av_2[counter]) * math.cos(
                                                phi_av_2[counter]),
                                            y_RIS + a_c[counter] * math.cos(theta_av_2[counter]) * math.sin(
                                                phi_av_2[counter]),
                                            z_RIS + a_c[counter] * math.sin(theta_av_2[counter])]
                                    if Coordinates[counter, 2] < 0:
                                        isCoordinatesDone = True
                                        self.G_c = Coordinates
                                        break
                                    else:
                                        if loop > 10:
                                            isCoordinatesDone = False
                                            break
                                        a_c[counter] = 0.8 * a_c[counter]
                                        # Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                                        # in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
                                        loop += 1

                        else:
                            Coordinates=self.G_c
                            isCoordinatesDone=True

                    # Plot Scatterers
                    a_c_rep_2 = []
                    for counter in range(C_2):
                        a_c_rep_2.extend(np.full((1, S_2[counter]), a_c[counter]).flatten())
                    if len(self.G_s) == 0:
                        Coordinates2 = np.zeros(shape=[sum(S_2), 3])  # for Scatterers
                        for counter2 in range(sum(S_2)):
                            if self.Scenario == 1:
                                Coordinates2[counter2, :] = [
                                    x_RIS - a_c_rep_2[counter2] * math.cos(theta_Tx_2[counter2]) * math.sin(
                                        phi_Tx_2[counter2]),
                                    y_RIS - a_c_rep_2[counter2] * math.cos(theta_Tx_2[counter2]) * math.cos(
                                        phi_Tx_2[counter2]),
                                    z_RIS + a_c_rep_2[counter2] * math.sin(theta_Tx_2[counter2])]
                            elif self.Scenario == 2:
                                Coordinates2[counter2, :] = [
                                    x_RIS - a_c_rep_2[counter2] * math.cos(theta_Tx_2[counter2]) * math.cos(
                                        phi_Tx_2[counter2]),
                                    y_RIS + a_c_rep_2[counter2] * math.cos(theta_Tx_2[counter2]) * math.sin(
                                        phi_Tx_2[counter2]),
                                    z_RIS + a_c_rep_2[counter2] * math.sin(theta_Tx_2[counter2])]
                        self.G_s = Coordinates2
                    else:
                        Coordinates2 = self.G_s

                    if len(self.G_i)==0:
                        ignore = []
                        for counter2 in range(sum(S_2)):
                            if Coordinates2[counter2, 2] < 0:  # only underground scatterers
                                ignore.extend([int(counter2)])

                        # updated indices
                        indices_2 = np.setdiff1d(np.arange(1, sum(S_2), 1), ignore)
                        self.G_i=indices_2
                    else:
                        indices_2=self.G_i
                    M_new_2 = len(indices_2)

                    # if you want to revert back Version 1.1 from Version 1.2, simply set
                    # indices_2=1:sum(S_2);
                    # M_new=sum(S_2)

                    # Necessary Loop to have at least one scatterer
                    if M_new_2 > 0:
                        break


                # Calculate Array Response (Step 4)
                # Array Response Calculation
                array_2 = np.zeros([sum(S_2), self.N], dtype=complex)
                for counter in indices_2:
                    counter2 = 0
                    for x in range(round(math.sqrt(self.N))):
                        for y in range(round(math.sqrt(self.N))):
                            array_2[counter, counter2] = cmath.exp(1j * self.k * self.dis * (
                                    x * math.sin(theta_Tx_2[counter]) + y * math.sin(phi_Tx_2[counter]) * math.cos(
                                theta_Tx_2[counter])))
                            counter2 += 1

                # Calculate Link Lengths (Step 3) and Path Loss (Step 5)
                phi_Rx_cs = np.zeros(sum(S_2))
                theta_Rx_cs = np.zeros(sum(S_2))
                b_cs_2 = np.zeros(sum(S_2))
                d_cs_2 = np.zeros(sum(S_2))

                for counter in indices_2:
                    b_cs_2[counter] = np.linalg.norm(
                        self.Rx_xyz - Coordinates2[counter, :])  # Distance between Scatterer and RIS
                    d_cs_2[counter] = a_c_rep_2[counter] + b_cs_2[counter]  # Total distance Tx-Scatterer-RIS

                    # AoA angles of Rx for NLOS RIS-Rx channel in an Outdoor
                    phi_av_Rx_NLOS = np.random.rand() * 180 - 90
                    theta_av_Rx_NLOS = np.random.rand() * 180 - 90

                    phi_Rx_cs[counter] = math.log(np.random.rand() / np.random.rand()) * math.sqrt(
                        25 / 2) + phi_av_Rx_NLOS
                    theta_Rx_cs[counter] = math.log(np.random.rand() / np.random.rand()) * math.sqrt(
                        25 / 2) + theta_av_Rx_NLOS

                array_Rx_cs = np.zeros([sum(S_2), self.Nr], dtype=complex)
                if self.ArrayType == 1:
                    for counter in indices_2:
                        counter2 = 0
                        for x in range(self.Nr):
                            array_Rx_cs[counter, counter2] = cmath.exp(
                                1j * self.k * self.dis * (
                                        x * math.sin(phi_Rx_cs[counter]) * math.cos(theta_Rx_cs[counter])))
                            counter2 += 1

                elif self.ArrayType == 2:
                    for counter in indices_2:
                        counter2 = 0
                        for x in range(round(math.sqrt(self.Nr))):
                            for y in range(round(math.sqrt(self.Nr))):
                                array_Rx_cs[counter, counter2] = cmath.exp(1j * self.k * self.dis * (
                                        x * math.sin(phi_Rx_cs[counter]) * math.cos(
                                    theta_Rx_cs[counter]) + y * math.sin(
                                    theta_Rx_cs[counter])))
                                counter2 += 1
                # Generate g
                g_NLOS = np.zeros([self.N, self.Nr], dtype=complex)
                for counter in indices_2:
                    X_sigma_2 = np.random.rand() * self.sigma_NLOS
                    Lcs_dB_2 = -20 * math.log10(4 * math.pi / self.waveLength) - 10 * self.n_NLOS * (
                            1 + self.b_NLOS * ((self.Frequency - self.f0) / self.f0)) * math.log10(
                        d_cs_2[counter]) - X_sigma_2
                    Lcs_2 = 10 ** (Lcs_dB_2 / 10)
                    beta_2 = ((np.random.randn() + 1j * np.random.randn()) / math.sqrt(2))
                    g_G = beta_2 * math.sqrt(Lcs_2) * (
                        cmath.sqrt(self.Gain * (math.cos(theta_Tx_2[counter])) ** (2 * self.q)))
                    g_t = g_G * np.transpose([array_2[counter, :]]) * [array_Rx_cs[counter, :]]

                    g_NLOS = g_NLOS + g_t
                g_NLOS = g_NLOS * math.sqrt(1 / M_new_2)
                g = g_NLOS + g_LOS
            G_new.append(g)
            self.G = G_new
            return self.G


if __name__ == '__main__':
    Tx_pos = np.array([0.0, 25.0, 2.0])
    Rx_pos = np.array([38.0, 48.0, 1.0])
    RIS_pos = np.array([40.0, 50.0, 2.0])
    sim = SimRIS(2, 2, 28, 2, 100, 4, 4, 10, Tx_pos, Rx_pos, RIS_pos)
    H = sim.generateH()
    for i in range(100):
        G = sim.generateG()
        print(i, sum(sum(G)))

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
