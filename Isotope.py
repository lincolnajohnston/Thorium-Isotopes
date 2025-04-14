import numpy as np
import math
import csv

class Isotope:
    # hardcoded for now, read in from file later
    # MT = [16,17,18,22,24,28,41,102,103,107,112] these are the reactions that we will consider, not all isotopes will have cross sections for all of these reactions though
    '''XS = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
    XS_Egrid = [[0.01,0.1,1],[0.01,0.1,1],[0.01,0.1,1],[0.01,0.1,1]]
    reaction_isotopes = [['Th-231'],['Th-230'],['Th-231'], ['fission], ['Th-233']] # cchange these to ZAIDs instead of strings
    ify_isotopes_thermal = []
    ify_probs_thermal = []
    ify_isotopes_fast = []
    ify_probs_fast = []
    lambda_t = 1.569E-18
    decay_isotopes = [88228]
    chi_d = [1]'''

    def __init__(self, __Z__, __A__, phi, phi_Egrid):
        self.Z = __Z__
        self.A = __A__
        self.ZAID = 1000*self.Z + self.A
        print("read in " + str(self.ZAID) + " nuclear data")


        xs_datafile = "XS_data/" + str(self.ZAID) + "-xs.csv"
        with open(xs_datafile, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            i = 0
            for row in reader:
                if i == 0:
                    self.MT = []
                    self.XS = []
                    self.XS_Egrid = []
                    self.numMT = len(row) - 1
                    self.XS = [[] for _ in range(self.numMT)]
                    self.XS_Egrid = [[] for _ in range(self.numMT)]
                    for mt in row[:-1]:
                        cur_mt = int(mt.split('=')[1].split(':')[0].strip())
                        self.MT.append(cur_mt)
                        if(cur_mt == 18):
                            self.read_fission_yields("XS_data/" + str(self.ZAID) + "-ify.csv")
                if i > 2:
                    XS_vals = row[0].split(';')
                    for i in range(self.numMT):
                        if XS_vals[i+1].strip():
                            self.XS_Egrid[i].append(float(XS_vals[0].strip()))
                            self.XS[i].append(10E-24 * float(XS_vals[i+1].strip()))
                i+= 1
        self.get_outgoing_reaction_isotopes()
        self.read_decay_data("XS_data/" + str(self.ZAID) + "-decay.csv")
        self.find_RRA(phi, phi_Egrid)

    def read_fission_yields(self, fission_yield_file):
        self.ify_isotopes_thermal = []
        self.ify_probs_thermal = []
        self.ify_isotopes_fast = []
        self.ify_probs_fast = []
        with open(fission_yield_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            i = 0
            for row in reader:
                if i > 0 and row:
                    if(row[7]):
                        self.ify_isotopes_thermal.append(int(row[0]) * 1000 + int(row[1]))
                        self.ify_probs_thermal.append(float(row[7].strip()))
                    if(row[9]):
                        self.ify_isotopes_fast.append(int(row[0]) * 1000 + int(row[1]))
                        self.ify_probs_fast.append(float(row[9].strip()))
                i += 1


    def get_outgoing_reaction_isotopes(self):
        self.reaction_isotopes = [[] for i in range(self.numMT)]
        for i,mt in enumerate(self.MT):
            if mt == 16: # (n,2n)
                self.reaction_isotopes[i].append(self.Z * 1000 + (self.A - 1))
            elif mt == 17: # (n,3n)
                self.reaction_isotopes[i].append(self.Z * 1000 + (self.A - 2))
            elif mt == 18: # fission
                self.reaction_isotopes[i].append('fission')
            elif mt == 22: # (n,n+alpha)
                self.reaction_isotopes[i].append((self.Z-2) * 1000 + (self.A - 4))
                self.reaction_isotopes[i].append(2 * 1000 + 4) # alpha particle
            elif mt == 24: # (n,2n+alpha)
                self.reaction_isotopes[i].append((self.Z-2) * 1000 + (self.A - 5))
                self.reaction_isotopes[i].append(2 * 1000 + 4) # alpha particle
            elif mt == 28: # (n,n+p)
                self.reaction_isotopes[i].append((self.Z-1) * 1000 + (self.A - 1))
                self.reaction_isotopes[i].append(1 * 1000 + 1) # H-1
            elif mt == 41: # (n,2n+p)
                self.reaction_isotopes[i].append((self.Z-1) * 1000 + (self.A - 2))
                self.reaction_isotopes[i].append(1 * 1000 + 1) # H-1
            elif mt == 102: # (n,gamma)
                self.reaction_isotopes[i].append((self.Z) * 1000 + (self.A + 1))
            elif mt == 103: # (n,p)
                self.reaction_isotopes[i].append((self.Z-1) * 1000 + (self.A))
                self.reaction_isotopes[i].append(1 * 1000 + 1) # H-1
            elif mt == 107: # (n,alpha)
                self.reaction_isotopes[i].append((self.Z-2) * 1000 + (self.A-3))
                self.reaction_isotopes[i].append(2 * 1000 + 4) # alpha particle
            elif mt == 112: # (n,p+alpha)
                self.reaction_isotopes[i].append((self.Z-3) * 1000 + (self.A-4))
                self.reaction_isotopes[i].append(2 * 1000 + 4) # alpha particle
                self.reaction_isotopes[i].append(1 * 1000 + 1) # H-1

    # create the lambda_t, decay_isotopes, and chi_d variables
    def read_decay_data(self, decay_file):
        self.decay_isotopes = []
        self.chi_d = []
        with open(decay_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            i = 0
            for row in reader:
                if i == 1:
                    half_life = float(row[14].strip())
                    self.lambda_t = math.log(2) / half_life
                    j = 0
                    while row[int(16 + 3*j)] and j < 3:
                        decay_type = row[int(16 + 3*j)]
                        self.decay_isotopes.append(self.get_new_isotope_from_reaction(decay_type))
                        self.chi_d.append(float(row[int(16 + 3*j + 1)]) / 100)
                        j += 1
                    break
                i += 1

    def get_new_isotope_from_reaction(self, reaction):
        if reaction == "A":
            dZ = -2
            dA = -4
        elif reaction == "B-":
            dZ = 1
            dA = 0
        elif reaction == "B+":
            dZ = -1
            dA = 0
        else:
            dZ = 0
            dA = 0
        return self.ZAID + dZ * 1000 + dA
    
    # binary search to get first value in ordered list over a threshold
    # got from here: https://stackoverflow.com/questions/71876716/get-index-of-first-value-above-threshold-in-ordered-list
    def first_over_ind(self, a, t, i=0) -> int:
        if a[0] > t:
            return i
        mid = int(len(a) / 2)
        if a[mid] > t:
            if a[mid-1] <= t:
                return i + mid
            return self.first_over_ind(a[:mid], t, i)
        return self.first_over_ind(a[mid:], t, i + mid)
    
    def interpolate_phi(self, phi, phi_Egrid, E):
        i_plus = self.first_over_ind(phi_Egrid, E)
        i_minus = i_plus - 1
        return phi[i_minus] + (E - phi_Egrid[i_minus]) * (phi[i_plus] - phi[i_minus]) / (phi_Egrid[i_plus] - phi_Egrid[i_minus]) # linear interpolation of phi values
    

    # given a flux distribution in energy, find the reaction rates per atom (RRA) for each reaction of this Isotope, independent of energy
    def find_RRA(self, phi, phi_Egrid):
        self.RRA = [0] * len(self.MT)
        for r in range(self.numMT):
            if self.MT[r] == 18:
                # find the indices of the cross section grid that are within the range of the flux distribution energy grid
                min_XS_Egrid_index_thermal = self.first_over_ind(self.XS_Egrid[r], phi_Egrid[0])
                max_XS_Egrid_index_thermal = self.first_over_ind(self.XS_Egrid[r], 1) # max out energy at 1 eV
                min_XS_Egrid_index_fast = max_XS_Egrid_index_thermal
                max_XS_Egrid_index_fast = self.first_over_ind(self.XS_Egrid[r], phi_Egrid[-1]) - 1

                # integrate over the thermal flux distribution * the microscopic cross section distribution, i is the index on the left side of the integration region
                integration_sum_thermal = 0
                for i in range(min_XS_Egrid_index_thermal, max_XS_Egrid_index_thermal):
                    El = self.XS_Egrid[r][i]
                    Er = self.XS_Egrid[r][i+1]
                    dE = Er - El
                    XS_avg = (self.XS[r][i] + self.XS[r][i+1]) / 2
                    phi_E = self.interpolate_phi(phi, phi_Egrid, (El + Er) / 2)
                    integration_sum_thermal += phi_E * XS_avg * dE
                # integrate over the fast flux distribution * the microscopic cross section distribution, i is the index on the left side of the integration region
                integration_sum_fast = 0
                for i in range(min_XS_Egrid_index_fast, max_XS_Egrid_index_fast):
                    El = self.XS_Egrid[r][i]
                    Er = self.XS_Egrid[r][i+1]
                    dE = Er - El
                    XS_avg = (self.XS[r][i] + self.XS[r][i+1]) / 2
                    phi_E = self.interpolate_phi(phi, phi_Egrid, (El + Er) / 2)
                    integration_sum_fast += phi_E * XS_avg * dE
                self.RRA[r] = [integration_sum_thermal, integration_sum_fast]
            else:
                # find the indices of the cross section grid that are within the range of the flux distribution energy grid
                min_XS_Egrid_index = self.first_over_ind(self.XS_Egrid[r], phi_Egrid[0])
                max_XS_Egrid_index = self.first_over_ind(self.XS_Egrid[r], phi_Egrid[-1]) - 1

                # integrate over the flux distribution * the microscopic cross section distribution, i is the index on the left side of the integration region
                integration_sum = 0
                for i in range(min_XS_Egrid_index, max_XS_Egrid_index):
                    El = self.XS_Egrid[r][i]
                    Er = self.XS_Egrid[r][i+1]
                    dE = Er - El
                    XS_avg = (self.XS[r][i] + self.XS[r][i+1]) / 2
                    phi_E = self.interpolate_phi(phi, phi_Egrid, (El + Er) / 2)
                    integration_sum += phi_E * XS_avg * dE
                self.RRA[r] = integration_sum