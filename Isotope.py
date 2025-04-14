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

    def __init__(self, __Z__, __A__):
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
                            self.XS[i].append(float(XS_vals[i+1].strip()))
                i+= 1
        self.get_outgoing_reaction_isotopes()
        self.read_decay_data("XS_data/" + str(self.ZAID) + "-decay.csv")

    def read_fission_yields(self, fission_yield_file):
        print("read fission yields")
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
                        self.ify_isotopes_thermal.append(float(row[0]) * 1000 + float(row[1]))
                        self.ify_probs_thermal.append(float(row[7].strip()))
                    if(row[9]):
                        self.ify_isotopes_fast.append(float(row[0]) * 1000 + float(row[1]))
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
        print("read decay data")
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

    # given a flux distribution in energy, find the reaction rates per atom for each reaction of this Isotope, independent of energy
    def find_RRA(self, phi, phi_Egrid):
        print("integrate phi with the reaction microscopic cross section distribution")
        self.RRA = [0] * len(self.MT)