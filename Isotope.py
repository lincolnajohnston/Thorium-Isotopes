import numpy as np
import math
import csv

class Isotope:
    # hardcoded for now, read in from file later
    MT = [16,17,18,22,24,28,41,102,103,107,112]
    '''XS = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
    XS_Egrid = [[0.01,0.1,1],[0.01,0.1,1],[0.01,0.1,1],[0.01,0.1,1]]
    reaction_isotopes = [['Th-231'],['Th-230'],['Th-231'], ['fission], ['Th-233']]
    ify_isotopes_thermal = []
    ify_probs_thermal = []
    ify_isotopes_fast = []
    ify_probs_fast = []
    chi_r = [[1],[1],[0.3,0.3,0.4],[1]]
    lambda_t = 1.569E-18
    decay_isotopes = ['Ra-228']
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
        reaction_isotopes = []
        # MT = 16
        reaction_isotopes.append[]
        
    # given a flux distribution in energy, find the reaction rates per atom for each reaction of this Isotope, independent of energy
    def find_RRA(self, phi, phi_Egrid):
        print("integrate phi with the reaction microscopic cross section distribution")
        self.RRA = [0] * len(self.MT)