import numpy as np
import math
import matplotlib.pyplot as plt
import os
import Isotope

#TODO:
# load in a bunch of cross sections and decay data:
#       -select MT 16,17,18,22,24,28,41,102,103,107,112 in JANIS data viewer, download csv, rename to ZZZAAA-xs.csv
#       -download fission yields from IAEA livechart, rename ZZZAAA-ify.csv
#       -download decay data by typing http://nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=AAA{element} in the url bar, rename to ZZAAA-decay.csv
# change isotope lists to dictionaries where appropriate


# units:
# Energy: eV, length: cm, time: s

class IsotopeEvolution:

    # using GEANT4 Mix 3 (Th-232 and U-233 mixture) from https://www.sciencedirect.com/science/article/pii/S0029549321005811
    phi = [2E11, 4E11, 1.05E11, 5E10, 1E10, 7E9, 2E9, 8E8, 4E8, 2E7, 2E5, 8E3] # units of neutrons per cm^2 per second per eV
    phi_Egrid = [3E-2, 8E-2, 1, 2, 1E2, 2E3, 1E5, 3E5, 1E6, 3E6, 1E7, 3E7] # units of eV

    #isotope_ZAIDs = [90232, 92233, 91233, 90233, 88228, 89228]
    isotope_ZAIDs = [90232, 90233]
    isotope_concs = np.zeros(len(isotope_ZAIDs))
    isotope_concs[0] = 1
    isotope_data = []

    untracked_isotopes = []
    untracked_concs = []
    untracked_concs_history = np.array([])

    T_max = 100000
    dt = 1000
    Nt = int(T_max/dt)

    def __init__(self):
        ######### DATA LOADING #########
        for iso_id,iso in enumerate(self.isotope_ZAIDs):
            self.isotope_data.append(Isotope.Isotope(math.floor(iso/1000), iso % 1000, self.phi, self.phi_Egrid))

        # plot the cross sections for each isotope read in
        '''for i in range(len(isotope_data[iso_id].XS)):
            plt.loglog(isotope_data[iso_id].XS_Egrid[i], isotope_data[iso_id].XS[i])
        plt.title("Cross section data for " + str(iso))
        plt.xlabel("E (eV)")
        plt.ylabel("Cross section (barns)")
        plt.figure()'''

    def isotopeDataExists(self, ZAID):
        decay_datafile = "XS_data/" + str(ZAID) + "-decay.csv"
        if os.path.exists(decay_datafile):
            return True

    # add isotope to list of isotopes being simulated, return True of done successfully
    def addIsotope(self, ZAID):
        if not self.isotopeDataExists(ZAID):
            if ZAID not in self.untracked_isotopes:
                self.untracked_isotopes.append(ZAID)
                self.untracked_concs.append(0)
                if len(self.untracked_concs_history):
                    self.untracked_concs_history = np.append(self.untracked_concs_history, np.zeros((self.Nt+1, 1)), axis = 1)
                else:
                    self.untracked_concs_history = np.zeros((self.Nt+1, 1))
            return False
        
        self.isotope_ZAIDs.append(ZAID)
        self.isotope_concs = np.append(self.isotope_concs, 0)
        self.isotope_concs_history = np.append(self.isotope_concs_history, np.zeros((self.Nt+1, 1)), axis = 1)
        self.isotope_data.append(Isotope.Isotope(math.floor(ZAID/1000), ZAID % 1000, self.phi, self.phi_Egrid))

        return True


######### CONCENTRATION SIMULATION IN TIME #########
    def doEvolution(self):
        # do time evolution (in units of seconds)
        T_grid = np.array(range(self.Nt + 1)) * self.dt
        self.isotope_concs_history = np.zeros((self.Nt + 1, len(self.isotope_ZAIDs)))
        self.isotope_concs_history[0,:] = self.isotope_concs # set the initial concentrations in the istope concentration history
        for t in T_grid[:-1]:
            ti = int(t/self.dt)
            conc_additions = [0] * len(self.isotope_concs) # tally change in isotope concentrations in this timestep
            # find the losses for each isotope
            for i, isotope in enumerate(self.isotope_data):
                if self.isotope_concs[i] == 0: # ignore any isotopes that have 0 concentration
                    continue
                total_N_loss = 0

                # reactions
                if isotope.hasReactions:
                    for ri, mt in enumerate(isotope.MT):
                        if mt == 18 and isotope.doesFission: # special case for fission reaction
                            RR_thermal_fis = self.isotope_concs[i] * isotope.RRA[ri][0] #  reaction rate for thermal fission
                            RR_fast_fis = self.isotope_concs[i] * isotope.RRA[ri][1] #  reaction rate for fast fission

                            outgoing_isotopes_ZAID_thermal = isotope.ify_isotopes_thermal
                            outgoing_isotopes_multiplicity_thermal = isotope.ify_probs_thermal
                            outgoing_isotopes_ZAID_fast = isotope.ify_isotopes_fast
                            outgoing_isotopes_multiplicity_fast = isotope.ify_probs_fast

                            for new_iso_i, new_iso_ZAID in enumerate(outgoing_isotopes_ZAID_thermal):
                                if new_iso_ZAID in self.isotope_ZAIDs:
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += RR_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i] * self.dt
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += RR_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i] * self.dt
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(RR_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i] * self.dt)
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                                total_N_loss += RR_thermal_fis * self.dt
                            for new_iso_i, new_iso_ZAID in enumerate(outgoing_isotopes_ZAID_fast):
                                if new_iso_ZAID in self.isotope_ZAIDs:
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += RR_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i] * self.dt
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += RR_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i] * self.dt
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(RR_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i] * self.dt)
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                                total_N_loss += RR_fast_fis * self.dt
                        else: # every other reaction
                            RR = self.isotope_concs[i] * isotope.RRA[ri] #  reaction rate for isotope i and reaction ri
                            outgoing_isotopes_ZAID = isotope.reaction_isotopes[ri]
                            for new_iso_ZAID in outgoing_isotopes_ZAID:
                                if new_iso_ZAID in self.isotope_ZAIDs: # if isotope is already being tracked in simulation
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += RR * self.dt
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += RR * self.dt
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(RR * self.dt)
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                                total_N_loss += RR * self.dt



                # decays
                N_loss_decay = self.isotope_concs[i] * isotope.lambda_t * self.dt
                total_N_loss += N_loss_decay
                for new_iso_index, new_iso_ZAID in enumerate(isotope.decay_isotopes):
                    if new_iso_ZAID in self.isotope_ZAIDs: # if isotope is already being tracked in simulation
                        conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += N_loss_decay * isotope.chi_d[new_iso_index]
                    else:
                        if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                            self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += N_loss_decay * isotope.chi_d[new_iso_index]
                        else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                            conc_additions.append(N_loss_decay * isotope.chi_d[new_iso_index])
                        #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                conc_additions[i] -= total_N_loss
            
            # update the isotope concentrations after all of the reactions/decays have been calculated
            for i in range(len(self.isotope_concs)):
                self.isotope_concs[i] += conc_additions[i]
            self.isotope_concs_history[ti + 1,:] = self.isotope_concs
            if self.untracked_concs:
                self.untracked_concs_history[ti,:] = self.untracked_concs

                
            print("isotope concentrations: ")
            print(self.isotope_concs)

        plt.loglog(self.phi_Egrid, self.phi)
        plt.title("flux plot")
        plt.xlabel("Energy (eV)")
        plt.ylabel("flux in neutrons/cm^2/s/MeV")
        plt.figure()

        plt.loglog(T_grid, self.isotope_concs_history)
        plt.legend(self.isotope_ZAIDs)
        plt.figure()

        plt.loglog(T_grid, self.untracked_concs_history)
        plt.legend(self.untracked_isotopes)
        plt.show()

        untracked_sorted_indices = np.argsort(self.untracked_concs)
        for i in untracked_sorted_indices:
            print("Isotope " + str(self.untracked_isotopes[i]) + ", conc = " + str(self.untracked_concs[i]))

sim = IsotopeEvolution()
sim.doEvolution()