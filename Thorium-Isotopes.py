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
# somehow 'fission' makes it into the isotope list, fix that


# units:
# Energy: eV, length: cm, time: s

class IsotopeEvolution:

    # using GEANT4 Mix 3 (Th-232 and U-233 mixture) from https://www.sciencedirect.com/science/article/pii/S0029549321005811
    phi = [2E11, 4E11, 1.05E11, 5E10, 1E10, 7E9, 2E9, 8E8, 4E8, 2E7, 2E5, 8E3] # units of neutrons per cm^2 per second per eV
    phi_scale = 0.1 # amount to scale flux by
    for i in range(len(phi)):
        phi[i] *= phi_scale
    phi_Egrid = [3E-2, 8E-2, 1, 2, 1E2, 2E3, 1E5, 3E5, 1E6, 3E6, 1E7, 3E7] # units of eV

    #isotope_ZAIDs = [90232, 92233, 91233, 90233, 88228, 89228]
    isotope_ZAIDs = [90232, 92233]
    isotope_concs = np.zeros(len(isotope_ZAIDs))
    isotope_concs[0] = 3.2
    isotope_concs[1] = 0.85
    isotope_data = []

    untracked_isotopes = []
    untracked_concs = []
    untracked_concs_history = np.array([])

    T_max = 100000 * 3600.0 + 3600 * 24 * 365 * 10
    dt = 1000000.0
    Nt = int(T_max/dt)
    power_cutoff_time = 100000 * 3600.0
    #power_cutoff_time = 1E30
    power_cutoff_index = int(min(T_max, power_cutoff_time) / dt)

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
            print("Time: ", t)
            ti = int(t/self.dt)
            conc_additions = [0] * len(self.isotope_concs) # tally change in isotope concentrations in this timestep
            # find the losses for each isotope
            for i, isotope in enumerate(self.isotope_data):
                if self.isotope_concs[i] == 0: # ignore any isotopes that have 0 concentration
                    continue
                total_N_loss = 0

                # reactions
                if isotope.hasReactions and t < self.power_cutoff_time:
                    for ri, mt in enumerate(isotope.MT):
                        if mt == 18 and isotope.doesFission: # special case for fission reaction
                            reactions_thermal_fis = min(self.isotope_concs[i] * isotope.RRA[ri][0] * self.dt, self.isotope_concs[i] - total_N_loss) #  reactions for thermal fission
                            total_N_loss += reactions_thermal_fis
                            reactions_fast_fis = min(self.isotope_concs[i] * isotope.RRA[ri][1] * self.dt, self.isotope_concs[i] - total_N_loss) #  reactions for fast fission
                            total_N_loss += reactions_fast_fis

                            outgoing_isotopes_ZAID_thermal = isotope.ify_isotopes_thermal
                            outgoing_isotopes_multiplicity_thermal = isotope.ify_probs_thermal
                            outgoing_isotopes_ZAID_fast = isotope.ify_isotopes_fast
                            outgoing_isotopes_multiplicity_fast = isotope.ify_probs_fast

                            for new_iso_i, new_iso_ZAID in enumerate(outgoing_isotopes_ZAID_thermal):
                                if new_iso_ZAID in self.isotope_ZAIDs:
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += reactions_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i]
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += reactions_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i]
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(reactions_thermal_fis * outgoing_isotopes_multiplicity_thermal[new_iso_i])
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                            for new_iso_i, new_iso_ZAID in enumerate(outgoing_isotopes_ZAID_fast):
                                if new_iso_ZAID in self.isotope_ZAIDs:
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += reactions_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i]
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += reactions_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i]
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(reactions_fast_fis * outgoing_isotopes_multiplicity_fast[new_iso_i])
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                        else: # every other reaction
                            reactions = min(self.isotope_concs[i] * isotope.RRA[ri] * self.dt, self.isotope_concs[i] - total_N_loss) #  reaction rate for isotope i and reaction ri
                            total_N_loss += reactions
                            outgoing_isotopes_ZAID = isotope.reaction_isotopes[ri]
                            if outgoing_isotopes_ZAID == 'fission':
                                print("WHAT!?!?!")
                            for new_iso_ZAID in outgoing_isotopes_ZAID:
                                if new_iso_ZAID in self.isotope_ZAIDs: # if isotope is already being tracked in simulation
                                    conc_additions[self.isotope_ZAIDs.index(new_iso_ZAID)] += reactions
                                else:
                                    if(new_iso_ZAID in self.untracked_isotopes or not self.addIsotope(new_iso_ZAID)): # isotope does not have data on file
                                        self.untracked_concs[self.untracked_isotopes.index(new_iso_ZAID)] += reactions
                                    else: # isotope data can be loaded, is added to the simulation and concentrattion updated
                                        conc_additions.append(reactions)
                                    #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")



                # decays
                N_loss_decay = min(self.isotope_concs[i] * isotope.lambda_t * self.dt, self.isotope_concs[i] - total_N_loss) # equation assumes the concentration does not change much within timestep, ensures N_loss doesn't exceed current concentration
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

                
            #print("isotope concentrations: ")
            #print(self.isotope_concs)

        plt.loglog(self.phi_Egrid, self.phi)
        plt.title("flux plot")
        plt.xlabel("Energy (eV)")
        plt.ylabel("flux in neutrons/cm^2/s/MeV")
        plt.figure()

        # compile concentrations of isotopes above threshold concentration in a separate list
        thresholded_concs_history = np.array([])
        thresholded_ZAIDs = []
        threshold_conc = 2E-3
        for i, ZAID in enumerate(self.isotope_ZAIDs):
            if self.isotope_concs_history[self.power_cutoff_index,i] > threshold_conc:
                thresholded_ZAIDs.append(ZAID)
                if len(thresholded_concs_history):
                    thresholded_concs_history = np.append(thresholded_concs_history, self.isotope_concs_history[:,i].reshape((self.Nt+1,1)), axis = 1)
                else:
                    thresholded_concs_history = self.isotope_concs_history[:,i].reshape((self.Nt+1,1))

        # compile concentrations of actinides in a separate list
        actinide_concs_history = np.array([])
        actinide_ZAIDs = []
        actinide_threshold_conc = 2E-3
        for i, ZAID in enumerate(self.isotope_ZAIDs):
            if ZAID > 80000 and self.isotope_concs_history[self.power_cutoff_index,i] > actinide_threshold_conc:
                actinide_ZAIDs.append(ZAID)
                if len(actinide_concs_history):
                    actinide_concs_history = np.append(actinide_concs_history, self.isotope_concs_history[:,i].reshape((self.Nt+1,1)), axis = 1)
                else:
                    actinide_concs_history = self.isotope_concs_history[:,i].reshape((self.Nt+1,1))

        # compile concentrations for notable isotopes
        notable_isotopes = [90232, 92233, 92232, 91231, 91233]
        N_notable_scaling = [1,1,1000,100,100] # scale the number densities for plotting
        notable_concs_history = np.array([])
        notable_ZAIDs = []
        for i, ZAID in enumerate(self.isotope_ZAIDs):
            if ZAID in notable_isotopes:
                notable_ZAIDs.append(ZAID)
                if len(notable_concs_history):
                    notable_concs_history = np.append(notable_concs_history, self.isotope_concs_history[:,i].reshape((self.Nt+1,1)), axis = 1)
                else:
                    notable_concs_history = self.isotope_concs_history[:,i].reshape((self.Nt+1,1))

        isotope_activity_history = self.isotope_concs_history  * [iso.lambda_t for iso in self.isotope_data]
        isotope_activities = isotope_activity_history[-1]

        # compile concentrations of isotopes with activity above threshold concentration in a separate list
        thresholded_activity_history = np.array([])
        thresholded_activity_ZAIDs = []
        threshold_activity = 1E-16
        for i, ZAID in enumerate(self.isotope_ZAIDs):
            if isotope_activity_history[-1,i] > threshold_activity:
                thresholded_activity_ZAIDs.append(ZAID)
                if len(thresholded_activity_history):
                    thresholded_activity_history = np.append(thresholded_activity_history, isotope_activity_history[:,i].reshape((self.Nt+1,1)), axis = 1)
                else:
                    thresholded_activity_history = isotope_activity_history[:,i].reshape((self.Nt+1,1))

        T_min = 0#.9 * self.power_cutoff_time
        conc_min = 1E-6
        conc_max = 1
        # plot all isotopes
        plt.semilogy(T_grid, self.isotope_concs_history)
        plt.legend(self.isotope_ZAIDs, bbox_to_anchor=(0.95, 1.00))
        plt.title("Concentrations of all tracked isotopes")
        plt.xlim([T_min,self.T_max])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Concentration")
        plt.figure()

        # plot all thresholded isotopes
        plt.semilogy(T_grid, thresholded_concs_history)
        plt.legend(thresholded_ZAIDs, bbox_to_anchor=(0.95, 0.5))
        plt.title("Concentrations of all tracked isotopes (with threshold concentration)")
        plt.xlim([T_min,self.T_max])
        plt.ylim([conc_min, conc_max])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Concentration")
        plt.figure()

        # plot all actinides
        plt.semilogy(T_grid, actinide_concs_history)
        plt.legend(actinide_ZAIDs, bbox_to_anchor=(0.95, 0.5))
        plt.title("Concentrations of isotopes with Z>80")
        plt.xlim([T_min,self.T_max])
        plt.ylim([conc_min, conc_max])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Concentration")
        plt.figure()

        # plot just notable (listed) isotopes
        color_list = ['purple', 'blue', 'black', 'green', 'red']
        for i in range(len(N_notable_scaling)):
            plt.plot(T_grid, notable_concs_history[:,i] * N_notable_scaling[i], color=color_list[i % len(color_list)])
        #plt.legend(notable_ZAIDs)
        plt.legend(["Th-232", "U-233", "U-232 (scaled by 1000)", "Pa-231 (scaled by 100)", "Pa-233 (scaled by 100)"])
        plt.title("Concentrations of notable isotopes")
        plt.xlim(0)
        plt.xlabel("Time (seconds)")
        plt.ylabel("Concentration")
        plt.figure()

        # plot all activity-thresholded isotopes
        plt.semilogy(T_grid, thresholded_activity_history)
        plt.legend(thresholded_activity_ZAIDs)
        plt.title("Activities of all tracked isotopes (with threshold activity)")
        plt.xlim([T_min,self.T_max])
        plt.ylim([1E-20, 1E-9])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Activity")
        plt.figure()

        # plot untracked isotopes
        #plt.loglog(T_grid, self.untracked_concs_history)
        #plt.legend(self.untracked_isotopes)
        #plt.title("Concentrations of all isotopes, tracked and untracked")
        plt.show()

        max_list_length = 10

        # print out isotopes with the most concentration at the end of the simulation
        print("==========Most prevalent isotopes at simulation completion==========")
        isotopes_sorted_indices = np.flip(np.argsort(self.isotope_concs))
        l = 0
        for i in isotopes_sorted_indices:
            if l >= max_list_length:
                break
            if self.isotope_ZAIDs[i] != 'fission' and int(self.isotope_ZAIDs[i]) > 0000:
                print("Isotope " + str(self.isotope_ZAIDs[i]) + ", conc = " + str(self.isotope_concs[i]))
                l+=1
        
        # print out actinides with the most concentration at the end of the simulation
        print("==========Most prevalent actinides at simulation completion==========")
        l=0
        for i in isotopes_sorted_indices:
            if l >= max_list_length:
                break
            if self.isotope_ZAIDs[i] != 'fission' and int(self.isotope_ZAIDs[i]) > 80000:
                print("Isotope " + str(self.isotope_ZAIDs[i]) + ", conc = " + str(self.isotope_concs[i]))
                l+=1

        
        # print out isotopes with the most activity at the end of the simulation
        print("==========Isotopes with most activity at simulation completion==========")
        isotopes_sorted_indices = np.flip(np.argsort(isotope_activities))
        l=0
        for i in isotopes_sorted_indices:
            if l >= max_list_length:
                break
            if self.isotope_ZAIDs[i] != 'fission' and int(self.isotope_ZAIDs[i]) > 0000:
                print("Isotope " + str(self.isotope_ZAIDs[i]) + ", activity = " + str(isotope_activities[i]))
                l+=1
        
        # print out actinides with the most activity at the end of the simulation
        print("==========Actinides with most activity at simulation completion==========")
        l=0
        for i in isotopes_sorted_indices:
            if l >= max_list_length:
                break
            if self.isotope_ZAIDs[i] != 'fission' and int(self.isotope_ZAIDs[i]) > 80000:
                print("Isotope " + str(self.isotope_ZAIDs[i]) + ", conc = " + str(isotope_activities[i]))
                l+=1

        # print out untracked isotopes with the highest concentrations (priority for adding in nuclear data)
        '''untracked_sorted_indices = np.argsort(self.untracked_concs)
        for i in untracked_sorted_indices:
            if self.untracked_isotopes[i] != 'fission' and int(self.untracked_isotopes[i]) > 0000:
                print("Isotope " + str(self.untracked_isotopes[i]) + ", conc = " + str(self.untracked_concs[i]))'''

sim = IsotopeEvolution()
sim.doEvolution()