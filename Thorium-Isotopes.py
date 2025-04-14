import numpy as np
import math
import matplotlib.pyplot as plt
import Isotope

# units:
# Energy: eV, length: cm, time: s

######### DATA LOADING #########

# using GEANT4 Mix 3 (Th-232 and U-233 mixture)
phi = [2E11, 4E11, 1.05E11, 5E10, 1E10, 7E9, 2E9, 8E8, 4E8, 2E7, 2E5, 8E3] # units of neutrons per cm^2 per second per eV
phi_Egrid = [3E-2, 8E-2, 1, 2, 1E2, 2E3, 1E5, 3E5, 1E6, 3E6, 1E7, 3E7] # units of eV

isotope_ZAIDs = [90232]
isotope_concs = [0 for _ in range(len(isotope_ZAIDs))]
isotope_concs[0] = 1
isotope_data = []

for iso_id,iso in enumerate(isotope_ZAIDs):
    isotope_data.append(Isotope.Isotope(math.floor(iso/1000), iso % 1000, phi, phi_Egrid))

    # plot the cross sections for each isotope read in
    '''for i in range(len(isotope_data[iso_id].XS)):
        plt.loglog(isotope_data[iso_id].XS_Egrid[i], isotope_data[iso_id].XS[i])
    plt.title("Cross section data for " + str(iso))
    plt.xlabel("E (eV)")
    plt.ylabel("Cross section (barns)")
    plt.figure()'''

######### CONCENTRATION SIMULATION IN TIME #########

# do time evolution (in units of seconds)
T_max = 100
dt = 0.5
for t in np.array(range(int(T_max/dt))) * dt:
    # find the losses for each isotope
    for i, isotope in enumerate(isotope_data):
        if isotope_concs[i] == 0: # ignore any isotopes that have 0 concentration
            continue
        total_N_loss = 0

        # reactions
        for ri, mt in enumerate(isotope.MT):
            if mt == 18: # special case for fission reaction
                RR_thermal_fis = isotope_concs[i] * isotope.RRA[ri][0] #  reaction rate for thermal fission
                RR_fast_fis = isotope_concs[i] * isotope.RRA[ri][1] #  reaction rate for fast fission

                outgoing_isotopes_ZAID_thermal = isotope.ify_isotopes_thermal
                outgoing_isotopes_multiplicity_thermal = isotope.ify_probs_thermal
                outgoing_isotopes_ZAID_fast = isotope.ify_isotopes_fast
                outgoing_isotopes_multiplicity_fast = isotope.ify_probs_fast

                for new_iso_ZAID in outgoing_isotopes_ZAID_thermal:
                    if new_iso_ZAID in isotope_ZAIDs:
                        isotope_concs[isotope_ZAIDs.index(new_iso_ZAID)] += RR_thermal_fis * outgoing_isotopes_multiplicity_thermal[i] * dt
                    else:
                        True
                        #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                    total_N_loss += RR_thermal_fis * dt
                for new_iso_ZAID in outgoing_isotopes_ZAID_fast:
                    if new_iso_ZAID in isotope_ZAIDs:
                        isotope_concs[isotope_ZAIDs.index(new_iso_ZAID)] += RR_fast_fis * outgoing_isotopes_multiplicity_fast[i] * dt
                    else:
                        True
                        #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                    total_N_loss += RR_fast_fis * dt
            else: # every other reaction
                RR = isotope_concs[i] * isotope.RRA[ri] #  reaction rate for isotope i and reaction ri
                outgoing_isotopes_ZAID = isotope.reaction_isotopes[ri]
                for new_iso_ZAID in outgoing_isotopes_ZAID:
                    if new_iso_ZAID in isotope_ZAIDs:
                        isotope_concs[isotope_ZAIDs.index(new_iso_ZAID)] += RR * dt
                    else:
                        True
                        #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
                    total_N_loss += RR * dt


                    
        # decays
        N_loss_decay = isotope_concs[i] * isotope.lambda_t * dt
        total_N_loss -= N_loss_decay
        for new_iso_index, new_iso_ZAID in enumerate(isotope.decay_isotopes):
            if new_iso_ZAID in isotope_ZAIDs:
                isotope_concs[isotope_ZAIDs.index(new_iso_ZAID)] += N_loss_decay * isotope.chi_d[new_iso_index]
            else:
                True
                #print("No data for isotope " + str(new_iso_ZAID) + ". It is ignored")
        isotope_concs[i] -= total_N_loss
            
    print("isotope concentrations: ")
    print(isotope_concs)

plt.loglog(phi_Egrid, phi)
plt.title("flux plot")
plt.xlabel("Energy (eV)")
plt.ylabel("flux in neutrons/cm^2/s/MeV")
plt.show()