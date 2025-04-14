import numpy as np
import math
import matplotlib.pyplot as plt
import Isotope

# using GEANT4 Mix 3 (Th-232 and U-233 mixture)
phi = [2E17, 4E17, 1.05E17, 5E16, 1E16, 7E15, 2E15, 8E14, 4E14, 2E13, 2E11, 8E9] # units of neutrons per cm^2 per second per MeV
phi_Egrid = [3E-8, 8E-8, 1E-6, 2E-6, 1E-4, 2E-3, 0.1, 0.3, 1, 3, 10, 30] # units of MeV

Th232 = Isotope.Isotope(90, 232)

for i in range(len(Th232.XS)):
    plt.loglog(Th232.XS_Egrid[i], Th232.XS[i])
plt.figure()

plt.loglog(phi_Egrid, phi)
plt.title("flux plot")
plt.xlabel("Energy (MeV)")
plt.ylabel("flux in neutrons/cm^2/s/MeV")
plt.show()