import numpy as np
import matplotlib.pyplot as plt

from sixseven.timestep.timestep import dyn_timestep
from sixseven.eos.eos_functions import *
from sixseven.nuclear.nuc_burn import *
from sixseven.transport.transport_simple import transport_step


def main():

    n_shells = 100
    dt = 1.0e14  # seconds 

    temps = np.linspace(1.5e7, 4e6, n_shells)
    density = np.linspace(160, 100, n_shells)
    # ben bug fix here: this array represents mass enclosed within a given radius
    # if the first element of this array is zero, then EoS breaks which makes burning break -> because the mass is zero! 
    m = np.linspace(1e-5, 1, n_shells) * 1e32 # multipled by this by a total mass in grams 

    structure = {
        "m": m,
        "Hp": np.ones(n_shells) * 1e9,
        "v_mlt": np.ones(n_shells) * 1e5,
        "is_convective": np.zeros(n_shells, dtype=bool),
        "grad_rad": np.ones(n_shells) * 0.4,
        "grad_ad": np.ones(n_shells) * 0.3,
        "grad_mu": np.ones(n_shells) * 0.01,
        "K": np.ones(n_shells) * 1e7,
        "Cp": np.ones(n_shells) * 1e8,
        "rho": density,
        "T": temps,
    }
    results = burn(temps=temps, rhos=density, time=dt, comps=None)

    comps = results 
    X_new = transport_step(comps, structure, dt)


    print(X_new)


main()
