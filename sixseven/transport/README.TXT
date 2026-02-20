Transport Module 

transport_simple.py is the main script for compositional transport. 
The only function you should need from this file is:

    transport_step(X, structure, dt)

Import it at the top of your script with:

    from sixseven.transport.test_transport_new import transport_step

transport_step() takes the current mass fractions X, the current stellar structure, and a timestep
dt (in seconds). It returns the updated mass fractions after one timestep of mixing.

This function includes:
- Convective mixing (MLT-style diffusion)
- Semiconvection (Which I'm pretty sure won't happen anyway since we are doing all or nothing as Emily mentioned)
- Thermohaline mixing

X should be a NumPy array of shape (n_shells, n_species).

The structure dictionary must contain:
- m              (mass coordinate)
- Hp             (pressure scale height)
- v_mlt          (convective velocity)
- is_convective  (boolean array)
- grad_rad
- grad_ad
- grad_mu
- K              (thermal conductivity)
- Cp             (specific heat)
- rho            (density)
- T              (temperature)

transport_step() ONLY updates composition.
It does NOT:
- call nuclear burning
- compute temperature
- compute density
- update internal energy


The diffusion scheme is explicit. If the model becomes unstable, the timestep may need to be reduced.
