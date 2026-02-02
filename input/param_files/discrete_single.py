# Build a single-star population with a Salpeter IMF
population = zams_single(IMF=imf.Salpeter, mass=1e5*u.Msun, metals=0, rot=0*u.km/u.s)

# Select what evolution parameters we want
evo_params = {}

# Select what output parameters we want
output_params = {'sn':['E','m','m_Z']}


# Select the times to emit an output (1000 years to 100 Myr, in 100 log-spaced bins)
times = np.logspace(-3,2,100)*u.Myr

# Evolve the population to build the yield tables
yields = evolve_population(population, evolution=evo_params, outputs=output_params, interval=times)
