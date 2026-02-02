imf = dist_funcs.imf.Salpeter(0.08 * u.Msun, 100 * u.Msun, alpha=2.3)

Mtot = 1e6 * u.Msun

population = DiscreteStellarPopulation(imf, Mtot, metallicity=0.0, fbin=0.0)

out_times = np.logspace(1, 8, 200) * u.yr

out_qtys = [population.sn_count]
