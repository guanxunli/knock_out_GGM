import numpy as np
import pandas as pd
from sergio import sergio
np.random.seed(1)

sim = sergio(number_genes=15, number_bins = 1, number_sc = 3000, 
             noise_params = 1, decays=18, sampling_state=50,
             noise_type='dpd')
sim

sim.build_graph('targetFile_2.py', 'regFile_2.py', shared_coop_state=2)
np.random.seed(1)
sim.simulate()

np.random.seed(1)
expr = sim.getExpressions()
expr = np.concatenate(expr, axis = 1)
expr

np.random.seed(1)
count_matrix = sim.convert_to_UMIcounts(expr)
np.savetxt('simulationOutput_2.csv',count_matrix, delimiter=',')