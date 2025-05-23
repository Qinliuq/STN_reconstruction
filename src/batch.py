from netpyne import batch
import numpy as np

params = {
    'gHCN_scale': np.linspace(0.5, 1.5, 20),
    'gCaT_scale': np.linspace(0.65, 2, 20),
}

b = batch.Batch('src/cfg.py', 'src/netParams.py', params=params, initCfg={})
b.batchLabel = 'gCaT-gHCN' # for Kir search it was 'kir-gCaT-gHCN'
b.saveFolder = b.batchLabel
b.method = 'grid'

b.runCfg = {'type': 'mpi_bulletin',
            'script': 'src/init.py',
            'skip': False}
b.run()