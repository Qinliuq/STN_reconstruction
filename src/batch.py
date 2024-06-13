from netpyne import batch
import numpy as np

params = {
    'gCaT_scale': np.linspace(1.0, 2.5, 3), # CaT coonductance increase
    'gSK_scale': np.linspace(1.0, 0.5, 3), # sKCa conductance decrease
}

b = batch.Batch('src/cfg.py', 'src/netParams.py', params=params, initCfg={})
b.batchLabel = 'gCaT-up_gSK-down'
b.saveFolder = b.batchLabel
b.method = 'grid'

b.runCfg = {'type': 'mpi_bulletin',
            'script': 'src/init.py',
            'skip': False}
b.run()