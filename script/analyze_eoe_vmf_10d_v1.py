# We use tensorflow version of scPhere: https://github.com/klarman-cell-observatory/scPhere
# but you can use the pytorch version: https://github.com/Ding-Group/scPhere
#

import pandas as pd
import numpy as np
from scphere.util.util import read_mtx
from scphere.util.trainer_dp import Trainer
from scphere.util.plot import plot_trace
from scphere.model.vae_dp import SCPHERE
import matplotlib.pyplot as plt

H_DIM = 128
Z_DIM = 10

x = read_mtx('/Users/jding/work/EoE/result/v2v3/eoe_var_gene.mtx')
x = x.transpose().todense()

model = SCPHERE(n_gene=x.shape[1], n_batch=[22, 3, 3],
                z_dim=Z_DIM, latent_dist='vmf', batch_invariant=False,
                observation_dist='nb', seed=0)

batch_p = pd.read_csv('/Users/jding/work/eoe/result/v2v3/eoe_batch_patient.tsv', header=None)
batch_h = pd.read_csv('/Users/jding/work/eoe/result/v2v3/eoe_batch_disease.tsv', header=None)
batch_l = pd.read_csv('/Users/jding/work/eoe/result/v2v3/eoe_batch_location.tsv', header=None)

batch = pd.concat([batch_p, batch_h, batch_l], axis=1)

trainer = Trainer(model=model, x=x, batch_id=batch.values, max_epoch=500,
                  mb_size=128, learning_rate=0.001)
trainer.train()

df = pd.DataFrame(list(zip(trainer.status['log_likelihood'],
                           trainer.status['kl_divergence'])),
                  columns=['log_likelihood', 'kl_divergence'])

df.to_csv('/Users/jding/work/eoe/output/trace_vmf_10d_500_epoch_v1.tsv',
          sep=' ', index=False)


z_mean = model.encode(x, batch.values)

save_path = '/Users/jding/work/eoe/output/'
np.savetxt(save_path +
           'latent_vmf_10d_500_epoch_v1.tsv',
           z_mean)

ll = model.get_log_likelihood(x, batch.values)
np.savetxt(save_path +
           'll_vmf_10d_500_epoch_v1.tsv',
           ll)

model_name = save_path + 'model_vmf_10d_500_epoch_v1'
model.save_sess(model_name)


