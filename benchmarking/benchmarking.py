import pandas as pd
import matplotlib.pyplot as plt
import os
import fnmatch

mnt = "/home/eleonora.renz/hpc4ds-project/benchmarking/"
experiments = ["light/", "heavy/", "cpu-increase"]

for experiment in experiments:

  df = []

  for filename in os.listdir(mnt+experiment):
      f = os.path.join(mnt+experiment, filename)
      if fnmatch.fnmatch(filename, '*sh.o*'):
        df_part = pd.read_csv(f, sep=',', names=['nnodes', 'ncpus', 'nproc', 'bytes', 'time'])
        df.append(df_part)

  df = pd.concat(df, sort=False)

  grouped_df = df.groupby(['nnodes', 'ncpus', 'nproc'])
  df = grouped_df[['time']].mean()
  df.reset_index(inplace = True)

  # plot
  xlabels = ['1', '2', '4', '8', '16', '32', '64', '128']
  plt.title(experiment.replace('/', ''))
  suptitle = "nodes:" + str(df.nnodes.unique()) + ", cpus:" + str(df.ncpus.unique())
  plt.suptitle(suptitle, y=0.05)
  plt.plot(df.nproc, df.time)
  plt.xticks(df.nproc, xlabels)
  plt.savefig(mnt+'plots/'+experiment.replace('/', ''))
  plt.show()