import re
import pandas as pd
from collections import defaultdict

project_tb = pd.read_csv("./tables/barseq_fastq_path.csv")

dfs = pd.DataFrame()

for exp in project_tb['experiment']:
  bt_csv_path="./bartender/"+exp+".clustered_cluster.csv"
  bc_freq_tb=pd.read_csv(bt_csv_path)[['Center', 'time_point_1']]
  bc_freq_tb.columns = ['barcode', exp]
  dfs=pd.merge(dfs, bc_freq_tb.set_index('barcode'), how='outer', left_index=True, right_index=True)

dfs.to_csv("tables/bartender_matrix.csv")
