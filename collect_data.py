import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing

def process_csv(namelist, batch_number):
    df_paths = pd.DataFrame(columns=["RESID", "LABEL1", "LABEL2", "CUTOFFS", "PATH"])

    for filename in tqdm(namelist):
        df = pd.read_csv(filename,
                         sep=" ", header=None,
                         names=["CUTOFF", "RESNAME", "RESID", "LABEL1", "LABEL2", "LEC"])
        resids = df.RESID.unique()

        for resid in resids:
            df_resid = df[df.RESID == resid]
            path = df_resid.LEC.values

            df_paths.loc[len(df_paths)] = [resid, df_resid.LABEL1.values[0], df_resid.LABEL2.values[0],
                                           df_resid.CUTOFF.values, path]

    # Save the result as a pickle file
    pickle_filename = f"all_paths_batch_{batch_number}.pkl"
    df_paths.to_pickle(pickle_filename)

import glob
def main():
    namelist = glob.glob('FIXED/*/*.label')
    print(namelist)
    # Calculate the number of batches
    num_batches = 14
    # Split the namelist into 14 equal-sized sublists
    namelist_batches = [namelist[i:i + len(namelist) // num_batches] for i in range(0, len(namelist), len(namelist) // num_batches)]
    # Process files in parallel using multiprocessing
    with multiprocessing.Pool(processes=num_batches) as pool:
        for i, namelist_batch in enumerate(namelist_batches):
            pool.apply_async(process_csv, (namelist_batch, i))
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
