from shortest_path_kernel import floyd_warshall
from shortest_path_kernel import sp_kernel
import scipy.io
import numpy as np
import os
import sys
import argparse
import math



if __name__ == '__main__':

    # Set up the parsing of command-line arguments
    parser = argparse.ArgumentParser(
        description="Compute shortest path kernel"
    )
    parser.add_argument(
        "--datadir",
        required=True,
        help="Path to input directory containing file MUTAG.mat"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Path to directory where graphs_output.txt will be created"
    )

    args = parser.parse_args()

    # Set the paths
    data_dir = args.datadir
    out_dir = args.outdir

    os.makedirs(args.outdir, exist_ok=True)

    # Read the file
    mat = scipy.io.loadmat("{}/{}".format(args.datadir, 'MUTAG.mat'))
    label = np.reshape(mat['lmutag'], (len(mat['lmutag'], )))
    data = np.reshape(mat['MUTAG']['am'], (len(label), )

    # Create the output file
    try:
        file_name = "{}/graphs_output.txt".format(args.outdir)
        f_out = open(file_name, 'w')
    except IOError:
        print("Output file {} cannot be created".format(file_name))
        sys.exit(1)

    cdict = {}
    cdict['mutagenic'] = 1
    cdict['non-mutagenic'] = -1
    lst_group = ['mutagenic', 'non-mutagenic']

    # Write header for output file
    f_out.write('{}\t{}\t{}\n'.format(
        'Pair of classes',
        'SP'))


    # Iterate through all combinations of pairs
    for idx_g1 in range(len(lst_group)):
        for idx_g2 in range(idx_g1, len(lst_group)):
            # Get the group data
            #group1 = data[data[:, 0] == cdict[lst_group[idx_g1]]]
            group1 = data[np.where(label==cdict[lst_group[idx_g1]])]
            #group2 = data[data[:, 0] == cdict[lst_group[idx_g2]]]
            group2 = data[np.where(label==cdict[lst_group[idx_g2]])]

            # Get average similarity
            count = 0
            vec_sim = 0
            for x in group1:
                for y in group2:
                    # Skip redundant calculations
                    if np.any(x!=):
                        continue

                    # Compute shortest path matrix with floyd_marshall
                    S1=floyd_warshall(x)
                    S2=floyd_warshall(y)
                    # Compute shortes path kernel
                    vec_sim += sp_kernel(S1,S2)

                    count += 1
            vec_sim /= count

            # Transform the vector of distances to a string
            str_sim = '\t'.join('{0:.2f}'.format(x) for x in vec_sim)

            # Save the output
            f_out.write(
                '{}:{}\t{}\n'.format(
                    lst_group[idx_g1], lst_group[idx_g2], str_sim)
            )
    f_out.close()

