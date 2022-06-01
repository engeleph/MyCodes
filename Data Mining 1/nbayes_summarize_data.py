import argparse
import os
import sys
import numpy as np
import math

def bayes(x,value):
    index=np.where(x==0)
    index=index[0][0]
    x=np.delete(x, index)
    count=0
    for i in range(len(x)):
        if x[i]==value:
            count=count+1
    prob=count/len(x)
    return prob

if __name__ == '__main__':

    # Set up the parsing of command-line arguments
    #
    parser = argparse.ArgumentParser(
        description="Compute distance functions on time-series"
    )
    parser.add_argument(
        "--traindir",
        required=True,
        help="Path to input directory containing training sest"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Path to directory where output_knn.txt will be created"
    )

    args = parser.parse_args()

    # Set the paths
    train_dir=args.traindir
    out_dir = args.outdir
    c_val=[0,2]

    os.makedirs(args.outdir, exist_ok=True)
    for k in c_val:
        # Read the file
        # If the output directory does not exist, then create it
        #
        try:
            if k==2:
                file_name = "{}/output_summary_class_2.txt".format(out_dir)
                f_out = open(file_name, 'w')
            else:
                file_name = "{}/output_summary_class_4.txt".format(out_dir)
                f_out = open(file_name, 'w')
        except IOError:
            print("Output file {} cannot be created".format(file_name))
            sys.exit(1)


        # Read the training and test data. For each dataset, get also the true labels.
        # Use the function load_data().
        # Important: Match the patientId between data points and class labels
        #
        matrix_train = np.loadtxt("{}/{}".format(train_dir,'tumor_info.txt'),delimiter='\t', converters = {i: lambda s: float(s.strip() or 0) for i in range(5)})
        #print(matrix_train)
        if k==2:
            rm=4
        else:
            rm=2
        class_vec=matrix_train[:,4]
        index=np.where(class_vec==rm)
        #print(index)
        index=index[0]
        print(index)
        matrix_train = np.delete(matrix_train, index, axis=0)
        #print(matrix_train)


        # Create the output file & write the header as specified in the homework sheet
        #
        f_out.write('{}\t{}\t{}\t{}\t{}\n'.format('Value', 'clump', 'uniformity', 'marginal', 'mitosis'))

        ############################## KNN algorithm ####################################

        # Create the k-NN object. (Hint about how to do it in the homework sheet)

        # Iterate through all possible values of k:
        # HINT: remember to set the number of neighbors for the KNN object through: knn_obj.set_k(k)

        # 1. Perform KNN training and classify all the test points. In this step, you will
        # obtain a prediction for each test point.
        #print(matrix_train[0:][0])
        #print(matrix_train[:,0][0])
        for j in range(10):
            value=j+1
            clump=float("{:.3f}".format(bayes(matrix_train[:,0],j+1)))
            uniformity=float("{:.3f}".format(bayes(matrix_train[:,1],j+1)))
            marginal=float("{:.3f}".format(bayes(matrix_train[:,2],j+1)))
            mitosis=float("{:.3f}".format(bayes(matrix_train[:,3],j+1)))

                # 3. Write performance results in the output file, as indicated the in homework
                # sheet. Close the file.
            f_out.write('{}\t{}\t{}\t{}\t{}\n'.format(value, clump, uniformity, marginal, mitosis))
        f_out.close()
