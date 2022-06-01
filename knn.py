"""
Homework : k-Nearest Neighbor and Naive Bayes
Course   : Data Mining (636-0018-00L)

Main program for k-NN.
Predicts the labels of the test data using the training data.
The k-NN algorithm is executed for different values of k (user-entered parameter)


Original author: Damian Roqueiro <damian.roqueiro@bsse.ethz.ch>
Extended by: Bastian Rieck <bastian.rieck@bsse.ethz.ch>
"""

import argparse
import os
import sys
import numpy as np
import math

# Import the file with the performance metrics 
import evaluation as eval

# Class imports
from knn_classifier import KNNClassifier


# Constants
# 1. Files with the datapoints and class labels
DATA_FILE  = "matrix_mirna_input.txt"
PHENO_FILE = "phenotype.txt"

# 2. Classification performance metrics to compute
PERF_METRICS = ["accuracy", "precision", "recall"]


def load_data(dir_path): # (TO DO)
    """
    Function for loading the data.
    Receives the path to a directory that will contain the DATA_FILE and PHENO_FILE.
    Loads both files into memory as numpy arrays. Matches the patientId to make
    sure the class labels are correctly assigned.

    Returns
     X : a matrix with the data points
     y : a vector with the class labels
    """

    x=np.loadtxt("{}/{}".format(dir_path,DATA_FILE),delimiter='\t', skiprows=1, usecols=np.arange(1,488))
    y=np.loadtxt("{}/{}".format(dir_path,PHENO_FILE),delimiter='\t', skiprows=1, usecols=(1), dtype=list)
    for i in range(len(y)):
         if y[i]=='+':
                 y[i]=1
         else:
                 y[i]=0
    return x,y


def obtain_performance_metrics(y_true, y_pred): # (TO DO)
    """
    Function obtain_performance_metrics
    Receives two numpy arrays with the true and predicted labels.
    Computes all classification performance metrics.
    
    In this function you might call the functions:
    compute_accuracy(), compute_precision(), compute_recall()
    from the evaluation.py file. You can call them by writing:
    evaluation.compute_accuracy, and similarly.

    Returns a vector with one value per metric. The positions in the
    vector match the metric names in PERF_METRICS.
    """

    # CODE GOES HERE
    accuracy=compute_accuracy(y_true, y_pred)
    precision=compute_precision(y_true, y_pred)
    recall=compute_recall(y_true, y_pred)
    perf_metrics=[accuracy,precision,recall]
    return perf_metrics



#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------

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
        "--testdir",
        required=True,
        help="Path to input directory containing testing set"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Path to directory where output_knn.txt will be created"
    )
    parser.add_argument(
        "--mink",
        required=True,
        help="Minimal k value"
    )
    parser.add_argument(
        "--maxk",
        required=True,
        help="Maximal k value"
    )

    args = parser.parse_args()

    #maximal and minimal k
    min_k=int(args.mink)
    max_k=int(args.maxk)
    # Set the paths
    train_dir=args.traindir
    test_dir=args.testdir
    out_dir = args.outdir

    os.makedirs(args.outdir, exist_ok=True)

    # Read the file
    # If the output directory does not exist, then create it
    #
    try:
        file_name = "{}/output_knn.txt".format(args.outdir)
        f_out = open(file_name, 'w')
    except IOError:
        print("Output file {} cannot be created".format(file_name))
        sys.exit(1)

    
    # Read the training and test data. For each dataset, get also the true labels.
    # Use the function load_data().
    # Important: Match the patientId between data points and class labels
    #
    result_train= load_data("{}".format(train_dir))
    matrix_train = result_train[0]
    label_train = result_train[1]
    result_test=load_data("{}".format(test_dir))
    matrix_test = result_test[0]
    label_test = result_test[1]


    # Create the output file & write the header as specified in the homework sheet
    #
    f_out.write('{}\t{}\t{}\t{}\n'.format(
        'Value of k',
        'Accuracy',
        'Precision',
        'Recall'))

    ############################## KNN algorithm ####################################

    # Create the k-NN object. (Hint about how to do it in the homework sheet)

    # Iterate through all possible values of k:
    # HINT: remember to set the number of neighbors for the KNN object through: knn_obj.set_k(k)

    # 1. Perform KNN training and classify all the test points. In this step, you will
    # obtain a prediction for each test point. 

    y_pred=np.zeros(((max_k-min_k+1),len(label_test)))
    knn = KNNClassifier(matrix_train,label_train,metric='euclidean')
    for i in range((max_k-min_k)+1):
        k=min_k+i
        knn.set_k(k)
        for j in range(len(label_test)):
            y_pred[i,j]=knn.predict(matrix_test[j,])





    # 2. Compute performance metrics given the true-labels vector and the predicted-
    # labels vector (you might consider to use obtain_performance_metrics() function)

    #precision=len(y_pred)*[0]
    #recall=len(y_pred)*[0]
    #accuracy=len(y_pred)*[0]
    #k_values=np.arange(min_k,max_k)
    kvalues=min_k
    for i in range(len(y_pred)):
        precision=float("{:.2f}".format(eval.compute_precision(label_test, y_pred[i,])))
        recall=float("{:.2f}".format(eval.compute_recall(label_test, y_pred[i,])))
        accuracy=float("{:.2f}".format(eval.compute_accuracy(label_test, y_pred[i,])))

        # 3. Write performance results in the output file, as indicated the in homework
        # sheet. Close the file.
        f_out.write('{}\t{}\t{}\t{}\n'.format(kvalues, accuracy, precision, recall))
        kvalues=kvalues+1
    # (TO DO)

    f_out.close()
