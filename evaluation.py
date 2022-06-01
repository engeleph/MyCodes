"""
Homework : k-Nearest Neighbor and Naive Bayes
Course   : Data Mining (636-0018-00L)

Auxiliary functions.

This file implements the metrics that are invoked from the main program.

Author: Damian Roqueiro <damian.roqueiro@bsse.ethz.ch>
Extended by: Bastian Rieck <bastian.rieck@bsse.ethz.ch>
"""

import numpy as np


def confusion_matrix(y_true, y_pred):
    '''
    Function for calculating TP, FP, TN, and FN.
    The input includes the vector of true labels
    and the vector of predicted labels
    '''
    # (TO DO)
    pos=(y_true==1)
    neg=(y_true==0)
    tp=sum(y_pred[pos]==1)
    fp=sum(y_pred[neg]==1)
    tn=sum(y_pred[neg]==0)
    fn=sum(y_pred[pos]==0)
    confusion_matrix = [[tp, fp], [fn, tn]]

    return confusion_matrix


def compute_precision(y_true, y_pred):
    """
    Function: compute_precision
    Invoke confusion_matrix() to obtain the counts
    """
    # (TO DO)
    pos=(y_true==1)
    neg=(y_true==0)
    tp=sum(y_pred[pos]==1)
    fp=sum(y_pred[neg]==1)
    precision=tp/(tp+fp)
    return precision


def compute_recall(y_true, y_pred):
    """
    Function: compute_recall
    Invoke confusion_matrix() to obtain the counts
    """
    # (TO DO)
    pos=(y_true==1)
    tp=sum(y_pred[pos]==1)
    fn=sum(y_pred[pos]==0)
    recall=tp/(tp+fn)
    return recall


def compute_accuracy(y_true, y_pred):
    """
    Function: compute_accuracy
    Invoke the confusion_matrix() to obtain the counts
    """
    # (TO DO)
    pos=(y_true==1)
    neg=(y_true==0)
    tp=sum(y_pred[pos]==1)
    fp=sum(y_pred[neg]==1)
    tn=sum(y_pred[neg]==0)
    fn=sum(y_pred[pos]==0)
    accuracy=(tp+tn)/(tp+tn+fp+fn)
    return accuracy

