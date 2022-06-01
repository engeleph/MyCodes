'''
Skeleton for Homework 4: Logistic Regression and Decision Trees
Part 1: Logistic Regression

Authors: Anja Gumpinger, Dean Bodenham, Bastian Rieck
'''

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import math

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import StandardScaler


def compute_metrics(y_true, y_pred):
    '''
    Computes several quality metrics of the predicted labels and prints
    them to `stdout`.

    :param y_true: true class labels
    :param y_pred: predicted class labels
    '''

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    print('Exercise 1.a')
    print('------------')
    print('TP: {0:d}'.format(tp))
    print('FP: {0:d}'.format(fp))
    print('TN: {0:d}'.format(tn))
    print('FN: {0:d}'.format(fn))
    print('Accuracy: {0:.3f}'.format(accuracy_score(y_true, y_pred)))


if __name__ == "__main__":

    #load training data
    train_file = "data/diabetes_train.csv"
    #read data from file using pandas
    df = pd.read_csv(train_file)
    # extract first 7 columns to data matrix X (actually, a numpy array)
    X_train = df.iloc[:, 0:7].values
    #scale data in X
    transform=StandardScaler()
    X_train=transform.fit_transform(X_train)
    # extract 8th column (labels) to numpy array
    y_train = df.iloc[:, 7].values

    #load test data
    test_file = "data/diabetes_test.csv"
    #read data from file using pandas
    df = pd.read_csv(test_file)
    # extract first 7 columns to data matrix X (actually, a numpy ndarray)
    X_test = df.iloc[:, 0:7].values
    #scale data in X
    X_test=transform.fit_transform(X_test)
    # extract 8th column (labels) to numpy array
    y_test = df.iloc[:, 7].values

    #perform logaristic regression with trainings set
    logReg = LogisticRegression()
    log_reg = logReg.fit(X_train, y_train)

    #predict y for test set
    y_predict=logReg.predict_proba(X_test)
    y_predict=y_predict[:,1]

    #convert probabilities to 0 and 1
    y_predict[y_predict>=0.5]=1
    y_predict[y_predict<0.5]=0

    compute_metrics(y_test, y_predict)
    print(' ')

    print('Exercise 1.b')
    print('------------')
    print('For this diabetes dataset I would choose Logistic Regression because the accuracy is a bit higher than in LDA')
    print('More important might be if we aim for a low false negative rate (LDA is better choice) or low false positive rate (LR is better choice)')
    print(' ')
    print('Exercise 1.c')
    print('------------')
    print('For another dataset I would again test both methods. Then I would choose again the method with the better accuracy')
    print('As metioned in 1 b) false positive or false negative rate is more important to focus on')
    print(' ')

    print('Exercise 1.d')
    print('------------')
    infl_coef=logReg.coef_[0]
    age=round(infl_coef[6],2)
    odds_age=round(math.exp(age),2)
    print('-The biggest influence on the evaluation have glucose concentration and diabetes pedigree function.')
    print(' ')
    print('-The coefficient for age is',age, '. Calculating the exponential function results in ',odds_age,',')
    print(' ')
    print(' which amounts to an increase in diabetes risk of', odds_age*100-100,' percent per additional year.')
    print(' ')

    train_file = "data/diabetes_train.csv"
    #read data from file using pandas
    df = pd.read_csv(train_file)
    # extract first 7 columns to data matrix X (actually, a numpy array)
    X_train = df.iloc[:, [0,1,2,4,5,6]].values
    #scale data in X
    transform=StandardScaler()
    X_train=transform.fit_transform(X_train)
    # extract 8th column (labels) to numpy array
    y_train = df.iloc[:, 7].values

    #load test data
    test_file = "data/diabetes_test.csv"
    #read data from file using pandas
    df = pd.read_csv(test_file)
    # extract first 7 columns to data matrix X (actually, a numpy ndarray)
    X_test = df.iloc[:, [0,1,2,4,5,6]].values
    #scale data in X
    X_test=transform.fit_transform(X_test)
    # extract 8th column (labels) to numpy array
    y_test = df.iloc[:, 7].values

    #perform logaristic regression with trainings set
    logReg = LogisticRegression()
    log_reg = logReg.fit(X_train, y_train)

    #predict y for test set
    y_predict=logReg.predict_proba(X_test)
    y_predict=y_predict[:,1]

    #convert probabilities to 0 and 1
    y_predict[y_predict>=0.5]=1
    y_predict[y_predict<0.5]=0
    skin=round(infl_coef[3],4)

    tn, fp, fn, tp = confusion_matrix(y_test, y_predict).ravel()
    print('TP: {0:d}'.format(tp))
    print('FP: {0:d}'.format(fp))
    print('TN: {0:d}'.format(tn))
    print('FN: {0:d}'.format(fn))
    print('Accuracy: {0:.3f}'.format(accuracy_score(y_test, y_predict)))
    print(' ')

    print('By comparing the performance and the coefficients obtained on the reduced')
    print('dataset with the ones on the model including all the attributes, I observe')
    print('that the evaluation does not change')
    print('My explanation is that the coeffiecient of skin, which is ', skin,' has a')
    print('tiny influence on the evaluation')


