#!/usr/bin/env python3

'''
Skeleton for Homework 4: Logistic Regression and Decision Trees
Part 2: Decision Trees

Authors: Anja Gumpinger, Bastian Rieck
'''

import numpy as np
import sklearn.datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
import statistics


if __name__ == '__main__':

    iris = sklearn.datasets.load_iris()
    X = iris.data
    y = iris.target
    print(X.shape)

    feature_names = iris.feature_names
    num_features = len(set(feature_names))

    def split_data(X, y, attribute_index, theta):
        bigger=X[:,attribute_index]>=theta
        smaller=X[:,attribute_index]<theta
        X1=X[bigger,:]
        y1=y[bigger]
        X2=X[smaller,:]
        y2=y[smaller]
        list=[X1,y1,X2,y2]
        return list

    def compute_information_content(y):
        y0=0
        y1=0
        y2=0
        n=len(y)
        for i in range(len(y)):
            if y[i]==0:
                y0=y0+1
            elif y[i]==1:
                y1=y1+1
            elif y[i]==2:
                y2=y2+1
        info=(-1)*((y0/n)*np.log2(y0/n)+(y1/n)*np.log2(y1/n)+(y2/n)*np.log2(y2/n))
        return info

    def compute_information_a(X, y, attribute_index, theta):
        splitData=split_data(X, y, attribute_index, theta)
        X1=splitData[0]
        y1=splitData[1]
        X2=splitData[2]
        y2=splitData[3]
        info1=compute_information_content(y1)
        info2=compute_information_content(y2)
        info_A=(len(y1)/len(y))*info1+(len(y2)/len(y))*info2
        return info_A

    def compute_information_gain(X, y, attribute_index, theta):
        info=compute_information_content(y)
        info_A=compute_information_a(X, y, attribute_index, theta)
        gain=info-info_A
        return gain


    print('Exercise 2.b')
    print('------------')
    sepal_length=round(compute_information_gain(X, y, 0, 5),2)
    sepal_width=round(compute_information_gain(X, y, 1, 3),2)
    petal_length=round(compute_information_gain(X, y, 2, 2.5),2)
    petal_width=round(compute_information_gain(X, y, 3, 1.5),2)
    print('Split ( sepal length (cm) < 5.0): information gain = ', sepal_length)
    print('Split ( sepal width (cm) < 3.0): information gain = ', sepal_width)
    print('Split ( petal length (cm) < 2.5): information gain = ', petal_length)
    print('Split ( petal width (cm) < 1.5): information gain = ', petal_width)

    print('')

    print('Exercise 2.c')
    print('------------')
    print('I would select ( sepal width (cm) < 3.0 ) to be the first split,')
    print('because it results in the greatest gain of information')

    print('')

    ####################################################################
    # Exercise 2.d
    ####################################################################

    # Do _not_ remove this line because you will get different splits
    # which make your results different from the expected ones...
    np.random.seed(42)
    decTree = DecisionTreeClassifier()
    cv=ShuffleSplit(n_splits=5)
    crossVal=cross_val_score(decTree, X, y, cv=cv)
    mean_accuray=statistics.mean(crossVal)
    mean_accuray=round(mean_accuray*100,2)
    print('Accuracy score using cross-validation')
    print('-------------------------------------\n')
    print('The mean accuracy is',mean_accuray,'%')


    print('')
    print('Feature importances for _original_ data set')
    print('-------------------------------------------\n')
    feature_matrix=np.zeros((5,4))
    count=0
    for train_index, test_index in cv.split(X, y):
        decTree = DecisionTreeClassifier()
        X_train, y_train = X[train_index], y[train_index]
        X_test, y_test = X[test_index], y[test_index]
        decTree.fit(X,y)
        feature_importances=decTree.feature_importances_
        feature_matrix[count,:]=feature_importances
        count=count+1
    mean_sepal_length=round(statistics.mean(feature_matrix[:,0]),3)
    mean_sepal_width=round(statistics.mean(feature_matrix[:,1]),3)
    mean_petal_length=round(statistics.mean(feature_matrix[:,2]),3)
    mean_petal_width=round(statistics.mean(feature_matrix[:,3]),3)
    print('Mean Feature Importance of sepal length:', mean_sepal_length)
    print('Mean Feature Importance of sepal width:', mean_sepal_width)
    print('Mean Feature Importance of petal length:', mean_petal_length)
    print('Mean Feature Importance of petal width:', mean_petal_width)
    print(' ')
    print('For the original data, the two most important features are:')
    print('-petal length')
    print('-petal width')


    print('')
    print('Feature importances for _reduced_ data set')
    print('------------------------------------------\n')
    X = X[y != 2]
    y = y[y != 2]
    feature_matrix=np.zeros((5,4))
    count=0
    for train_index, test_index in cv.split(X, y):
        decTree = DecisionTreeClassifier()
        X_train, y_train = X[train_index], y[train_index]
        X_test, y_test = X[test_index], y[test_index]
        decTree.fit(X,y)
        feature_importances=decTree.feature_importances_
        feature_matrix[count,:]=feature_importances
        count=count+1
    mean_sepal_length=round(statistics.mean(feature_matrix[:,0]),3)
    mean_sepal_width=round(statistics.mean(feature_matrix[:,1]),3)
    mean_petal_length=round(statistics.mean(feature_matrix[:,2]),3)
    mean_petal_width=round(statistics.mean(feature_matrix[:,3]),3)
    print('Mean Feature Importance of sepal length:', mean_sepal_length)
    print('Mean Feature Importance of sepal width:', mean_sepal_width)
    print('Mean Feature Importance of petal length:', mean_petal_length)
    print('Mean Feature Importance of petal width:', mean_petal_width)
    print(' ')
    print('For the reduced data, the most important feature is:')
    print('-petal length')
    print(' ')
    print('This means that petal width is especially important to classify class 2')
