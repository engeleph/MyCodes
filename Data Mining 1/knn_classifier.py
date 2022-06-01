import os
import numpy as np
import math


class KNNClassifier:
    '''
    A class object that implements the methods of a k-Nearest Neighbor classifier
    The class assumes there are only two labels, namely POS and NEG

    Attributes of the class
    -----------------------
    k : Number of neighbors
    X : A matrix containing the data points (train set)
    y : A vector with the labels
    dist : Distance metric used. Possible values are: 'euclidean', 'hamming', 'minkowski', and others
           For a full list of possible metrics have a look at:
           http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html

    HINT: for using the attributes of the class in the class' methods, you can use: self.attribute_name
          (e.g. self.X for accessing the value of the attribute named X)
    '''     
    def __init__(self, X, y, metric):
        '''
        Constructor when X and Y are given.
        
        Parameters
        ----------
        X : Matrix with data points
        Y : Vector with class labels
        metric : Name of the distance metric to use
        '''
        # Default values
        self.verbose = False
        self.k = 1

        # Parameters
        self.X = X
        self.y = y
        self.metric = metric


    def debug(self, switch):
        '''
        Method to set the debug mode.
        
        Parameters
        ----------
        switch : String with value 'on' or 'off'
        '''
        self.verbose = True if switch == "on" else False


    def set_k(self, k):
        '''
        Method to set the value of k.
        
        Parameters
        ----------
        k : Number of nearest neighbors
        '''
        self.k = k


    def _compute_distances(self, X, x):
        '''
        Private function to compute distances. 
        Compute the distance between x and all points in X
    
        Parameters
        ----------
        x : a vector (data point)
        '''

        # (TO DO)
        # CODE GOES HERE!
        distances=np.array(len(X)*[0.0])
        x=np.array(x)
        for i in range(len(X)):
            distances[i]=float(math.sqrt(np.sum((np.array(X[i,])-x)**2)))

        return distances


    def predict(self, x):
        '''
        Method to predict the label of one data point.
        Here you actually code the KNN algorithm.
       
        Hint: for calling the method _compute_distance 
              (which is private), you can use: 
              self._compute_distances(self.X, x) 
        
        Parameters
        ----------
        x : Vector from the test data.
        '''
        
        # (TO DO)
        if self.k%2==0:
            k=self.k-1
            dist=self._compute_distances(self.X, x)
            patient=k*[0]
            for i in range(k):
                index=np.where(dist==min(dist))
                patient[i]=index[0][0]
                dist=np.delete(dist, patient[i])
            pos=0
            neg=0
            for i in range(len(patient)):
                if self.y[patient[i]]==1:
                    pos=pos+1
                else:
                    neg=neg+1
            if pos>neg:
                y_predicted=1
            else:
                y_predicted=0
        else:
            dist=self._compute_distances(self.X, x)
            patient=self.k*[0]
            for i in range(self.k):
                index=np.where(dist==min(dist))
                patient[i]=index[0][0]
                dist=np.delete(dist, patient[i])
            pos=0
            neg=0
            for i in range(len(patient)):
                if self.y[patient[i]]==1:
                    pos=pos+1
                else:
                    neg=neg+1
            if pos>neg:
                y_predicted=1
            else:
                y_predicted=0
        
        return y_predicted
