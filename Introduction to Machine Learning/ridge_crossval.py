from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
import numpy as np
import pandas as pd

np.random.seed(0)

data=pd.read_csv('./train.csv')
data=data.sample(frac=1)
data=data.to_numpy()
Y=data[:,0]
X=data[:,1:14]

def rmse(Y_pred, Y):
    return np.sqrt(((Y.reshape(-1) - Y_pred.reshape(-1)) ** 2).mean())

values=[]
lambdas=[0.1,1,10,100,200]
for lam in lambdas:
    rm=[]
    kf = KFold(n_splits=10)
    for train, test in kf.split(X):

        X_train=X[train,:]
        Y_train=Y[train]
        X_test=X[test,:]
        Y_test=Y[test]

        ridge_reg=Ridge(alpha=lam, fit_intercept=True, solver="auto")
        ridge_reg.fit(X_train,Y_train)

        Y_pred=pd.DataFrame(ridge_reg.predict(X_test))
        Y_pred=Y_pred.to_numpy()
        rm.append(rmse(Y_pred,Y_test))

    rm=sum(rm)/len(rm)
    values.append(rm)
values=pd.DataFrame(values)
values.to_csv('values.csv', index=False, header=False)

