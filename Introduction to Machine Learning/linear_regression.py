import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

data=pd.read_csv('./train.csv')
Y=data['y']
Y=Y.to_numpy()
X=data.drop('y', axis=1)
X=X.to_numpy()
X=np.delete(X, 0, 1)
X2=X**2
X3=np.exp(X)
X4=np.cos(X)
X_comb=np.concatenate((X,X2,X3,X4),axis=1)
X_comb=np.c_[X_comb, np.ones(len(X))]
print(X_comb.shape)

reg = LinearRegression(fit_intercept=False).fit(X_comb, Y)
coef=pd.DataFrame(reg.coef_)
#coef.to_csv('values.csv', index=False, header=False)
