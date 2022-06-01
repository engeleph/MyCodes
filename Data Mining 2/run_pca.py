"""
Homework: Principal Component Analysis
Course  : Data Mining II (636-0019-00L)
"""

# Import all necessary functions
from utils import *
from pca import *
'''
Main Function
'''
if __name__ in "__main__":
    # Initialise plotting defaults
    initPlotLib()

    ####################################
    # Exercise 2:
    
    # Simulate Data
    data = simulateData()
    X = data['data']
    Y = data['target']
    # Perform a PCA
    # 1. Compute covariance matrix
    cov_matrix = computeCov(X)
    # 2. Compute PCA by computing eigen values and eigen vectors
    pca, eigen_val = computePCA(cov_matrix)
    # 3. Transform your input data onto a 2-dimensional subspace using the first two PCs
    transformed_matrix = transformData(pca,X)
    # 4. Plot your transformed data and highlight the four different sample classes
    plotTransformedData(transformed_matrix,Y,filename="exercise1.pdf")
    # 5. How much variance can be explained with each principle component?
    var = computeVarianceExplained(eigen_val)# Compute Variance Explained
    np.set_printoptions(precision=2)
    print("Variance Explained Exercise 2.1: ")
    # Uncomment the following 3 lines!
    for i in range(15):
        print("PC %d: %.2f"%(i+1,var[i]))
    # 6. Plot cumulative variance explained per PC
    plotCumSumVariance(var,filename="cumsum.pdf")
    
    ####################################
    # Exercise 2 Part 2:
    
    # 1. Normalise data
    X_std=dataNormalisation(X)
    # 2. Compute covariance matrix
    cov_matrix_std = computeCov(X_std)
    # 3. Compute PCA
    pca_std, eigen_val_std = computePCA(cov_matrix_std)
    # 4. Transform your input data inot a 2-dimensional subspace using the first two PCs
    transformed_matrix_std = transformData(pca_std,X_std)
    # 5. Plot your transformed data
    plotTransformedData(transformed_matrix_std,Y,filename="exercise2.pdf")
    # 6. Compute Variance Explained
    var = computeVarianceExplained(eigen_val_std) # Compute Variance Explained
    np.set_printoptions(precision=2)
    print("Variance Explained Exercise 2.2: ")
    # Uncomment the following 3 lines!
    for i in range(15):
        print("PC %d: %.2f"%(i+1,var[i]))
    # 7. Plot Cumulative Variance
    plotCumSumVariance(var,filename="cumsum_2.pdf")

    #Extra: Implement sklearn to compare results
    #from sklearn.decomposition import PCA as sklearnPCA
    #sklearn_pca = sklearnPCA(n_components=2)
    #Y_sklearn = sklearn_pca.fit_transform(X_std)
    #plot transformed data
    #plotTransformedData(Y_sklearn,Y,filename="control_plot.pdf")
    #print(sklearn_pca.explained_variance_ratio_)


