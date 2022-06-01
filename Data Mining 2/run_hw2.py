#import all necessary functions
from utils import *
from pca import *
from scipy import linalg
from pinv import *
from noise_reduction import *
import imageio

'''
Main Function
'''
if __name__ in "__main__":
    # Initialise plotting defaults
    initPlotLib()

    ####################################
    # Exercise 1:

    # Get Iris Data
    data = loadIrisData()
    X=data['data']
    Y=data['target']
    
    #Perform a PCA using covariance matrix and eigen-value decomposition
    # 1. Compute covariance matrix
    cov_matrix = computeCov(X)
    # 2. Compute PCA by computing eigen values and eigen vectors
    pca = computePCA(cov_matrix)
    # 3. Transform your input data onto a 2-dimensional subspace using the first two PCs
    transformed_mat = transformData(pca[1],X)[:,[0,1]]
    # 4. Plot your transformed data and highlight the three different sample classes
    plotTransformedData(transformed_mat,Y,filename="PCA.png")
    # 5. How much variance can be explained with each principle component?
    var = computeVarianceExplained(pca[0])
    plotCumSumVariance(var,filename="Variance.png")
    print("Variance Explained PCA: ")
    # Uncomment the following 2 lines:
    for i in range(var.shape[0]):
         print("PC %d: %.2f"%(i+1,var[i]))


    # Perform a PCA using SVD
    # 1. Normalise data by substracting the mean
    X_zm = zeroMean(X)
    # 2. Compute PCA by computing eigen values and eigen vectors
    U,S,V=linalg.svd(X_zm)
    # 3. Transform your input data onto a 2-dimensional subspace using the first two PCs
    transformed_mat=U[:, :2]*S[:2]
    # 4. Plot your transformed data and highlight the three different sample classes
    plotTransformedData(transformed_mat,Y,filename="PCA_normalized.png")
    # 5. How much variance can be explained with each principle component?
    var = computeVarianceExplained(S)
    plotCumSumVariance(var,filename="Variance_svd.png")
    print("Variance Explained SVD: ")
    # Uncomment the following 2 lines:
    for i in range(var.shape[0]):
        print("PC %d: %.2f"%(i+1,var[i]))

    ####################################
    # Exercise 2:
    # 1. Compute the Moore-Penrose Pseudo-Inverse on the Iris data
    X_plus = compute_pinv(X)
    X_plus2 = linalg.pinv(X)

    # 2. Check Properties

    print("\nChecking status exercise 3:")
    status = False
    if np.dot(X,np.dot(X_plus,X)).all()==X.all():
        status = True
    print(f"X X^+ X = X is {status}")
    status = False
    if np.dot(X_plus,np.dot(X,X_plus)).all()==X_plus.all():
        status = True
    print(f"X^+ X X^+ = X^+ is {status}")



    # Exercise 3
    ####################################
    ##treshold=1000
    # 1. Loading the images
    img_sp = imageio.imread('images/greece_s&p.jpg', as_gray = True)
    
    # 2.  Perform noise reduction via SVD, i.e. call the function NoiseReduction
    img = NoiseReduction(img_sp, 1000)

    # 3. Save the denoised images using the command: 
    imageio.imwrite('images/noise_reduced_1000.jpg', img)

    ##treshold=500
    # 1. Loading the images
    img_sp = imageio.imread('images/greece_s&p.jpg', as_gray = True)

    # 2.  Perform noise reduction via SVD, i.e. call the function NoiseReduction
    img = NoiseReduction(img_sp, 500)

    # 3. Save the denoised images using the command:
    imageio.imwrite('images/noise_reduced_500.jpg', img)

    #treshold=300
    # 1. Loading the images
    img_sp = imageio.imread('images/greece_s&p.jpg', as_gray = True)

    # 2.  Perform noise reduction via SVD, i.e. call the function NoiseReduction
    img = NoiseReduction(img_sp, 3000)

    # 3. Save the denoised images using the command:
    imageio.imwrite('images/noise_reduced_3000.jpg', img)
