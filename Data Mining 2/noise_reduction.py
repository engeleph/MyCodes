import imageio
import scipy.linalg as linalg
import numpy as np


'''
Function for performing the noise reduction through SVD.
This is a simple implementation where the reduction of
the noise is performed via the so-called nullification 
of singular values (and corresponding singular vectors 
after the multiplication) below a fixed threshold.

Input
----------------
img:	    noisy image
threshold:  singular values threshold

Output
---------------- 
img_dn:		denoised image
'''
def NoiseReduction(img, threshold):
	# 1. Perform SVD on the noisy image, i.e., img
	U,S,V = linalg.svd(img)
	# 2. Nullify the singular values lower than the 
	#    threshold, i.e., threshold
	S[S<=threshold]=0
	Sig = np.zeros((img.shape[0], img.shape[1]))
	Sig[0:len(S), 0:len(S)] = np.diag(S)
	# 3. Reconstruct the image using the modified 
	#    singular values
	return np.dot(U,np.dot(Sig,V))

