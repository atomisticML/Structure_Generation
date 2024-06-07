import numpy as np

avdesc = False
ndescs = 11
A = np.zeros((0,ndescs))
print (A.shape)
print (A)
for i in range(1,11):
    fname = 'gsqs_%d_all_descriptors.npy' % i
    ai = np.load(fname)
    if avdesc:
        ai =np.average(ai,axis=0)
        ai=np.array([ai])
    print ('i %d,'%i,ai.shape)
    A =np.concatenate( (A,ai),axis=0)

print (A.shape)

#A = A[:,[iii for iii in range(1,np.shape(A)[1] -1)  ]]
A = A[:,[iii for iii in range(0,np.shape(A)[1] )  ]]
from sklearn import preprocessing
At = preprocessing.normalize(A)
#print(At)
print('amat shape',np.shape(At))
from sklearn.decomposition import PCA
pca = PCA(n_components=10,whiten=False)
pca.fit(At)
print ('explained variance',pca.explained_variance_ratio_)


ts = pca.transform(At)
