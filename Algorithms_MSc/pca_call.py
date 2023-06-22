# %%
from AggelosKozwnakis import conventionalPCA, PCAplots

# %% To generate the plots i used combination of the variables below. I changed variance, std and method each time.
from sklearn.datasets import make_blobs
X, y = make_blobs(n_samples=1000,n_features=20,centers=3,cluster_std=1.0,random_state=69)  

proj, load, var_ex, eigval, eigvec = conventionalPCA(data=X,variance=0.9, method = 'eig')

#PCAplots(proj,load, var_ex, plot='all')
#PCAplots(proj,load, var_ex, dims=3, pcs=[1,2,3], plot='all') #For 3D. pcs list can be the preferred PCs, given they are included in the dimensions kept.


# %% Barplot for loadings of PC1 + PC2 -- Just the sample code is provided since multiple plots were produced
import matplotlib.pyplot as plt
import numpy as np

#Generated with the help of a smart machine.
features = np.arange(len(load[0])) + 1
bar_width = 0.35
plt.bar(features, load[0], width=bar_width, label='PC1')
plt.bar(features + bar_width, load[1], width=bar_width, label='PC2')
plt.xlabel('Features')
plt.ylabel('Correlation Value')
plt.title('Loadings of PC1 and PC2')
plt.xticks(features + bar_width/2, features)
plt.legend()
plt.show()

# %% Heatmap for loadings of PC1 + PC2 -- Just the sample code is provided since multiple plots were produced
import matplotlib.pyplot as plt
import numpy as np

plt.imshow(load[:2, :], cmap= 'coolwarm', aspect = 'auto')
plt.colorbar()
plt.xlabel('Features')
plt.ylabel('Principal Components')
plt.title('Loadings Heatmap')
plt.yticks([0, 1], ['PC1', 'PC2'])
plt.xticks(np.arange(load.shape[1]), np.arange(1, load.shape[1]+1))
plt.show()


