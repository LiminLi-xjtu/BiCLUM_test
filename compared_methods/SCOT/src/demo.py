import sys
import utils as ut
import evals as evals
import scot2 as sc
import numpy as np

X=np.genfromtxt("../data/scGEM_methylation.csv", delimiter=",")
y=np.genfromtxt("../data/scGEM_expression.csv", delimiter=",")
print("Dimensions of input datasets are: ", "X= ", X.shape, " y= ", y.shape)

# initialize SCOT object
scot=sc.SCOT(X, y)
# call the alignment with z-score normalization
X_new, y_new = scot.align( k=35, e=5e-3,  normalize=True, norm="l2")

fracs=evals.calc_domainAveraged_FOSCTTM(X_new, y_new)
print("Average FOSCTTM score for this alignment is: ", np.mean(fracs))

import matplotlib.pyplot as plt
legend_label="SCOT alignment FOSCTTM \n average value: "+str(np.mean(fracs))
plt.plot(np.arange(len(fracs)), np.sort(fracs), "r--", label=legend_label)
plt.legend()
plt.xlabel("Cells")
plt.ylabel("Sorted FOSCTTM")
plt.show()

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

pca=PCA(n_components=2)
Xy_pca=pca.fit_transform(np.concatenate((X_new, y_new), axis=0))
X_pca=Xy_pca[0: 177,]
y_pca=Xy_pca[177:,]

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(X_pca[:,0], X_pca[:,1], c="k", s=15, label="Gene Expression")
ax1.scatter(y_pca[:,0], y_pca[:,1], c="r", s=15, label="DNA Methylation")
ax1.legend()
ax1.set_title("Colored based on domains")

ax2.scatter(X_pca[:,0], X_pca[:,1], cmap="rainbow")# , c=Xlabels, s=15)
ax2.scatter(y_pca[:,0], y_pca[:,1], cmap="rainbow")# , c=ylabels, s=15)
ax2.set_title("Colored based on cell type")
plt.show()