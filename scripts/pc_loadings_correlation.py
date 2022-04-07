from scipy.stats import spearmanr
import pandas as pd
import matplotlib.pyplot as plt

dis=pd.read_table('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings.txt',sep=' ', header=0)
val=pd.read_table('analysis/pca/pca_glass_gbm_idhwt_rt_tmz_local_log2fc_glassFALSE_PC1_loadings.txt',sep=' ', header=0)

data=pd.merge(dis,val,how='inner', left_index=True, right_index=True)

rho, p = spearmanr(data['x_x'], data['x_y'])

print(rho)
print(p)

plt.figure()
plt.scatter(data['x_x'], data['x_y'], c='gray',edgecolors='none',s=8, alpha=0.5)
plt.xlabel('Discovery PC1 loadings')
plt.ylabel('Validation PC1 loadings')
plt.text(-0.1,0.04,"Spearman's rank pval="+str(round(p,sigfigs = 2)))
plt.show()

plt.savefig('analysis/pca/pc1_correlation.pdf')
