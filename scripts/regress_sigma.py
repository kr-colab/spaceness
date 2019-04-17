#extra trees regression on mate selection radius sd
import pandas, numpy as np, sklearn, os
from sklearn import model_selection
from sklearn import ensemble
import pickle
import matplotlib.pyplot as plt
#os.chdir("/projects/kernlab/cbattey2/pytheas/nwf/")
os.chdir("/Users/cj/spaceness")

############################ train extra trees classifier ###########################
#load summary stats
ss=np.loadtxt("sumstats/ss_msp_grid.txt",skiprows=1)
ss=pandas.DataFrame(ss)
ss=ss.drop(6,1)
ss=np.array(ss)
statnames=open("sumstats/ss_msp_grid.txt","r").readline().split(" ")

rand_order=np.random.choice(np.arange(0,len(ss),1),len(ss),replace=False)
ss=ss[rand_order,:]
train_stats=ss[0:int(len(ss)*.9),1:]
test_stats=ss[int(len(ss)*.9)+1:len(ss),1:]
train_lab=ss[0:int(len(ss)*.9),0]
test_lab=ss[int(len(ss)*.9)+1:len(ss),0]

#train extra trees regressor
param_grid_forest = {"max_depth": [3, 5, 10,50],
                     "min_samples_split": [2, 3, 10, 50],
                     "min_samples_leaf": [1, 3, 10, 50]}
clf = sklearn.ensemble.ExtraTreesRegressor(n_estimators=100,bootstrap=False)
clf_grid = sklearn.model_selection.GridSearchCV(clf,
                                            param_grid=param_grid_forest,
                                            cv=3,n_jobs=10)
clf_grid.fit(train_stats, train_lab)
clf_grid.best_params_
pytheas=clf_grid.best_estimator_

r2=pytheas.score(test_stats,test_lab)
pred=pytheas.predict(test_stats)

fig = plt.figure(figsize=(2,4),dpi=200)
plt.rcParams.update({'font.size': 7})
ax1=fig.add_axes([0,0.55,1,.45])
ax1.plot(test_lab,test_lab,"-",color="k")
ax1.plot(test_lab,pred,"o",color="green")
ax1.set_ylabel("Predicted Sigma")
ax1.set_xlabel("Simulated Sigma")
ax1.annotate(s="$R^2=$"+str(np.round(r2,3)),xy=(0.02,0.09))

ax2=fig.add_axes([0,0,1,.45])
ax2.barh(y=statnames[1:19],width=pytheas.feature_importances_[0:18],color="green")
ax2.set_xlabel("Feature Importance")
fig.savefig("/Users/cj/spaceness/figures/regression_summary.pdf",bbox_inches='tight')
