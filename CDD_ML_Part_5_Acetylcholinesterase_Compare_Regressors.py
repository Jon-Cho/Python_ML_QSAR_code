import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
import lazypredict
from lazypredict.Supervised import LazyRegressor

# load data
df = pd.read_csv('acetylcholinesterase_06_bioactivity_data_3class_pIC50_substructure_fp.csv')

X = df.drop('pIC50', axis=1)
Y= df.pIC50

# data pre-processing
X.shape
Y.shape

# remove low variance features
from sklearn.feature_selection import VarianceThreshold
selection = VarianceThreshold(threshold=(.8*(1-.8)))
X=selection.fit_transform(X)
X.shape

# Perform data splitting using 80/20 ratio
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# compare ML algorithms
# defines and builds the lazyclassifier
clf = LazyRegressor(verbose=0, ignore_warnings='True', custom_metric=None)
models_train,predictions_train = clf.fit(X_train, X_train, Y_train, Y_train)
models_test,predictions_test = clf.fit(X_train, X_test, Y_train, Y_test)

# performance table of the training set (80% subset)
predictions_train

# performance table of the test set (20% subset)
predictions_test

# data visulization of model performance
# barplot of R-square values
import matplotlib.pyplot as plt
import seaborn as sns

# train['R-squared'] = [0 if i<0, else i for i in train.iloc[:,0]]
plt.figure(figsize=(5,10))
sns.set_theme(style='whitegrid')
ax = sns.barplot(y=predictions_train.index, x='R-Squared', data=predictions_train)
ax.set(xlim=(0,0.45))
plt.tight_layout()
plt.savefig('predictions_train(model performance).pdf')

# Bar plot of RMSE values
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=predictions_train.index, x="RMSE", data=predictions_train)
ax.set(xlim=(0, 6.5))
plt.tight_layout()
plt.savefig('predictions_train(RMSE_model performance).pdf')

# Bar plot of calculation time
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=predictions_train.index, x="Time Taken", data=predictions_train)
ax.set(xlim=(0, 4))
plt.tight_layout()
plt.savefig('predictions_train(Time_model performance).pdf')
