import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

df = pd.read_csv('acetylcholinesterase_06_bioactivity_data_3class_pIC50_substructure_fp.csv')

# input features
X = df.drop('pIC50', axis=1)
X

# output features
Y = df.pIC50
Y

# examine the data dimension
X.shape
Y.shape

# remove low variance features
from sklearn.feature_selection import VarianceThreshold
selection = VarianceThreshold(threshold=(.8*(1-.8)))
X = selection.fit_transform(X)

X.shape

# data split (80/20 ratio)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

X_train.shape

Y_train.shape

# building a regression model using random forest
model = RandomForestRegressor(n_estimators=100)
model.fit(X_train, Y_train)
r2 = model.score(X_test, Y_test)
r2

Y_pred = model.predict(X_test)

# scatter plot of experimental vs predicted pIC50 values
import matplotlib.pyplot as plt

sns.set(color_codes=True)    
sns.set_style('white')

ax = sns.regplot(Y_test, Y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large',fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)

plt.savefig('randomforest_regression.pdf')