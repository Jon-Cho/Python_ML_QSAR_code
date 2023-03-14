import pandas as pd

dataset = pd.read_csv('delaney_solubility_with_descriptors.csv')
dataset

# data set
X = dataset.drop(['logS'], axis=1)
Y = dataset.iloc[:,-1]

Y.sample(frac=1, replace=False, random_state=0)

# model from original X-Y pairs
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# data set
X = dataset.drop(['logS'], axis=1)
Y = dataset.iloc[:,-1]
# data split
X_train, X_test, Y_train, Y_test = train_test_split(X,Y, test_size=0.2, random_state=0)

# model building
model = LinearRegression()
model.fit(X_train, Y_train)
# make prediction
Y_train_pre = model.predict(X_train)
# compute R2 score
r2a = r2_score(Y_train, Y_train_pre)
r2a

# plot of predicted vs actual
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(5,5))
plt.scatter(x=Y_train, y=Y_train_pre, c='#7CAE00', alpha=0.3)

# Add trendline
# https://stackoverflow.com/questions/26447191/how-to-add-trendline-in-python-matplotlib-dot-scatter-graphs
z = np.polyfit(Y_train, Y_train_pre, 1)
p = np.poly1d(z)

plt.plot(Y,p(Y), '#F8766D')
plt.ylabel('Predicted LogS')
plt.xlabel('Experimental LogS')

plt.savefig('Y-scrambled model_no shuffle.pdf')

# building Y-scrambled models
# model building
r2b_train_list = []
for i in range(1000):
    # Y-shuffling
    Y_train = Y_train.sample(frac=1, replace=False, random_state=i)
    #model building
    model = LinearRegression()
    model.fit(X_train, Y_train)
    # make prediction
    Y_train_pre = model.predict(X_train)
    # compute R2 score
    r2b = r2_score(Y_train, Y_train_pre)
    r2b_train_list.append(r2b)

# Histogram plot
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
plt.figure(figsize = (20, 2)) 

ax = sns.histplot(data=r2b_train_list, color="skyblue", kde=True, bins = 10)
plt.axvline(r2a)
ax.set(xlim=(0, 1))

plt.savefig('Y-scrambled_vs_before.pdf')
