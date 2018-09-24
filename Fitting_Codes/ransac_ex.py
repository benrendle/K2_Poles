''' RANdom SAmple Consensus least squares algorithm example script adpated from
    http://scikit-learn.org/stable/auto_examples/linear_model/plot_ransac.html
    23/08/2018
    '''

import numpy as np
from matplotlib import pyplot as plt
from sklearn import linear_model, datasets
import pandas as pd
import sys

n_samples = 1000
n_outliers = 50

''' Generation of artificial regression set and outliers '''
# X, y, coef = datasets.make_regression(n_samples=n_samples, n_features=1, n_informative=1, noise=10, coef=True, random_state=0)

# Add outlier data
# np.random.seed(0)
# X[:n_outliers] = 3 + 0.5 * np.random.normal(size=(n_outliers, 1))
# y[:n_outliers] = -3 + 10 * np.random.normal(size=n_outliers)

''' Data in '''
df = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Normal/Aug_2018/AS_03082018')
df['Av_sig'] = (df['Av_68U']-df['Av_68L'])/2.
df['L'] = df['rad']**2 * (df['Teff']/5777)**4
X = (df['L'])
X = X.values.reshape(-1,1)
y = (df['Av'])

''' Fit line using all data '''
lr = linear_model.LinearRegression()
lr.fit(X, y)

''' Robustly fit linear model with RANSAC algorithm '''
a = np.zeros((100,2))
b = np.zeros((100,422)) # Size of inlier mask
for i in range(100):
    ransac = linear_model.RANSACRegressor(min_samples=5)
    ransac.fit(X, y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)

    ''' Predict data of estimated models '''
    line_X = np.arange(X.min(), X.max(), 0.01)#[:, np.newaxis]
    line_X = line_X.reshape(-1,1)
    line_y = lr.predict(line_X)
    line_y_ransac = ransac.predict(line_X)
    ''' Compare estimated coefficients '''
    # print("Estimated coefficients (linear regression, RANSAC):")
    # print(lr.coef_, ransac.estimator_.coef_)
    # print(ransac.estimator_.coef_,ransac.estimator_.intercept_)
    # print(dir(ransac))
    a[i,0], a[i,1] = ransac.estimator_.coef_,ransac.estimator_.intercept_
    b[i,:] = inlier_mask

    lw = 2
    # plt.errorbar(X[inlier_mask], y[inlier_mask], yerr=df['Av_sig'][inlier_mask], color='yellowgreen', marker='.', label='Inliers',fmt='o')
    # plt.errorbar(X[outlier_mask], y[outlier_mask], yerr=df['Av_sig'][outlier_mask], color='gold', marker='.', label='Outliers',fmt='o')
    # plt.plot(line_X, line_y, color='navy', linewidth=lw, label='Linear regressor')
    plt.plot(line_X, line_y_ransac, color='cornflowerblue', linewidth=lw, label='_nolegend_', alpha=0.1)
    # plt.legend(loc='lower right')
    plt.xlabel("Luminosity")
    plt.ylabel("Extinction")

b = b.astype('bool')
c = np.median(a[:,0])
medIdx = (np.abs(a[:,0] - c)).argmin()
# print(np.median(a[:,0]),a[medIdx,0],b[medIdx])

plt.errorbar(X, y, yerr=df['Av_sig'], color='gold', marker='.', label='Outliers',fmt='o')
plt.errorbar(X[b[medIdx]], y[b[medIdx]], yerr=df['Av_sig'][b[medIdx]], color='yellowgreen', marker='.', label='Inliers',fmt='o')
plt.plot(line_X, a[medIdx,0]*line_X + a[medIdx,1], color='k', linewidth=3, label='_nolegend_', alpha=0.75)

plt.show()
