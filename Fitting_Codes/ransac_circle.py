''' RANdom SAmple Consensus least squares algorithm example script for solving
    circular problems with RANSAC.
    19/09/2018
    '''

from __future__ import division
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from matplotlib import pyplot as plt
from sklearn import linear_model, datasets
from skimage.measure import CircleModel, ransac
import pandas as pd
import sys
import random

class circle():

    def y_intercept(self, p, m):
        C = p[1] - m*p[0]
        return C

    def intersection(self, m1, c1, m2, c2):
        if m1 == m2:
            return None, None
        x = (c2 - c1)/(m1 - m2)
        y = m1*x + c1
        return x,y

    def estimate(self,points):
        ''' Determine gradients of perpendicular lines from random points and
            calculate the expected circle radius from this. '''
        p1 = random.choice(points)
        p2 = random.choice(points)
        p3 = random.choice(points)
        p4 = [(p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.] # mid-point of p1p2
        p5 = [(p1[0]+p3[0])/2.,(p1[1]+p3[1])/2.] # mid-point of p1p3

        m1 = (p1[1]-p2[1])/(p1[0]-p2[0])
        m2 = (p1[1]-p3[1])/(p1[0]-p3[0])
        m3 = -1./m1 # grad. with p4
        m4 = -1./m2 # grad. with p5

        c1 = self.y_intercept(p4,m3)
        c2 = self.y_intercept(p5,m4)

        x,y = self.intersection(m3,c1,m4,c2)
        if x == None:
            x = 0.
            y = 0.
        r = np.sqrt((p1[0]-x)**2 + (p1[1]-y)**2)
        self.params = (x,y,r)

        return True

    def residuals(self, points):
        xc, yc, r = self.params # central coords + radius
        x = points[:,0]
        y = points[:,1]

        return r - np.sqrt((x-xc)**2 + (y-yc)**2)


    def predict_xy(self, t, params=None):
        """Taken from source code:
        Predict x- and y-coordinates using the estimated model.

        Parameters
        ----------
        t : array
            Angles in circle in radians. Angles start to count from positive
            x-axis to positive y-axis in a right-handed system.
        params : (3, ) array, optional
            Optional custom parameter set.

        Returns
        -------
        xy : (..., 2) array
            Predicted x- and y-coordinates.

        """
        if params is None:
            params = self.params
        xc, yc, r = params
        x = xc + r * np.cos(t)
        y = yc + r * np.sin(t)

        return np.concatenate((x[..., None], y[..., None]), axis=t.ndim)


nsamples = 1000
n_outliers = 150

''' Generation of artificial regression set and outliers '''
X = np.linspace(0,1.,nsamples)
np.random.seed(0)
X[:n_outliers] = 1 * np.random.normal(size=(n_outliers))
X1 = np.linspace(0,-1,nsamples)
X1[:n_outliers] = 1 * np.random.normal(size=(n_outliers))
X2 = np.concatenate((X,X,X1,X1),axis=0)
X2 = X2+5
r = 1
y1 = np.sqrt(r - X**2)
y2 = -np.sqrt(r - X**2)
Y = np.concatenate((y1,y2,y1,y2),axis=0)

''' Data to single array '''
points = np.column_stack([X2,Y])

''' Fitting of 'circle' model with RANSAC '''
model_robust, inliers = ransac(points, circle, min_samples=3, residual_threshold=0.05, max_trials=1000)

''' Extract circle centroid coords and radius '''
fit = (model_robust.params)
print(fit)

''' Generate best fit data '''
t = np.linspace(0,2*np.pi,nsamples)
predict = model_robust.predict_xy(t,fit)

''' Plotting '''
fig, ax = plt.subplots()
ax.scatter(X2,Y,alpha=0.2)
ax.scatter(points[inliers,0],points[inliers,1],color='m')
plt.scatter(predict[:,0],predict[:,1])
plt.show()
