#!/usr/bin/env python
"""
Created on Fri Dec  9 11:28:23 2016

Modified on Tue Sep 5 17:11:00 2017

@author: hamedkhakzad
"""

import pandas
from sklearn import cross_validation
from sklearn.ensemble import BaggingClassifier
from sklearn.tree import DecisionTreeClassifier

def TX_Learning_Train(url_train, num_folds, num_trees):
    names_train = ['lightmonomz', 'heavymonomz', 'rtapexlight', 'rtapexheavy', 'masslight', 'massheavy', 'intlight', 'intheavy', 'intscore', 'rtscore', 'massscore', 'scoresum', 'z', 'xlink', 'train']
    trainframe = pandas.read_csv(url_train, sep=',', names=names_train, header=0)
    
    trainframe = trainframe.fillna(0)
    array_train = trainframe.values
    
    X_train = array_train[:,0:12]
    Y_train = trainframe['train'].fillna(0)

    num_instances = len(X_train)
    seed = 7
    kfold = cross_validation.KFold(n=num_instances, n_folds=num_folds, random_state=seed)
    cart = DecisionTreeClassifier()

    model = BaggingClassifier(base_estimator=cart, n_estimators=num_trees, random_state=seed)
    results = cross_validation.cross_val_score(model, X_train, Y_train, cv=kfold)
    print('Accuracy on k-fold validation is: ',results.mean())
    model.fit(X_train, Y_train, sample_weight=None)
    
    return model