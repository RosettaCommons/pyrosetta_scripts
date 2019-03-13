#!/usr/bin/env python
"""
Created on Fri Dec  9 11:28:23 2016

Modified on Tue Sep 5 17:11:00 2017

@author: hamedkhakzad
"""

import pandas
import numpy as np
from csv2db import csv2db
from sklearn import cross_validation
from sklearn.ensemble import BaggingClassifier
from sklearn.tree import DecisionTreeClassifier

def TX_Learning_Test(model, url_test, num):


    names_test = ['lightmonomz', 'heavymonomz', 'rtapexlight', 'rtapexheavy', 'masslight', 'massheavy', 'intlight', 'intheavy', 'intscore', 'rtscore', 'massscore', 'scoresum', 'z','xlink']
    names_out =  ['lightmonomz', 'heavymonomz', 'rtapexlight', 'rtapexheavy', 'masslight', 'massheavy', 'intlight', 'intheavy', 'intscore', 'rtscore', 'massscore', 'scoresum', 'z','xlink', 'model_tag']    
    testframe = pandas.read_csv(url_test, sep='\t', names=names_test, header=0)
    
    array_test = testframe.values

    X_test = array_test[:,0:12]
    xlink_list = array_test[:,13]

    model_tag = model.predict(X_test)

    print('Number of rows with tag 1: ', sum(model_tag))

    result_array = np.c_[array_test, model_tag]
    result = pandas.DataFrame(data=result_array[0:,0:], columns=names_out)

    top_XL_list = csv2db(result, num)

    top_XL_file_name = 'top_XL'+url_test+'.txt'
    top_XL_file = open(top_XL_file_name, 'w')
    for item in top_XL_list:
    	top_XL_file.write("%s\n" % item)