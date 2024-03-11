# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:02:01 2023

@author: zmzhai
"""

import numpy as np
import rc_multi
from bayes_opt import BayesianOptimization
import time
import scipy
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import pickle
import warnings
from joblib import Parallel, delayed
import random
from scipy.ndimage import gaussian_filter1d
import netCDF4
import datetime

warnings.filterwarnings("ignore")

start = time.time()

tt = []
AMOC0, AMOC1, AMOC2 = [], [], []
GM = []

with open('../data.txt') as f:
    for line in f:
        data = line.split()
        tt.append(data[0])
        AMOC0.append(data[1])
        AMOC1.append(data[2])
        AMOC2.append(data[3])
        GM.append(data[4])
        # print(data)
        
tt = [float(i) for i in tt[1:]]
AMOC0 = [float(i) for i in AMOC0[1:]]
AMOC1 = [float(i) for i in AMOC1[1:]]
AMOC2 = [float(i) for i in AMOC2[1:]]
GM = [float(i) for i in GM[1:]]

# sigma = 2.0  # Adjust the sigma value as needed (controls the width of the Gaussian distribution)
# y_smoothed = gaussian_filter1d(AMOC1, sigma)
# y_smoothed = gaussian_filter1d(y_smoothed, sigma)

# AMOC1 = y_smoothed

tt_add = list(np.linspace(2021, 2121, 100*12+1))

tt = tt + tt_add[1:]

# only input one dimensional data: AMOC1
data = np.array(AMOC1)
data = np.expand_dims(data, axis=1)

coupling = False
vali_length = 24
forecast_iteration = True

def target_amoc(d, rho, gamma, alpha, beta, bias, iter_time=5, proportion=1):
    forecast_iteration=True
    
    dim = 1
    # vali_length_all = int(np.floor(450 / vali_length))
    vali_length_all = 300

    config = {}
    config['n'] = 500
    config['d'] = d
    config['alpha'] = alpha
    config['beta'] = beta
    config['gamma'] = gamma
    config['rho'] = rho
    config['bias'] = bias
    
    config['train_length'] = 900 - np.random.randint(20)
    config['wash_length'] = 0
    config['vali_length'] = vali_length
    config['vali_length_all'] = vali_length_all
    config['input_dim'] = dim
    config['output_dim'] = dim
    
    rmse_all = []
    for i in range(iter_time):
        rc_model = rc_multi.Reservoir(data=data, config=config, Win_type=1, forecast=True, forecast_iteration=forecast_iteration)
        rc_model.data_preprocessing() # normalization
        rc_model.initialize_rc()
        train_preditions, train_x = rc_model.train()
        
        rmse_vali = []
        
        rmse, vali_real, vali_pred, _ = rc_model.validation(r_index=0, u_update=False)
        rmse_vali.append(rmse)
        vali_real_all, vali_pred_all = vali_real, vali_pred
        
        for pred_i in range(vali_length_all-1):
            rmse, vali_real, vali_pred, _ = rc_model.validation(r_index=0, u_update=True)
            if pred_i % vali_length == 0:
                vali_real_all = np.concatenate((vali_real_all, vali_real))
                vali_pred_all = np.concatenate((vali_pred_all, vali_pred))
            rmse_vali.append(rmse)
        
        rmse_all.append(np.mean(rmse))
    
    rmse_mean = np.average(sorted(rmse_all)[:int(proportion * iter_time)])
    
    print(rmse_mean)

    return 1 / rmse_mean

def target_amoc_noiter(d, rho, gamma, alpha, beta, bias, iter_time=5, proportion=1):
    
    forecast_iteration=False
    
    dim = 1
    # vali_length_all = int(np.floor(450 / vali_length))
    vali_length_all = 300

    config = {}
    config['n'] = 500
    config['d'] = d
    config['alpha'] = alpha
    config['beta'] = beta
    config['gamma'] = gamma
    config['rho'] = rho
    config['bias'] = bias
    
    config['train_length'] = 900 - np.random.randint(20)
    config['wash_length'] = 0
    config['vali_length'] = vali_length
    config['vali_length_all'] = vali_length_all
    config['input_dim'] = dim
    config['output_dim'] = dim * vali_length
    
    rmse_all = []
    for i in range(iter_time):
        rc_model = rc_multi.Reservoir(data=data, config=config, Win_type=1, forecast=True, forecast_iteration=forecast_iteration)
        rc_model.data_preprocessing() # normalization
        rc_model.initialize_rc()
        train_preditions, train_x = rc_model.train()
        
        rmse, vali_real, vali_pred = rc_model.validation_noiteration(r_index=0, u_update=False)
        rmse_all.append(np.mean(rmse))
    
    rmse_mean = np.average(sorted(rmse_all)[:int(proportion * iter_time)])
    
    print(rmse_mean)

    return 1 / rmse_mean

for i in range(5):
    optimizer = BayesianOptimization(target_amoc,
                                      {'d': (0.01, 1), 'rho': (0.01, 5), 'gamma': (0.01, 5), 'alpha': (0.01, 1), 'beta': (-7, -1), 'bias': (-5, 5)},)

    optimizer.maximize(n_iter=200)
    print('rapid')
    print(optimizer.max)
    
    pkl_file = open('./save_opt/rc_opt_amoc_{}_{}'.format(vali_length, i) + '.pkl', 'wb')
    pickle.dump(optimizer.max, pkl_file)
    pkl_file.close()
    
for i in range(5):
    optimizer = BayesianOptimization(target_amoc_noiter,
                                      {'d': (0.01, 1), 'rho': (0.01, 5), 'gamma': (0.01, 5), 'alpha': (0.01, 1), 'beta': (-7, -1), 'bias': (-5, 5)},)

    optimizer.maximize(n_iter=200)
    print('rapid')
    print(optimizer.max)
    
    pkl_file = open('./save_opt/rc_opt_amoc_noiter_{}_{}'.format(vali_length, i) + '.pkl', 'wb')
    pickle.dump(optimizer.max, pkl_file)
    pkl_file.close()


end = time.time()
print(end - start)


















































