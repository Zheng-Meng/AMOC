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

file2read = netCDF4.Dataset('../RAPID_data/moc_transports.nc', 'r')

print(file2read.variables.keys())
print(file2read)
# read t
time = file2read.variables['time']
print(time)
t = time[:]
# read layer transports
therm10 = file2read.variables['t_therm10']
aiw = file2read.variables['t_aiw10']
ud = file2read.variables['t_ud10']
ld = file2read.variables['t_ld10']
bw = file2read.variables['t_bw10']
therm10 = therm10[:]
aiw = aiw[:]
ud = ud[:]
ld = ld[:]
bw = bw[:]

# read moc
moc = file2read.variables['moc_mar_hc10']
umo = file2read.variables['t_umo10']
gs = file2read.variables['t_gs10']
ek = file2read.variables['t_ek10']
moc = moc[:]
umo = umo[:]
gs = gs[:]
ek = ek[:]
# set time
start_date = datetime.datetime(2004, 4, 1)
num_intervals = len(t)
interval = datetime.timedelta(days=0.5)
date_list = [start_date + interval * i for i in range(num_intervals)]

moc = moc[10:12193]
umo = umo[10:12193]
gs = gs[10:12193]
ek = ek[10:12193]
t = t[10:12193]
date_list = date_list[10:12193]

moc = np.expand_dims(moc, axis=1)
umo = np.expand_dims(umo, axis=1)
gs = np.expand_dims(gs, axis=1)
ek = np.expand_dims(ek, axis=1)

data = np.concatenate((moc, umo, gs, ek), axis=1)

vali_length = 90

def target_amoc(d, rho, gamma, alpha, beta, bias, iter_time=5, proportion=1):
    forecast_iteration=True
    
    dim = 4
    # vali_length_all = int(np.floor(450 / vali_length))
    vali_length_all = 2000

    config = {}
    config['n'] = 500
    config['d'] = d
    config['alpha'] = alpha
    config['beta'] = beta
    config['gamma'] = gamma
    config['rho'] = rho
    config['bias'] = bias
    
    config['train_length'] = 6000 + np.random.randint(100) - 50
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
    
    dim = 4
    # vali_length_all = int(np.floor(450 / vali_length))
    vali_length_all = 2000

    config = {}
    config['n'] = 500
    config['d'] = d
    config['alpha'] = alpha
    config['beta'] = beta
    config['gamma'] = gamma
    config['rho'] = rho
    config['bias'] = bias
    
    config['train_length'] = 6000 + np.random.randint(100) - 50
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
    
    pkl_file = open('./save_opt/rc_opt_rapid_{}_{}'.format(vali_length, i) + '.pkl', 'wb')
    pickle.dump(optimizer.max, pkl_file)
    pkl_file.close()
    
for i in range(5):
    optimizer = BayesianOptimization(target_amoc_noiter,
                                      {'d': (0.01, 1), 'rho': (0.01, 5), 'gamma': (0.01, 5), 'alpha': (0.01, 1), 'beta': (-7, -1), 'bias': (-5, 5)},)

    optimizer.maximize(n_iter=200)
    print('rapid')
    print(optimizer.max)
    
    pkl_file = open('./save_opt/rc_opt_rapid_noiter_{}_{}'.format(vali_length, i) + '.pkl', 'wb')
    pickle.dump(optimizer.max, pkl_file)
    pkl_file.close()


end = time.time()
print(end - start)


















































