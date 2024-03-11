# -*- coding: utf-8 -*-
"""
Created on Tue May 16 15:38:38 2023

@author: zmzhai
"""

import numpy as np
import random
import scipy.stats as stats
import scipy.sparse as sparse
import networkx as nx
import scipy
import os.path
import datetime
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import pickle
from scipy.ndimage import gaussian_filter1d
import netCDF4
import datetime
import pandas as pd
import scipy.io

class Reservoir:
    def __init__(self, data=None, config=None, Win_type=1, forecast=True, forecast_iteration=False):
        self.data = data
        self.config = config
        self.Win_type = Win_type
        self.forecast = forecast
        self.forecast_iteration = forecast_iteration
        
        # reservoir setting
        self.n = config['n']
        self.d = config['d']
        self.rho = config['rho']
        self.gamma = config['gamma']
        self.alpha = config['alpha']
        self.beta = 10 ** config['beta']
        self.bias = config['bias']
        self.train_length = config['train_length']
        self.wash_length = config['wash_length']
        self.vali_length = config['vali_length']
        self.vali_length_all = config['vali_length_all']
        self.input_dim = config['input_dim']
        self.output_dim = config['output_dim']
        
        self.vali_strat = self.train_length + 1

    def data_preprocessing(self, noise_level=0.0):
        # random_start = np.random.randint(10000, 50000)
        random_start = 0
        self.data = self.data[random_start:, :]
        
        # normalization
        scaler = MinMaxScaler()
        self.data = scaler.fit_transform(self.data)

        noise = np.random.normal(0, noise_level, self.data[:self.train_length, :].shape)
        self.data[:self.train_length, :] = self.data[:self.train_length, :] + noise

    def initialize_rc(self):
        if self.Win_type == 1:
            self.Win = np.random.uniform(-self.gamma, 
                                         self.gamma, (self.n, self.input_dim))
        elif self.Win_type == 2:
            self.Win = np.zeros((self.n, self.input_dim))
            n_win = self.n - self.n % self.input_dim
            index = list(range(n_win))
            random.shuffle(index)
            index = np.reshape(index, (int(n_win/self.input_dim), self.input_dim))
            for di in range(self.input_dim):
                self.Win[index[:, di], di] = np.random.uniform(-self.gamma, 
                                             self.gamma, (int(n_win/self.input_dim), 1)).reshape(1, -1)
        
        graph = nx.erdos_renyi_graph(self.n, self.d, 42, False)
        for (u, v) in graph.edges():
            graph.edges[u, v]['weight'] = np.random.normal(0.0, 1.0)
        self.A = nx.adjacency_matrix(graph).todense()
        rho = max(np.linalg.eig(self.A)[0])
        self.A = (self.rho / abs(rho)) * self.A

    
    def train(self):
        r_train = np.zeros((self.n, self.train_length - self.wash_length))
        y_train = np.zeros((self.output_dim, self.train_length - self.wash_length))
        self.r_end = np.zeros((self.n, 1))
        
        train_x = - np.ones((self.train_length, self.input_dim))
        train_y = np.zeros((self.train_length, self.output_dim))
        
        if self.forecast:
            if self.forecast_iteration:
                train_y[:, :] = self.data[1:self.train_length+1, :]
            else:
                for d_i in range(self.train_length):
                    train_y[d_i, :] = self.data[1+d_i:1+d_i+self.vali_length, :].flatten()
        else:
            train_y[:, :] = self.data[:self.train_length, :]
        
        train_x = self.data[:self.train_length, :]
        
        train_x = np.transpose(train_x)
        train_y = np.transpose(train_y)
        
        r_all = np.zeros((self.n, self.train_length+1))
        for ti in range(self.train_length):
            r_all[:, ti+1] = (1 - self.alpha) * r_all[:, ti] + \
                self.alpha * np.tanh( np.dot(self.A, r_all[:, ti]) + np.dot(self.Win, train_x[:, ti]) +  + self.bias * np.ones((1, self.n))  )
        
        r_out = r_all[:, self.wash_length+1:]
        self.r_end[:] = r_all[:, -1].reshape(-1, 1)
        
        r_train[:, :] = r_out
        y_train[:, :] = train_y[:self.output_dim, self.wash_length:]
        
        self.Wout = np.dot(np.dot(y_train, np.transpose(r_train)), np.linalg.inv(np.dot(r_train, np.transpose(r_train)) + self.beta * np.eye(self.n)) )
        
        train_preditions = np.dot(self.Wout, r_out)
        return train_preditions, train_x
    
    def validation(self, r_index=0, u_update=True):
        # obs_type: will the machine get the input of the obs or not.
        vali_pred = np.zeros((self.vali_length, self.output_dim))
        vali_real = np.zeros((self.vali_length, self.output_dim))
        if self.vali_strat + self.vali_length <= np.shape(self.data)[0]:
            vali_real[:, :] = self.data[self.vali_strat:self.vali_strat + self.vali_length, :]
        else:
            vali_real[:np.shape(self.data[self.vali_strat:, :])[0], :] = self.data[self.vali_strat:, :]
                
        r = self.r_end
        r_record = self.r_end
        
        u = -np.ones((self.input_dim, 1))
        u_record = -np.ones((self.input_dim, 1))
        if u_update:
            u[:] = self.u_record[:]
        else:
            u[:] = self.data[self.vali_strat-1, :].reshape(-1, 1)
        
        for ti in range(self.vali_length):
            r = (1 - self.alpha) * r + self.alpha * np.tanh(np.dot(self.A, r) + np.dot(self.Win, u) + self.bias * np.ones((self.n, 1)))
            r_out = r
                
            pred = np.dot(self.Wout, r_out)
            vali_pred[ti, :] = pred.reshape(1, -1)
            
            u = -np.ones((self.input_dim, 1))
            u[:] = pred[:].reshape(-1, 1)
            
            if ti == r_index:
                r_record = r_out
                u_record = self.data[self.vali_strat+ti, :].reshape(-1, 1)
        
        rmse = rmse_calculation(vali_pred, vali_real)
        # self.vali_strat = self.vali_strat + self.vali_length
        self.vali_strat = self.vali_strat + r_index + 1
        self.r_end[:] = r_record
        self.u_record = u_record
        
        return rmse, vali_real, vali_pred, self.u_record
    
    def validation_noiteration(self, r_index=0, u_update=True, coupling=False):
        vali_pred = np.zeros((self.vali_length_all, self.output_dim))
        vali_real = np.zeros((self.vali_length_all, self.output_dim))
        if self.vali_strat + self.vali_length_all < np.shape(self.data)[0]:
            for d_i in range(self.vali_length_all):
                vali_real[d_i, :] = self.data[self.vali_strat+d_i:self.vali_strat+d_i + self.vali_length, :].flatten()
        else:
            pass
                
        r = self.r_end
        r_record = self.r_end
        
        u = -np.ones((self.input_dim, 1))
        u_record = -np.ones((self.input_dim, 1))
        if u_update:
            u[:] = self.u_record[:]
        else:
            u[:] = self.data[self.vali_strat-1, :].reshape(-1, 1)
            
        for ti in range(self.vali_length_all):
            r = (1 - self.alpha) * r + self.alpha * np.tanh(np.dot(self.A, r) + np.dot(self.Win, u) + self.bias * np.ones((self.n, 1)))
            r_out = r
                
            pred = np.dot(self.Wout, r_out)
            vali_pred[ti, :] = pred.reshape(1, -1)
    
            u = self.data[self.vali_strat+ti, :].reshape(-1, 1)
            
        rmse = rmse_calculation(vali_pred, vali_real)
        self.vali_strat = self.vali_strat + self.vali_length_all
        # self.vali_strat = self.vali_strat + r_index + 1
        self.r_end[:] = r_record
        self.u_record = u_record
        
        return rmse, vali_real, vali_pred


def rmse_calculation(A, B):
    # calculate rmse
    return (np.sqrt(np.square(np.subtract(A, B)).mean()))

def mse_calculation(A, B):
    # calculate mse
    return (np.square(np.subtract(A, B)).mean())

def slope_calculation(A, B):
    x = range(len(A))
    
    m_A, _ = np.polyfit(x, A, 1)
    m_B, _ = np.polyfit(x, B, 1)
    
    return np.abs(m_A - m_B)


if __name__ == '__main__':
    print('------rc------')
    
    ###################
    # sst dataset
    
    # time = []
    # AMOC0, AMOC1, AMOC2 = [], [], []
    # GM = []

    # with open('./data.txt') as f:
    #     for line in f:
    #         data = line.split()
    #         time.append(data[0])
    #         AMOC0.append(data[1])
    #         AMOC1.append(data[2])
    #         AMOC2.append(data[3])
    #         GM.append(data[4])
    #         # print(data)
            
    # time = [float(i) for i in time[1:]]
    # AMOC0 = [float(i) for i in AMOC0[1:]]
    # AMOC1 = [float(i) for i in AMOC1[1:]]
    # AMOC2 = [float(i) for i in AMOC2[1:]]
    # GM = [float(i) for i in GM[1:]]
    
    # # sigma = 2.0  # Adjust the sigma value as needed (controls the width of the Gaussian distribution)
    # # y_smoothed = gaussian_filter1d(AMOC1, sigma)
    # # y_smoothed = gaussian_filter1d(y_smoothed, sigma)
    
    # # AMOC1 = y_smoothed
    
    # time_add = list(np.linspace(2021, 2121, 100*12+1))
    
    # time = time + time_add[1:]
    # date_list = time
    
    # # only input one dimensional data: AMOC1
    # data = np.array(AMOC1)
    # data = np.expand_dims(data, axis=1)
    
    # pkl_file = open('./save_opt/rc_opt_amoc_24_4' + '.pkl', 'rb') # 2 is bad
    
    # # pkl_file = open('./save_opt/rc_opt_amoc_noiter_24_4' + '.pkl', 'rb') # 2 is bad
    # opt_results = pickle.load(pkl_file)
    # pkl_file.close()
    # opt_params = opt_results['params']
    
    # config = {}

    # train_length = 900
    # vali_length = 36
    # # vali_length_all = int(np.floor(800 / vali_length))
    # vali_length_all = 870
    
    # config['n'] = 500
    # config['d'] = opt_params['d']
    # config['alpha'] = opt_params['alpha']
    # config['beta'] = opt_params['beta']
    # config['gamma'] = opt_params['gamma']
    # config['rho'] = opt_params['rho']
    # config['bias'] = opt_params['bias']
    
    # # config['alpha'] = 1
    # # config['rho'] = 0
    
    # config['train_length'] = train_length
    # config['wash_length'] = 0
    # config['vali_length'] = vali_length
    # config['vali_length_all'] = vali_length_all
    # config['input_dim'] = 1
    # config['output_dim'] = 1
    
    # config['output_dim'] = vali_length
    

    ###################
    # rapid dataset
    
    file2read = netCDF4.Dataset('./RAPID_data/moc_transports.nc', 'r')

    print(file2read.variables.keys())
    print(file2read)
    # read t
    time = file2read.variables['time']
    print(time)
    t = time[:]
    # read layer transports

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
    
    pkl_file = open('./save_opt/rc_opt_rapid_90_2' + '.pkl', 'rb')
    # pkl_file = open('./save_opt/rc_opt_rapid_noiter_90_3' + '.pkl', 'rb')
    opt_results = pickle.load(pkl_file)
    pkl_file.close()
    opt_params = opt_results['params']
    
    config = {}

    train_length = 6000
    vali_length = 40
    # vali_length_all = int(np.floor(800 / vali_length))
    vali_length_all = 6000
    
    config['n'] = 500
    config['d'] = opt_params['d']
    config['alpha'] = opt_params['alpha']
    config['beta'] = opt_params['beta']
    config['gamma'] = opt_params['gamma']
    config['rho'] = opt_params['rho']
    config['bias'] = opt_params['bias']
    
    config['train_length'] = train_length
    config['wash_length'] = 0
    config['vali_length'] = vali_length
    config['vali_length_all'] = vali_length_all
    config['input_dim'] = 4
    config['output_dim'] = 4
    
    # config['output_dim'] = vali_length * 4
    
    ###################
    # move dataset
    
    # file2read = netCDF4.Dataset('./recover_data/os_MOVE_TRANSPORTS.nc', 'r')

    # print(file2read.variables.keys())
    # print(file2read)
    # # read t
    # time = file2read.variables['TIME']
    # print(time)
    # t = time[:]
    
    # total = file2read.variables['TRANSPORT_TOTAL']
    # total = total[:]
    # total = total[32:6742]
    # missing_indices = np.where(np.isnan(total))[0]
    # # Interpolate missing values using linear interpolation
    # total[missing_indices] = np.interp(missing_indices, 
    #                                           np.where(~np.isnan(total))[0], 
    #                                           total[~np.isnan(total)])
    
    # # set time
    # start_date = datetime.datetime(2000, 1, 31)
    # num_intervals = len(t)
    # interval = datetime.timedelta(days=1)
    # date_list = [start_date + interval * i for i in range(num_intervals)]
    
    # t = t[32:6742]
    
    # data = np.expand_dims(total, axis=1)
    
    # pkl_file = open('./save_opt/rc_opt_move_60_2' + '.pkl', 'rb')
    # # pkl_file = open('./save_opt/rc_opt_move_noiter_60_4' + '.pkl', 'rb')
    # opt_results = pickle.load(pkl_file)
    # pkl_file.close()
    # opt_params = opt_results['params']
    
    # config = {}

    # train_length = 3000
    # vali_length = 20
    # # vali_length_all = int(np.floor(800 / vali_length))
    # vali_length_all = 3600
    
    # config['n'] = 500
    # config['d'] = opt_params['d']
    # config['alpha'] = opt_params['alpha']
    # config['beta'] = opt_params['beta']
    # config['gamma'] = opt_params['gamma']
    # config['rho'] = opt_params['rho']
    # config['bias'] = opt_params['bias']
    
    # config['train_length'] = train_length
    # config['wash_length'] = 0
    # config['vali_length'] = vali_length
    # config['vali_length_all'] = vali_length_all
    # config['input_dim'] = 1
    # config['output_dim'] = 1
    
    
    # config['output_dim'] = vali_length
    
    ###################
    # amoc model dataset
    # file2read = scipy.io.loadmat('./recover_data/amoc_model_data.mat')
    # # file2read = netCDF4.Dataset('./recover_data/os_MOVE_TRANSPORTS.nc', 'r')
    
    # data = file2read['ts_train']
    # data = data[500:]
    # t = range(len(data))
    
    # # set time
    # start_date = datetime.datetime(2000, 1, 31)
    # num_intervals = len(t)
    # interval = datetime.timedelta(days=1)
    # date_list = [start_date + interval * i for i in range(num_intervals)]
    
    # # t = t[32:6742]
    
    # # data = np.expand_dims(total, axis=1)
    
    # pkl_file = open('./save_opt/rc_opt_model_40_0' + '.pkl', 'rb')
    # opt_results = pickle.load(pkl_file)
    # pkl_file.close()
    # opt_params = opt_results['params']
    
    # config = {}

    # train_length = 5000
    # vali_length = 100
    # # vali_length_all = int(np.floor(800 / vali_length))
    # vali_length_all = 5600
    
    # config['n'] = 500
    # config['d'] = opt_params['d']
    # config['alpha'] = opt_params['alpha']
    # config['beta'] = opt_params['beta']
    # config['gamma'] = opt_params['gamma']
    # config['rho'] = opt_params['rho']
    # config['bias'] = opt_params['bias']
    
    # config['train_length'] = train_length
    # config['wash_length'] = 0
    # config['vali_length'] = vali_length
    # config['vali_length_all'] = vali_length_all
    # config['input_dim'] = 1
    # config['output_dim'] = 1
    
    
    # config['output_dim'] = vali_length
    
    
    ###########################################
    # in common
    # autoregressive or direct prediction
    
    # autoregressive prediction
    forecast_iteration = True
    
    rc = Reservoir(data=data, config=config, Win_type=1, forecast=True, forecast_iteration=forecast_iteration)
    rc.data_preprocessing() # normalization
    rc.initialize_rc()
    train_preditions, train_x = rc.train()
    
    rmse_all = []
    slope_all = []
    
    rmse, vali_real, vali_pred, u_record = rc.validation(r_index=0, u_update=False)
    rmse_all.append(rmse)
    slope_all.append(slope_calculation(vali_real, vali_pred))
    vali_real_all, vali_pred_all, u_record_all = vali_real, vali_pred, u_record
    vali_real_matrix, vali_pred_matrix = [vali_real], [vali_pred]
    
    for pred_i in range(vali_length_all-1):
        rmse, vali_real, vali_pred, u_record = rc.validation(r_index=0, u_update=True)
        u_record_all = np.concatenate((u_record_all, u_record))
        vali_real_matrix.append(vali_real)
        vali_pred_matrix.append(vali_pred)
        if pred_i % vali_length == 0:
            vali_real_all = np.concatenate((vali_real_all, vali_real))
            vali_pred_all = np.concatenate((vali_pred_all, vali_pred))
        rmse_all.append(rmse)
        slope_all.append(slope_calculation(vali_real, vali_pred))
        
    rmse_mean = np.mean(rmse_all)
    vali_real_matrix = np.array(vali_real_matrix)
    vali_pred_matrix = np.array(vali_pred_matrix)
    
    # save file for move model
    # pkl_file = open('./save_matrix/rc_pred_rapid_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_real_matrix, pkl_file)
    # pickle.dump(vali_pred_matrix, pkl_file)
    # pickle.dump(vali_real_all, pkl_file)
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real, pkl_file)
    # pickle.dump(vali_pred, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()
    
    # save file for move model
    # pkl_file = open('./save_matrix/rc_pred_move_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_real_matrix, pkl_file)
    # pickle.dump(vali_pred_matrix, pkl_file)
    # pickle.dump(vali_real_all, pkl_file)
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real, pkl_file)
    # pickle.dump(vali_pred, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()
    
    # save file for amoc
    # pkl_file = open('./save_matrix/rc_pred_amoc_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_real_matrix, pkl_file)
    # pickle.dump(vali_pred_matrix, pkl_file)
    # pickle.dump(vali_real_all, pkl_file)
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real, pkl_file)
    # pickle.dump(vali_pred, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()
    
    # direct prediction
    # forecast_iteration = False
    
    # rc = Reservoir(data=data, config=config, Win_type=1, forecast=True, forecast_iteration=forecast_iteration)
    # rc.data_preprocessing() # normalization
    # rc.initialize_rc()
    # train_preditions, train_x = rc.train()
    
    # rmse_all = []
    
    # rmse, vali_real, vali_pred = rc.validation_noiteration(r_index=0, u_update=False)
    # rmse_all.append(rmse)
    
    # rmse_mean = np.mean(rmse_all)
    
    # vali_pred_all = []
    # for i in range(vali_length_all):
    #     if i % vali_length == 0:
    #         vali_pred_all.append(vali_pred[i, :])
 
    # vali_pred_all = np.array(vali_pred_all).flatten()
    
    # pkl_file = open('./save_matrix/rc_timeseries_rapid_noiter_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real, pkl_file)
    # pickle.dump(vali_pred, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()

    # fig, ax = plt.subplots(1, 1, figsize=(8,6), constrained_layout=True)
 
    # ax.plot(date_list[:train_length], train_x[0, :], color='blue', label='real')
    # ax.plot(date_list[:train_length], train_preditions[0, :], color='orange', label='train')
    # plot_length = len(vali_pred[:, 0])
    # # ax.plot(time[train_length:train_length+len(u_record_all[:, 0])], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    # ax.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real[:plot_length, 0]), color='blue', alpha=0.8)
    # ax.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length], color='green', label='pred', linewidth=3)
    
    # ax.set_xlabel('time')
    # ax.set_ylabel('moc')
    # ax.legend()
    
    # plt.show()
    
    # # plot zoom in figures
    # fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
    
    # plot_length = len(vali_pred_all[:])
    # # ax.plot(time[train_length+vali_length:train_length+len(u_record_all[:, 0])+vali_length], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    # ax.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real[:plot_length, 0]), color='blue', alpha=0.8)
    # ax.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length], color='green', label='pred', linewidth=3)
    
    # ax.set_xlabel('time')
    # ax.set_ylabel('moc')
    # ax.legend()
    
    # plt.show()
    ###############################################
    # plot figures 
    # plot figures
    
    # pkl_file = open('./save_matrix/rc_timeseries_move_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real_all, pkl_file)
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()
    
    fig, ax = plt.subplots(1, 1, figsize=(8,6), constrained_layout=True)
    
    ax.plot(date_list[:train_length], train_x[0, :], color='blue', label='real')
    ax.plot(date_list[:train_length], train_preditions[0, :], color='orange', label='train')
    plot_length = len(vali_pred_all[:, 0])
    # ax.plot(time[train_length:train_length+len(u_record_all[:, 0])], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    ax.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real_all[:plot_length, 0]), color='blue', alpha=0.8)
    ax.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 0], color='green', label='pred', linewidth=3)
    
    ax.set_xlabel('time')
    ax.set_ylabel('moc')
    ax.legend()
    
    plt.show()
    
    # # plot zoom in figures
    # fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
    
    # plot_length = len(vali_pred_all[:, 0])
    # # ax.plot(time[train_length+vali_length:train_length+len(u_record_all[:, 0])+vali_length], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    # ax.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real_all[:plot_length, 0]), color='blue', alpha=0.8)
    # ax.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 0], color='green', label='pred', linewidth=3)
    
    # ax.set_xlabel('time')
    # ax.set_ylabel('moc')
    # ax.legend()
    
    # plt.show()

    # # plot 4 dimensional systems.
    # # plot figures
    
    # pkl_file = open('./save_matrix/rc_timeseries_rapid_valilen_{}'.format(vali_length) + '.pkl', 'wb')
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(train_x, pkl_file)
    # pickle.dump(train_preditions, pkl_file)
    # pickle.dump(vali_real_all, pkl_file)
    # pickle.dump(vali_pred_all, pkl_file)
    # pickle.dump(date_list, pkl_file)
    # pickle.dump(train_length, pkl_file)
    # pkl_file.close()
    
    fig, ax = plt.subplots(4, 1, figsize=(10,10), constrained_layout=True)
    ax0, ax1, ax2, ax3 = ax[0], ax[1], ax[2], ax[3]
    
    ax0.plot(date_list[:train_length], train_x[0, :], color='blue', label='real')
    ax0.plot(date_list[:train_length], train_preditions[0, :], color='orange', label='train')
    plot_length = len(vali_pred_all[:, 0])
    # ax.plot(time[train_length:train_length+len(u_record_all[:, 0])], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    ax0.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real_all[:plot_length, 0]), color='blue', alpha=0.8)
    ax0.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 0], color='green', label='pred', linewidth=3)
    
    ax0.set_xlabel('time')
    ax0.set_ylabel('moc')
    ax0.legend()
    
    ax1.plot(date_list[:train_length], train_x[1, :], color='blue', label='real')
    ax1.plot(date_list[:train_length], train_preditions[1, :], color='orange', label='train')
    plot_length = len(vali_pred_all[:, 0])
    ax1.plot(date_list[train_length-1:train_length+plot_length], [train_x[1, :][-1]]+list(vali_real_all[:plot_length,1]), color='blue', alpha=0.8)
    ax1.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 1], color='green', label='pred', linewidth=3)
    
    ax1.set_xlabel('time')
    ax1.set_ylabel('upper mid ocean')
    ax1.legend()
    
    ax2.plot(date_list[:train_length], train_x[2, :], color='blue', label='real')
    ax2.plot(date_list[:train_length], train_preditions[2, :], color='orange', label='train')
    plot_length = len(vali_pred_all[:, 0])
    ax2.plot(date_list[train_length-1:train_length+plot_length], [train_x[2, :][-1]]+list(vali_real_all[:plot_length, 2]), color='blue', alpha=0.8)
    ax2.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 2], color='green', label='pred', linewidth=3)
    
    ax2.set_xlabel('time')
    ax2.set_ylabel('gulf stream')
    ax2.legend()
    
    ax3.plot(date_list[:train_length], train_x[3, :], color='blue', label='real')
    ax3.plot(date_list[:train_length], train_preditions[3, :], color='orange', label='train')
    plot_length = len(vali_pred_all[:, 0])
    ax3.plot(date_list[train_length-1:train_length+plot_length], [train_x[3, :][-1]]+list(vali_real_all[:plot_length, 3]), color='blue', alpha=0.8)
    ax3.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 3], color='green', label='pred', linewidth=3)
    
    ax3.set_xlabel('time')
    ax3.set_ylabel('ekman')
    ax3.legend()

    plt.show()
    
    # # plot zoom in figures
    # fig, ax = plt.subplots(4, 1, figsize=(10, 10), constrained_layout=True)
    # ax0, ax1, ax2, ax3 = ax[0], ax[1], ax[2], ax[3]
    
    # plot_length = len(vali_pred_all[:, 0])
    # # ax.plot(time[train_length+vali_length:train_length+len(u_record_all[:, 0])+vali_length], u_record_all[:, 0], color='lightpink', label='trace', alpha=0.6)
    # ax0.plot(date_list[train_length-1:train_length+plot_length], [train_x[0, :][-1]]+list(vali_real_all[:plot_length, 0]), color='blue', alpha=0.8)
    # ax0.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 0], color='green', label='pred', linewidth=3)
    
    # ax0.set_xlabel('time')
    # ax0.set_ylabel('moc')
    # ax0.legend()
    
    # ax1.plot(date_list[train_length-1:train_length+plot_length], [train_x[1, :][-1]]+list(vali_real_all[:plot_length,1]), color='blue', alpha=0.8)
    # ax1.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 1], color='green', label='pred', linewidth=3)
    
    # ax1.set_xlabel('time')
    # ax1.set_ylabel('upper mid ocean')
    # ax1.legend()
    
    # ax2.plot(date_list[train_length-1:train_length+plot_length], [train_x[2, :][-1]]+list(vali_real_all[:plot_length, 2]), color='blue', alpha=0.8)
    # ax2.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 2], color='green', label='pred', linewidth=3)
    
    # ax2.set_xlabel('time')
    # ax2.set_ylabel('gulf stream')
    # ax2.legend()
    
    # ax3.plot(date_list[train_length-1:train_length+plot_length], [train_x[3, :][-1]]+list(vali_real_all[:plot_length, 3]), color='blue', alpha=0.8)
    # ax3.plot(date_list[train_length:train_length+plot_length], vali_pred_all[:plot_length, 3], color='green', label='pred', linewidth=3)
    
    # ax3.set_xlabel('time')
    # ax3.set_ylabel('ekman')
    # ax3.legend()
    
    # plt.show()
    
    # smooth data
    # window = 10
    # data_smooth = data

    # df = pd.DataFrame({'dd': data[:train_length, :].flatten()})
    # smoothed_data = df['dd'].rolling(window=window).mean()
    
    
    # data_smooth[window:train_length, :] = np.expand_dims(smoothed_data[window:], axis=1)

    # fig, ax = plt.subplots(2, 1, figsize=(8, 10))
    # ax0, ax1 = ax[0], ax[1]
    
    # ax0.plot(data[:train_length, :])
    # # ax1.plot(data_smooth[:train_length, :])
    
    # plt.show()


































