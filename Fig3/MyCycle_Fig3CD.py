#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:51:19 2021

@author: sku

take LFP  of each trial, filter to theta range [4, 12], cycle-by-cycle detect features and plot amplitude distribution 
comparing distCA1 vs. proxCA1
accumulate all old trials together
"""


import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
import statsmodels.api as sm
import pandas as pd
from collections import defaultdict
#%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns
#import neurodsp
from neurodsp.filt import filter_signal
from bycycle.features import compute_features
import os
import h5py
sns.set_style('whitegrid')
plt.style.use('/home/sku/MyPython/MyJNotebook/Cole_2018.mplstyle')

from config import config_dict
from neurodsp.plts.time_series import plot_time_series
from neurodsp.plts.spectral import plot_power_spectra
from neurodsp.utils import create_times
# Import spectral power functions
from neurodsp.spectral import compute_spectrum, rotate_powerlaw
import openpyxl
#%%
# Make directory for saving figures
if not os.path.exists('figs/1'):
    os.makedirs('figs/1')
    
#%%
#define Cohen's d
from numpy import std, mean, sqrt

#correct if the population S.D. is expected to be equal for the two groups.
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)
# define z score.
def z_score2(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (mean(x) - mean(y)) / sqrt(((std(x, ddof=1) ** 2)/nx + (std(y, ddof=1) ** 2) / ny))
    

##%%
# Load LFP
##Nlx

animal='LE84'
lfp_filename_distCA1 ='/home/sku/MyPython/My_testdata/LE84/LE84_tet2.mat'
lfp_filename_proxCA1 ='/home/sku/MyPython/My_testdata/LE84/LE84_tet12.mat'
lfp_filename_distCA1_intrial ='/home/sku/MyPython/My_testdata/LE84/LE84_tet2_test.mat'

fs =4000

#%%
f_d=h5py.File(lfp_filename_distCA1,'r')
lfp_d = np.array(f_d['datalfp'])
lfp_study_d = np.array(f_d['datalfp_study'])
lfp_bsl_d = np.array(f_d['datalfp_bsl'])

f_d_in=h5py.File(lfp_filename_distCA1_intrial,'r')
lfp_d_in = np.array(f_d_in['datalfp_intrial'])
lfp_d_pre = np.array(f_d_in['datalfp_pre'])

f_p=h5py.File(lfp_filename_proxCA1,'r')
lfp_p = np.array(f_p['datalfp'])

lfp_study_p = np.array(f_p['datalfp_study'])
lfp_bsl_p = np.array(f_p['datalfp_bsl'])
 #%%
lfp_d=lfp_d.transpose()
lfp_bsl_d=lfp_bsl_d.transpose()
lfp_study_d=lfp_study_d.transpose()
lfp_p=lfp_p.transpose()
lfp_bsl_p=lfp_bsl_p.transpose()
lfp_study_p=lfp_study_p.transpose()
lfp_d_in=lfp_d_in.transpose()    
lfp_d_pre=lfp_d_pre.transpose()    

oldnew = np.array(f_d['label']) #old: repeated trials; new: non-repeated trials

old = np.where(oldnew == 0) #trial for old are old[0]
new = np.where(oldnew == 1) #trial for new are new[0]
#%%
# start the loop of going through all old trials
mycount=np.arange(0,10,1)
data1=[]
data2=[]
data3=[]
data4=[]
data5=[]
data6=[]
data7=[]
data8=[]
data9=[]
data10=[]
data11=[]
data12=[]

data1_amp=[]
data2_amp=[]
data3_amp=[]
data4_amp=[]
data5_amp=[]
data6_amp=[]
data7_amp=[]
data8_amp=[]
data9_amp=[]
data10_amp=[]
data11_amp=[]
data12_amp=[]
data1_amp_c=[]
data2_amp_c=[]
data3_amp_c=[]
data4_amp_c=[]
data5_amp_c=[]
data6_amp_c=[]
data7_amp_c=[]
data8_amp_c=[]
data9_amp_c=[]
data10_amp_c=[]
data1_period=[]
data2_period=[]
data3_period=[]
data4_period=[]
data5_period=[]
data6_period=[]
data7_period=[]
data8_period=[]
data9_period=[]
data10_period=[]

for n in mycount:
 
    trial=old[0][n]
    trial_new=new[0][n]
    
    #%%
    lfp1=lfp_d[:,1]
    times = create_times(len(lfp1)/fs, fs)

    #%%
    lfp1  = filter_signal(lfp_d[:,n], fs, 'bandpass', (6,12))
    lfp2  = filter_signal(lfp_p[:,n], fs, 'bandpass', (6,12))
    lfp3  = filter_signal(lfp_bsl_d[:,n], fs, 'bandpass', (6,12))
    lfp4  = filter_signal(lfp_bsl_p[:,n], fs, 'bandpass', (6,12))
    lfp5  = filter_signal(lfp_study_d[:,n], fs, 'bandpass', (6,12))
    lfp6  = filter_signal(lfp_study_p[:,n], fs, 'bandpass', (6,12))
    lfp7  = filter_signal(lfp_d[:,n+10], fs, 'bandpass', (6,12))
    lfp8  = filter_signal(lfp_p[:,n+10], fs, 'bandpass', (6,12))
    lfp9  = filter_signal(lfp_d_in[:,old[0][n]], fs, 'bandpass', (6,12))
    lfp10  = filter_signal(lfp_d_in[:,new[0][n]], fs, 'bandpass', (6,12))
    lfp11  = filter_signal(lfp_d_pre[:,old[0][n]], fs, 'bandpass', (6,12))
    lfp12  = filter_signal(lfp_d_pre[:,new[0][n]], fs, 'bandpass', (6,12))
    times = create_times(len(lfp1)/fs, fs)
    #%%
    res=np.where(np.isnan(lfp1)== False )
    t = np.arange(0, len(lfp1)/fs, 1/fs)
    tlim = (res[0][0]/fs, res[0][-1]/fs)
    tidx = np.logical_and(t >= tlim[0], t < tlim[1])
    
    #%%
    # Compute cycle features
    osc_kwargs = {'amplitude_fraction_threshold':0,
                  'amplitude_consistency_threshold':.6,
                  'period_consistency_threshold':.75,
                  'monotonicity_threshold':.8,
                  'N_cycles_min':3}
    # test old dist
    df1 = compute_features(lfp1, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    #  test old prox
    df2 = compute_features(lfp2, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # bsl dist
    df3 = compute_features(lfp3, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # bsl prox
    df4 = compute_features(lfp4, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # study dist
    df5 = compute_features(lfp5, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # study prox
    df6 = compute_features(lfp6, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # test new dist
    df7 = compute_features(lfp7, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    #  test new prox
    df8 = compute_features(lfp8, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # test old dist intrial
    df9 = compute_features(lfp9, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    #  test new dist intrial
    df10 = compute_features(lfp10, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    # pre test old dist intrial
    df11 = compute_features(lfp11, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    #  pre test new dist intrial
    df12 = compute_features(lfp12, fs, [6,12], center_extrema='P',burst_detection_kwargs = osc_kwargs)
    
    data1=np.append(data1,df1['period_consistency'])
    data2=np.append(data2,df2['period_consistency'])
    data3=np.append(data3,df3['period_consistency'])
    data4=np.append(data4,df4['period_consistency'])
    data5=np.append(data5,df5['period_consistency'])
    data6=np.append(data6,df6['period_consistency'])
    
    data1_amp_c=np.append(data1_amp_c,df1['amp_consistency'])
    data2_amp_c=np.append(data2_amp_c,df2['amp_consistency'])
    data3_amp_c=np.append(data3_amp_c,df3['amp_consistency'])
    data4_amp_c=np.append(data4_amp_c,df4['amp_consistency'])
    data5_amp_c=np.append(data5_amp_c,df5['amp_consistency'])
    data6_amp_c=np.append(data6_amp_c,df6['amp_consistency'])

    data1_period=np.append(data1_period,df1['period'])
    data2_period=np.append(data2_period,df2['period'])
    data3_period=np.append(data3_period,df3['period'])
    data4_period=np.append(data4_period,df4['period'])
    data5_period=np.append(data5_period,df5['period'])
    data6_period=np.append(data6_period,df6['period'])
    
    data1_amp=np.append(data1_amp,df1['volt_amp'])
    data2_amp=np.append(data2_amp,df2['volt_amp'])
    data3_amp=np.append(data3_amp,df3['volt_amp'])
    data4_amp=np.append(data4_amp,df4['volt_amp'])
    data5_amp=np.append(data5_amp,df5['volt_amp'])
    data6_amp=np.append(data6_amp,df6['volt_amp'])
    data7_amp=np.append(data7_amp,df7['volt_amp'])
    data8_amp=np.append(data8_amp,df8['volt_amp'])
    data9_amp=np.append(data9_amp,df9['volt_amp'])
    data10_amp=np.append(data10_amp,df10['volt_amp'])
    data11_amp=np.append(data11_amp,df11['volt_amp'])
    data12_amp=np.append(data12_amp,df12['volt_amp'])

     
#%%
# statistics

# plot amp histogram
plt.figure(figsize=(5,5))
plt.hist(df3['volt_amp'], weights=[1/len(df3)]*len(df3),
         bins=np.arange(0, 10000, 50), color='k', alpha=.5, label='baseline')
plt.hist(df2['volt_amp'], weights=[1/len(df2)]*len(df2),
         bins=np.arange(0, 10000, 50), color='g', alpha=.5, label='pCA1')
plt.hist(df1['volt_amp'], weights=[1/len(df1)]*len(df1),
         bins=np.arange(0, 10000, 50), color='r', alpha=.5, label='dCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('Cycle amplitude (a.u.)')
plt.ylabel('fraction of cycles')

#%%
# plot rise-decay asymmetry histogram
plt.figure(figsize=(5,5))
plt.hist(df3['time_rdsym'], weights=[1/len(df3)]*len(df3),
         bins=np.arange(0, 1, .005), color='k', alpha=.5, label='baseline')
plt.hist(df2['time_rdsym'], weights=[1/len(df2)]*len(df2),
         bins=np.arange(0, 1, .005), color='g', alpha=.5, label='pCA1')
plt.hist(df1['time_rdsym'], weights=[1/len(df1)]*len(df1),
         bins=np.arange(0, 1, .005), color='r', alpha=.5, label='dCA1')

plt.legend(fontsize=15)
plt.xlim((0.4,0.6))
plt.xlabel('rise-decay asymmetry')
plt.ylabel('fraction of cycles')

#%%
# plot period consistency histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data2, weights=[1/len(data2)]*len(data2),
         bins=np.arange(0, 1, .025), color='k', alpha=.5, label='pCA1')
plt.hist(data1, weights=[1/len(data1)]*len(data1),
         bins=np.arange(0, 1, .025), color='r', alpha=.5, label='dCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('period_consistency')
plt.ylabel('fraction of cycles')

#%%
# plot amplitude histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data2_amp, weights=[1/len(data2_amp)]*len(data2_amp),
         bins=np.arange(0, 4000,100), color='k', alpha=.5, label='proxCA1')
plt.hist(data1_amp, weights=[1/len(data1_amp)]*len(data1_amp),
         bins=np.arange(0, 4000, 100), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude old')
plt.ylabel('fraction of cycles')

#%%
# plot amplitude histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data4_amp, weights=[1/len(data4_amp)]*len(data4_amp),
         bins=np.arange(0, 4000,100), color='k', alpha=.5, label='proxCA1')
plt.hist(data3_amp, weights=[1/len(data3_amp)]*len(data3_amp),
         bins=np.arange(0, 4000, 100), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude baseline')
plt.ylabel('fraction of cycles')

#%%
# plot amplitude histogram  diff between disCA1 and proxCA1 study phases
plt.figure(figsize=(5,5))

plt.hist(data6_amp, weights=[1/len(data6_amp)]*len(data6_amp),
         bins=np.arange(0, 4000,100), color='k', alpha=.5, label='proxCA1')
plt.hist(data5_amp, weights=[1/len(data5_amp)]*len(data5_amp),
         bins=np.arange(0, 4000, 100), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude study')
plt.ylabel('fraction of cycles')


#%%
# plot amplitude consistency histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data2_amp_c, weights=[1/len(data2_amp_c)]*len(data2_amp_c),
         bins=np.arange(0, 1, .03), color='k', alpha=.5, label='proxCA1')
plt.hist(data1_amp_c, weights=[1/len(data1_amp_c)]*len(data1_amp_c),
         bins=np.arange(0, 1, .03), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude_consistency old')
plt.ylabel('fraction of cycles')


#%%
# plot amplitude consistency histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data6_amp_c, weights=[1/len(data6_amp_c)]*len(data6_amp_c),
         bins=np.arange(0, 1, .03), color='k', alpha=.5, label='proxCA1')
plt.hist(data5_amp_c, weights=[1/len(data5_amp_c)]*len(data5_amp_c),
         bins=np.arange(0, 1, .03), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude_consistency study')
plt.ylabel('fraction of cycles')

#%%
# plot amplitude consistency histogram  slightly diff between disCA1 and proxCA1
plt.figure(figsize=(5,5))

plt.hist(data4_amp_c, weights=[1/len(data4_amp_c)]*len(data4_amp_c),
         bins=np.arange(0, 1, .03), color='k', alpha=.5, label='proxCA1')
plt.hist(data3_amp_c, weights=[1/len(data3_amp_c)]*len(data3_amp_c),
         bins=np.arange(0, 1, .03), color='r', alpha=.5, label='distCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('amplitude_consistency bsl')
plt.ylabel('fraction of cycles')


#%%
# plot period histogram  disCA1 is slightly shorter than proxCA1
plt.figure(figsize=(5,5))
plt.hist(data2_period, weights=[1/len(data2_period)]*len(data2_period),
         bins=np.arange(300, 700, 10), color='k', alpha=.5, label='pCA1')

plt.hist(data1_period, weights=[1/len(data1_period)]*len(data1_period),
         bins=np.arange(300, 700, 10), color='r', alpha=.5, label='dCA1')
#plt.xticks(np.arange(7))
plt.legend(fontsize=15)
#plt.xlim((0,6))
plt.xlabel('period')
plt.ylabel('fraction of cycles')




#%% stats Liliefors test (normality) and Mann-Whitney U-test(non-parametrix 2 independent sample test)
stat, p = mannwhitneyu(data1_amp, data2_amp)
# Use Wilcoxon signed-rank test for dependent data
stat_bsl, p_bsl=wilcoxon(data1_amp,data2_amp)
# concatenate the old and new trials during retrieval
np.mean(np.concatenate((data9_amp,data10_amp)))
sm.stats.diagnostic.lilliefors(data6_amp,dist='norm')


# stats for proximosital paper sup. table 5
p_ret_stu=mannwhitneyu(np.concatenate((data9_amp, data10_amp)),data5_amp)
p_ret_bsl=mannwhitneyu(np.concatenate((data9_amp, data10_amp)),data3_amp)
p_stu_bsl=mannwhitneyu(data5_amp,data3_amp)
p_old_new=mannwhitneyu(data9_amp,data10_amp)
p_pre_int=mannwhitneyu(np.concatenate((data9_amp, data10_amp)),np.concatenate((data11_amp, data12_amp)))

amp_bsl=np.mean(data3_amp)
amp_stu=np.mean(data5_amp)
amp_old=np.mean(data9_amp)
amp_new=np.mean(data10_amp)
amp_pre=np.mean(np.concatenate((data11_amp,data12_amp)))
amp_int=np.mean(np.concatenate((data9_amp,data10_amp)))

