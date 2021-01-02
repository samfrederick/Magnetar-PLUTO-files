#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 12:45:39 2020

@author: Sam Frederick
"""
import os
os.chdir('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/Analysis')

import pandas as pd
import shutil as sh
import matplotlib.pyplot as plt
import seaborn as sns
import math as mt
from datetime import datetime
from differentiate import TimeDerivative # Must set path to Analysis folder
sns.set_style('darkgrid')

data_path = '/media/sam/BE48068348063B23/Simulation_Results/201206/'
extern_drive_path = '/media/sam/ASTRO_DATA/201206/'

folder = 'Ellipticity Analysis'

def MOI_Import(filename):
    """
    Parameters
    ----------
    filename : string
        Filename of dataset. Can be relative, but need to set cwd to folder
        where data are placed.

    Returns
    -------
    dataset : Pandas DataFrame
        Dataset with moment of inertia tensor component values.
        Ellipticity is calculated for each recorded timestamp in datafile.
    """
    # Import .txt file with principal MOI values
    if filename.endswith('.txt'):
        dataset = pd.read_csv(filename, delimiter=r'\s+', header=0)

        # Round the time column, set as index
        dataset['t'] = round(dataset['t'], 4)
        dataset = dataset.set_index('t', drop=True)

    # Import comma delimited inertia tensor csv files
    if filename.endswith('.csv'):
        dataset = pd.read_csv(filename, index_col=('t'))

    # Compute ellipticity using Izz(t=0) for I_0 denominator term
    dataset['ellip'] = (dataset.Izz - dataset.Ixx)/(dataset.loc[0.0, 'Izz'])

    # Compute time derivative for ellipticity
    dataset = TimeDerivative(dataset, column='ellip', time_column='t')
    dataset = TimeDerivative(dataset, column='d_ellip', time_column='t')

    # Compute first time derivative of principal MOIs
    for moi in ['Ixx', 'Iyy', 'Izz']:
        dataset = TimeDerivative(dataset, column=moi, time_column='t')

    # Compute second time derivative of principal MOIs
    for d_moi in ['d_Ixx', 'd_Iyy', 'd_Izz']:
        dataset = TimeDerivative(dataset, column=d_moi, time_column='t')

    return dataset

def Plot_Ellip_Timeseries(df, savefig=False):
    """
    """
    fig, axs = plt.subplots(3, 1, figsize=(9, 7))
    plt.subplots_adjust(hspace=.35)

    time_lims = -0.5, mt.ceil(df.index[-1])

    axs[0].set_title(r'Stellar Ellipticity ($\epsilon$)')
    axs[0].plot(df.index, df.ellip, c='#aa4ab0')
    axs[0].legend(labels=[r'$\epsilon$'], loc='lower right')
    axs[0].axhline(y=0, color='#949494', linestyle='--')
    axs[0].set_ylim(-0.15, 0.15)
    axs[0].set_xlim(time_lims)

    axs[1].set_title('First Derivative')
    axs[1].plot(df.index, df.d_ellip, c='#d6004c')
    axs[1].axhline(y=0, color='#949494', linestyle='--')
    axs[1].legend(labels=[r'$\partial_t$ $\epsilon$'], loc='lower right')
    axs[1].set_ylim(-0.15, 0.15)
    axs[1].set_xlim(time_lims)

    axs[2].set_title('Second Derivative')
    axs[2].plot(df.index, df.d_d_ellip, c='#465aea')
    axs[2].axhline(y=0, color='#949494', linestyle='--')
    axs[2].legend(labels=[r'$\partial^2_t$ $\epsilon$'], loc='lower right')
    axs[2].set_ylim(-1.5, 1.5)
    axs[2].set_xlim(time_lims)
    axs[2].set_xlabel('Time (s)')

    if savefig:
        today = datetime.now().strftime('%Y%m%d_%H%M%S')
        plt.savefig(folder + '/201206_MOI_ellip_'+today+'.png', dpi=300)
        #plt.close()


def Plot_MOI_Timeseries(df, savefig=False):
    """
    """
    fig, axs = plt.subplots(3, 1, figsize=(9, 8))
    plt.subplots_adjust(hspace=.35)
    time_lims = -0.5, mt.ceil(df.index[-1])

    axs[0].set_title('Principal Moments of Inertia I$_{xx}$ and I$_{zz}$')
    axs[0].plot(df.index, df.Ixx, df.index, df.Izz)
    axs[0].legend(labels=['I$_{xx}$', 'I$_{zz}$'])
    axs[0].set_xlim(time_lims)

    axs[1].set_title('First Derivative')
    axs[1].plot(df.index, df.d_Ixx, df.index, df.d_Izz)
    axs[1].legend(labels=[r'$\partial_t$ I$_{xx}$',
                          r'$\partial_t$ I$_{zz}$'])
    axs[1].axhline(y=0, color='#949494', linestyle='--')
    axs[1].set_ylim(-2.5e44, 2.5e44)
    axs[1].set_xlim(time_lims)

    axs[2].set_title('Second Derivative')
    axs[2].plot(df.index, df.d_d_Ixx, df.index, df.d_d_Izz)
    axs[2].legend(labels=[r'$\partial^2_t$ I$_{xx}$',
                          r'$\partial^2_t$ I$_{zz}$'])
    axs[2].axhline(y=0, color='#949494', linestyle='--')
    axs[2].set_ylim(-2.5e45, 2.5e45)
    axs[2].set_xlabel('Time (s)')
    axs[2].set_xlim(time_lims)

    if savefig:
        today = datetime.now().strftime('%Y%m%d_%H%M%S')
        plt.savefig(folder + '/201206_MOI_'+today+'.png', dpi=300)
        #plt.close()


def Copy_Data(path, destination):
    """
    """
    # Copy every tenth data file from data directory to external drive
    for cwd, folders, files in os.walk(data_path):
        for filename in files:
            if 'data' in filename:
                file_n = int(filename.replace('data.', '').replace('.vtk', ''))
                if file_n % 50 == 0:
                    # Only write if file not in folder
                    if not os.path.exists(extern_drive_path+filename):
                        sh.copyfile(data_path+filename,
                                    extern_drive_path+filename)


# -----------------------------------------------------------------------------

if __name__ == '__main__':
    #Copy_Data(path=data_path, destination=extern_drive_path)
     
    df = MOI_Import(folder +'/201206_InertiaTensor.csv')
    
    Plot_Ellip_Timeseries(df, savefig=True)
    
    Plot_MOI_Timeseries(df, savefig=True)
