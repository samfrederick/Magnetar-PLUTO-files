#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 12:45:39 2020

@author: Sam Frederick
"""
import pandas as pd
import os
import shutil as sh
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
sns.set_style('darkgrid')


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

    # Compute first time derivative of principal MOIs
    for moi in ['Ixx', 'Iyy', 'Izz']:
        dataset = Differentiate(dataset, column=moi, time_column='t')

    # Compute second time derivative of principal MOIs
    for d_moi in ['d_Ixx', 'd_Iyy', 'd_Izz']:
        dataset = Differentiate(dataset, column=d_moi, time_column='t')

    return dataset


def Differentiate(df, column, time_column):
    """
    Compute the time derivative for pandas column data using the central
    difference method and foward and backward schemes for endpoints.
    """
    idx_name = df.index.name
    df = df.reset_index()
    end_idx = len(df) - 1

    # Backward diff for column tail
    tail_data = df[[time_column, column]].tail(n=2)
    tail_dt = (tail_data.loc[end_idx, time_column]
               - tail_data.loc[end_idx-1, time_column])
    tail_deriv = (1/tail_dt)*(tail_data.loc[end_idx, column]
                              - tail_data.loc[end_idx-1, column])

    # Forward diff for column head
    head_data = df[[time_column, column]].head(n=2)
    head_dt = (head_data.loc[1, time_column]
               - head_data.loc[0, time_column])
    head_deriv = (1/head_dt)*(head_data.loc[1, column]
                              - head_data.loc[0, column])

    # Central diff for column body
    central_data = df.loc[1:end_idx-1, [time_column, column]]
    central_dt = pd.Series(central_data[1:][time_column].values
                           - central_data[:-1][time_column].values)
    central_diff = pd.Series(central_data[1:-1][column].values
                             - central_data[0:-2][column].values)
    central_deriv = (1/central_dt)*central_diff
    central_deriv.index = central_data[0:-1][time_column]

    # Merge derivative values to staging dataframe
    deriv_df = pd.DataFrame(index=df[idx_name], columns=['d' + column])
    deriv_df.iloc[0, 0] = head_deriv
    deriv_df.iloc[-1, 0] = tail_deriv
    deriv_df.iloc[1:-1, 0] = central_deriv[:]

    # Join derivative column with passed dataframe
    df = df.set_index(idx_name, drop=True)
    df['d_'+column] = deriv_df['d' + column]

    return df


def Plot_MOI_Ellip(df, savefig=False):
    """
    """
    fig, axs = plt.subplots(2, 1, figsize=(9, 7))
    plt.subplots_adjust(hspace=.35)

    axs[0].set_title('Principal Moments of Inertia')
    axs[0].plot(df.index, df.Ixx, df.index, df.Iyy, df.index, df.Izz)
    axs[0].legend(labels=['Ixx', 'Iyy', 'Izz'])

    axs[1].set_title('Stellar Ellipticity')
    axs[1].plot(df.index, df.ellip, c='purple')
    axs[1].legend(labels=['Ellipticity'])
    axs[1].axhline(y=0, color='#949494', linestyle='--')
    axs[1].set_ylim(-0.1, 0.1)

    if savefig:
        today = datetime.now().strftime('%Y%m%d_%H%M%S')
        plt.savefig('201206_MOI_ellip_'+today+'.png', dpi=300)
        plt.close()


df = MOI_Import('201206_InertiaTensor.csv')
Plot_MOI_Ellip(df, savefig=True)
df[['Ixx', 'Iyy', 'Izz']].plot()
df[['d_Ixx', 'd_Iyy', 'd_Izz']].plot()
df[['d_d_Ixx', 'd_d_Iyy', 'd_d_Izz']].plot()

# -----------------------------------------------------------------------------
data_path = '/media/sam/BE48068348063B23/Simulation_Results/201206/'
extern_drive_path = '/media/sam/ASTRO_DATA/201206/'

# Copy every tenth data file from data directory to external drive
for cwd, folders, files in os.walk(data_path):
    for filename in files:
        if 'data' in filename:
            file_n = int(filename.replace('data.', '').replace('.vtk', ''))
            if file_n % 10 == 0:
                # Only write if file not in folder
                if not os.path.exists(extern_drive_path+filename):
                    sh.copyfile(data_path+filename, extern_drive_path+filename)
