#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 14:54:37 2021

@author: samfrederick
"""
import pandas as pd

def TimeDerivative(df, column, time_column):
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

    try:
        deriv_df.iloc[0, 0] = head_deriv
    except ValueError as e:
        print('head:', e)
        pass

    try:
        deriv_df.iloc[-1, 0] = tail_deriv
    except ValueError as e:
        print('tail:', e)
        pass

    try:
        # Remove issue at t = 0.45 with duplicated index
        central_deriv = central_deriv[~central_deriv.index.duplicated()]
        deriv_df.iloc[1:-1, 0] = central_deriv[:]
    except ValueError as e:
        print('body:', e)
        pass

    # Join derivative column with passed dataframe
    df = df.set_index(idx_name, drop=True)
    df['d_'+column] = deriv_df['d' + column]

    return df