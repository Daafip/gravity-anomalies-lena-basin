# import warnings
# warnings.simplefilter("ignore")
import os
os.environ['USE_PYGEOS'] = '0'

import glob
import re
import time
import regionmask
import datetime
# to fix MJD
import astropy
from astropy.time import Time

# data management
import numpy as np
import pandas as pd
import xarray as xr

# imported from provided py code: long functions
from calc_geoid_change_alt import readstrokescoefficients, plm
from multiprocessing import Pool

def run_initialisation_raw():
    grace_files = glob.glob('Data\\Grace\*.gfc')
    love_numbers_kl = np.loadtxt('Data\\loadLoveNumbers_60.txt')[:,1]
    l = 60
    m = 60

    fname_stoke_coeff1 = f'Data\\Grace\\degree1_stokes_coeff.txt'
    df_stokes_1 = pd.read_csv(fname_stoke_coeff1, skiprows=116, delimiter=" ",
                              names=['GRCOF2', "_1", "_2", "_3", "l", "_4", "_5", "m", "Clm", "Slm", "sd_Clm", "sd_Slm",
                                     "begin_date", "end_date"])
    ### drop unwanted
    df_stokes_1.drop(columns=["_1", "_2", "_3", "_4", "_5", "GRCOF2", "sd_Clm", "sd_Slm"], inplace=True)
    ### Reformat dates
    df_stokes_1["end_date"] = df_stokes_1.apply(
        lambda x: pd.Timestamp(f'{str(x.end_date)[0:4]}-{str(x.end_date)[4:6]}-{str(x.end_date)[6:8]}'), axis=1)
    ### dirty fix to make more fit
    df_stokes_1["begin_date"] = df_stokes_1.apply(
        lambda x: pd.Timestamp(f'{str(x.begin_date)[0:4]}-{str(x.begin_date)[4:6]}-01'), axis=1)

    df_1_0 = df_stokes_1[df_stokes_1['m'] == 0].set_index('begin_date')
    df_1_1 = df_stokes_1[df_stokes_1['m'] == 1].set_index('begin_date')

    # easiest way to fix the indexing in matchting the correct indexes
    df_index_replace = df_1_1.index.to_numpy().copy()
    for i, index in enumerate(df_1_1.index):
        if i == 0:
            pass
        if (df_1_1.index[i - 1] - index) == pd.Timedelta(0):
            replace = pd.Timestamp(f'{index.year}-{index.month + 1}-{index.day}')
            df_index_replace[i] = replace
    df_1_1.index = df_index_replace
    df_1_0.index = df_index_replace

    fname_stoke_coeff2 = f'Data\\Grace\\degree2_stokes_coeff.txt'
    col_names = ['MJD begin', "Year fraction begin", "C20", "C20 - C20_mean (1.0E-10)", "sig_C20 (1.0E-10)",
                 "C30", "C30 - C30_mean (1.0E-10)", "sig_C30 (1.0E-10)", 'MJD end', "Year fraction end"]
    df_stokes_2_3 = pd.read_csv(fname_stoke_coeff2, skiprows=37, delimiter="\s+", names=col_names)
    # fix date format
    df_stokes_2_3["begin_date"] = df_stokes_2_3.apply(lambda x: MJD_to_ts(x['MJD begin']), axis=1)
    df_stokes_2_3["end_date"] = df_stokes_2_3.apply(lambda x: MJD_to_ts(x['MJD begin']), axis=1)
    df_stokes_2_3 = df_stokes_2_3[["begin_date", "C20", "C30", "end_date"]].set_index("begin_date")

    # allign indexes and replace like in C_1_1..
    df_stokes_2_3 = df_stokes_2_3.iloc[:-2]  # remove last two months to make same length
    df_stokes_2_3.index = df_index_replace

    # add trend for Glacial isostatic Adjustment (GIA)
    fname_gia = f'Data\\GIA\\GIA_stoke_coeff_trend.gz'
    df_1 = pd.read_csv(fname_gia, compression='gzip', skiprows=1, nrows=6, names=["l", "m", "Clm", "Slm"],
                       delimiter="\s+")
    df_2 = pd.read_csv(fname_gia, compression='gzip', skiprows=9, nrows=6, names=["l", "m", "Clm", "Slm"],
                       delimiter="\s+")
    df_3 = pd.read_csv(fname_gia, compression='gzip', skiprows=17, names=["l", "m", "Clm", "Slm"], delimiter="\s+")

    df_combined = pd.concat([df_2, df_3])
    df_combined["l"] = df_combined["l"].astype(int)
    df_combined["m"] = df_combined["m"].astype(int)

    GIA_C = np.zeros((l + 1, m + 1))
    GIA_S = np.zeros((l + 1, m + 1))
    for index, row in df_combined.iterrows():
        if row.m <= m and row.l <= l:
            GIA_C[int(row.m), int(row.l)] = row["Clm"]
            GIA_S[int(row.m), int(row.l)] = row["Slm"]

    #############################################################
    # Names and dates needed, obtained from file name
    grace_names = [file[-11:-4] for file in grace_files]
    times = [pd.Timestamp(grace_names_i) for grace_names_i in grace_names]
    # load in all coefficient
    C, S, = [], []
    for i, file in enumerate(grace_files):
        C_i, S_i, R, GM = readstrokescoefficients(file)

        # replace C_1_0,C_1_1, S_1_1, C_2_0, C_3_0
        try:
            test = df_1_0.loc[times[i], "Clm"]
            new_time = times[i]
        except KeyError:  # issue with finding correct value, this is easiest fix
            new_time = pd.Timestamp(f'{times[i].year}-{times[i].month - 1}-{times[i].day}')
            # print(new_time)

        C_i[1, 0] = df_1_0.loc[new_time, "Clm"]
        C_i[1, 1] = df_1_1.loc[new_time, "Clm"]
        S_i[1, 1] = df_1_1.loc[new_time, "Slm"]
        C_i[2, 0] = df_stokes_2_3.loc[new_time, "C20"]
        C_30 = df_stokes_2_3.loc[new_time, "C30"]
        if np.isnan(C_30):
            pass
        else:
            C_i[3, 0] = C_30

        # add glacia isostatic adjustment
        C_i += GIA_C
        S_i += GIA_S

        C.append(C_i), S.append(S_i)
    #############################################################
    # calculate means of coefficients
    C = np.array(C)
    S = np.array(S)
    C_mean = C.sum(axis=0) / len(C)
    S_mean = S.sum(axis=0) / len(S)
    # remove mean coefficients
    dc1_store = []
    ds1_store = []
    for i, c in enumerate(C):
        dc1 = C[i] - C_mean
        ds1 = S[i] - S_mean
        dc1_store.append(dc1)
        ds1_store.append(ds1)

    ##############################################################
    rho_av__rho_w = 5.5  # aprox
    # convert to dh
    dh_c1_store = []
    dh_s1_store = []

    for k, dc1 in enumerate(dc1_store):
        C_dhw_i_1 = np.zeros((l + 1, l + 1))
        S_dhw_i_1 = np.zeros((l + 1, l + 1))
        for i in range(l + 1):
            for j in range(i + 1):
                multiplication_factor = (R * (2 * i + 1) * rho_av__rho_w) / (3 * (1 + love_numbers_kl[i]))
                C_dhw_i_1[i, j] = (dc1_store[k][i, j] * multiplication_factor)
                S_dhw_i_1[i, j] = (ds1_store[k][i, j] * multiplication_factor)
        dh_c1_store.append(C_dhw_i_1)
        dh_s1_store.append(S_dhw_i_1)

    _lambda = np.pi / 180 * np.arange(270, 330, 1) - np.pi  # 90 - 150 # deg lon
    theta = np.pi - np.pi / 180 * np.arange(180 - 40, 180 - 10, 1)  # 80 - 50 deg lat

    return _lambda, theta, dh_c1_store, dh_s1_store, times

def plm(theta, degree):
    p_lm = np.zeros((degree + 1, degree + 1))
    u = np.sqrt(1 - np.cos(theta)**2)
    for l in range(degree +1):
        if l == 0:
            p_lm[l, l] = 1
        elif l == 1:
            p_lm[l, l] = u * np.sqrt(3)
        elif l > 1:
            p_lm[l, l] = u * np.sqrt((2*l + 1) / (2*l)) * p_lm[l-1, l-1]

    for l in range(degree + 1):
        for m in range(0, l+1):
            if m == l-1:
                if m == 0:
                    delta = 1
                else:
                    delta = 0
                a_lm = (2/np.sqrt(1+delta)) * ((m+1)/np.sqrt((l-m)*(l+m+1)))
                p_lm[l, m] = a_lm * (np.cos(theta)/u) * p_lm[l, m+1]

    for l in range(degree + 1):
        for m in range(l-2, 0 - 1, -1):
            if m == 0:
                delta = 1
            else:
                delta = 0
            a_lm = (2/np.sqrt(1+delta)) * ((m+1)/np.sqrt((l-m)*(l+m+1)))
            b_lm = (1 / np.sqrt(1 + delta)) * ((np.sqrt((l + m + 2) * (l - m - 1))) / (np.sqrt((l - m) * (l + m + 1))))
            p_lm[l, m] = a_lm * (np.cos(theta)/u) * p_lm[l, m + 1] - b_lm * p_lm[l, m+2]
    return p_lm


def run_loop(input):
    _lambda , theta, dh_c1_store, dh_s1_store = input
    ewh_i_1 = np.zeros((len(theta), len(_lambda)))
    m, l = 60, 60
    for i in range(len(theta)):  # loop over all thetas
        # print(f'{i=} (out of {len(theta)})',end='\r')
        P_lm = plm(theta[i], l)  # all Legendre Functions for one theta
        for j in range(len(_lambda)):  # loop over all lambdas
            for k in range(l + 1):  # loop over all degrees
                for t in range(k + 1):  # loop over negative orders
                    sin_t_lambda = np.sin(t * _lambda[j])  # negative orders
                    cos_t_lambda = np.cos(t * _lambda[j])  # non-negative orders
                    # compute here equivalent water heights
                    ewh_i_1[i, j] = ewh_i_1[i, j] + (dh_s1_store[k, t] * P_lm[k, t] * sin_t_lambda)
                    ewh_i_1[i, j] = ewh_i_1[i, j] + (dh_c1_store[k, t] * P_lm[k, t] * cos_t_lambda)

    return ewh_i_1

def MJD_to_ts(mjd):
    ## thanks to https://stackoverflow.com/questions/72597699/how-do-i-convert-mjdjulian-date-to-utc-using-astropy
    # Start with some time in modified julian date (MJD)
    # Convert to Julian Date
    mjd = float(mjd)
    jd = mjd + 2400000.5
    # Convert to astropy Time object
    t = astropy.time.Time(jd, format='jd')
    # Convert to datetime
    str = t.to_datetime()
    return str


# for monthly data differences
if __name__ == "__main__":
    #### last full run  514.2109310626984 second - ~10 min
    #### last run  240. second - ~4 min, n=80
    print("start")
    # first way, using multiprocessing
    _lambda, theta, dh_c1_store, dh_s1_store, times = run_initialisation_raw()
    debug = False
    if debug:
        n = 5
        args = [[_lambda, theta, dh_c1_store[i], dh_s1_store[i]]
                for i in range(len(dh_s1_store[:n]))]
        times = times[:(n)]  # debug
        fname = "Data\\anomally_waterhead_raw_test_multicore.nc" # debug
    else:
        args = [[_lambda, theta, dh_c1_store[i], dh_s1_store[i]]
                for i in range(len(dh_s1_store))]
        fname = "Data\\anomally_waterhead_raw_multicore.nc"
    print("initialisation done")
    start_time = time.time()
    print(datetime.datetime.now())
    with Pool() as pool:
        store_ewh = pool.map(run_loop, args)

    finish_time = time.time()
    print("Program finished in {} min - using multiprocessing".format((finish_time - start_time)/60))

    ds = xr.DataArray(store_ewh, dims=("time", "lat", "lon"), coords={"lon": _lambda / np.pi * 180,
                                                                      "lat": 90 - theta/np.pi*180,
                                                                      "time": times}, name="dh(m)")


    ds.to_netcdf(fname)
