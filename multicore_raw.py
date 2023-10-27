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

def run_initialisation():
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

    # pd.read_csv didnt work so must read 'manually'
    reach_end_of_head = False
    list_lines = []
    with open(fname_stoke_coeff2) as fin:
        for line in fin:
            line = line.strip()
            if not reach_end_of_head:
                if line.startswith("Product:"):
                    reach_end_of_head = True
            else:
                line = line.split(' ')
                line = np.array(line)
                line = line[line != " "]
                line = line[line != ""]
                list_lines.append(line)
    arr_lines = np.array(list_lines)

    # convert array into df
    df_stokes_2_3 = pd.DataFrame(arr_lines,
                                 columns=['MJD begin', "Year fraction begin", "C20", "C20 - C20_mean (1.0E-10)",
                                          "sig_C20 (1.0E-10)",
                                          "C30", "C30 - C30_mean (1.0E-10)", "sig_C30 (1.0E-10)",
                                          'MJD end', "Year fraction end"
                                          ])
    # fix date format
    df_stokes_2_3["begin_date"] = df_stokes_2_3.apply(lambda x: MJD_to_ts(x['MJD begin']), axis=1)
    df_stokes_2_3["end_date"] = df_stokes_2_3.apply(lambda x: MJD_to_ts(x['MJD begin']), axis=1)
    df_stokes_2_3 = df_stokes_2_3[["begin_date", "C20", "C30", "end_date"]].set_index("begin_date")

    # allign indexes and replace like in C_1_1..
    df_stokes_2_3 = df_stokes_2_3.iloc[:-2]  # remove last two months to make same length
    df_stokes_2_3.index = df_index_replace

    # Names and dates needed, obtained from file name
    grace_names = [file[-11:-4] for file in grace_files]
    times = [pd.Timestamp(grace_names_i) for grace_names_i in grace_names]
    #############################################################
    # load in all coefficient
    C, S = [], []
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
        if C_30 == "NaN":
            pass
        else:
            C_i[3, 0] = C_30

        C.append(C_i), S.append(S_i)
    #############################################################
    # get differences for each pair
    dc1_store = []
    dc2_store = []
    ds1_store = []
    ds2_store = []
    for i, c in enumerate(C):
        if i == 0:
            pass  # first has no prev
        else:
            dc1 = C[i] - ((C[i] + C[i - 1]) / 2)
            dc2 = C[i - 1] - ((C[i - 1] + C[i]) / 2)
            ds1 = S[i] - ((S[i] + S[i - 1]) / 2)
            ds2 = S[i - 1] - ((S[i - 1] + S[i]) / 2)
            dc1_store.append(dc1), dc2_store.append(dc2)
            ds1_store.append(ds1), ds2_store.append(ds2)

    ##############################################################
    rho_av__rho_w = 5.5  # aprox
    # convert to dh
    dh_c1_store = []
    dh_c2_store = []
    dh_s1_store = []
    dh_s2_store = []

    for k, dc1 in enumerate(dc1_store):
        C_dhw_i_1 = np.zeros((l + 1, l + 1))
        S_dhw_i_1 = np.zeros((l + 1, l + 1))
        C_dhw_i_2 = np.zeros((l + 1, l + 1))
        S_dhw_i_2 = np.zeros((l + 1, l + 1))
        for i in range(l + 1):
            for j in range(i + 1):
                multiplication_factor = (R * (2 * i + 1) * rho_av__rho_w) / (3 * (1 + love_numbers_kl[i]))
                C_dhw_i_1[i, j] = (dc1_store[k][i, j] * multiplication_factor)
                S_dhw_i_1[i, j] = (ds1_store[k][i, j] * multiplication_factor)
                C_dhw_i_2[i, j] = (dc2_store[k][i, j] * multiplication_factor)
                S_dhw_i_2[i, j] = (ds2_store[k][i, j] * multiplication_factor)
        dh_c1_store.append(C_dhw_i_1), dh_c2_store.append(C_dhw_i_2)
        dh_s1_store.append(S_dhw_i_1), dh_s2_store.append(S_dhw_i_2)

    _lambda = np.pi / 180 * np.arange(270, 330, 1) - np.pi  # 90 - 150 # deg lon
    theta = np.pi - np.pi / 180 * np.arange(180 - 40, 180 - 10, 1)  # 80 - 50 deg lat



    return _lambda, theta, dh_c1_store, dh_c2_store, dh_s1_store, dh_s2_store, times

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
    _lambda , theta, dh_c1_store, dh_c2_store, dh_s1_store, dh_s2_store = input
    ewh_i_1 = np.zeros((len(theta), len(_lambda)))
    ewh_i_2 = np.zeros((len(theta), len(_lambda)))
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
                    ewh_i_2[i, j] = ewh_i_2[i, j] + (dh_s2_store[k, t] * P_lm[k, t] * sin_t_lambda)
                    ewh_i_1[i, j] = ewh_i_1[i, j] + (dh_c1_store[k, t] * P_lm[k, t] * cos_t_lambda)
                    ewh_i_2[i, j] = ewh_i_2[i, j] + (dh_c2_store[k, t] * P_lm[k, t] * cos_t_lambda)

    return ewh_i_1 - ewh_i_2


def readstrokescoefficients(filename):
    with open(filename) as f:
        reach_end_of_head = 0
        for line in f:
            line = line.strip()
            if reach_end_of_head == 0:
                if line.startswith("earth_gravity_constant"):
                    GM = float(line[line.index(" ") + 1:])
                elif "radius" in line:
                    line = re.sub(' +', ' ', line)
                    R = float(line.split(" ")[1])
                elif line.startswith("max_degree"):
                    line = re.sub(' +', ' ', line)
                    max_degree = int(line.split(" ")[1])
                    C = np.zeros((max_degree + 1, max_degree + 1))
                    S = np.zeros((max_degree + 1, max_degree + 1))
                else:
                    if line.startswith("end_of_head"):
                        reach_end_of_head = 1
            else:
                line = re.sub(' +', ' ', line)
                line = line.split()
                L = int(line[1])
                M = int(line[2])
                C[L, M] = float(line[3])
                S[L, M] = float(line[4])
    return C, S, R, GM


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



if __name__ == "__main__":
    #### last full run  514.2109310626984 second - ~10 min
    #### last run  240. second - ~4 min, n=80
    print("start")
    # first way, using multiprocessing
    _lambda, theta, dh_c1_store, dh_c2_store, dh_s1_store, dh_s2_store, times = run_initialisation()
    debug = False
    if debug:
        n = 10
        args = [[_lambda, theta, dh_c1_store[i], dh_c2_store[i], dh_s1_store[i], dh_s2_store[i]]
                for i in range(len(dh_s2_store[:n]))]
        times = times[:(n + 1)]  # debug
        fname = "Data\\anomally_waterhead_raw_test_multicore2.nc" # debug
    else:
        args = [[_lambda, theta, dh_c1_store[i], dh_c2_store[i], dh_s1_store[i], dh_s2_store[i]]
                for i in range(len(dh_s2_store))]
        fname = "Data\\anomally_waterhead_raw_multicore.nc"
    print("initialisation done")
    start_time = time.time()
    print(datetime.datetime.now())
    with Pool() as pool:
        store_ewh = pool.map(run_loop, args)

    finish_time = time.time()
    print("Program finished in {} seconds - using multiprocessing".format(finish_time - start_time))

    ds = xr.DataArray(store_ewh, dims=("time", "lat", "lon"), coords={"lon": _lambda / np.pi * 180,
                                                                      "lat": 90 - theta/np.pi*180,
                                                                      "time": times[1:]}, name="dh(m)")


    ds.to_netcdf(fname)