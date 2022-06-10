import numpy as np
import pandas as pd
import csv
import wget
from os import mkdir
from os.path import basename, splitext, isdir, exists
from datetime import date

data_dir = 'original_data/'
file_link = 'http://www.iceandclimate.nbi.ku.dk/data/NGRIP_d18O_and_dust_5cm.xls'
table_link = 'https://www.iceandclimate.nbi.ku.dk/data/Rasmussen_et_al_2014_QSR_Table_2.xlsx'
file_name = basename(file_link)
table_name = basename(table_link)
data_file = data_dir + file_name


def download():
    if not exists(data_dir + '/accessed.txt'):
        print('downloading data from:', file_link)
        opath = data_dir + '/' + file_name
        if not isdir(data_dir):
            mkdir(data_dir)
        r = wget.download(file_link,
                          out=data_dir + '/' + file_name)
        p = wget.download(table_link,
                          out=data_dir + '/' + table_name)

        today = date.today().strftime('%Y/%m/%d')
        with open(data_dir + '/accessed.txt', 'w') as comment:
            comment.write('data downloaded from ' + file_link + ' on ' + today)

        print('data stored into:', data_dir)

    else:
        accessed = np.genfromtxt(data_dir + '/accessed.txt',
                                 dtype='str',
                                 delimiter=',')
        print(accessed)
        print('reading existing data')

    colnames = ['depth', 'd18o', 'dust', 'age', 'age_err']
    data = pd.read_excel(data_file,
                         sheet_name='NGRIP-2 d18O and Dust',
                         header=None,
                         skiprows=[0],
                         usecols=[0, 1, 2, 3, 4],
                         names=colnames)
    return data


def binning(t_ax, data, bins=None):
    '''
    INPUT
    -----
    data:= data of a time series
    t_ax:= time axis of a time series
    bins:= array like, equally spaced. 
        centers of bins will define new t_axis


    Output

    binned_data:= the i-th entry of binned_data is the 
        the mean of all data points which are located 
        between the i-th and i+1th elements of bins on the 
        t_ax. 
    binned_t_ax := center value of the bins
    '''

    if bins is None:
        res = np.max(np.diff(t_ax))
        bins = np.arange(t_ax[0]-res // 2,
                         t_ax[-1]+res,
                         res)

    binned_data = np.array([np.mean(data[(t_ax > bins[i]) &
                                         (t_ax < bins[i+1])])
                            for i in range(len(bins)-1)])

    binned_t_ax = bins[:-1]+0.5*(bins[1]-bins[0])

    idx = np.argwhere(np.isnan(binned_data))
    diff = np.diff(idx.flatten())
    mask = (np.append(diff, 0) == 1) + (np.append(0, diff) == 1)

    if any(np.diff(idx.flatten()) == 1):
        print('adjacent empty data')

    single_skip = idx[~mask]

    for i in single_skip:
        binned_data[i] = (binned_data[i-1] + binned_data[i+1]) / 2

    return binned_t_ax, np.array(binned_data)


def NGRIP_stadial_mask(age):
    '''
    this function takes an age axis as an input, that is compatibles 
    with the NGRIP GICC05modelext chronology in terms of covered time. 

    it returen a mask of the same length as the input age axis, which 
    is true, where the age of input corresponds a Greenland stadial. 
    '''

    # load dataset on GI onsets and GS onsets
    stratigraphic = pd.read_excel(data_dir + table_name,
                                  header=None,
                                  skiprows=range(23),
                                  names=['event', 'age'],
                                  usecols='A,C')

    # create a mask from the above data which is True for all
    # GI events
    stadials = np.array(['GS' in x for x in stratigraphic['event']])

    # from that mask, a mask is derived that selects only transitions
    # between GI and GS from the stratigraphic data set (the set includes
    # minor GI events that follow major events => GI to GI to GS for example)
    # if transitions[i] is True, there is a transition between the i-th and the
    # i+1-th event from the stratigraphic data set. Since time goes in the
    # opposite direction than age, the age corresponding to the i-th event
    # is the correct age of the transition (that is the point in time where
    # the preceeding phase ended).

    transitions = np.append(np.array(stadials[:-1]) != np.array(stadials[1:]),
                            False)
    transition_ages = stratigraphic['age'][transitions].values

    max_age = np.max(age)
    min_age = np.min(age)

    start_idx = 0
    while transition_ages[start_idx] < min_age:
        start_idx += 1

    end_idx = len(transition_ages)-1
    while transition_ages[end_idx] > max_age:
        end_idx -= 1

    if stadials[start_idx]:
        GS = age < transition_ages[start_idx]

        for i in range(start_idx + 1, end_idx, 2):
            GS_mask = ((transition_ages[i] < age)
                       & (age < transition_ages[i+1]))
            GS = GS | GS_mask
    else:
        GS = np.full(len(age), False)
        for i in range(start_idx, end_idx, 2):
            GS_mask = ((transition_ages[i] < age)
                       & (age < transition_ages[i+1]))
            GS = GS | GS_mask

    return GS, transition_ages[start_idx: end_idx+1]
