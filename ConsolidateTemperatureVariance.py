import cantrips as can
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import datetime

def getPlotColors():
    channel_colors = {'Ref Resistors':'k', 'x-axis':'r','y-axis':'g', 'z-axis':'b', 'V+':'violet', 'V-':'orange', 'Ref Therm/Resist':'cyan'}
    return channel_colors

def readInDataFromTempVarFile(file, data_dir, date_str):
    data = can.readInColumnsToList(file, file_dir = data_dir, delimiter = ',', verbose = 0)
    header = [col[0] for col in data]
    data_cols = [col[1:] for col in data]
    time = [int(elem) for elem in header[0].split(' ')[-1].split(':')]
    start_time_elems = [int(elem) for elem in date_str.split('_')] + time #Y, M, d, h, m, s
    start_time = datetime.datetime(*start_time_elems).timestamp()
    obs_times = [float(elem) + start_time for elem in data_cols[0]]
    data_dict = {header[i].strip():[float(elem) for elem in data_cols[i]] for i in range(1, len(header))}
    data_dict['t'] = obs_times

    return data_dict

if __name__ == "__main__":
    args = sys.argv[1:]
    add_dates = 1
    labelsize = 16
    ticklabelsize = 14
    figsize = (10, 5)
    date_strs = []
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/TempProbe/data/'
    save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/TempProbe/results/'
    file_prefix = 'TempVar_Arduino_7Channel'

    while add_dates:
        date_str = can.getUserInputWithDefault('Enter the next date_str that will be added (just hit [RETURN] when you are done): ', '-1')
        if date_str == '-1':
            add_dates = 0

        elif os.path.exists(data_dir + date_str + '/'):
            date_strs = date_strs + [date_str]
        else:
            print ('I could not find the directory corresponding to that date (should look like): ')
            print (data_dir + date_str + '/')
            print ('Double check formatting? ')
    if len(date_strs) == 0:
        print ('No usable dates requested.  We will not do anything.' )
        sys.exit()
    print ('I am going to read in all temperature variation data from the following directories:')
    for date_str in date_strs:
        print (data_dir + date_str + '/')
    data_files = [[] for date_str in date_strs]
    for i in range(len(date_strs)):
        date_dir = data_dir + date_strs[i] + '/'
        date_files = [f for f in os.listdir(date_dir) if file_prefix in f]
        date_indeces = [int(f.split('_')[-1].split('.')[0]) for f in date_files]
        date_files, date_indeces = can.safeSortOneListByAnother(date_indeces, [date_files, date_indeces])
        print ('Here are all the Temperature Variation files for the date ' + str(date_str) + ': ')
        print (date_files)
        first_index = int(can.getUserInputWithDefault('Which is the first index from this directory that I should use (default ' + str(date_indeces[0]) + ') : ', str(date_indeces[0])))
        if first_index in date_indeces:
            first_index = date_indeces.index(first_index)
        else:
            print ('Index ' + str(first_index) + ' not an valid index for these data files. Going with default.')
            first_index = date_indeces[0]
        last_index = int(can.getUserInputWithDefault('Which is the last index from this directory that I should use (default ' + str(date_indeces[-1]) + ') : ', str(date_indeces[-1])))
        if last_index in date_indeces:
            last_index = date_indeces.index(last_index)
        else:
            print ('Index ' + str(last_index) + ' not an valid index for these data files. Going with default.')
            last_index = date_indeces[-1]
        date_files, date_indeces = [date_files[first_index:last_index + 1], date_indeces[first_index:last_index + 1]]
        still_indeces_to_ignore = 1
        indeces_to_ignore = []
        while still_indeces_to_ignore:
            index_to_ignore = (can.getUserInputWithDefault('Enter indeces that I should ignore from this date, one at a time (just hit [RETURN] when done): ', '-1'))
            if index_to_ignore.isdigit():
                index_to_ignore = int(index_to_ignore)
            if index_to_ignore in date_indeces:
                indeces_to_ignore = indeces_to_ignore + [index_to_ignore]
            elif str(index_to_ignore) == '-1':
                still_indeces_to_ignore = 0
            else:
                print ('Index ' + str(index_to_ignore) + ' not in list of indeces.  Remember, just the trailing number, no whitespace.')
        for index_to_ignore in indeces_to_ignore:
            index = date_indeces.index(index_to_ignore)
            date_files, date_indeces = [can.removeListElement(date_files, index), can.removeListElement(date_indeces, index)]
        data_files[i] = date_files

    all_data_files = can.flattenListOfLists(data_files)
    print ('Okay.  I am going to plot data from the following Temperature Variation data files: ')
    print (all_data_files)
    all_data_dicts = can.flattenListOfLists([[readInDataFromTempVarFile(data_file, data_dir + date_strs[i] + '/', date_strs[i]) for data_file in data_files[i]] for i in range(len(date_strs))])
    shared_header_elems = []
    for key in all_data_dicts[0]:
        if np.all([key in data_dict for data_dict in all_data_dicts]):
            shared_header_elems = shared_header_elems + [key]
    full_data_dict = {}
    for key in shared_header_elems:
        full_data_dict[key] = can.flattenListOfLists([data_dict[key] for data_dict in all_data_dicts])

    save_data_file_root = (can.getUserInputWithDefault('What is the root of the file to save (not including file type; default: TempVariance_' + 'and'.join(date_strs) + '): ', 'TempVariance_' + 'and'.join(date_strs)))
    save_data_file = save_data_file_root + '.txt'
    save_data_plot = save_data_file_root + '.pdf'

    can.saveListsToColumns([full_data_dict[key] for key in shared_header_elems], save_data_file, save_dir, header = shared_header_elems)
    print ('full_data_dict.keys() = ' + str(full_data_dict.keys()))
    Delta_ts = np.array(full_data_dict['t']) - np.min(full_data_dict['t'])
    xs = full_data_dict['x-axis TVar (K^2)']
    ys = full_data_dict['y-axis TVar (K^2)']
    zs = full_data_dict['z-axis TVar (K^2)']
    ref_rr = full_data_dict['Ref Resistors TVar (K^2)']
    xs_err = full_data_dict['x-axis TVarErr (K^2)']
    ys_err = full_data_dict['y-axis TVarErr (K^2)']
    zs_err = full_data_dict['z-axis TVarErr (K^2)']
    ref_rr_err = full_data_dict['Ref Resistors TVarErr (K^2)']
    legend_strs = []
    scats = []
    f, ax = plt.subplots(1,1, figsize = figsize)
    x_scat = ax.scatter(Delta_ts, xs, c = getPlotColors()['x-axis'])
    #ax.errorbar(Delta_ts, xs, yerr = xs_err, fmt = 'none', ecolor = getPlotColors()['x-axis'])
    y_scat = ax.scatter(Delta_ts, ys, c = getPlotColors()['y-axis'])
    #ax.errorbar(Delta_ts, ys, yerr = ys_err, fmt = 'none', ecolor = getPlotColors()['y-axis'])
    z_scat = ax.scatter(Delta_ts, zs, c = getPlotColors()['z-axis'])
    #ax.errorbar(Delta_ts, zs, yerr = zs_err, fmt = 'none', ecolor = getPlotColors()['z-axis'])
    ref_rr_scat = ax.scatter(Delta_ts, ref_rr, c = getPlotColors()['Ref Resistors'])
    #ax.errorbar(Delta_ts, ref_rr, yerr = ref_rr_err, fmt = 'none', ecolor = getPlotColors()['Ref Resistors'])

    ax.set_xlabel(r'$\Delta t$ (s)', fontsize= labelsize)
    ax.set_ylabel(r'Temp variance over 10 s (K$^2$)', fontsize= labelsize)
    ax.set_yscale('log')
    plt.tick_params(axis='both', which='major', labelsize=ticklabelsize)
    plt.tick_params(axis='both', which='minor', labelsize=ticklabelsize)
    ax.legend([x_scat, y_scat, z_scat, ref_rr_scat], ['x-axis', 'y-axis', 'z-axis', 'Ref Resistors'])
    plt.tight_layout()
    plt.savefig(save_dir + save_data_file_root + '.pdf')
    plt.show()


    print ('Done. ')
