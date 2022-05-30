#Read in data from Jim McArthur's arduino

import serial
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
import cantrips as can
import os
#from datetime import date
import datetime
import time

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def takeFourierTransform(t_sep, ys, missing_data_edges, data_gap_lengths ):
    ys = list(ys)
    new_ys = ys[missing_data_edges[0]:missing_data_edges[1]]
    for j in range(1, len(missing_data_edges) - 1):
        fill_val = (ys[missing_data_edges[j]] + ys[missing_data_edges[j] + 1]) / 2
        new_ys = new_ys + [fill_val for i in range(int(data_gap_lengths[j-1] / t_sep))]  + ys[missing_data_edges[j]:missing_data_edges[j+1]]
    y_fft = fft(np.array(ys) - np.mean(ys))
    x_fft = fftfreq(len(ys), t_sep)

    return (x_fft[1:len(y_fft)//2]), (y_fft[1:len(y_fft)//2])


def convertVoltagesToRRatios(V_plus, V_minus, V_out, gain = 1, average_V_plus = 1,  average_V_minus = 1, ):
    """
    If V+ and V- are very stable, then user should average over them.

    """
    if average_V_plus:
        V_plus = np.mean(V_plus)
    if average_V_minus:
        V_minus = np.mean(V_minus)
    V_out = V_out / gain
    R_top_div_R_bottom = (V_out - V_plus) / (V_minus - V_out)
    return R_top_div_R_bottom

def convertRRatioToTempDiff(RR, alpha, A1 = 1, A2 = 1):
    """
    Assumes that R is related to T according to:
    R_i(T) = A_i exp(alpha_i T)
    And assumes that all alpha_i's are the same
    """
    DeltaT = 1.0 / alpha * np.log( A2 / A1 * RR )
    return DeltaT

def getPlotColors():
    channel_colors = {'Ref Resistors':'k', 'x-axis':'r','y-axis':'g', 'z-axis':'b', 'V+':'violet', 'V-':'orange', 'Ref Therm/Resist':'cyan'}
    return channel_colors

def getChannelSequenceMapping():
    channel_mapping = {'V+':5, 'V-':6,'Ref Resistors':0, 'x-axis':1, 'y-axis':2, 'z-axis':3, 'Ref Therm/Resist':4}
    return channel_mapping

def computeVarianceOverSectionOfArray(delta_ts, data_to_average, n_points_to_average, n_sig_clip = np.inf):
    n_data_points = len(data_to_average) // n_points_to_average
    print ('n_data_points = ' + str(n_data_points))
    #averaged_data = [can.sigClipMean(data_to_average[i*n_points_to_average:(i+1)*n_points_to_average], sig_clip = n_sig_clip) for i in range(n_data_points)]
    #averaged_stdErrs = [can.sigClipStd(data_to_average[i*n_points_to_average:(i+1)*n_points_to_average], sig_clip = n_sig_clip, standard_error = 1) for i in range(n_data_points)]
    averaged_data = [np.mean(data_to_average[i*n_points_to_average:(i+1)*n_points_to_average]) for i in range(n_data_points)]
    data_variance = [ np.mean( (np.array(data_to_average[i*n_points_to_average:(i+1)*n_points_to_average]) - averaged_data[i]) ** 2.0 ) for i in range(n_data_points) ]
    variance_uncertainty = [np.std( (np.array(data_to_average[i*n_points_to_average:(i+1)*n_points_to_average]) - averaged_data[i]) ** 2.0 ) / np.sqrt(n_points_to_average) for i in range(n_data_points)]
    delta_t_means = [np.mean(delta_ts[i*n_points_to_average:(i+1)*n_points_to_average]) for i in range(n_data_points)]
    return [delta_t_means, data_variance, variance_uncertainty]

def measureTemperatureVariance(all_data, average_time_s, save_root,
                               plot_file_ext = '.pdf',   save_file_ext = '.csv',   save_plot_dir = '',   save_data_dir = '', save_plot = 1, save_data = 1,
                               key_strs = ['Ref Resistors', 'x-axis','y-axis', 'z-axis'],
                               CT_label = 'TempVar_', fig_height = 12, fig_horiz_elem_size = 16, labelsize = 14,
                               legendsize = 10, thermistor_percent_R_per_deg_C = 0.04, gain = 1, delta_t_tolerance = 0.1, update_tally = 0, sep = ', ', start_time = [0, 0, 0]):

    thermistor_exp_term = np.log(1.0 / (1.0 + thermistor_percent_R_per_deg_C))
    delta_ts = all_data[0]
    period = np.mean([delta_ts[i+1] - delta_ts[i] for i in range(len(delta_ts)-1)])
    cols = all_data[1:]
    n_V_cols = len(cols)
    channel_mapping = getChannelSequenceMapping()
    channel_colors = getPlotColors()

    delta_Vs = [col - np.mean(col) for col in cols]
    V_plus, V_minus = [np.mean(cols[channel_mapping['V+']]), np.mean(cols[channel_mapping['V-']])]  #We should read these from data stream, once we have more channels
    RR_cols = [np.array(convertVoltagesToRRatios(V_plus, V_minus, np.array(V_out), gain = gain)) for V_out in cols ]
    DeltaT_cols = [convertRRatioToTempDiff(RR, thermistor_exp_term) for RR in RR_cols ]
    n_T_cols = len(DeltaT_cols)

    n_points_to_average = int(average_time_s / period)
    mid_t_times, T_variance_cols, T_variance_stdErr = [[[] for key in key_strs], [[] for key in key_strs], [[] for key in key_strs]]


    for i in range(len(key_strs)):
        channel_mapping
        col = DeltaT_cols[channel_mapping[key_strs[i]]]
        mid_t_times[i], T_variance_cols[i], T_variance_stdErr[i] = computeVarianceOverSectionOfArray(delta_ts, col, n_points_to_average)
    mid_t_times = mid_t_times[0]

    f, axarr = plt.subplots(1,1, figsize = (fig_horiz_elem_size * 1, fig_height ) )
    VT_plots = []
    legend_strs = []
    if save_plot:
        for to_plot_str in key_strs:
            data_col = channel_mapping[to_plot_str]
            if data_col < n_T_cols:
                VT_plots = VT_plots + [axarr.scatter(mid_t_times, T_variance_cols[data_col] , c = channel_colors[to_plot_str])]
                axarr.errorbar(mid_t_times, T_variance_cols[data_col], yerr = T_variance_stdErr[data_col], fmt = 'none', ecolor = channel_colors[to_plot_str])
                legend_strs = legend_strs + [to_plot_str ]
        axarr.legend(VT_plots, legend_strs, ncol = int(np.ceil(len(VT_plots) / 2)), loc = 'upper center', fontsize = legendsize)
        axarr.set_xlabel(r'Average time since sequence start (s)', fontsize = labelsize)
        axarr.set_ylabel(r'$\Delta T$ Variance (K$^2$)', fontsize = labelsize)

        plt.tight_layout()
        if not(os.path.isdir(save_plot_dir)):
            os.mkdir(save_plot_dir)
        tally_num = getTally(save_plot_dir, update_tally = update_tally)
        axarr.set_yscale('log')
        plt.savefig(save_plot_dir + CT_label + save_root + '_' + str(tally_num) + plot_file_ext)
        #plt.show()
        plt.close('all')
    if save_data:
        col_nums = [channel_mapping[key] for key in key_strs]
        ordered_keys = can.safeSortOneListByAnother(col_nums, [key_strs])[0]
        header = 'Delta t ' + ':'.join([str(elem) for elem in start_time]) + sep + (' TVar (K^2)' + sep).join(ordered_keys) + ' TVar (K^2)' + sep + (' TVarErr (K^2)' + sep).join(ordered_keys) + ' TVarErr (K^2)'
        if not(os.path.isdir(save_data_dir)):
            os.mkdir(save_data_dir)
        tally_num = getTally(save_data_dir, update_tally = update_tally)
        can.saveListsToColumns([mid_t_times] + T_variance_cols + T_variance_stdErr, CT_label + save_root + '_' + str(tally_num) + save_file_ext, save_data_dir, sep = sep, header = header)

    return 1



def plotDataStream(all_data, plot_to_save_root, #freq,
                   average_time_s = None,
                   V_keys_to_plot = ['V+', 'V-', 'Ref Therm/Resist', 'Ref Resistors', 'x-axis','y-axis', 'z-axis'  ],
                   T_keys_to_plot = ['Ref Resistors', 'x-axis','y-axis', 'z-axis'], save_dir = '', plot_file_ext = '.pdf',
                   fft_label = 'fft_', full_stream_label = 'Full_', full_plot_prefix = 'Full_Feed',
                   timestamp_key = 't', smoothing_box = 20, full_plot_alpha = 0.2, fig_height = 12, fig_horiz_elem_size = 16, labelsize = 14,
                   legendsize = 10, thermistor_percent_R_per_deg_C = 0.04, gain = 1, delta_t_tolerance = 0.1, max_n_data_segments_to_plot = 5,
                   update_tally = 1):

    thermistor_exp_term = np.log(1.0 / (1.0 + thermistor_percent_R_per_deg_C))

    delta_ts = all_data[0]
    N = len(delta_ts)
    freq = 1.0 / np.mean([delta_ts[i+1] - delta_ts[i] for i in range(len(delta_ts)-1)])
    period = 1.0 / freq
    cols = all_data[1:]
    n_V_cols = len(cols)
    channel_mapping = getChannelSequenceMapping()
    channel_colors = getPlotColors()
    if average_time_s != None:
        n_t_bins = int(np.ceil((delta_ts[-1] - delta_ts[0]) / average_time_s))
        n_points_to_time_average = int(average_time_s / period )
    else:
        n_t_bins = None
        n_points_to_time_average = None

    delta_Vs = [col - np.mean(col) for col in cols]
    V_plus, V_minus = [np.mean(cols[channel_mapping['V+']]), np.mean(cols[channel_mapping['V-']])]  #We should read these from data stream, once we have more channels
    #print ('[V_plus, V_minus] = ' + str([V_plus, V_minus]))
    RR_cols = [np.array(convertVoltagesToRRatios(V_plus, V_minus, np.array(cols[channel_mapping[key]]), gain = gain)) for key in T_keys_to_plot]
    #print ('RR_cols = ' + str(RR_cols))
    DeltaT_cols = [convertRRatioToTempDiff(RR, thermistor_exp_term) for RR in RR_cols ]
    fft_cols = [takeFourierTransform(1 / freq, list(DeltaT_col - np.mean( DeltaT_col)), [0, len(DeltaT_col)], []) for DeltaT_col in DeltaT_cols]
    n_T_cols = len(DeltaT_cols)

    delta_V_plots = []
    legend_V_strs = []
    legend_T_strs = []
    full_V_plots = []
    fft_plots = []
    delta_T_plots = []
    print ('Starting to plot data')
    f, axarr = plt.subplots(4,1, figsize = (fig_horiz_elem_size * 1, fig_height ))
    for to_plot_str in V_keys_to_plot:
        data_col = channel_mapping[to_plot_str]
        if data_col < n_V_cols:
            [axarr[0].plot(delta_ts, cols[data_col], c = channel_colors[to_plot_str], alpha = full_plot_alpha)[0]]
            full_V_plots = full_V_plots + [axarr[0].plot(delta_ts[smoothing_box:-smoothing_box], smooth(cols[data_col], smoothing_box)[smoothing_box:-smoothing_box], c = channel_colors[to_plot_str])[0]]
            legend_V_strs = legend_V_strs + [to_plot_str + ' (Smooth x' + str(smoothing_box) + ')']
            [axarr[1].plot(delta_ts, delta_Vs[data_col], c = channel_colors[to_plot_str], alpha = full_plot_alpha)[0]]
            delta_V_plots = delta_V_plots + [axarr[1].plot(delta_ts[smoothing_box:-smoothing_box], smooth(delta_Vs[data_col], smoothing_box)[smoothing_box:-smoothing_box], c = channel_colors[to_plot_str])[0]]

    for to_plot_str in T_keys_to_plot:
        data_col = channel_mapping[to_plot_str]
        if data_col < n_T_cols:
            [axarr[2].plot(delta_ts, DeltaT_cols[data_col] - np.mean(DeltaT_cols[data_col] ), c = channel_colors[to_plot_str], alpha = full_plot_alpha)[0]]
            delta_T_plots = delta_T_plots + [axarr[2].plot(delta_ts[smoothing_box:-smoothing_box], smooth(DeltaT_cols[data_col] - np.mean(DeltaT_cols[data_col]), smoothing_box)[smoothing_box:-smoothing_box] , c = channel_colors[to_plot_str])[0]]
            if n_t_bins != None:
                #for i in range(len(n_points_to_time_average) + 1):
                #    start_i = time.time()
                #    t = delta_ts[min(i * n_points_to_time_average, len(delta_ts) - 1)]
                #    axarr[2].axvline(t, c = 'gray', linestyle = '--', alpha = 0.5)
                #    end_i = time.time()
                #    print ('Took ' + str(end_i - start_i) + 's for i = ' + str(i))
                [axarr[2].axvline(delta_ts[min(i * n_points_to_time_average, len(delta_ts) - 1)], c = 'gray', linestyle = '--', alpha = 0.5) for i in range(0, n_t_bins + 1)]
            legend_T_strs = legend_T_strs + [to_plot_str + ' (Smooth x' + str(smoothing_box) + ')']
            fft_plots = fft_plots + [axarr[3].plot(fft_cols[data_col][0][1:], 2.0/N * np.abs(fft_cols[data_col][1][1:]), c = channel_colors[to_plot_str], alpha = full_plot_alpha)[0]]
            #fft_plots = fft_plots + [axarr[3].plot(fft_cols[data_col][0][smoothing_box:-smoothing_box], smooth(2.0/N * np.abs(fft_cols[data_col][1]), smoothing_box)[smoothing_box:-smoothing_box], c = channel_colors[to_plot_str]) [0]]
    axarr[3].set_yscale('log')
    for i in range(4):
        y_lim = axarr[i].get_ylim()
        axarr[i].set_ylim([y_lim[0], y_lim[1] * 1.5])
    axarr[0].legend(full_V_plots, legend_V_strs, ncol = int(np.ceil(len(full_V_plots) / 2)), loc = 'upper center', fontsize = legendsize)
    axarr[1].legend(delta_V_plots, legend_V_strs, ncol = int(np.ceil(len(delta_V_plots) / 2)), loc = 'upper center', fontsize = legendsize)
    axarr[2].legend(delta_T_plots, legend_T_strs, ncol = int(np.ceil(len(delta_T_plots) / 2)), loc = 'upper center', fontsize = legendsize)
    axarr[3].legend(fft_plots, legend_T_strs, ncol = int(np.ceil(len(fft_plots) / 2)), loc = 'upper center', fontsize = legendsize)
    axarr[0].set_xlabel(r'$\Delta t$ (s)', fontsize = labelsize)
    axarr[1].set_xlabel(r'$\Delta t$ (s)', fontsize = labelsize)
    axarr[2].set_xlabel(r'$\Delta t$ (s)', fontsize = labelsize)
    axarr[3].set_xlabel(r'$f$ (Hz)', fontsize = labelsize)
    axarr[0].set_ylabel(r'$V_{Out}$ (V)', fontsize = labelsize)
    axarr[1].set_ylabel(r'$V_{Out} - <V_{Out}>$ (V)', fontsize = labelsize)
    axarr[2].set_ylabel(r'$\Delta T - $ ' +    r'$ <\Delta T>$ (K)', fontsize = labelsize)
    axarr[3].set_ylabel(r'Fourier Amplitude' + '\n' +  'of $\Delta T$ (K)', fontsize = labelsize)
    print ('Trying to put in tight layout')
    plt.tight_layout()
    print ('Done plotting')

    if not(os.path.isdir(save_dir)):
        os.mkdir(save_dir)
    tally_num = getTally(save_dir, update_tally = update_tally)
    plt.savefig(save_dir + full_plot_prefix + plot_to_save_root + '_' + str(tally_num) + plot_file_ext)
    plt.close('all')
    #plt.show()

    return 1

def getTally(data_dir, tally_file_name = 'Tally.txt', update_tally = 1):
     if not(os.path.isfile(data_dir + tally_file_name)):
         with open(data_dir + tally_file_name, 'w') as f:
             f.write('0')
     with open(data_dir + tally_file_name, 'r') as f:
          if update_tally:
              tally_num = int(f.readline()) + 1
          else:
              tally_num = int(f.readline())
     with open(data_dir + tally_file_name, 'w') as f:
         f.write(str(tally_num ))
     return tally_num

def saveDataStreamToFile(all_data, save_root, #freq,
                   V_keys_to_save = ['V+', 'V-', 'Ref Therm/Resist', 'Ref Resistors', 'x-axis','y-axis', 'z-axis'  ],
                   sep = ', ', save_file_ext = '.csv', save_dir = '', update_tally = 1, start_time = [0, 0, 0]):
    channel_mapping = getChannelSequenceMapping()
    col_nums = [channel_mapping[key] for key in V_keys_to_save]
    header = 'Delta t ' +':'.join([str(t) for t in start_time]) + sep + sep.join(can.safeSortOneListByAnother(col_nums, [V_keys_to_save])[0])
    if not(os.path.isdir(save_dir)):
        os.mkdir(save_dir)
    tally_num = getTally(save_dir, update_tally = update_tally)
    can.saveListsToColumns(all_data, save_root + '_' + str(tally_num) + save_file_ext, save_dir, sep = sep, header = header)

    return 1


def ADUToVoltageConversion(bit_vals,
                           ground_bit = 33536.8, bit_to_voltage_slope = 8.03e-5,
                           ground_bit_err = 1, bit_to_voltage_slope_err = 1e-7 ):
    #These are measured and/or known from arduino
    voltage_vals = bit_to_voltage_slope * (bit_vals - ground_bit)
    voltage_errs = np.sqrt( (bit_vals - ground_bit) ** 2.0 * bit_to_voltage_slope ** 2.0
                                      + bit_to_voltage_slope ** 2.0 * ground_bit_err ** 2.0 )
    return voltage_vals, voltage_errs

if __name__=="__main__":
    #Find out port by:
    #   $ ls /tty/USBACM*
    #You may need to make that port publicly accessible:
    #   $ sudo chmod a+rw /dev/ttyACM0
    #
    # Arduino provides 16 bit A/D converter, so full well in bits is 65535
    # Bit to voltage mapping is: 0V => 33536.8 bits, 0.999V => 45988.5 bits, -0.998V => 21111 bits
    # This yields bit to voltage relation of V = (8.03 +/- 0.01) e-5 V / (ADU - 33536.8 +/- 1)
    # \Delta Bits = XXXXX (38411.5 - 26041.5)=> \Delta V = (3.000 - 2.003) volt
    # Ground is 1231 bits
    #Execute as:
    #   $ python3 ReadDataFromJimArduinoCircuit.py "/dev/ttyACM0" 60
    args = sys.argv[1:]
    port = args[0]
    sample_time_s = int(args[1])
    average_time_s = int(args[2])
    save_full = int(args[3])

    save_file_root = 'Arduino_7Channel'

    today = datetime.date.today()
    date_str = today.strftime("%Y_%m_%d")
    now = datetime.datetime.now()
    h, m, s = [now.hour, now.minute, now.second]
    data_save_dir = '../data/' + date_str + '/'
    plot_save_dir = '../plots/' + date_str + '/'
    sample_rate_Hz = 1000
    #The data array consistently has an inconsistency at element 160
    #  So we burn in that many samples before our "real" data stream begins.
    n_burn_in = 0
    n_samples = sample_time_s * sample_rate_Hz + n_burn_in
    baud_rate = 9600 #Higher values appear to increase fidelity, but reduce precision
    timeout = 2
    arduino = serial.Serial(port, baud_rate, timeout = timeout)
    arduino.reset_input_buffer()
    arduino.reset_output_buffer()

    full_well = 2 ** 15 - 1
    full_well_V = 4

    data = ['' for i in range(n_samples)]
    print ('Beginning to sample from arduino')
    for i in range(n_samples):
        data[i] = arduino.readline()
    print ('Done sampling')
    rows = [(line.decode('UTF-8')[0:-2].split(',')) for line in data]
    n_cols = int(np.median([len(row) for row in rows]))
    print ('I see ' + str(n_cols) + ' channels off of arduino.')
    bad_rows = [i for i in range(len(rows)) if len(rows[i]) != n_cols]
    print ('Data points with obvious glitch: ' + str(bad_rows))
    n_burn_in = max([n_burn_in] + [row + 1 for row in bad_rows])
    print ('We will ignore the first ' + str(n_burn_in)+ ' points. ')
    cols = [[ADUToVoltageConversion(int(line.decode('UTF-8')[0:-2].split(',')[i]))[0]  for line in data[n_burn_in+1:]] for i in range(n_cols)]
    #cols = [[int(line.decode('UTF-8')[0:-2].split(',')[i])  for line in data[n_burn_in+1:]] for i in range(n_cols)]
    delta_ts = [i * 1 / sample_rate_Hz for i in range(len(cols[0]))]
    if save_full:
        print ('save_full = ' + str(save_full))
        saveDataStreamToFile([delta_ts] + cols, save_file_root + '_' + date_str,  save_file_ext = '.csv',   sep = ', ', save_dir = data_save_dir, update_tally = 1, start_time = [h, m, s])
        plotDataStream([delta_ts] + cols, save_file_root  + '_' + date_str,  average_time_s = average_time_s, plot_file_ext = '.pdf',   save_dir = plot_save_dir, update_tally = 1)
    measureTemperatureVariance([delta_ts] + cols, average_time_s, save_file_root  + '_' + date_str,  plot_file_ext = '.pdf',   save_file_ext = '.csv',   save_plot_dir = plot_save_dir,   save_data_dir = data_save_dir, save_plot = 1, save_data = 1, update_tally = (0 if save_full else 1), start_time = [h, m, s])

    arduino.close()
    print ('Done')
