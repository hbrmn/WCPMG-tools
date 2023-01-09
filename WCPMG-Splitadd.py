# -*- coding: utf-8 -*-
""" Code for the

Created on Tue Oct 25 10:33:41 2022

@author: Henrik Bradtmüller - mail@bradtmueller.net - https://hbrmn.github.io/

"""
import csv
import numpy as np
import nmrglue as ng
from matplotlib import pyplot as plt

##### User Input
name = r'LNS45_echo0-3'
path = r'C:\Users\edwu5ea1\data_work\600MHz SC\nmr\93Nb-SiLiNb\58\pdata\1'
# path = r'C:\Users\edwu5ea1\data_work\Projects\1_SiO2-Li2O-Nb2O5\93Nb_NMR\Static\240MHz\221014-93Nb-LNS45-UFWCPMG.fid'

first_echo = 1 # First acquisition period to consider to the addition
number_of_echoes_to_add = 3 # Number of subsequent acquisition periods to add
vendor = 'bruker' # 'varian' or 'bruker'
zero_fill = 18 * 1024
normalize_spectrum = 'y' # 'y' or anything else
export_data = 'n' # 'y' or anything else

##### Built-in functions

##### Export the spectrum in DMFit format for easy use in Origin and ssNake

def export(xaxis, yaxis, xlab, ylab, result=None):

    exp_var = zip(xaxis, yaxis)
    with open(path + '\\' + name
              + '.txt', 'w') as file:
        writer = csv.writer(file, delimiter=' ', lineterminator='\n')
        writer.writerow(('ti:' + name, '', ''))
        writer.writerow((xlab, ylab, result))
        writer = csv.writer(file, delimiter='\t', lineterminator='\n')
        for word in exp_var:
            writer.writerows([word])

##### Load experimental parameters

if vendor == 'bruker':
    dic, rawdata = ng.bruker.read(path)
    rawdata = ng.bruker.remove_digital_filter(dic, rawdata)

    spectrometer_frequency = dic['procs']['SF']

    pulse_length = dic['acqus']['P'][1]
    echo_length = dic['acqus']['D'][6] * 1e6
    blank_length = dic['acqus']['D'][3] * 1e6
    spectral_width = dic['acqus']['SW_h']
    dwell_time = np.round(1/(2*spectral_width) * 1e6, decimals=4)

    carrier_freq = dic['acqus']['SFO1']
    spectral_width = dic['acqus']['SW_h']
    relative_offset_frequency = (dic['acqus']['O1']
                                 - (dic['procs']['SF']
                                    - dic['acqus']['BF1']) * 1e6)

    ##### Calculate delays and pulse durations in points
    pulse_length_points = (pulse_length / dwell_time) / 2
    blank_length_points = (blank_length / dwell_time) / 2
    echo_length_points  = (echo_length  / dwell_time) / 2
    blank_block  = pulse_length_points + 2 * blank_length_points
    echo_to_echo = blank_block + echo_length_points

    data = rawdata

elif vendor == 'varian':
    dic, rawdata = ng.varian.read(path)

    spectrometer_frequency = float(dic['procpar']['sfrq']['values'][0])

    carrier_freq = np.float64(
        dic['procpar']['reffrq1']['values'][0])
    spectral_width = np.float64(
        dic['procpar']['sw']['values'][0])
    relative_offset_frequency = (
        np.float64(dic['procpar']['sfrq']['values'][0])
        - carrier_freq) * 1e6
    dwell_time = np.round(1/spectral_width * 1e6, decimals=4)

    # rof1 is the delay before receiver unblank
    rof1 = float(dic['procpar']['rof1']['values'][0])
    # rof2 is (probably) the delay ...?
    rof2 = float(dic['procpar']['rof2']['values'][0])
    # rof3 is (probably) the delay after receiver blank
    rof3 = float(dic['procpar']['rof3']['values'][0])

    ##### Blanking and Pulsing block
    loop_length = float(dic['procpar']['tauXcpmg']['values'][0])
    pulse_length = float(dic['procpar']['pwXcpmg']['values'][0])
    pre_pulse_blank = float(dic['procpar']['r2Xcpmg']['values'][0])
    post_pulse_blank = float(dic['procpar']['r3Xcpmg']['values'][0])
    blank_length = (pulse_length + pre_pulse_blank + post_pulse_blank +
                    rof2 + 0.1) + 3.3 # 3.3 unknown

    echo_length = loop_length - (pulse_length + pre_pulse_blank +
                                 post_pulse_blank + rof2
                                 + 0.1) + rof1 + 2 * rof3 # origin of 0.1 unknown

    ##### Calculate delays and pulse durations in points
    blank_length_points = (blank_length / dwell_time)
    echo_length_points = (echo_length/dwell_time)
    blank_block = blank_length_points
    echo_to_echo = blank_block + echo_length_points

    data = rawdata

    ##### Initial zero filling, here achieved by rolling the data,
    ##### to be able to use the same code as for Bruker
    zero_fill_length = np.round(echo_length_points -
                                np.where(rawdata == 0)[0][0] + blank_block)
    data = ng.process.proc_base.roll(data, pts=zero_fill_length, neg=False)

##### Initialize fid variable to store the added FIDs
fid = np.array([0])
spec = np.array([0])

##### Adding of FIDs
for i in range(number_of_echoes_to_add):

    echo_start = (int(blank_block * (i+first_echo)
                      + echo_length_points * (i+first_echo-1)))
    echo_end = int(echo_start + echo_length_points)

    fid = fid + data[echo_start:echo_end]

##### Zero filling and Fourier transformation
fid = ng.proc_base.zf_size(fid, zero_fill)
spec = ng.proc_base.fft(fid)

##### Creates ppm axis with correct referencing for Bruker data
if vendor == 'bruker':
    spec = ng.proc_base.rev(spec)
length = spec.shape[data.ndim - 1]

freq_scale = (
    (np.arange((spectral_width/2) + relative_offset_frequency,
               - spectral_width/2 + relative_offset_frequency,
               - spectral_width/length)))

# Create ppm axis
ppm_scale = freq_scale/carrier_freq

# Transmitter offset frequency in points
# Selecting data size out of array
sfo_point = (int(np.round(length/2)+ (relative_offset_frequency /
                                      (spectral_width/length))))

##### Show spikelett spectrum for debugging
spec_compare = ng.proc_base.zf_size(rawdata, zero_fill)
spec_compare = ng.proc_base.fft(spec_compare)

spec_compare = np.abs(spec_compare)
spec = np.abs(spec)

if vendor == 'bruker':
    spec_compare = ng.proc_base.rev(spec_compare)

if normalize_spectrum == 'y':
    spec_compare = spec_compare/np.max(spec_compare)
    spec = spec/np.max(spec)


##### Plot spectrum from added FIDs
figure = fig = plt.figure()
plt.plot(ppm_scale, spec_compare, color='gray')
plt.plot(ppm_scale, spec, color='b')
plt.gca().invert_xaxis()

if export_data == 'y':
    export(freq_scale, abs(spec), '##freq',
           np.round(spectrometer_frequency, decimals=5))










