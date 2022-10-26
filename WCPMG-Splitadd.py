# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 10:33:41 2022

@author: edwu5ea1
Idea to circumvent FID addition if the total echo length in points is an odd
number: Slice the FID, FFT and then add the resulting spectra
"""
import csv
import numpy as np
import nmrglue as ng
from matplotlib import pyplot as plt

##### User Input
name = r'Test'
path = r'C:\Users\HB\data_work\600MHz SC\nmr\93Nb-SiLiNb\58\pdata\1'
number_of_echoes_to_add = 20
first_echo = 1
#####

dic, rawdata = ng.bruker.read(path)
rawdata = ng.bruker.remove_digital_filter(dic, rawdata)
data = rawdata

##### Phasing of the FID for debugging
# data = ng.proc_autophase.autops(rawdata, "acme", p0=0, p1=0)
# data = np.abs(rawdata)

##### Get parameters from dictionary
num_echoes = dic['acqus']['L'][22]
pulse_length = dic['acqus']['P'][1]
echo_length = dic['acqus']['D'][6] * 1e6
blank_length = dic['acqus']['D'][3] * 1e6
spectral_width = dic['acqus']['SW_h']
dwell_time = np.round(1/(2*spectral_width) * 1e6, decimals=4)

##### Calculate delays and pulse durations in points
pulse_length_points = (pulse_length / dwell_time) / 2
blank_length_points = (blank_length / dwell_time) / 2
echo_length_points = (echo_length/dwell_time) / 2
blank_block = pulse_length_points + 2 * blank_length_points
echo_to_echo = blank_block + echo_length_points

##### Initialize an empty figure
figure = fig = plt.figure()

##### Initialize fid variable to store the added FIDs
fid = np.array([0])

for i in range(number_of_echoes_to_add):

    echo_start = (int(blank_block * (i+first_echo) 
                      + echo_length_points * (i+first_echo-1)))
    echo_end = int(echo_start + echo_length_points)
    # Plot all individual FIDs for debugging
    # plt.plot(rawdata[echo_start:echo_end])
    fid = fid + data[echo_start:echo_end]

##### Zero filling and Fourier transformation
fid = ng.proc_base.zf_size(fid, 18*1024)
spec = ng.proc_base.fft(fid)

##### Creates ppm axis with correct referencing for Bruker data

spec = ng.proc_base.rev(spec)
length = spec.shape[data.ndim - 1]
carrier_freq = dic['acqus']['SFO1']
spectral_width = dic['acqus']['SW_h']
relative_offset_frequency = (dic['acqus']['O1']- (dic['procs']['SF']
                                                  - dic['acqus']['BF1']) * 1e6)

freq_scale = (
    (np.arange((spectral_width/2) + relative_offset_frequency,
               - spectral_width/2 + relative_offset_frequency,
               - spectral_width/length)))

# Create ppm axis
ppm_scale = freq_scale/carrier_freq

# Transmitter offset frequency in points
# Selecting data size out of array
sfo_point = (int(np.round(length/2)+ (relative_offset_frequency / 
                                      (spectral_width/length)
                         )))

##### Show spikelett spectrum for debugging
spec_compare = ng.proc_base.zf_size(rawdata, 18*1024)
spec_compare = ng.proc_base.fft(spec_compare)
spec_compare = ng.proc_base.rev(spec_compare)
plt.plot(ppm_scale, np.abs(spec_compare))

##### Plot spectrum from added FIDs
plt.plot(ppm_scale, np.abs(spec))
plt.gca().invert_xaxis()

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

export(freq_scale, abs(spec), '##freq', np.round(dic['procs']['SF'], decimals=5))










