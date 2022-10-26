import os
from subprocess import run
from cpmg_tools import processing
import numpy as np
import lmfit
from collections import defaultdict
import sys


def constant_factory(value):
    """This function allows the use of custom values for the initialization of defaultdicts
    Use:
    dict = defaultdict(lambda: defaultdict(constant_factory(1.0)))
    to initialize each new entry as 1.0.

    Args:
        value (int, float, str): Default value for new entries in the defaultdict

    Returns:
        Given value
    """
    return lambda: value


def create_simpson(output_path, output_name, input_dict=None, proc_dict=None, params_scaling={}):  # Write simpson input files
    """This generates a custom Simpson inputfile, which can be used from the terminal with: 'simpson <output_name>'

    Args:
        output_path (str): Path to save the generated file
        output_name (str): Name of generated inputfile
        input_dict (dict): Used to specify simulation parameter. Not used parameter use the following defaults:
            "nuclei":'1H',
            "cs_iso":0.0,
            "csa":0.0,
            "csa_eta":0.0,
            "alpha":0.0,
            "beta":0.0,
            "gamma":0.0,
            "pulse":'',
            "spin_rate":0.0,
            "proton_frequency":500.0e6,
            "start_operator":'Inx',
            "detect_operator":'Inp',
            "crystal_file":'rep100',
            "gamma_angles":10,
            "method":'direct',
            "sw":2e6,
            "np":16384,
            "lb":1000,
            "lb_ratio":1.0,
            "si":'np*2',
            "scaling_factor":1.0
    """

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    simpson_variables = {}
    if input_dict is not None:
        if any(isinstance(i,dict) for i in input_dict.values()):
            for keys in input_dict:
                simpson_variables[keys] = {}

            for keys in input_dict:
                simpson_variables[keys] = {
                    "nuclei":'1H',
                    "cs_iso":0.0,
                    "csa":0.0,
                    "csa_eta":0.0,
                    "alpha":0.0,
                    "beta":0.0,
                    "gamma":0.0,
                    "pulse":'',
                    "spin_rate":0.0,
                    "proton_frequency":500.0e6,
                    "start_operator":'Inx',
                    "detect_operator":'Inp',
                    "crystal_file":'rep100',
                    "gamma_angles":10,
                    "method":'direct',
                    "sw":2e6,
                    "np":16384,
                    "lb":1000,
                    "lb_ratio":1.0,
                    "si":'np*2',
                    "scaling_factor":1.0,
                    "output_name":output_name}
            for keys in input_dict:
                for subkeys in input_dict[keys]:
                    simpson_variables[keys][subkeys] = input_dict[keys][subkeys]
                    if keys in params_scaling:
                        if subkeys in params_scaling[keys]:
                            simpson_variables[keys][subkeys] = input_dict[keys][subkeys]*params_scaling[keys][subkeys]

        else:
            simpson_variables = {
                "nuclei":'1H',
                "cs_iso":0.0,
                "csa":0.0,
                "csa_eta":0.0,
                "alpha":0.0,
                "beta":0.0,
                "gamma":0.0,
                "pulse":'',
                "spin_rate":0.0,
                "proton_frequency":500.0e6,
                "start_operator":'Inx',
                "detect_operator":'Inp',
                "crystal_file":'rep100',
                "gamma_angles":10,
                "method":'direct',
                "sw":2e6,
                "np":16384,
                "lb":1000,
                "lb_ratio":1.0,
                "si":'np*2',
                "scaling_factor":1.0,
                "output_name":output_name}
            for keys in input_dict:
                if keys in params_scaling:
                    simpson_variables[keys] = input_dict[keys]*params_scaling[keys]
                else:
                    simpson_variables[keys] = input_dict[keys]
    else:
        simpson_variables = {
                    "nuclei":'1H',
                    "cs_iso":0.0,
                    "csa":0.0,
                    "csa_eta":0.0,
                    "alpha":0.0,
                    "beta":0.0,
                    "gamma":0.0,
                    "pulse":'',
                    "spin_rate":0.0,
                    "proton_frequency":500.0e6,
                    "start_operator":'Inx',
                    "detect_operator":'Inp',
                    "crystal_file":'rep100',
                    "gamma_angles":10,
                    "method":'direct',
                    "sw":2e6,
                    "np":16384,
                    "lb":1000,
                    "lb_ratio":1.0,
                    "si":'np*2',
                    "scaling_factor":1.0,
                    "output_name":output_name}

    # TODO: REMOVE INDENTATION OF MULTILINE STRINGS IN FINAL VERSION
    simpson = {'spinsys': '', 'par': '', 'pulseq': '', 'main': ''}

    simpson['spinsys'] = """
    spinsys {{
        channels {nuclei}
        nuclei {nuclei}
        shift 1 {cs_iso}p {csa}p {csa_eta} {alpha} {beta} {gamma}
    }}
    """

    simpson['par'] = """
    par {{
        spin_rate        {spin_rate}
        proton_frequency {proton_frequency}
        start_operator   {start_operator}
        detect_operator  {detect_operator}
        method           {method}
        crystal_file     {crystal_file}
        gamma_angles     {gamma_angles}
        sw               {sw}
        variable tsw     1e6/sw
        verbose          0000
        np               {np}
        variable si      {si}
    }}
    """

    simpson['pulseq'] = """
    proc pulseq {{}} {{
        global par
        {pulse}
        acq_block {{
            delay $par(tsw)
        }}
    }}
    """

    simpson['main'] = """
    proc main {{}} {{
        global par

        set f [fsimpson]

        faddlb $f {lb} {lb_ratio}

        fsave $f {output_name}.xy -xreim
    }}
    """

    if proc_dict is not None:
        for keys in proc_dict:
            simpson[keys] = proc_dict[keys]

    if any(isinstance(i,dict) for i in simpson_variables.values()):
        counter = 0
        for keys in simpson_variables:
            with  open(output_path + output_name,'w') as myfile:
                myfile.write(''.join(simpson.values()).format(**simpson_variables[keys]))

            os.chdir(output_path)
            run(['simpson', output_name])
            data, timescale = processing.read_ascii_fid(output_path + output_name + '.xy')
            if counter == 0:
                data = data*simpson_variables[keys]['scaling_factor']
            else:
                data = data*simpson_variables[keys]['scaling_factor'] + data_temp
            data_temp = data
            counter += 1
    else:
        with  open(output_path + output_name,'w') as myfile:
                myfile.write(''.join(simpson.values()).format(**simpson_variables))
        os.chdir(output_path)
        run(['simpson', output_name])
        data, timescale = processing.read_ascii_fid(output_path + output_name + '.xy')
    
    return (data, timescale)


def create_simpson_fit(output_path, output_name, input_dict=None, proc_dict=None, params_scaling={}):  # Write simpson input files
    """This generates a custom Simpson inputfile, which can be used from the terminal with: 'simpson <output_name>'

    Args:
        output_path (str): Path to save the generated file
        output_name (str): Name of generated inputfile
        input_dict (dict): Used to specify simulation parameter. Not used parameter use the following defaults:
            "nuclei":'1H',
            "cs_iso":0.0,
            "csa":0.0,
            "csa_eta":0.0,
            "alpha":0.0,
            "beta":0.0,
            "gamma":0.0,
            "pulse":'',
            "spin_rate":0.0,
            "proton_frequency":500.0e6,
            "start_operator":'Inx',
            "detect_operator":'Inp',
            "crystal_file":'rep100',
            "gamma_angles":10,
            "method":'direct',
            "sw":2e6,
            "np":16384,
            "lb":1000,
            "lb_ratio":1.0,
            "si":'np*2',
            "scaling_factor":1.0
    """

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    simpson_variables = {}
    if input_dict is not None:
        simpson_variables = {
            "nuclei":'1H',
            "cs_iso":0.0,
            "csa":0.0,
            "csa_eta":0.0,
            "alpha":0.0,
            "beta":0.0,
            "gamma":0.0,
            "pulse":'',
            "spin_rate":0.0,
            "proton_frequency":500.0e6,
            "start_operator":'Inx',
            "detect_operator":'Inp',
            "crystal_file":'rep100',
            "gamma_angles":10,
            "method":'direct',
            "sw":2e6,
            "np":16384,
            "lb":1000,
            "lb_ratio":1.0,
            "si":'np*2',
            "scaling_factor":1.0,
            "output_name":output_name}
        for keys in input_dict:
            if keys in params_scaling:
                simpson_variables[keys] = input_dict[keys]*params_scaling[keys]
            else:
                simpson_variables[keys] = input_dict[keys]
    else:
        simpson_variables = {
                    "nuclei":'1H',
                    "cs_iso":0.0,
                    "csa":0.0,
                    "csa_eta":0.0,
                    "alpha":0.0,
                    "beta":0.0,
                    "gamma":0.0,
                    "pulse":'',
                    "spin_rate":0.0,
                    "proton_frequency":500.0e6,
                    "start_operator":'Inx',
                    "detect_operator":'Inp',
                    "crystal_file":'rep100',
                    "gamma_angles":10,
                    "method":'direct',
                    "sw":2e6,
                    "np":16384,
                    "lb":1000,
                    "lb_ratio":1.0,
                    "si":'np*2',
                    "scaling_factor":1.0,
                    "output_name":output_name}

    # TODO: REMOVE INDENTATION OF MULTILINE STRINGS IN FINAL VERSION
    simpson = {'spinsys': '', 'par': '', 'pulseq': '', 'main': '', 'rms': ''}

    simpson['spinsys'] = """
    spinsys {{
        channels {channels}
        nuclei {nuclei}
        shift 1 -218p 0p 0.0 0 0 0
        shift 2 128p 0p 0.0 0 0 0
    }}
    """

    simpson['par'] = """
    par {{
        spin_rate        {spin_rate}
        proton_frequency {proton_frequency}
        start_operator   {start_operator}
        detect_operator  {detect_operator}
        method           {method}
        crystal_file     {crystal_file}
        gamma_angles     {gamma_angles}
        sw               {sw}
        variable tsw     1e6/sw
        verbose          0000
        np               {np}
        variable si      {si}
    }}
    """

    simpson['pulseq'] = """
    proc pulseq {{}} {{
        global par
        {pulse}
        acq_block {{
            delay $par(tsw)
        }}
    }}
    """

    # TODO: Insert correct load file for comparison
    simpson['main'] = """
    proc main {{}} {{
        global par

        lappend ::auto_path /usr/local/tcl
        lappend ::auto_path ./opt 1.1
        if {{![namespace exists opt]}} {{
            package require opt 1.1
            package require simpson-utils
            namespace import simpson-utils::*
        }}
        set par(verb) 1

        set par(exp) [fload {comp_file}]

        opt::function rms

        opt::newpar csa1 -200 1 -700 -100
        opt::newpar eta1 0.1 0.1 0 0.5
        opt::newpar csa2 -1100 1 -1200 -700
        opt::newpar eta2 0.25 0.1 0 0.5

        opt::minimize 1.0e-5

        rms 1
    }}
    """

    simpson['rms'] = """
    proc rms {{{{save 0}}}} {{
        global par

        set f [fsimpson [list \\
                [list shift_1_aniso ${{opt::csa1}}p] \\
                [list shift_1_eta $opt::eta1] \\
                [list shift_2_aniso ${{opt::csa2}}p] \\
                [list shift_2_eta $opt::eta2] \\
                ]]

        faddlb $f {lb} {lb_ratio}
        fzerofill $f $par(si)
        fft $f

        fautoscale $f $par(exp) -re
        set rms [frms $f $par(exp) -re ]
        if {{$save == 1}} {{
            set fileId [open {output_name}.results "w"]
            puts $fileId [format " \[%s\] %10.3f %10.3f %10.3f %10.3f %10.3f" \\
            FINAL $opt::csa1 $opt::csa2 $opt::eta1 $opt::eta2 $rms]
            fsave $f $par(name)_final.spe
        }}
        flush stdout
        funload $f
        return $rms
    }}
    """

    if proc_dict is not None:
        for keys in proc_dict:
            simpson[keys] = proc_dict[keys]

    with  open(output_path + output_name,'w') as myfile:
            myfile.write(''.join(simpson.values()).format(**simpson_variables))
    os.chdir(output_path)
    run(['simpson', output_name])

    return ()


def run_simpson(input_file, working_dir, *args):
    """This function executes a given simpson inputfile in the given directory. *args can be used in simpson via argc and argv.

    Args:
        input_file (str): Absolute path to the simpson input file
        working_dir (str): Absolute path to the working directory. Results are saved in this directory
    """
    os.chdir(working_dir)
    run(['simpson', input_file, *args])

    return()


def fit_helper(params, data, output_path, output_name, input_dict, proc_dict, params_scaling, si, verb):
    """This function is called by the simpson_fit wrapper. It generates the simpson input file, executes simpson and calculates the residual.

    Args:
        params (tuple() or tuple(tuple())): Parameter to be fitted
        data (ndarray): Comparison data for fitting
        output_path (str): Path to save and execute simpson files
        output_name (str): Name of simpson file
        input_dict (dict or nested dict): Dict to specify simulation parameter
        proc_dict (dict): Dict to specify custom pulseprogs
        params_scaling (dict): Scaling of each parameter inside the fit. Scales them to the range 0 to 1
        si (int): Value for zerofilling of simulation
        verb (bool): Print debugging information
    """
    if input_dict is not None:
        if any(isinstance(i,dict) for i in input_dict.values()):
            for keys in params.valuesdict():
                input_dict[str(keys.rsplit('_', 1)[-1])][str(keys.rsplit('_', 1)[0])] = params.valuesdict()[keys]
        else:
            for keys in params.valuesdict():
                input_dict[keys] = params.valuesdict()[keys]

    if input_dict is not None:
        if any(isinstance(i,dict) for i in input_dict.values()):
            params_scaling_input = defaultdict(lambda: defaultdict(constant_factory(1.0)))
            for keys in params_scaling:
                params_scaling_input[str(keys.rsplit('_', 1)[-1])][str(keys.rsplit('_', 1)[0])] = params_scaling[keys]
        else:
            params_scaling_input = defaultdict(lambda: defaultdict(constant_factory(1.0)))
            for keys in params_scaling:
                params_scaling_input[keys] = params_scaling[keys]

    # print(input_dict)
    data_model, timescale_model  = create_simpson(output_path, output_name, input_dict=input_dict, proc_dict=proc_dict, params_scaling=params_scaling_input)
    if(len(data_model)>si):
        si = len(data_model)
        print('SI value to low. Value was corrected to: ' + str(si))
    data_model_fft, _ = processing.asciifft(data_model, timescale_model, si=si)
    
    # Scaling
    data_model_fft = data_model_fft/max(data_model_fft.real)
    data = data/max(data.real)

    # Check SI compliance
    if(len(data_model_fft)!=len(data)):
        sys.exit('Model and comparison data are not of same length. Check SI value. ' + str(len(data_model_fft)) + ' vs. ' + str(len(data)))

    residual = processing.calc_logcosh(data_model_fft.real, data.real)

    if verb:
        print(residual)
    # print(processing.calc_mse(data_model_fft.real, data.real))
    # return(data_model_fft.real - data.real)
    return(residual)


def fit_simpson(output_path, output_name, params_input, data, si, input_dict=None, proc_dict=None, verb=True, method='powell', **fit_kws):
    """This function is a wrapper to combine the simpson-python pipeline with lmfit. Most methods from lmfit can be used,
       but the easiest method seems to be 'powell'.
       Both the experimental data and the simpson simulation have to use the same number of points.
    
    Args:
        output_path (str): Path to save and execute simpson files
        output_name (str): Name of simpson file
        params_input (tuple() or tuple(tuple())): Parameter to be fitted
        data (ndarray): Comparison data for fitting
        si (int): Value for zerofilling of simulation
        input_dict (dict): Used to specify simulation parameter. Not used parameter use the following defaults:
            "nuclei":'1H',
            "cs_iso":0.0,
            "csa":0.0,
            "csa_eta":0.0,
            "alpha":0.0,
            "beta":0.0,
            "gamma":0.0,
            "pulse":'',
            "spin_rate":0.0,
            "proton_frequency":500.0e6,
            "start_operator":'Inx',
            "detect_operator":'Inp',
            "crystal_file":'rep100',
            "gamma_angles":10,
            "method":'direct',
            "sw":2e6,
            "np":16384,
            "lb":1000,
            "lb_ratio":1.0,
            "si":'np*2',
            "scaling_factor":1.0
        For multiple species these Parameter have to be given as nested dicts following this scheme:
            input_dict{'1':{"cs_iso":10.0,}, '2':{"cs_iso":20.0,}}
        proc_dict (dict, optional): Dict to specify custom pulseprogs

            proc_dict = {}

            proc_dict['spinsys'] = " " " (Remove spaces between " ")

            spinsys {{

                channels 1H

                nuclei 1H 1H

                shift 1 {cs_iso}p {csa}p {csa_eta} {alpha} {beta} {gamma}

                shift 2 50p 50p 0.5 {alpha} {beta} {gamma}

                dipole 1 2 -10000 0 0 0

            }}

            " " " (Remove spaces between " ")
        verb (bool, optional): [description]. Defaults to True.
        method (str, optional): [description]. Defaults to 'powell'.

    Returns:
        lmfit fit report
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.chdir(output_path)

    params = lmfit.Parameters()
    params_scaling = {}
    if any(isinstance(i,tuple) for i in params_input):
        for elem in params_input:
            name, value, vary, valmin, valmax, expr, brute_step = elem
            scalefactor = 1
            if value:
                if(abs(value)>scalefactor):
                    scalefactor = abs(value)
            if valmin != -np.inf and valmin != None:
                if(abs(valmin)>scalefactor):
                    scalefactor = abs(valmin)
            if valmax != np.inf and valmax != None:
                if(abs(valmax)>scalefactor):
                    scalefactor = abs(valmax)
            if value:
                value = value/scalefactor
            if valmin != -np.inf and valmin != None:
                valmin = valmin/scalefactor
            if valmax != np.inf and valmax != None:
                valmax = valmax/scalefactor
            if brute_step:
                brute_step = brute_step/scalefactor
            params.add(*(name, value, vary, valmin, valmax, expr, brute_step))
            params_scaling[name] = scalefactor
    else:
        name, value, vary, valmin, valmax, expr, brute_step = params_input
        # print(name, value, vary, valmin, valmax, expr, brute_step)
        scalefactor = 1
        if value:
            if(abs(value)>scalefactor):
                scalefactor = abs(value)
        if valmin != -np.inf and valmin != None:
            if(abs(valmin)>scalefactor):
                scalefactor = abs(valmin)
        if valmax != np.inf and valmax != None:
            if(abs(valmax)>scalefactor):
                scalefactor = abs(valmax)
        if value:
            value = value/scalefactor
        if valmin != -np.inf and valmin != None:
            valmin = valmin/scalefactor
        if valmax != np.inf and valmax != None:
            valmax = valmax/scalefactor
        if brute_step:
            brute_step = brute_step/scalefactor
        # print(name, value, vary, valmin, valmax, expr, brute_step)
        params.add(*(name, value, vary, valmin, valmax, expr, brute_step))
        params_scaling[name] = scalefactor

    print('Starting Fit.')
    print('This are the parameters you passed for optimization, after scaling:')
    params.pretty_print()
    print('Scale parameters are:')
    print(params_scaling)

    out = lmfit.minimize(fit_helper, params, args=(data, output_path, output_name, input_dict, proc_dict, params_scaling, si, verb), method=method, **fit_kws)

    print('--------------------------------------')
    print('These are the resulting fit parameter:')
    print('--------------------------------------')
    print('Parameter    Value    Min    Max    Stderr')
    for name, param in out.params.items():
        print('{:7s} {:11.5f} {:11.5f} {:11.5f} {}'.format(name, param.value*params_scaling[name], param.min*params_scaling[name], param.max*params_scaling[name], param.stderr))

    return(out)
