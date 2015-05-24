#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""LigFit

Usage:
  ligfit.py fit <input> [<output>] (--lsq | --odr) [--precision=<dps>]
  ligfit.py makeinput [<filename>]

Options:
  -l --lsq             fit using least squares (error only in Y)
  -o --odr             fit using orthogonal distance regression, error in all
                       axes
  --precision=<dps>    specify the decimal place precision for calculations
  -h --help            show this help message and exit
  -v --version         show version and exit
  -q --quiet           report only file names

"""

from docopt import docopt  # Used for command-line argument parsing
from scipy.optimize import curve_fit
import numpy  # Used for arrays
import json  # Used to store/retrieve input information
import matplotlib.pyplot as plt

#mpmath.mp.precision

#Configuration of input() for support in both Python 2 and Python 3
try:
    input = raw_input
except NameError:
    pass


def parse_input_file(file_path):
    with open(file_path, 'r') as input_file:
        input_lines = input_file.readlines()
        #Load the JSON string on the second line, parse it
        fitting_params = json.loads(input_lines[1])
        #Parse the space-delimited data
        total_ligand, y_obs, y_err, weight = [], [], [], []
        for data_line in input_lines[4:]:
            vals = data_line.split(' ')
            total_ligand.append(vals[0])
            y_obs.append(vals[1])
            y_err.append(vals[2])
            weight.append(vals[3])
    data = [numpy.array(i, numpy.float64) for i in [total_ligand, y_obs, y_err, weight]]
    return fitting_params, data


def make_input_file():
    #If the filename was not specified at command, ask for it
    if arguments['<filename>'] is None:
        arguments['<filename>'] = input('Name of input file?: ')

    print('''Please enter values in uM units, notation like \"1e3\" is \
valid.''')
    vals = {}
    #Get the protein concentration used, a constant value
    vals['prot_total'] = float(input('Total protein used in experiments?: '))
    #Get initial guesses for K_d and alpha (cooperativity constant)
    vals['init_kd'] = float(input('Initial guess for K_d?: '))
    vals['init_alpha'] = float(input('Inital guess for alpha: '))

    #Open the file
    with open(arguments['<filename>'], 'w') as out:
        out.write('The next line contains a JSON dump of fitting parameters\n')
        out.write(json.dumps(vals))
        out.write('''
The remaining lines contain the raw data for fitting in a space-separated \
format in the following order:
[L]_total (M), Y_obs (RFU or other), Y-error, weighting (0-1)
<Your Data Here>''')


def calc_abc(kd, alpha, p_total, l_total):
    '''
    Calculates the cubic polynomial constants a, b, and c from the given Kd,
    alpha, p_total, and l_total values.
    '''
    a = (2.0 * kd) / alpha + (2.0 * l_total) - p_total
    b = (kd + 2.0 * l_total - 2.0 * p_total) * (kd / alpha)
    c = (-1.0 * numpy.power(kd, 2.0) * p_total) / alpha
    return a, b, c


def calc_qr(a, b, c):
    '''
    Calculates Q and R (functions of a, b, and c) for the cubic solution. The
    value of Q^3 + R^2 will determine if the solution needs to be calculated in
    cartesian or polar coordinates.
    '''
    q = (3.0 * b - numpy.power(a, 2.0)) / 9.0
    r = (9.0 * a * b - 27.0 * c - 2.0 * numpy.power(a, 3.0)) / 54.0
    return q, r


def cartesian_cubic(a, q, r):
    '''
    The cartesian coordinate solution of the cubic function; for use when
    Q^3+R^2 > 0
    '''
    #Calculate the three terms of the cubic
    first = a / -3.0
    second = numpy.power(r + numpy.sqrt(numpy.power(q, 3) + numpy.power(r, 2)), 1.0 / 3.0)
    third = numpy.power(r - numpy.sqrt(numpy.power(q, 3) + numpy.power(r, 2)), 1.0 / 3.0)
    return first + second + third


def polar_cubic(a, q, r):
    '''
    The polar coordinate solution of the cubic function; for use when
    Q^3+R^2 < 0
    '''
    theta = numpy.arccos(r / numpy.sqrt(-1.0 * numpy.power(q, 3.0)))
    val = numpy.cos(theta / 3.0) * numpy.sqrt(-1.0 * q) * 2.0 - a / 3.0
    return val


def get_pl(kd, alpha, l_total, p):
    '''
    After solving for [P]-free, this function will return the concentration of
    monomer protein bound to ligand: [PL].
    '''
    numerator = 2.0 * kd * l_total * p
    denominator = numpy.power(kd, 2.0) + 2.0 * kd * p + alpha * numpy.power(p, 2.0)
    return numerator / denominator

def get_plp(kd, alpha, l_total, p):
    '''
    After solving for [P]-free, this function will return the concentration of
    dimer protein bound to ligand: [PLP].
    '''
    numerator = alpha * l_total * numpy.power(p, 2.0)
    denominator = numpy.power(kd, 2.0) + 2.0 * kd * p + alpha * numpy.power(p, 2.0)
    return numerator / denominator

def model_func(kd, alpha, p_total, l_total):
    a, b, c = calc_abc(kd, alpha, p_total, l_total)
    q, r = calc_qr(a, b, c)
    p = []
    for a_val, q_val, r_val in zip(a, q, r):
        if numpy.power(q_val, 3.0) + numpy.power(r_val, 2.0) > 0:
            p.append(cartesian_cubic(a_val, q_val, r_val))
        else:
            p.append(polar_cubic(a_val, q_val, r_val))
    p = numpy.array(p)
    pl = get_pl(kd, alpha, l_total, p)
    plp = get_plp(kd, alpha, l_total, p)
    print(p)
    return p, pl, plp


def model_fitting(l_total, kd, alpha, p_total):
    a, b, c = calc_abc(kd, alpha, p_total, l_total)
    q, r = calc_qr(a, b, c)
    p = []
    for a_val, q_val, r_val in zip(a, q, r):
        if numpy.power(q_val, 3.0) + numpy.power(r_val, 2.0) > 0:
            p.append(cartesian_cubic(a_val, q_val, r_val))
        else:
            p.append(polar_cubic(a_val, q_val, r_val))
    p = numpy.array(p)
    plp = get_plp(kd, alpha, l_total, p)
    return numpy.nan_to_num(plp)

def fit():
    #Read in the information from the input file
    parameters, data = parse_input_file(arguments['<input>'])
    #Unpack the data list into the components
    total_ligand, y_obs, y_err, weight = data
    print(total_ligand.dtype)

    #Unless working with units of concentration in the Y-observed (not my
    #present use-case) the Y-observed must be scaled and normalized in a way
    #that allows the model's predicted Y values to be relevant. This
    #consequently adds another dimension for fitting.

    #Let's pick some hypothetical values for an M8 experiment, units are in uM
    #kd = 2.5
    #alpha = 1000
    #p_total = 0.1

    if arguments['--lsq']:
        #First argument must be the independent argument, the others will be
        #fitted
        #p_total = 0.1
        y_obs = y_obs / y_obs.max()
        popt, pcov = curve_fit(lambda total_ligand, kd, alpha: model_fitting(total_ligand, kd, alpha, 0.1), total_ligand, y_obs, p0=(800, 100))
        print(popt)

    plt.plot(total_ligand, y_obs, 'x', label='rawdata')
    newx = numpy.logspace(0,5.1)
    plt.plot(newx, model_fitting(newx, popt[0], popt[1], 0.1), label='fit')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.grid()
    plt.show()

    #kd = .02
    #alpha = .05
    #p_total = 0.1
    #ligand_range = 0.00001 * numpy.power(10.0, numpy.linspace(1, 8, 71))
    #p, pl, plp = model_func(kd, alpha, p_total, ligand_range)
    #plt.plot(ligand_range, p / p_total, label='[P]')
    #plt.plot(ligand_range, pl / p_total, label='[PL]')
    #plt.plot(ligand_range, plp / p_total, label='[PLP]')
    #plt.legend(loc='center right')
    #plt.xscale('log')
    #plt.grid()
    #plt.show()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.0.1')
    if arguments['--precision']:
        mp.dps = int(arguments['--precision'])
    if arguments['makeinput']:
        make_input_file()
    elif arguments['fit']:
        fit()
