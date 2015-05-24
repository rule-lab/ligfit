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
import scipy  # Used for fitting algorithms
import mpmath
from mpmath import mp, mpf  # Used for arbitrary precision floats
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
            total_ligand.append(mpf(vals[0]))
            y_obs.append(mpf(vals[1]))
            y_err.append(mpf(vals[2]))
            weight.append(mpf(vals[3]))
    data = [numpy.array(i,dtype=object) for i in [total_ligand, y_obs, y_err, weight]]
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

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    a = mpf('2') * kd / alpha + mpf('2') * l_total - p_total
    b = (kd + mpf('2') * l_total - mpf('2') * p_total) * (kd / alpha)
    c = (mpf('-1') * numpy.power(kd, mpf('2')) * p_total) / alpha
    return a, b, c


def calc_qr(a, b, c):
    '''
    Calculates Q and R (functions of a, b, and c) for the cubic solution. The
    value of Q^3 + R^2 will determine if the solution needs to be calculated in
    cartesian or polar coordinates.

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    q = (mpf('3') * b - numpy.power(a, mpf('2'))) / mpf('9')
    r = (mpf('9') * a * b - mpf('27') * c - mpf('2') * numpy.power(a, mpf('3'))) / mpf('54')
    return q, r


def cartesian_cubic(a, q, r):
    '''
    The cartesian coordinate solution of the cubic function; for use when
    Q^3+R^2 > 0

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    #Calculate the three terms of the cubic
    first = a / mpf('-3')
    second = numpy.power(r + numpy.power(numpy.power(q, mpf('3')) + numpy.power(r, mpf('2')), mpf('0.5')), mpf('1/3'))
    third = numpy.power(r - numpy.power(numpy.power(q, mpf('3')) + numpy.power(r, mpf('2')), mpf('0.5')), mpf('1/3'))
    return first + second + third


def polar_cubic(a, q, r):
    '''
    The polar coordinate solution of the cubic function; for use when
    Q^3+R^2 < 0

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    theta = mpmath.acos(r / numpy.power(mpf('-1') * numpy.power(q, mpf('3')), mpf('0.5')))
    return mpmath.cos(theta / mpf('3')) * numpy.power(mpf('-1') * q, mpf('0.5')) * mpf('2') - a / mpf('3')


def get_pl(kd, alpha, l_total, p):
    '''
    After solving for [P]-free, this function will return the concentration of
    monomer protein bound to ligand: [PL].

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    numerator = mpf('2') * kd * l_total * p
    denominator = numpy.power(kd, mpf('2')) + mpf('2') * kd * p + alpha * numpy.power(p, mpf('2'))
    return numerator / denominator

def get_plp(kd, alpha, l_total, p):
    '''
    After solving for [P]-free, this function will return the concentration of
    dimer protein bound to ligand: [PLP].

    This version formulated to work with mpmath for arbitrary floating point
    precision and numpy arrays.
    '''
    numerator = alpha * l_total * numpy.power(p, mpf('2'))
    denominator = numpy.power(kd, mpf('2')) + mpf('2') * kd * p + alpha * numpy.power(p, mpf('2'))
    return numerator / denominator

def model_func(kd, alpha, p_total, l_total):
    a, b, c = calc_abc(kd, alpha, p_total, l_total)
    q, r = calc_qr(a, b, c)
    p = []
    for a_val, q_val, r_val in zip(a, q, r):
        if numpy.power(q_val, mpf('3')) + numpy.power(r_val, mpf('2')) > mpf('0'):
            p.append(cartesian_cubic(a_val, q_val, r_val))
        else:
            p.append(polar_cubic(a_val, q_val, r_val))
    p = numpy.array(p)
    pl = get_pl(kd, alpha, l_total, p)
    plp = get_plp(kd, alpha, l_total, p)
    return p, pl, plp


def fit():
    #Read in the information from the input file
    parameters, data = parse_input_file(arguments['<input>'])
    #Unpack the data list into the components
    total_ligand, y_obs, y_err, weight = data
    print(total_ligand)

    #Unless working with units of concentration in the Y-observed (not my
    #present use-case) the Y-observed must be scaled and normalized in a way
    #that allows the model's predicted Y values to be relevant. This
    #consequently adds another dimension for fitting.

    if arguments['--lsq']:
        from scipy.optimize import curve_fit
        #First argument must be the independent argument, the others will be
        #fitted
        p_total = mpf('0.1')

        def func(l_total, kd, alpha, scaling):
            kd = mpf(kd)
            alpha = mpf(alpha)
            scaling = mpf(scaling)
            p, pl, plp = model_func(kd, alpha, p_total, total_ligand)
            return scaling * plp
        popt, pcov = curve_fit(f=func, xdata=total_ligand, ydata=y_obs)
        print(popt)



    #kd = mpf('0.02')
    #alpha = mpf('0.05')
    #p_total = mpf('0.1')
    #ligand_range = mpf('0.00001') * numpy.power(mpf('10'), numpy.linspace(mpf('1'), mpf('8'), mpf('71')))

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

