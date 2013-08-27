#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The code has the straightforward aim of producing figures to illustrate the
kinetic model for protein homo-dimerization around a bifunctional ligand. It
will generate a 3D plot of the solution space in normalized units.
"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FormatStrFormatter


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
    #pl = get_pl(kd, alpha, l_total, p)
    plp = get_plp(kd, alpha, l_total, p)
    #All we care about here is PLP
    return plp

def max_plp(kd, p_total):
    '''
    Returns the ligand concentration at which the maximum concentration of PLP
    will be achieved for a given Kd and P_total.
    '''
    return (kd + p_total) / 2.0

#The following values are set for modelling a normalized system
p_total = 1.0
kd = 1.0

#The following values are our independent variables
l_total = numpy.logspace(-2, 3.5, 100)
alpha = numpy.logspace(-2, 4, 100)

#I considered using meshgrid here, but was confounded by the "polar or cubic"
#operation, so I went for a bit less elegance
plp_concentration = []
for alpha_value in alpha:
    plp = model_func(kd, alpha_value, p_total, l_total)
    plp_concentration.append(plp)
plp_soln = numpy.array(plp_concentration)


lm, am = numpy.meshgrid(l_total, alpha)
lm = numpy.log10(lm)
am = numpy.log10(am)

#Create the figure
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(lm, am, plp_soln, rstride=1, cstride=1,
                       cmap=cm.coolwarm, linewidth=0, antialiased=False, alpha=0.6)
#cset = ax.contour(lm, am, plp_soln, zdir='z', offset=0, cmap=cm.winter)
cset = ax.contour(lm, am, plp_soln, zdir='x', offset=-2, cmap=cm.winter)
cset = ax.contour(lm, am, plp_soln, zdir='y', offset=4, cmap=cm.winter)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('Log-Total Ligand [L] (conc.)')
ax.set_ylabel('Log-Cooperativity')
ax.set_zlabel('Dimerized Protein [PLP](conc.)')

plt.show()
