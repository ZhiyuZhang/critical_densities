# critical_densities

Calculation for critical densities of molecular lines.

Based on the collisional rate data from LAMDA database

http://home.strw.leidenuniv.nl/~moldata/


Usage:

    from critical_density import *
    ncrit(molecule,J_low,Tkin,o/p_ratio,verbose)

Example:

    from critical_density import *
    ncrit('HCN',3,100,3,True)

Written by Zhi-Yu Zhang

pmozhang@gmail.com

Last update:
    23 Jan. 2018


Next to do:

    Interpolate Cij using values from nearby temperatures
    Currently, only the temperature value from the datafile can be used
