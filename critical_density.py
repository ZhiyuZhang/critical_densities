'''
Calculation for critical density of molecular lines. 

Based on the collisional rate data from Lamda database 

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

'''


import numpy as np
from astropy import constants as const

h=const.h.cgs.value   # 6.626075E-27          #[erg/s]
k=const.k_B.cgs.value # 1.380658E-16          #[erg/K]

def input(molecule, ): 
    if molecule == '12CO':
        T_row_skip  = 100
        Cp_row_beg  = 102 
        Cp_row_read = 820 
        Co_row_beg  = 931 
        Co_row_read = 820 
        hd_row_beg  = 51 
        hd_row_read = 40 
        filename    = '12co.dat' 
    elif molecule == '13CO': 
        T_row_skip  = 100
        Cp_row_beg  = 102 
        Cp_row_read = 820 
        Co_row_beg  = 931 
        Co_row_read = 820 
        hd_row_beg  = 51 
        hd_row_read = 40 
        filename    = '13co.dat'
    elif molecule == 'C18O': 
        T_row_skip  = 100
        Cp_row_beg  = 102 
        Cp_row_read = 820 
        Co_row_beg  = 931 
        Co_row_read = 820 
        hd_row_beg  = 51 
        hd_row_read = 40 
        filename    = 'c18o.dat'
    elif molecule == 'HCN': 
        T_row_skip  = 70
        Cp_row_beg  = 72 
        Cp_row_read = 325
        Co_row_beg  = 72 
        Co_row_read = 325
        hd_row_beg  = 36 
        hd_row_read = 25 
        filename    = 'hcn.dat'
    else:
        print("input error")

    T_list = np.genfromtxt(filename, skip_header = T_row_skip,max_rows    = 1)
#   print(T_list)
#   collisional coefficient for para-H2 
    Cp     = np.genfromtxt(filename, skip_header = Cp_row_beg, max_rows= Cp_row_read)
#   print(Cp.shape)
#   collisional coefficient for ortho-H2 
    Co     = np.genfromtxt(filename, skip_header = Co_row_beg,  max_rows= Co_row_read) 
#   print(Co.shape)
    header = np.genfromtxt(filename, skip_header = hd_row_beg,  max_rows= hd_row_read) 
#   print(header.shape)
    A      = header[:,3]
    freq   = header[:,4] * 1E9 # Hz
    Eu     = header[:,5]
    return T_list, Cp, Co, A, freq, Eu 



def g(j): 
    g=2*j+1
    return g


# Cr = > collisional rate 
def Cr(J_low, T, T_list, Co, Cp, ortho_para): 
    if ortho_para == 'p':
        i    = np.where(Cp[:,2]     == J_low+1)[0]
        k    = np.where(T_list      == T      )[0]+3  
        cij  = Cp[i,k]
        ju   = Cp[i,1]
    elif ortho_para == 'o':
        i    = np.where(Co[:,2]     == J_low+1)[0]
        k    = np.where(T_list      == T      )[0]+3  
        cij  = Co[i,k]
        ju   = Co[i,1]
    else:
        print("error in ortho para setup")
        return 
    return cij


def gamma(J_low,T, T_list, Co, Cp, freq, op): 
    gamma_ij= g(J_low+1) / g(J_low)  * Cr(J_low,T, T_list,Co,Cp,op) * np.exp( -1. * h * freq[J_low] / (k * T))
    return gamma_ij 


def Sigma_gamma_pre(J_low,T, T_list, Co, Cp, freq, op):
        i     = np.where(Cp[:,1]     == J_low+1)[0]
        k     = np.where(T_list      == T      )[0]+3  
        cij   = Cp[i,k]
        print(i)
        print("           ")
        print(k)
        print("           ")
        print(cij)
        S     = np.sum(cij)
#       S     = 0
        return S


def ncrit(molecule, J_low, T, op_ratio, verbose): 
    T_list, Cp, Co, A, freq, Eu = input(molecule) 
    ncrit_o = A[J_low] / (np.sum(gamma(J_low,T,T_list,Co,Cp,freq,'o')) + Sigma_gamma_pre(J_low,T,T_list,Co,Cp,freq,'o'))
    ncrit_p = A[J_low] / (np.sum(gamma(J_low,T,T_list,Co,Cp,freq,'p')) + Sigma_gamma_pre(J_low,T,T_list,Co,Cp,freq,'p'))
    ncrit   = (op_ratio * ncrit_o + 1 * ncrit_p )/(op_ratio+1)  
    trans     = "J="+str(J_low+1)+"-"+str(J_low)
    transition = trans.replace(" ", "")
    if verbose: 
        print(molecule , transition)
        print("Tkin = ", T, "K")
        print("ortho-para ratio            :","%.2E" % op_ratio,       )
        print("collisional partner ortho-H2:","%.2E" % ncrit_o , "cm-3")
        print("collisional partner para-H2 :","%.2E" % ncrit_p , "cm-3")
        print("Final critical density      :","%.2E" % ncrit   , "cm-3")
    return(ncrit)


