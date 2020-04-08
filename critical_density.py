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
    09 Apr. 2020 

updated the bugs in calculation the upper-word coefficient. Now it gives similar results as Shirley's results
updated interpolation for any given temperature


Next to do: 

1. Be more automatic in reading collision files. Now we are still specifying each of them.

2. recognise ortho- and para- parameters automatically.

3. try more complicated molecular transitions, not be limited to rotational transitions.

'''


import numpy as np
from scipy.interpolate import interp1d
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
    elif molecule == 'HCO+': 
        T_row_skip  = 80
        Cp_row_beg  = 82 
        Cp_row_read = 465 
        Co_row_beg  = 82 
        Co_row_read = 465 
        hd_row_beg  = 41 
        hd_row_read = 30 
        filename    = 'hco+.dat'

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
    Eu=np.insert(Eu,0,0)
    return T_list, Cp, Co, A, freq, Eu 



def g(j): 
    g=2*j+1
    return g


def Cij(J_low, J_up, T, T_list, Cr): 
    i    = np.where( (Cr[:,2]  == J_low+1) & (Cr[:,1] == J_up+1))[0]
    x    = T_list
    y    = Cr[i,3:] 
    f    = interp1d(x, y)
    cij  = f(T)
    return cij




def gamma(J_low, J_up,T, T_list, Cr, Eu): 
    cij     = Cij(J_low, J_up, T, T_list, Cr) 
    Eu_up = Eu[J_up]
    Eu_low= Eu[J_low]
    gamma_ij= g(J_up) / g(J_low)  * cij * np.exp( -1. * (Eu_up - Eu_low) /  T)
    return gamma_ij 

def Sigma_gamma(J_low,T, T_list, Cr, Eu):
    S_gamma=[]
    # all transitions that go upwards by collisions. 
    for J_up in  np.arange(J_low+2, 24):
        S_gamma.append(gamma(J_low+1,J_up,T,T_list,Cr,Eu))
    Sigma_gamma=np.sum(S_gamma)
#   print(S_gamma)
    return Sigma_gamma

def Sigma_gamma_pre(J_low,T, T_list, Cr, Eu,):
    i     = np.where(Cr[:,1]     == J_low+2)[0]
    # all transitions that go downside by collisions. 
    k     = np.where(T_list      == T      )[0]+3  
    cij   = Cr[i,k]
    S   = np.sum(cij)
    print(S)
    return S


def ncrit(molecule, J_low, T, op_ratio, verbose): 
    T_list, Cp, Co, A, freq, Eu = input(molecule) 
    ncrit_o = A[J_low] / (Sigma_gamma(J_low,T,T_list,Co,Eu,) + Sigma_gamma_pre(J_low,T,T_list,Co,Eu,))
    ncrit_p = A[J_low] / (Sigma_gamma(J_low,T,T_list,Cp,Eu,) + Sigma_gamma_pre(J_low,T,T_list,Cp,Eu,))
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


