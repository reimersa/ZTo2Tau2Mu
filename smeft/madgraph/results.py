#! /usr/bin/env python

from argparse import ArgumentParser
import ROOT
from tdrstyle_all import *
from array import array
from utils import *
from parallelize import *
from printing_utils import *
from collections import defaultdict, OrderedDict
import parse
from tqdm import tqdm
import numpy as np
import sympy as sp
from sympy.abc import x as xvar
from sympy.abc import y as yvar
from scipy.optimize import fsolve, root

import os, sys, math, json
import subprocess
from yaml import safe_load
from limits_r_smp import get_limits_r
import copy


workarea = '/work/areimers/ZTo2Tau2Mu/smeft/madgraph'
plotfolder = os.path.join(workarea, 'plots/fixed_formula/fixed_r/fixed_acc')
filefolder = os.path.join(workarea, 'files/fixed_formula/fixed_r/fixed_acc')
# plotfolder = os.path.join(workarea, 'plots/xcheck')
ensureDirectory(plotfolder)
ensureDirectory(filefolder)


r_sm  = 0.9
br_sm = 0.1739
f_sm = 0.839

# operators =  ['cll2222', 'cll2233', 'cll2332', 'cle2222', 'cle2233', 'cle3322', 'cle2332', 'cee2222', 'cee2233']
operators =  ['cll2233', 'cll2332', 'cle2233', 'cle3322', 'cle2332', 'cee2233']

operatornames_pretty = {
    'cll2233': 'C_{LL}^{2233}',
    'cll2222': 'C_{LL}^{2222}', 
    'cll2332': 'C_{LL}^{2332}', 
    'cee2222': 'C_{RR}^{2222}', 
    'cee2233': 'C_{RR}^{2233}', 
    'cle2222': 'C_{LR}^{2222}', 
    'cle2233': 'C_{LR}^{2233}', 
    'cle3322': 'C_{LR}^{3322}', 
    'cle2332': 'C_{LR}^{2332}',
}

# operatornames_pretty = {
#     'cll2233': 'C_{LL}^{#mu#mu#tau#tau}',
#     'cee2233': 'C_{RR}^{#mu#mu#tau#tau}', 
#     'cle2233': 'C_{LR}^{#mu#mu#tau#tau}', 
#     'cle3322': 'C_{LR}^{#tau#tau#mu#mu}', 
#     'cle2332': 'C_{LR}^{#mu#tau#tau#mu}',
# }

def main():
    ROOT.gROOT.SetBatch(1)

    # always read in the json with coefficients
    with open(os.path.join(workarea, 'coefficients.json'), 'r') as j:
        coefficients = safe_load(j)
    with open(os.path.join(workarea, 'acceptances.json'), 'r') as j:
        acceptances = safe_load(j)


    # r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r(r_sm=r_sm)

        
    # get limits from known formula
    # limits_1d_per_operator = {}
    # limits_1d_per_operator_comb = {}
    as_per_operator = {}
    bs_per_operator = {}
    aincs_per_operator = {}
    bincs_per_operator = {}
    proc_affected_per_operator = {}
    acceptances_int_per_operator = {}
    acceptances_bsm_per_operator = {}
    acceptances_sm_per_operator = {}
    abrs_per_operator = {}
    bbrs_per_operator = {}

    for op in coefficients:
        coefficients_per_proc = coefficients[op]
        
        ### for single operators, only exactly one of the two processes can be affected in SMEFT. Makes the computation much easier and even analytic.
        # Find which process is affected
        proc_affected = find_procs_affected(coefficients_per_proc=coefficients_per_proc)
        if len(proc_affected) == 1:
            proc_affected = proc_affected[0]
            a = coefficients_per_proc[proc_affected]['interference']
            b = coefficients_per_proc[proc_affected]['purebsm']
            if op in acceptances:
                acc_op_int = acceptances[op][proc_affected]['interference'] if 'interference' in acceptances[op][proc_affected] else 0.
                acc_op_bsm = acceptances[op][proc_affected]['purebsm'] if 'purebsm' in acceptances[op][proc_affected] else 0.
                acc_sm = acceptances['sm'][proc_affected]['sm']
            # elif '+' in op:
            #     acc_op_bsm = acceptances[op][proc_affected]['purebsm'] if 'purebsm' in acceptances[op][proc_affected] else 0.
        elif len(proc_affected) == 2:
            continue
        else:
            print '-- operator:', op
            print 'procs affected:', proc_affected
            raise ValueError('Not 1 or 2 process affected, huh?')
        
        if affects_incttmm(coefficients_per_proc=coefficients_per_proc):
            ainc = coefficients_per_proc['incttmm']['interference']
            binc = coefficients_per_proc['incttmm']['purebsm']

        if affects_tmvv(coefficients_per_proc=coefficients_per_proc):
            abr = coefficients_per_proc['tmvv']['interference']
            bbr = coefficients_per_proc['tmvv']['purebsm']
        else:
            abr = bbr = 0.

        # limits_1d_per_operator[op]      = get_limits_at_cl(a=a, b=b, proc_affected=proc_affected, order='interference', r_sm=r_sm, cl=0.95)
        # limits_1d_per_operator_comb[op] = get_limits_at_cl(a=a, b=b, proc_affected=proc_affected, order='comb',         r_sm=r_sm, cl=0.95)
        as_per_operator[op] = a
        bs_per_operator[op] = b
        aincs_per_operator[op] = ainc
        bincs_per_operator[op] = binc
        abrs_per_operator[op] = abr
        bbrs_per_operator[op] = bbr
        proc_affected_per_operator[op] = proc_affected
        if op in acceptances:
            acceptances_int_per_operator[op] = acc_op_int
            acceptances_bsm_per_operator[op] = acc_op_bsm
            acceptances_sm_per_operator[op] = acc_sm
        # elif '+' in op:
        #     acceptances_bsm_per_operator[op] = acc_op_bsm
    
    for op in operators:
        a = as_per_operator[op]
        b = bs_per_operator[op]
        ainc = aincs_per_operator[op]
        binc = bincs_per_operator[op]
        # proc_affected = proc_affected_per_operator[op]
        acc_int = acceptances_int_per_operator[op]
        acc_bsm = acceptances_bsm_per_operator[op]
        acc_sm = acceptances_sm_per_operator[op]
        abr = abrs_per_operator[op]
        bbr = bbrs_per_operator[op]

        limits_r = get_limits_r_at_cl(cl=0.95)
        limits_comb = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0, ax=0., bx=0., ay=a, by=b, c=0., aincx=0., bincx=0., aincy=ainc, bincy=binc, cinc=0., acc_int_x=0., acc_int_y=acc_int, acc_bsm_x=0., acc_bsm_y=acc_bsm, acc_bsm_intpart=0., acc_sm_x=acc_sm, acc_sm_y=acc_sm, abrx=0., abry=abr, bbrx=0., bbry=bbr, limit=limits_r[0])
        limits_int = get_limit_y_for_limit_x_interference_withacc_withf_withbr_allsolutions(x=0., ax=0., ay=a, aincx=0., aincy=ainc, acc_int_x=0., acc_int_y=acc_int, acc_sm_x=acc_sm, acc_sm_y=acc_sm, abrx=0., abry=abr, limit=limits_r[0])
        print '\n=== Operator: %s' % (op)
        # print 'limit interference:', limits_int
        # print 'limit combined:    ', limits_comb

        try:
            # tau BR modification
            x = sp.Symbol('x')
            taubrresult = sp.solve_univariate_inequality((1 + abr*x) > 0, x, relational=False)
            print 'region, in which Gamma(tau->mu) > 0 [interference]:', taubrresult
            print 'region, in which only Gamma(tau->mu) > 0 [interference]:', sp.solve_univariate_inequality((1 + abr*x) > 0, x, relational=False)
            

            # N after correction > 0
            nresult = sp.solve_univariate_inequality(1+a*x > 0, x, relational=False, domain=taubrresult)
            print 'region, in which N > 0 [interference]:', nresult
            print 'region, in which only N > 0 [interference]:', sp.solve_univariate_inequality(1+a*x > 0, x, relational=False)
    
            # acceptance after correction > 0
            accresult = sp.solve_univariate_inequality((1+(a*x)*(acc_sm/acc_int))/(1+a*x) > 0, x, relational=False, domain=nresult)
            print 'region, in which acc-correction > 0 [interference]:', accresult
            print 'region, in which only acc-correction > 0 [interference]:', sp.solve_univariate_inequality((1+(a*x)*(acc_sm/acc_int))/(1+a*x) > 0, x, relational=False)
    
            # f after correction > 0
            fresult = sp.solve_univariate_inequality(1+ainc*x > 0, x, relational=False, domain=accresult)
            print 'region, in which inclusive pp -> ttmm xsec > 0 [interference]:', fresult
            print 'region, in which only inclusive pp -> ttmm xsec > 0 [interference]:', sp.solve_univariate_inequality(1+ainc*x > 0, x, relational=False)
    
            # f after correction > 0
            fbelowoneresult = sp.solve_univariate_inequality((1+a*x)/(1+ainc*x) < 1./f_sm, x, relational=False, domain=fresult)
            print 'region, in which f (after correction) < 1 [interference]:', fbelowoneresult
            print 'region, in which only f (after correction) < 1 [interference]:', sp.solve_univariate_inequality((1+a*x)/(1+ainc*x) < 1./f_sm, x, relational=False)
    
            # R positive in the first place
            rposresult = sp.solve_univariate_inequality((1+a*x)**2/(1+ainc*x) * ((1+(a*x)*(acc_sm/acc_int))/(1+a*x)) / ((1+abr*x)/(1+br_sm*(abr*x)))**2 > 0, x, relational=False, domain=fbelowoneresult)
            print 'region, in which R > 0 [interference]:', rposresult
            print 'region, in which only R > 0 [interference]:', sp.solve_univariate_inequality((1+a*x)**2/(1+ainc*x) * ((1+(a*x)*(acc_sm/acc_int))/(1+a*x)) / ((1+abr*x)/(1+br_sm*(abr*x)))**2 > 0, x, relational=False)
    
            # R below limit
            rlimresult = sp.solve_univariate_inequality((1+a*x)**2/(1+ainc*x) * ((1+(a*x)*(acc_sm/acc_int))/(1+a*x)) / ((1+abr*x)/(1+br_sm*(abr*x)))**2 < limits_r[0]/r_sm, x, relational=False, domain=rposresult)
            print 'region, in which R < limit [interference]:', rlimresult
            print 'region, in which only R < limit [interference]:', sp.solve_univariate_inequality((1+a*x)**2/(1+ainc*x) * ((1+(a*x)*(acc_sm/acc_int))/(1+a*x)) / ((1+abr*x)/(1+br_sm*(abr*x)))**2 < limits_r[0]/r_sm, x, relational=False)
    
            # Final result
            allallowedregions = rlimresult
            print 'all allowed regions [interference]:', allallowedregions
        except:
            pass

        try:
            # tau BR modification
            x = sp.Symbol('x')
            taubrresult = sp.solve_univariate_inequality((1 + abr*x + bbr*x**2) > 0, x, relational=False)
            # print 'region, in which Gamma(tau->mu) > 0 [combined]:', taubrresult

            # N after correction > 0
            nresult = sp.solve_univariate_inequality(1+a*x+b*x**2 > 0, x, relational=False, domain=taubrresult)
            # print 'region, in which N > 0 [combined]:', nresult
    
            # acceptance after correction > 0
            acc_int_safe = acc_int if acc_int != 0 else -1.
            accresult = sp.solve_univariate_inequality((1+(a*x)*(acc_sm/acc_int_safe)+(b*x**2)*(acc_sm/acc_bsm))/(1+a*x+b*x**2) > 0, x, relational=False, domain=nresult)
            # print 'region, in which acc > 0 [combined]:', accresult
    
            # f after correction > 0
            fresult = sp.solve_univariate_inequality((1+ainc*x+binc*x**2) > 0, x, relational=False, domain=accresult)
            # print 'region, in which inclusive pp -> ttmm xsec > 0 [combined]:', fresult
    
            # f after correction > 0
            fbelowoneresult = sp.solve_univariate_inequality((1+a*x+b*x**2)/(1+ainc*x+binc*x**2) < 1./f_sm, x, relational=False, domain=fresult)
            # print 'region, in which f (after correction) < 1 [interference]:', fbelowoneresult
    
            # R positive in the first place
            rposresult = sp.solve_univariate_inequality((1+a*x+b*x**2)**2/(1+ainc*x+binc*x**2) * ((1+(a*x)*(acc_sm/acc_int_safe)+(b*x**2)*(acc_sm/acc_bsm))/(1+a*x+b*x**2)) / ((1 + abr*x + bbr*x**2)/(1+br_sm*(abr*x + bbr*x**2)))**2 > 0, x, relational=False, domain=fbelowoneresult)
            # print 'region, in which R > 0 [combined]:', rposresult
    
            # R below limit
            rlimresult = sp.solve_univariate_inequality((1+a*x+b*x**2)**2/(1+ainc*x+binc*x**2) * ((1+(a*x)*(acc_sm/acc_int_safe)+(b*x**2)*(acc_sm/acc_bsm))/(1+a*x+b*x**2)) / ((1 + abr*x + bbr*x**2)/(1+br_sm*(abr*x + bbr*x**2)))**2 < limits_r[0]/r_sm, x, relational=False, domain=rposresult)
            # print 'region, in which R < limit [combined]:', rlimresult
    
            # Final result
            allallowedregions = rlimresult
            print 'all allowed regions [combined]:', allallowedregions
        except:
            pass















        # plot_r_vs_c(a=a, b=0., op=op, plotname='R_interference_vs_%s.pdf' % (op), acc_int=acc_int, acc_bsm=acc_sm, acc_sm=acc_sm, ainc=ainc, binc=0., abr=abr, bbr=0., options='afb')
    #     plot_r_vs_c(a=a, b=b, op=op, plotname='R_comb_vs_%s.pdf' % (op), acc_int=acc_int, acc_bsm=acc_bsm, acc_sm=acc_sm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options='afb')


        # for order in ['interference', 'comb']:
        #     # if order != 'interference': continue
        #     l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options='')
        #     plot_limits_vs_cl(op=op, order=order, l_obs=l_obs, l_exp=l_exp, l_68_low=l_68_low, l_68_high=l_68_high, l_95_low=l_95_low, l_95_high=l_95_high, plottag='')

        #     l_withacc_obs, l_withacc_exp, l_withacc_68_low, l_withacc_68_high, l_withacc_95_low, l_withacc_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options='a')
        #     plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_obs, l_exp=l_withacc_exp, l_68_low=l_withacc_68_low, l_68_high=l_withacc_68_high, l_95_low=l_withacc_95_low, l_95_high=l_withacc_95_high, plottag='withacc')

        #     l_withf_obs, l_withf_exp, l_withf_68_low, l_withf_68_high, l_withf_95_low, l_withf_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=ainc, binc=binc, options='f')
        #     plot_limits_vs_cl(op=op, order=order, l_obs=l_withf_obs, l_exp=l_withf_exp, l_68_low=l_withf_68_low, l_68_high=l_withf_68_high, l_95_low=l_withf_95_low, l_95_high=l_withf_95_high, plottag='withf')

        #     l_withacc_withf_obs, l_withacc_withf_exp, l_withacc_withf_68_low, l_withacc_withf_68_high, l_withacc_withf_95_low, l_withacc_withf_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options='af')
        #     plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_withf_obs, l_exp=l_withacc_withf_exp, l_68_low=l_withacc_withf_68_low, l_68_high=l_withacc_withf_68_high, l_95_low=l_withacc_withf_95_low, l_95_high=l_withacc_withf_95_high, plottag='withacc_withf')

        #     l_withacc_withf_withbr_obs, l_withacc_withf_withbr_exp, l_withacc_withf_withbr_68_low, l_withacc_withf_withbr_68_high, l_withacc_withf_withbr_95_low, l_withacc_withf_withbr_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options='afb')
        #     plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_withf_withbr_obs, l_exp=l_withacc_withf_withbr_exp, l_68_low=l_withacc_withf_withbr_68_low, l_68_high=l_withacc_withf_withbr_68_high, l_95_low=l_withacc_withf_withbr_95_low, l_95_high=l_withacc_withf_withbr_95_high, plottag='withacc_withf_withbr')

    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233')

    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')

    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')

    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    # plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')


    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233')

    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')

    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')

    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    # plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')













    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')



    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')

    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')




    
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cee2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')

    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2332', opy='cll2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cll2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
    # plot_limits_2d_both(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2332', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')







    # plot_limits_2d_summary(operators=[('cll2233', 'cle2233'), ('cll2233', 'cee2233'), ('cll2233', 'cle2332'), ('cee2233', 'cle2332')])
    # plot_limits_2d_summary(operators=[('cll2233', 'cle2233'), ('cee2233', 'cle2233'), ('cll2233', 'cle2332'), ('cee2233', 'cle2332')])
    # plot_limits_2d_summary(operators=[('cll2233', 'cle2233'), ('cee2233', 'cle2233'), ('cll2233', 'cee2233')])

    # only_plot_limits_2d_both(operators=[('cll2233', 'cle2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cll2233', 'cle3322')], options='afb')
    # only_plot_limits_2d_both(operators=[('cll2233', 'cee2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cll2332', 'cle2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cll2332', 'cle3322')], options='afb')
    # only_plot_limits_2d_both(operators=[('cll2332', 'cee2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2233', 'cle3322')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2233', 'cee2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle3322', 'cee2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2332', 'cee2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cee2233', 'cle2233')], options='afb')

    # only_plot_limits_2d_both(operators=[('cll2332', 'cll2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2332', 'cll2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2332', 'cle2233')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2332', 'cle3322')], options='afb')
    # only_plot_limits_2d_both(operators=[('cle2332', 'cll2332')], options='afb')
    opxopys = [
        ('cll2233', 'cee2233'),
        ('cll2332', 'cee2233'),
        ('cle2332', 'cle2233'),
        ('cle2332', 'cll2332'),
        ('cll2233', 'cle2233'),
        ('cll2233', 'cle3322'),
        ('cll2332', 'cle2233'),
        ('cll2332', 'cle3322'),
        ('cle2233', 'cle3322'),
        ('cle2233', 'cee2233'),
        ('cle3322', 'cee2233'),
        ('cle2332', 'cee2233'),
        ('cee2233', 'cle2233'),
        ('cll2332', 'cll2233'),
        ('cle2332', 'cll2233'),
        ('cle2332', 'cle3322'),
    ]
    # for (opx, opy) in opxopys:
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb', order='interference')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='af', order='interference')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='a', order='interference')
        # draw_limits_function_2d(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx=opx, opy=opy, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='', order='interference')
    draw_limits_function_2d_summary(operators=[('cle2332', 'cll2233'), ('cll2332', 'cle2233'), ('cle2332', 'cll2332'), ('cll2233', 'cle2233')], as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb', limittype='obs')
    draw_limits_function_2d_summary(operators=[('cle2332', 'cll2233'), ('cll2332', 'cle2233'), ('cle2332', 'cll2332'), ('cll2233', 'cle2233')], as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator, options='afb', limittype='exp')

def draw_limits_function_2d_summary(operators, as_per_operator, bs_per_operator, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, abrs_per_operator=None, bbrs_per_operator=None, options='', limittype='obs'):

    if limittype not in ['obs', 'exp']:
        raise ValueError('\'limittype\' must be either \'obs\' or \'exp\', but it is %s.' % (limittype))

    colors_per_opcomb = {
        # ('cll2233', 'cle2233'): ROOT.kRed+1, 
        # ('cll2233', 'cee2233'): ROOT.kOrange, 
        # ('cle2233', 'cee2233'): ROOT.kOrange, 
        # ('cee2233', 'cle2233'): ROOT.kOrange, 
        # ('cle2332', 'cee2233'): ROOT.kAzure-2, 
        # ('cee2233', 'cle2332'): ROOT.kAzure-2, 
        # ('cll2233', 'cle2332'): ROOT.kGreen-2, 
        # ('cll2233', 'cll2332'): ROOT.kGreen-2
        # ('cll2233', 'cle2233'): ROOT.TColor.GetColor('#d7191c'), 
        # ('cee2233', 'cle2233'): ROOT.TColor.GetColor('#fdae61'), 
        # ('cee2233', 'cle2332'): ROOT.TColor.GetColor('#abd9e9'), 
        # ('cll2233', 'cle2332'): ROOT.TColor.GetColor('#2c7bb6'), 
        # ('cll2233', 'cee2233'): ROOT.TColor.GetColor('#abd9e9'), 
        # ('cll2332', 'cll2233'): ROOT.kOrange, 
        ('cll2233', 'cle2233'): ROOT.TColor.GetColor('#abd9e9'),
        ('cll2332', 'cle2233'): ROOT.TColor.GetColor('#fdae61'),
        ('cle2332', 'cll2233'): ROOT.TColor.GetColor('#d7191c'),
        ('cle2332', 'cll2332'): ROOT.TColor.GetColor('#2c7bb6'),
    }


    limits_r = get_limits_r_at_cl(cl=0.95)

    funcs_obs = []
    graphs_obs = OrderedDict()
    funcs_obs_zoom = []
    graphs_obs_zoom = OrderedDict()
    # sf = 0.2
    sf = 5.
    scale = {}
    for (opx, opy) in operators:
        (ax, ay, bx, by, c, aincx, aincy, bincx, bincy, cinc, acrintx, acrinty, acrbsmx, acrbsmy, acrbsmintpart, abrx, abry, bbrx, bbry) = get_coefficients_for_operators(opx=opx, opy=opy, options=options, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator)

        if options == 'afb':
            def func_for_root(var):
                # result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) / (1./(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) + (ax*var[0])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrintx + (ay*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrinty + (bx*var[0]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmx + (by*var[1]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmy + (c*var[0]*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmintpart) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2

                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2
                return result
            
        elif options == 'af':
            def func_for_root(var):
                # result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) / (1./(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) + (ax*var[0])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrintx + (ay*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrinty + (bx*var[0]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmx + (by*var[1]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmy + (c*var[0]*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmintpart) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2

                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
                return result

        elif options == 'a':
            def func_for_root(var):
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
                return result

        elif options == '':
            def func_for_root(var):
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
                return result
        else:
            raise ValueError('Function for root undefined for options \'%s\'' % (options))

        xmin_orig = xmin = -1E4
        xmax_orig = xmax = +1E4
        ymin_orig = ymin = -1E4
        ymax_orig = ymax = +1E4
        if '2332' in opx or '2332' in opy:
            xmin_orig = xmin = -5E4
            xmax_orig = xmax = +5E4
            ymin_orig = ymin = -5E4
            ymax_orig = ymax = (+8./1.5)*1E4        

            if opx=='cll2332' and opy=='cll2233':
                ymax_orig = ymax = (12./1.5)*1E4/12*10
                ymin_orig = ymin = -5E4
                xmin_orig = xmin = -2.5E4
                xmax_orig = xmax = +2.5E4
        ymax_orig = ymax = ymax_orig*1.5
        funcs_obs.append(ROOT.TF2('func_obs', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig))
        np = 800
        # np = 100
        funcs_obs[-1].SetNpx(np)
        funcs_obs[-1].SetNpy(np)

        canvdummy = ROOT.TCanvas()
        if limittype == 'obs':
            funcs_obs[-1].SetContour(1, array('d', [limits_r[0]/r_sm]))
        elif limittype == 'exp':
            funcs_obs[-1].SetContour(1, array('d', [limits_r[1]/r_sm]))
        funcs_obs[-1].Draw('CONT LIST')
        canvdummy.Update()
        conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_obs[(opx, opy)] = [copy.deepcopy(x) for x in conts.At(0)]

        
        if not ('2332' in opx or '2332' in opy):
            scale[(opx, opy)] = True
        else:
            scale[(opx, opy)] = False
        if scale[(opx, opy)]:
            for g in graphs_obs[(opx, opy)]:
                for idx in range(g.GetN()):
                    g.GetX()[idx] *= sf
                    g.GetY()[idx] *= sf
    


    maxdigits = 3
    lumitag='138 fb^{-1} (13 TeV)'
    canv, axishist = tdrCanvas(canvName='canv', x_min=-3E4, x_max=+3E4, y_min=-3E4, y_max=+5E4, nameXaxis='C_{X} / #Lambda^{2} #scale[1.15]{[}#scale[0.90]{TeV^{-2}}#scale[1.15]{]}', nameYaxis='C_{Y} / #Lambda^{2} #scale[1.15]{[}#scale[0.90]{TeV^{-2}}#scale[1.15]{]}', square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.13, 0.13), maxdigits=(maxdigits, maxdigits), return_hist=True, ndivisions=(505, 505))
    axishist.GetXaxis().SetTitleOffset(0.90)
    axishist.GetYaxis().SetTitleOffset(1.00)

    # leg = tdrLeg(0.31,0.67,0.89,0.92, textSize=0.036)
    leg = tdrLeg(0.16,0.68,0.89,0.86, textSize=0.037)
    # leg = tdrLeg(0.31,0.73,0.86,0.92, textSize=0.034)
    if limittype == 'obs':
        leg.SetHeader('95% CL observed limits on C_{Y} vs. C_{X}')
    elif limittype == 'exp':
        leg.SetHeader('95% CL expected limits on C_{Y} vs. C_{X}')
    leg.SetNColumns(2)

    lstyle = 1 if limittype == 'obs' else 2
    for (opx, opy) in graphs_obs.keys():
        for g in graphs_obs[opx, opy]:
            tdrDraw(g, "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=lstyle, lwidth=-303)
        legentryname = '%s vs. %s' % (operatornames_pretty[opy], operatornames_pretty[opx])
        if scale[(opx, opy)]:
            if int(sf) == float(sf):
                legentryname += ' (#times%i)' % (int(sf))
            else:
                legentryname += ' (#times%1.1f)' % (sf)
        leg.AddEntry(graphs_obs[(opx, opy)][0], legentryname, 'L')
    
    for i, (opx, opy) in enumerate(graphs_obs.keys()):
        if i==1: continue
        for g in graphs_obs[opx, opy]:
            tdrDraw(g, "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=lstyle, lwidth=-303)



    # canv.RedrawAxis()
    # leg.Draw('SAME')

    width_smallpad  = 0.278
    height_smallpad = 0.278
    xfrac_left = (1. - canv.GetLeftMargin() - canv.GetRightMargin()) * (abs(xmin_orig))/(abs(xmax_orig - xmin_orig)) + canv.GetLeftMargin() - width_smallpad/2.
    xfrac_right = (1. - canv.GetLeftMargin() - canv.GetRightMargin()) * (abs(xmin_orig))/(abs(xmax_orig - xmin_orig)) + canv.GetLeftMargin() + width_smallpad/2.
    yfrac_left = (1. - canv.GetTopMargin() - canv.GetBottomMargin()) * (abs(ymin_orig))/(abs(ymax_orig - ymin_orig)) + canv.GetBottomMargin() - height_smallpad/2.
    yfrac_right = (1. - canv.GetTopMargin() - canv.GetBottomMargin()) * (abs(ymin_orig))/(abs(ymax_orig - ymin_orig)) + canv.GetBottomMargin() + height_smallpad/2.
    padsmall = ROOT.TPad('padsmall', 'padsmall', xfrac_left, yfrac_left, xfrac_right, yfrac_right)
    padsmall.SetFillColor(0)
    padsmall.SetBorderMode(1)
    padsmall.SetFrameFillStyle(0)
    padsmall.SetFrameBorderMode(0)
    padsmall.SetTopMargin(0.13)
    padsmall.SetRightMargin(0.15)
    padsmall.SetBottomMargin(0.13)
    padsmall.SetLeftMargin(0.15)
    padsmall.Draw()
    padsmall.cd()

    hsmall = padsmall.DrawFrame(-1E2, -1E4, +1E2, +1E4)
    hsmall.GetXaxis().SetNdivisions(502, False)
    hsmall.GetYaxis().SetNdivisions(502, False)
    hsmall.GetXaxis().SetMaxDigits(3)
    hsmall.GetYaxis().SetMaxDigits(3)
    hsmall.GetXaxis().SetLabelSize(0.10)
    hsmall.GetYaxis().SetLabelSize(0.10)
    hsmall.GetXaxis().SetLabelOffset(0.03)
    # hsmall.GetYaxis().SetLabelOffset(0.02)
    hsmall.Draw('AXIS')

    for (opx, opy) in operators:

        (ax, ay, bx, by, c, aincx, aincy, bincx, bincy, cinc, acrintx, acrinty, acrbsmx, acrbsmy, acrbsmintpart, abrx, abry, bbrx, bbry) = get_coefficients_for_operators(opx=opx, opy=opy, options=options, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator)

        def func_for_root(var):
            result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) / (1./(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) + (ax*var[0])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrintx + (ay*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrinty + (bx*var[0]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmx + (by*var[1]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmy + (c*var[0]*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmintpart) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2
            return result


        if not ('2332' in opx or '2332' in opy): continue    
        funcs_obs_zoom.append(ROOT.TF2('func_obs_zoom', func_for_root, -1E2, +1E2, -1E4, +1E4))
        if (opx, opy) == ('cle2332', 'cll2332'):
            funcs_obs_zoom[-1].SetRange(-6E1, -6E1, +6E1, +6E1)
        np = 200
        # if (opx == 'cll2332' and opy == 'cll2233'):
        #     funcs_obs_zoom[-1].SetRange(xmin*2, ymin*2, xmax*2, ymax*2)
        #     np = 350
        funcs_obs_zoom[-1].SetNpx(np)
        funcs_obs_zoom[-1].SetNpy(np)
    
        canvdummy2 = ROOT.TCanvas()
        if limittype == 'obs':
            funcs_obs_zoom[-1].SetContour(1, array('d', [limits_r[0]/r_sm]))
        elif limittype == 'exp':
            funcs_obs_zoom[-1].SetContour(1, array('d', [limits_r[1]/r_sm]))
        funcs_obs_zoom[-1].Draw('CONT LIST')
        canvdummy2.Update()
        conts_zoom = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_obs_zoom[(opx, opy)] = [copy.deepcopy(x) for x in conts_zoom.At(0)]
        if (opx, opy) == ('cle2332', 'cll2332'):
            for g in graphs_obs_zoom[(opx, opy)]:
                for idx in range(g.GetN()):
                    # g.GetX()[idx] *= 10
                    g.GetY()[idx] *= 100

    padsmall.cd()
    for (opx, opy) in graphs_obs_zoom:
        for g in graphs_obs_zoom[opx, opy]:
            tdrDraw(g, "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=lstyle, lwidth=-202)
    for i, (opx, opy) in enumerate(graphs_obs_zoom):
        if i==1: continue
        for g in graphs_obs_zoom[opx, opy]:
            tdrDraw(g, "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=lstyle, lwidth=-202)


    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.11)
    smtext.SetTextAlign(32)
    smtext.DrawLatex(0-0.15*abs(-1E2), 0., 'SM')


    scaletext = ROOT.TLatex()
    scaletext.SetTextFont(62)
    scaletext.SetTextSize(0.06)
    scaletext.SetTextAlign(32)
    scaletext.SetTextColor(colors_per_opcomb[('cle2332', 'cll2332')])
    scaletext.DrawLatex(0-0.21*abs(-1E2), 3E3, 'C_{LL}^{2332}#times100')

    canv.cd()
    for (opx, opy) in graphs_obs.keys():
        if ('2332' in opx and '2332' not in opy) or ('2332' in opy and '2332' not in opx): continue
        # if ('2332' in opx and '2332' in opy): continue
        for g in graphs_obs[opx, opy]:
            tdrDraw(g, "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=1, lwidth=-303)
        # legentryname = '%s vs. %s' % (operatornames_pretty[opy], operatornames_pretty[opx])
        # if scale[(opx, opy)]:
        #     if int(sf) == float(sf):
        #         legentryname += ' (#times%i)' % (int(sf))
        #     else:
        #         legentryname += ' (#times%1.1f)' % (sf)
        # leg.AddEntry(graphs_obs[opx, opy][0], legentryname, 'L')



    canv.RedrawAxis()
    leg.Draw('SAME')



    plotname = 'R_comb_2d_summary.pdf'
    if limittype == 'exp':
        plotname = plotname.replace('.pdf', '_expected.pdf')
    if options == 'a': plotname = plotname.replace('_comb_', '_comb_withacc_')
    if options == 'f': plotname = plotname.replace('_comb_', '_comb_withf_')
    if options == 'af': plotname = plotname.replace('_comb_', '_comb_withacc_withf_')
    if options == 'afb': plotname = plotname.replace('_comb_', '_comb_withacc_withf_withbr_')
    canv.SaveAs(os.path.join(plotfolder, plotname))



def draw_limits_function_2d(as_per_operator, bs_per_operator, opx, opy, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, abrs_per_operator=None, bbrs_per_operator=None, options='', order='comb'):
    limits_r = get_limits_r_at_cl(cl=0.95)
    if order not in ['interference', 'comb']: raise ValueError('Invalid order \'%s\'.')
    if order == 'interference' and 'cle2332' in [opx, opy]: return

    (ax, ay, bx, by, c, aincx, aincy, bincx, bincy, cinc, acrintx, acrinty, acrbsmx, acrbsmy, acrbsmintpart, abrx, abry, bbrx, bbry) = get_coefficients_for_operators(opx=opx, opy=opy, options=options, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, abrs_per_operator=abrs_per_operator, bbrs_per_operator=bbrs_per_operator)
    if order == 'interference':
        bx = by = c = bincx = bincy = cinc = acrbsmx = acrbsmy = acrbsmintpart = bbrx = bbry = 0.
    
    if options == 'afb':
        def func_for_root(var):
            # result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) / (1./(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) + (ax*var[0])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrintx + (ay*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrinty + (bx*var[0]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmx + (by*var[1]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmy + (c*var[0]*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmintpart) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])**2/(1 + aincx*var[0] + aincy*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty)/(1 + ax*var[0] + ay*var[1])) / ((1 + abrx*var[0] + abry*var[1])/(1+br_sm*(abrx*var[0] + abry*var[1])))**2
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2
            return result
            
        def func_tauwidth(var):
            if order == 'interference':
                # result = (1 + ax*var[0] + ay*var[1])
                result = (1 + abrx*var[0] + abry*var[1])
            else:
                result = 1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2
            return result
        
        def func_acc(var):
            acrintxsafe=acrintx if acrintx != 0 else -1.
            acrintysafe=acrinty if acrinty != 0 else -1.
            if order == 'interference':
                result = ((1 + ax*var[0]/acrintxsafe + ay*var[1]/acrintysafe)/(1 + ax*var[0] + ay*var[1]))
            else:
                result = ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
            return result
        
        def func_f(var):
            if order == 'interference':
                result = 1 + aincx*var[0] + aincy*var[1]
            else:
                result = 1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]
            return result
        
        def func_funphys(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])/(1 + aincx*var[0] + aincy*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1])
            return result
            
        def func_n(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
            return result
    elif options == 'af':
        def func_for_root(var):
            # result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) / (1./(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) + (ax*var[0])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrintx + (ay*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrinty + (bx*var[0]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmx + (by*var[1]**2)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmy + (c*var[0]*var[1])/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])*acrbsmintpart) / ((1 + abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)/(1+br_sm*(abrx*var[0] + abry*var[1] + bbrx*var[0]**2 + bbry*var[1]**2)))**2
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])**2/(1 + aincx*var[0] + aincy*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty)/(1 + ax*var[0] + ay*var[1]))
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])**2/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
            return result
            
        def func_n(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
            return result
        
        def func_f(var):
            if order == 'interference':
                result = 1 + aincx*var[0] + aincy*var[1]
            else:
                result = 1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1]
            return result
        
        def func_funphys(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])/(1 + aincx*var[0] + aincy*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])/(1 + aincx*var[0] + aincy*var[1] + bincx*var[0]**2 + bincy*var[1]**2 + cinc*var[0]*var[1])
            return result
        
        def func_acc(var):
            acrintxsafe=acrintx if acrintx != 0 else -1.
            acrintysafe=acrinty if acrinty != 0 else -1.
            if order == 'interference':
                result = ((1 + ax*var[0]/acrintxsafe + ay*var[1]/acrintysafe)/(1 + ax*var[0] + ay*var[1]))
            else:
                result = ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
            return result

    elif options == 'a':
        def func_for_root(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty)/(1 + ax*var[0] + ay*var[1]))
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]) * ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
            return result
            
        def func_n(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
            return result
        
        def func_acc(var):
            acrintxsafe=acrintx if acrintx != 0 else -1.
            acrintysafe=acrinty if acrinty != 0 else -1.
            if order == 'interference':
                result = ((1 + ax*var[0]/acrintxsafe + ay*var[1]/acrintysafe)/(1 + ax*var[0] + ay*var[1]))
            else:
                result = ((1 + ax*var[0]/acrintx + ay*var[1]/acrinty + bx*var[0]**2/acrbsmx + by*var[1]**2/acrbsmy + c*var[0]*var[1]/acrbsmintpart)/(1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1]))
            return result

    elif options == '':
        def func_for_root(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
            return result
            
        def func_n(var):
            if order == 'interference':
                result = (1 + ax*var[0] + ay*var[1])
            else:
                result = (1 + ax*var[0] + ay*var[1] + bx*var[0]**2 + by*var[1]**2 + c*var[0]*var[1])
            return result
    else:
        raise ValueError('Function for root undefined for options \'%s\'' % (options))

    xmin_orig = xmin = -1E4 if order != 'interference' else -1E5 # -5E4
    xmax_orig = xmax = +1E4 if order != 'interference' else 1E5 # 5E4
    ymin_orig = ymin = -1E4 if order != 'interference' else -2E5 # -5E4
    ymax_orig = ymax = +1E4 if order != 'interference' else 2E5 # 5E4
    if 'b' in options:
        if '2332' in opx or '2332' in opy:
            xmin_orig = xmin = -5E4           if order != 'interference' else -4E4
            xmax_orig = xmax = +5E4           if order != 'interference' else 4E4
            ymin_orig = ymin = -5E4           if order != 'interference' else -1E6
            ymax_orig = ymax = (+8./1.5)*1E4  if order != 'interference' else (4/1.5)*1E6

            if opx=='cll2332' and opy=='cll2233':
                xmin_orig = xmin = -2.5E4               if order != 'interference' else -4E4
                xmax_orig = xmax = +2.5E4               if order != 'interference' else 4E4
                ymin_orig = ymin = -5E4                 if order != 'interference' else -1E6
                ymax_orig = ymax = (12./1.5)*1E4/12*10  if order != 'interference' else (4/1.5)*1E6
            # if opx=='cll2332' and opy=='cee2233' and order=='interference':
            #     ymax_orig = ymax = 3E6
            #     ymin_orig = ymin = -3E6
            #     xmin_orig = xmin = -3E6
            #     xmax_orig = xmax = 3E6

    ymax_orig = ymax = ymax_orig*1.5
        



    func_95_hi = ROOT.TF2('func_95_hi', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    func_95_lo = ROOT.TF2('func_95_lo', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    func_68_hi = ROOT.TF2('func_68_hi', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    func_68_lo = ROOT.TF2('func_68_lo', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    func_exp   = ROOT.TF2('func_exp', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    func_obs   = ROOT.TF2('func_obs', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
    if order == 'interference':
        func_nneg            = ROOT.TF2('func_nneg', func_n, xmin_orig*2, xmax_orig*2, ymin_orig*2, ymax_orig*2)
        if 'b' in options:
            func_tauwidthneg = ROOT.TF2('func_tauwidthneg', func_tauwidth, xmin_orig*2, xmax_orig*2, ymin_orig*2, ymax_orig*2)
        if 'a' in options:
            func_accneg      = ROOT.TF2('func_accneg', func_acc, xmin_orig*2, xmax_orig*2, ymin_orig*2, ymax_orig*2)
        if 'f' in options:
            func_fneg      = ROOT.TF2('func_fneg', func_f, xmin_orig*2, xmax_orig*2, ymin_orig*2, ymax_orig*2)
            func_funphyshigh= ROOT.TF2('func_funphys', func_funphys, xmin_orig*2, xmax_orig*2, ymin_orig*2, ymax_orig*2)
    np = 600 if 'b' in options and ('2332' in opx or '2332' in opy) and order != 'interference' else 100
    if (opx == 'cll2332' and opy == 'cll2233') and order != 'interference':
        func_95_hi.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        func_95_lo.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        func_68_hi.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        func_68_lo.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        func_exp.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        func_obs.SetRange(xmin_orig*6.0, ymin_orig*6.0, xmax_orig*6.0, ymax_orig*6.0)
        np = 2100 if 'b' in options and ('2332' in opx or '2332' in opy) else 600
    # if order == 'interference' and not ('b' in options and opx=='cll2332' and opy=='cee2233'):
    if order == 'interference':
        factor_range = 10.
        if 'b' in options and opy == 'cee2233':
            factor = 2.
        func_95_hi.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        func_95_lo.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        func_68_hi.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        func_68_lo.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        func_exp.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        func_obs.SetRange(xmin_orig*factor_range, ymin_orig*factor_range, xmax_orig*factor_range, ymax_orig*factor_range)
        if 'b' in options and ('2332' in opx or '2332' in opy):
            func_95_hi.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
            func_95_lo.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
            func_68_hi.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
            func_68_lo.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
            func_exp.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
            func_obs.SetRange(xmin_orig*1.2, 0.001, -0.001, ymax_orig)
    # np = int(100 * xmax_orig / xmax)
    # if 'b' in options and opx=='cll2332' and opy=='cee2233' and order=='interference':
    # if order == 'interference' and ('2332' in opx or '2332' in opy):
    # if order == 'interference' and 'b' in options:
    if order == 'interference':
        np = 1500
        if 'b' in options and not ('2332' in opx or '2332' in opy):
            np = 200
            if opy == 'cee2233':
                np = 10
        # if '2332' in opx or '2332' in opy:
        #     np = 800
        # if 'b' in options:
            # np = 2000
    print 'number of grid points in x-direction: %i' % (np)
    func_95_hi.SetNpx(np)
    func_95_hi.SetNpy(np)
    func_95_lo.SetNpx(np)
    func_95_lo.SetNpy(np)
    func_68_hi.SetNpx(np)
    func_68_hi.SetNpy(np)
    func_68_lo.SetNpx(np)
    func_68_lo.SetNpy(np)
    func_exp.SetNpx(np)
    func_exp.SetNpy(np)
    func_obs.SetNpx(np)
    func_obs.SetNpy(np)

    # factor = 1.6
    # if '2332' in opx and opy == 'cee2233':
    #     factor = 3.
    if order == 'interference':
        func_nneg.SetNpx(2400)
        func_nneg.SetNpy(2400)
        if 'b' in options:
            func_tauwidthneg.SetNpx(2400)
            func_tauwidthneg.SetNpy(2400)
        if 'a' in options:
            func_accneg.SetNpx(2400)
            func_accneg.SetNpy(2400)
        if 'f' in options:
            func_fneg.SetNpx(2400)
            func_fneg.SetNpy(2400)
            func_funphyshigh.SetNpx(2400)
            func_funphyshigh.SetNpy(2400)

    # ROOT.gStyle.SetLineStyleString(11, '4 4')
    # ROOT.gStyle.SetLineStyleString(12, '4 48')

    canvdummy = ROOT.TCanvas()
    func_exp.SetContour(1, array('d', [limits_r[1]/r_sm]))
    func_exp.Draw('CONT LIST')
    canvdummy.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    graphs_exp = [copy.deepcopy(x) for x in conts.At(0)]

    canvdummyobs = ROOT.TCanvas()
    func_obs.SetContour(1, array('d', [limits_r[0]/r_sm]))
    func_obs.Draw('CONT LIST')
    canvdummyobs.Update()
    contsobs = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    graphs_obs = [copy.deepcopy(x) for x in contsobs.At(0)]

    if order == 'interference':
        canvdummy2 = ROOT.TCanvas()
        func_95_hi.SetContour(1, array('d', [limits_r[5]/r_sm]))
        func_95_hi.Draw('CONT LIST')
        canvdummy2.Update()
        conts2 = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_95_hi = [copy.deepcopy(x) for x in conts2.At(0)]


        canvdummy3 = ROOT.TCanvas()
        func_68_hi.SetContour(1, array('d', [limits_r[3]/r_sm]))
        func_68_hi.Draw('CONT LIST')
        canvdummy3.Update()
        conts3 = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_68_hi = [copy.deepcopy(x) for x in conts3.At(0)]


        canvdummy4 = ROOT.TCanvas()
        func_68_lo.SetContour(1, array('d', [limits_r[2]/r_sm]))
        func_68_lo.Draw('CONT LIST')
        canvdummy4.Update()
        conts4 = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_68_lo = [copy.deepcopy(x) for x in conts4.At(0)]


        canvdummy5 = ROOT.TCanvas()
        func_95_lo.SetContour(1, array('d', [limits_r[4]/r_sm]))
        func_95_lo.Draw('CONT LIST')
        canvdummy5.Update()
        conts5 = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_95_lo = [copy.deepcopy(x) for x in conts5.At(0)]

    maxdigits = 3
    lumitag='138 fb^{-1} (13 TeV)'
    canv = tdrCanvas(canvName='canv', x_min=xmin_orig, x_max=xmax_orig, y_min=ymin_orig, y_max=ymax_orig, nameXaxis='%s / #Lambda^{2} #scale[1.15]{[}#scale[0.90]{TeV^{-2}}#scale[1.15]{]}' % (operatornames_pretty[opx]), nameYaxis='%s / #Lambda^{2} #scale[1.15]{[}#scale[0.90]{TeV^{-2}}#scale[1.15]{]}' % (operatornames_pretty[opy]), square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), extraText='Supplementary', ndivisions=(505, 505))

    func_95_hi.SetContour(1, array('d', [limits_r[5]/r_sm]))
    func_68_hi.SetContour(1, array('d', [limits_r[3]/r_sm]))
    func_68_lo.SetContour(1, array('d', [limits_r[2]/r_sm]))
    func_95_lo.SetContour(1, array('d', [limits_r[4]/r_sm]))

    if order == 'interference':
        # tdrDraw(func_95_hi, 'CONT3', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=303)
        if not ('b' in options and ('2332' in opx or '2332' in opy)): # normal settings
            if 'b' in options:
                tdrDraw(graphs_95_lo[0], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
                tdrDraw(graphs_68_lo[0], 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1, fstyle=1001, lwidth=-9901)
                tdrDraw(graphs_68_hi[0], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
                tdrDraw(graphs_95_hi[0], 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite, fstyle=1001, lwidth=-9901)
                # tdrDraw(graphs_95_lo[0], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
                # tdrDraw(graphs_68_lo[0], 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
                # tdrDraw(graphs_68_hi[0], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
                # tdrDraw(graphs_95_hi[0], 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite)
            else:
                for g in graphs_95_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
                for g in graphs_68_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1, fstyle=1001, lwidth=-9901)
                for g in graphs_68_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
                for g in graphs_95_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite, fstyle=1001, lwidth=-9901)
        else: #special cases
            print len(graphs_95_lo), len(graphs_68_lo)
            # for g in graphs_95_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
            tdrDraw(graphs_95_lo[-1], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
            # for g in graphs_68_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1, fstyle=1001, lwidth=-9901)
            tdrDraw(graphs_68_lo[-1], 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1, fstyle=1001, lwidth=-9901)
            # for g in graphs_68_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
            tdrDraw(graphs_68_hi[-1], 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange, fstyle=1001, lwidth=-9901)
            # for g in graphs_95_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite, fstyle=1001, lwidth=-9901)
            tdrDraw(graphs_95_hi[-1], 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite, fstyle=1001, lwidth=-9901)
            # graphs_95_lo_forwhite = [copy.deepcopy(x) for x in graphs_95_lo]
            # for g in graphs_95_lo_forwhite: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kWhite, fcolor=ROOT.kWhite, fstyle=1001, lwidth=801)
            # for g in graphs_95_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        # for g in graphs_95_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        # for g in graphs_68_lo: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        # for g in graphs_68_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        # for g in graphs_95_hi: tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    else:
        if (opx == 'cll2332' and opy == 'cll2233'):
            tdrDraw(func_95_lo, 'CONT0', mcolor=0, fcolor=ROOT.kOrange, lcolor=ROOT.kOrange)
            tdrDraw(func_68_lo, 'CONT0', mcolor=0, fcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1)
            tdrDraw(func_68_hi, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
            tdrDraw(func_95_hi, 'CONT0', mcolor=0, fcolor=ROOT.kWhite)
            # ass

        else:
            tdrDraw(func_95_hi, 'CONT0', mcolor=0, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
            tdrDraw(func_68_hi, 'CONT0', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
            tdrDraw(func_68_lo, 'CONT0', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kOrange)
            tdrDraw(func_95_lo, 'CONT0', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kWhite)


    if order != 'interference' or not ('b' in options and ('2332' in opx or '2332' in opy)):
        for graph_exp in graphs_exp:
            tdrDraw(graph_exp, 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    else: #hack
        tdrDraw(graphs_exp[-1], 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    # func_obs.SetContour(1, array('d', [limits_r[0]/r_sm]))

    if order != 'interference' or not ('b' in options and ('2332' in opx or '2332' in opy)):
        for graph_obs in graphs_obs:
            tdrDraw(graph_obs, 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack)
    else:
        tdrDraw(graphs_obs[-1], 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack)
    # tdrDraw(func_obs, 'SURF', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack)

    if order == 'interference':
        if 'a' in options:
            func_accneg.SetContour(1, array('d', [0.]))
            tdrDraw(func_accneg, 'CONT0', mcolor=0, lcolor=ROOT.kGray, fcolor=ROOT.kGray, lstyle=1, lwidth=2)
            # tdrDraw(func_accneg, 'CONT3', mcolor=0, lcolor=ROOT.kRed, lstyle=1, lwidth=2)
        if 'f' in options:
            func_funphyshigh.SetContour(1, array('d', [1/f_sm]))
            tdrDraw(func_funphyshigh, 'CONT0', mcolor=0, lcolor=ROOT.kGray, fcolor=ROOT.kGray, lstyle=1, lwidth=2)
            # tdrDraw(func_funphyshigh, 'CONT3', mcolor=0, lcolor=ROOT.kBlue, lstyle=1, lwidth=2)

        canvdummy_n = ROOT.TCanvas()
        func_nneg.SetContour(1, array('d', [0.]))
        func_nneg.Draw('CONT LIST')
        canvdummy_n.Update()
        conts_n = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_nneg = [copy.deepcopy(x) for x in conts_n.At(0)]
        canv.cd()
        for g in graphs_nneg: 
            tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGray, fcolor=ROOT.kGray, fstyle=1001, lwidth=9901)
            # tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGreen-7, lwidth=2)

        if 'f' in options:
            canvdummy_fneg = ROOT.TCanvas()
            func_fneg.SetContour(1, array('d', [0.]))
            func_fneg.Draw('CONT LIST')
            canvdummy_fneg.Update()
            conts_fneg = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
            graphs_fneg = [copy.deepcopy(x) for x in conts_fneg.At(0)]
            canv.cd()
            for g in graphs_fneg: 
                tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGray, fcolor=ROOT.kGray, fstyle=1001, lwidth=9901)
                # tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kCyan, lwidth=2)

        if 'b' in options:
            canvdummy_tauwidthneg = ROOT.TCanvas()
            func_tauwidthneg.SetContour(1, array('d', [0.]))
            func_tauwidthneg.Draw('CONT LIST')
            canvdummy_tauwidthneg.Update()
            conts_tauwidthneg = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
            graphs_tauwidthneg = [copy.deepcopy(x) for x in conts_tauwidthneg.At(0)]
            canv.cd()
            for g in graphs_tauwidthneg: 
                tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kGray, fcolor=ROOT.kGray, fstyle=1001, lwidth=9901)
                # tdrDraw(g, 'L', mcolor=0, lcolor=ROOT.kOrange, lwidth=2)




    legcomb = tdrLeg(0.53,0.70,0.89,0.90, textSize=0.038)
    legcomb.SetHeader('95% CL limits')
    legcomb.AddEntry(graphs_obs[-1], 'Observed', 'L')
    legcomb.AddEntry(graphs_exp[-1], 'Median expected', 'L')
    if (opx == 'cll2332' and opy == 'cll2233'):
        if order == 'interference':
            legcomb.AddEntry(graphs_68_lo[-1], '68% expected', 'LF')
            legcomb.AddEntry(graphs_95_lo[-1], '95% expected', 'LF')
        else:
            legcomb.AddEntry(func_68_lo, '68% expected', 'LF')
            legcomb.AddEntry(func_95_lo, '95% expected', 'LF')
    else:
        if order == 'interference':
            if 'b' in options:
                legcomb.AddEntry(graphs_68_lo[0], '68% expected', 'LF')
                legcomb.AddEntry(graphs_95_lo[0], '95% expected', 'LF')
            else:
                legcomb.AddEntry(graphs_68_lo[-1], '68% expected', 'LF')
                legcomb.AddEntry(graphs_95_lo[-1], '95% expected', 'LF')
        else:
            legcomb.AddEntry(func_68_hi, '68% expected', 'LF')
            legcomb.AddEntry(func_95_hi, '95% expected', 'LF')
    legcomb.Draw('SAME')

    CMS_lumi(pad=canv, iPosX=11, lumitag=lumitag, extraText='Supplementary')



    canv.RedrawAxis()

    if ('2332' in opx or '2332' in opy) and 'b' in options and order != 'interference':
        xmin = -1E2
        xmax = +1E2
        ymin = -1E4
        ymax = +1E4
        if opx == 'cle2332' and opy == 'cll2332':
            xmin = -6E1
            xmax = +6E1
            ymin = -6E1
            ymax = +6E1

        width_smallpad = 0.15  
        height_smallpad = 0.15 
        if opx == 'cle2332' and opy == 'cll2332':
            width_smallpad = 0.14  
            height_smallpad = 0.14 
        xfrac_left = (1. - canv.GetLeftMargin() - canv.GetRightMargin()) * (abs(xmin_orig))/(abs(xmax_orig - xmin_orig)) + canv.GetLeftMargin() - width_smallpad/2.
        xfrac_right = (1. - canv.GetLeftMargin() - canv.GetRightMargin()) * (abs(xmin_orig))/(abs(xmax_orig - xmin_orig)) + canv.GetLeftMargin() + width_smallpad/2.
        yfrac_left = (1. - canv.GetTopMargin() - canv.GetBottomMargin()) * (abs(ymin_orig))/(abs(ymax_orig - ymin_orig)) + canv.GetBottomMargin() - height_smallpad/2.
        yfrac_right = (1. - canv.GetTopMargin() - canv.GetBottomMargin()) * (abs(ymin_orig))/(abs(ymax_orig - ymin_orig)) + canv.GetBottomMargin() + height_smallpad/2.
        padsmall = ROOT.TPad('padsmall', 'padsmall', xfrac_left, yfrac_left, xfrac_right, yfrac_right)
        padsmall.SetFillColor(0)
        padsmall.SetBorderMode(1)
        padsmall.SetFrameFillStyle(0)
        padsmall.SetFrameBorderMode(0)
        padsmall.SetTopMargin(0.13)
        padsmall.SetRightMargin(0.15)
        padsmall.SetBottomMargin(0.13)
        padsmall.SetLeftMargin(0.15)
        padsmall.Draw()
        padsmall.cd()

        hsmall = padsmall.DrawFrame(xmin, ymin, xmax, ymax)
        # if opx == 'cle2332' and opy == 'cll2332':
        #     hsmall = padsmall.DrawFrame(-1E2, -1E2, +1E2, +1E2)
        # else:
        #     hsmall = padsmall.DrawFrame(-2E2, -3E4, +2E2, +3E4)
        hsmall.GetXaxis().SetNdivisions(502, False)
        hsmall.GetYaxis().SetNdivisions(502, False)
        if opx == 'cle2332' and opy == 'cll2332':
            hsmall.GetXaxis().SetNdivisions(302, False)
            hsmall.GetYaxis().SetNdivisions(302, False)
        hsmall.GetXaxis().SetMaxDigits(3)
        hsmall.GetYaxis().SetMaxDigits(3)
        hsmall.GetXaxis().SetLabelSize(0.10)
        hsmall.GetYaxis().SetLabelSize(0.10)
        hsmall.Draw('AXIS')

        # func_95_hi_zoom = ROOT.TF2('func_95_hi_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        # func_95_lo_zoom = ROOT.TF2('func_95_lo_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        # func_68_hi_zoom = ROOT.TF2('func_68_hi_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        # func_68_lo_zoom = ROOT.TF2('func_68_lo_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        # func_exp_zoom   = ROOT.TF2('func_exp_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        # func_obs_zoom   = ROOT.TF2('func_obs_zoom', func_for_root, xmin_orig, xmax_orig, ymin_orig, ymax_orig)
        func_95_hi_zoom = ROOT.TF2('func_95_hi_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        func_95_lo_zoom = ROOT.TF2('func_95_lo_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        func_68_hi_zoom = ROOT.TF2('func_68_hi_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        func_68_lo_zoom = ROOT.TF2('func_68_lo_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        func_exp_zoom   = ROOT.TF2('func_exp_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        func_obs_zoom   = ROOT.TF2('func_obs_zoom', func_for_root, xmin*2, xmax*2, ymin*2, ymax*2)
        if opx == 'cll2332':
            func_95_hi_zoom.SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
            func_95_lo_zoom.SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
            func_68_hi_zoom.SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
            func_68_lo_zoom.SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
            func_exp_zoom  .SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
            func_obs_zoom  .SetRange(xmin*2, ymin*5, xmax*2, ymax*5)
        # np = int(100 * xmax_orig / xmax)
        np = 100
        print 'number of grid points in x-direction: %i' % (np)
        func_95_hi_zoom.SetNpx(np)
        func_95_hi_zoom.SetNpy(np)
        func_95_lo_zoom.SetNpx(np)
        func_95_lo_zoom.SetNpy(np)
        func_68_hi_zoom.SetNpx(np)
        func_68_hi_zoom.SetNpy(np)
        func_68_lo_zoom.SetNpx(np)
        func_68_lo_zoom.SetNpy(np)
        func_exp_zoom.SetNpx(np)
        func_exp_zoom.SetNpy(np)
        func_obs_zoom.SetNpx(np)
        func_obs_zoom.SetNpy(np)

        canvdummy2 = ROOT.TCanvas()
        func_exp_zoom.SetContour(1, array('d', [limits_r[1]/r_sm]))
        func_exp_zoom.Draw('CONT LIST')
        canvdummy2.Update()
        conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        graphs_exp_zoom = [copy.deepcopy(x) for x in conts.At(0)]
        # graph_exp_zoom = copy.deepcopy(graph_exp_zoom_tmp)
        padsmall.cd()
    
        func_95_hi_zoom.SetContour(1, array('d', [limits_r[5]/r_sm]))
        func_68_hi_zoom.SetContour(1, array('d', [limits_r[3]/r_sm]))
        # tdrDraw(func_68_hi_zoom, 'CONT0', mcolor=0, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        func_68_lo_zoom.SetContour(1, array('d', [limits_r[2]/r_sm]))
        # tdrDraw(func_68_lo_zoom, 'CONT0', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kOrange)
        func_95_lo_zoom.SetContour(1, array('d', [limits_r[4]/r_sm]))
        if opx == 'cle2332':
            tdrDraw(func_95_lo_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
            tdrDraw(func_68_lo_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kGreen+1)
            tdrDraw(func_68_hi_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
            tdrDraw(func_95_hi_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kWhite)
        elif opx == 'cll2332':
            tdrDraw(func_95_hi_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
            tdrDraw(func_68_hi_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kGreen+1)
            tdrDraw(func_68_lo_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
            tdrDraw(func_95_lo_zoom, 'CONT0', mcolor=0, fcolor=ROOT.kWhite)
    
        # func_exp_zoom.SetContour(1, array('d', [limits_r[1]/r_sm]))
        # tdrDraw(func_exp_zoom, 'CONT3', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=11, lwidth=1)
        for i in range(len(graphs_exp_zoom)):
            tdrDraw(graphs_exp_zoom[i], 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2, lwidth=1)
        
        func_obs_zoom.SetContour(1, array('d', [limits_r[0]/r_sm]))
        tdrDraw(func_obs_zoom, 'CONT3', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lwidth=1)

        padsmall.RedrawAxis()

        marker_sm = ROOT.TMarker(0., 0., 8)
        marker_sm.SetMarkerSize(0.8)
        marker_sm.Draw('SAME')
    
    
    
        smtext = ROOT.TLatex()
        smtext.SetTextFont(42)
        smtext.SetTextSize(0.15)
        smtext.SetTextAlign(32)
        smtext.DrawLatex(0-0.15*abs(xmin), 0., 'SM')


        # if (opx == 'cll2332' and opy == 'cll2233'):
        #     canv.cd()
        #     tdrDraw(func_95_lo, 'CONT0', mcolor=0, fcolor=ROOT.kOrange, lcolor=ROOT.kOrange)
        #     tdrDraw(func_68_lo, 'CONT0', mcolor=0, fcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1)
        #     tdrDraw(func_68_hi, 'CONT0', mcolor=0, fcolor=ROOT.kOrange)
        #     tdrDraw(func_95_hi, 'CONT0', mcolor=0, fcolor=ROOT.kWhite)
        #     for graph_exp in graphs_exp:
        #         tdrDraw(graph_exp, 'L', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        #     tdrDraw(func_obs, 'CONT3', mcolor=0, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack)
        #     canv.RedrawAxis()





    
    else:
        marker_sm = ROOT.TMarker(0., 0., 8)
        marker_sm.Draw('SAME')

        smtext = ROOT.TLatex()
        smtext.SetTextFont(42)
        smtext.SetTextSize(0.035)
        smtext.SetTextAlign(32)
        smtext.DrawLatex(0-0.03*abs(xmin), 0., 'SM')
    
        if order == 'interference':
            unphystext = ROOT.TLatex()
            unphystext.SetTextFont(52)
            unphystext.SetTextSize(0.031)
            unphystext.SetTextAlign(22)
            unphystext.DrawLatex(xmin_orig/2., ymin_orig/2., 'unphysical')


    plotname = 'R_comb_2d_%s_vs_%s.pdf' % (opx, opy)
    if options == 'a': plotname = plotname.replace('_comb_', '_comb_withacc_')
    if options == 'f': plotname = plotname.replace('_comb_', '_comb_withf_')
    if options == 'af': plotname = plotname.replace('_comb_', '_comb_withacc_withf_')
    if options == 'afb': plotname = plotname.replace('_comb_', '_comb_withacc_withf_withbr_')
    if order == 'interference': plotname = plotname.replace('_comb_', '_interference_')
    canv.SaveAs(os.path.join(plotfolder, plotname))

def only_plot_limits_2d_both(operators, options):
    
    limits_per_ops_comb = OrderedDict()
    limits_per_ops_interference = OrderedDict()
    for (opx, opy) in operators:
        infilename_comb = os.path.join(filefolder, 'limits_comb_%s_%s_vs_%s.root' % (options, opx, opy))
        infile_comb = ROOT.TFile(infilename_comb, 'READ')
        if not 'cle2332' in [opx, opy]:
            infilename_interference = os.path.join(filefolder, 'limits_interference_%s_%s_vs_%s.root' % (options, opx, opy))
            infile_interference = ROOT.TFile(infilename_interference, 'READ')

        limits_per_ops_comb[(opx, opy)] = (infile_comb.Get('g_obs_p'), infile_comb.Get('g_obs_n'), infile_comb.Get('g_exp_p'), infile_comb.Get('g_exp_n'), infile_comb.Get('g_68_p'), infile_comb.Get('g_95_p'), infile_comb.Get('g_68_n'), infile_comb.Get('g_95_n'))
        if not 'cle2332' in [opx, opy]:
            limits_per_ops_interference[(opx, opy)] = (infile_interference.Get('g_obs_p'), infile_interference.Get('g_obs_n'), infile_interference.Get('g_exp_p'), infile_interference.Get('g_exp_n'), infile_interference.Get('g_68_p'), infile_interference.Get('g_95_p'), infile_interference.Get('g_68_n'), infile_interference.Get('g_95_n'))

    maxdigits = 3
    lumitag='138 fb^{-1} (13 TeV)'
    xmin = -1E4
    xmax = +1E4
    ymin = -1E4
    ymax = +1E4
    if '2332' in opx or '2332' in opy:
        xmin = -5E4
        xmax = +5E4
        ymin = -5E4
        ymax = (+8./1.5)*1E4

        # xmin/=40.
        # xmax /=40.
        if opx=='cll2332' and opy=='cll2233':
            ymax = (12./1.5)*1E4
    ccomb = tdrCanvas(canvName='ccomb', x_min=xmin, x_max=xmax, y_min=ymin, y_max=ymax*1.5, nameXaxis=operatornames_pretty[opx], nameYaxis=operatornames_pretty[opy], square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), extraText='Supplementary', ndivisions=(505, 505))


    legcomb = tdrLeg(0.50,0.70,0.89,0.90, textSize=0.038)
    legcomb.SetHeader('95% CL limits')
    legcomb.AddEntry(limits_per_ops_comb[(opx, opy)][0], 'Observed', 'L')
    legcomb.AddEntry(limits_per_ops_comb[(opx, opy)][2], 'Median expected', 'L')
    legcomb.AddEntry(limits_per_ops_comb[(opx, opy)][4], '68% expected', 'LF')
    legcomb.AddEntry(limits_per_ops_comb[(opx, opy)][5], '95% expected', 'LF')
    legcomb.Draw('SAME')

    for (opx, opy) in limits_per_ops_comb.keys():
        tdrDraw(limits_per_ops_comb[(opx, opy)][5], "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        tdrDraw(limits_per_ops_comb[(opx, opy)][7], "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        tdrDraw(limits_per_ops_comb[(opx, opy)][4], "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        tdrDraw(limits_per_ops_comb[(opx, opy)][6], "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        tdrDraw(limits_per_ops_comb[(opx, opy)][2], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        tdrDraw(limits_per_ops_comb[(opx, opy)][3], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        # tdrDraw(limits_per_ops_comb[(opx, opy)][0], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=303)
        # tdrDraw(limits_per_ops_comb[(opx, opy)][1], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-303)
        tdrDraw(limits_per_ops_comb[(opx, opy)][0], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
        tdrDraw(limits_per_ops_comb[(opx, opy)][1], "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
    

    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)
    smtext.DrawLatex(0., 4.5E2, 'SM')

    ccomb.RedrawAxis()
    legcomb.Draw('SAME')

    outname = os.path.join(plotfolder, 'Limits2d_paperstyle_comb_%s_vs_%s.pdf' % (opx, opy))
    if options == 'a': outname = outname.replace('_comb_', '_comb_withacc_')
    if options == 'f': outname = outname.replace('_comb_', '_comb_withf_')
    if options == 'af': outname = outname.replace('_comb_', '_comb_withacc_withf_')
    if options == 'afb': outname = outname.replace('_comb_', '_comb_withacc_withf_withbr_')
    ccomb.SaveAs(outname)
    
def plot_r_vs_c(a, b, op, plotname, acc_int=None, acc_bsm=None, acc_sm=None, ainc=None, binc=None, abr=None, bbr=None, options=''):

    limits_r = get_limits_r_at_cl(cl=0.95)

    if op == 'cle2332': acc_int = 1.
    if not 'a' in options: 
        acc_int = acc_bsm = acc_sm = 1.
    if not 'f' in options:
        ainc = a
        binc = b
    if not 'b' in options: 
        abr = bbr = 0.

    def rvsc(x, par):
        return rvsc2(x[0], par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8])

    def rvsc2(x, a, b, acc_int, acc_bsm, acc_sm, ainc, binc, abr, bbr):
        acc_int_safe = acc_int if acc_int != 0 else -1.
        acc_bsm_safe = acc_bsm if acc_bsm != 0 else -1.
        return (1+a*x+b*x**2)**2/(1+ainc*x+binc*x**2) * ((1+a*x*acc_sm/acc_int+b*x**2*acc_sm/acc_bsm)/(1+a*x+b*x**2)) / ((1+abr*x+bbr*x**2)/(1+br_sm*(abr*x + bbr*x**2)))**2

    xmin = -25000.
    xmax = +25000.
    rfunc = ROOT.TF1('r', rvsc, xmin, xmax, 9)
    rfunc.FixParameter(0, a)
    rfunc.FixParameter(1, b)
    rfunc.FixParameter(2, acc_int)
    rfunc.FixParameter(3, acc_bsm)
    rfunc.FixParameter(4, acc_sm)
    rfunc.FixParameter(5, ainc)
    rfunc.FixParameter(6, binc)
    rfunc.FixParameter(7, abr)
    rfunc.FixParameter(8, bbr)

    c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=-1., y_max=10., nameXaxis=operatornames_pretty[op], nameYaxis='R / R^{SM}', square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(3,3), ndivisions=(505, 505))
    rfunc.SetLineWidth(2)
    rfunc.Draw('SAME')

    print limits_r[0], limits_r[0]/r_sm
    line_obs = ROOT.TLine(xmin, limits_r[0]/r_sm , xmax, limits_r[0]/r_sm)
    line_obs.SetLineWidth(2)
    line_obs.Draw('SAME')

    c.SaveAs(os.path.join(plotfolder, plotname))

    # zoom around 0

    czoom = tdrCanvas(canvName='c', x_min=-3210, x_max=-2400, y_min=-2., y_max=10., nameXaxis=operatornames_pretty[op], nameYaxis='R / R^{SM}', square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(2,2), ndivisions=(505, 505))
    rfunc.SetLineWidth(2)
    rfunc.Draw('SAME')

    print limits_r[0], limits_r[0]/r_sm
    line_obs_zoom = ROOT.TLine(-3210, limits_r[0]/r_sm , -2400, limits_r[0]/r_sm)
    line_obs_zoom.SetLineWidth(2)
    line_obs_zoom.Draw('SAME')

    czoom.SaveAs(os.path.join(plotfolder, plotname.replace('.pdf', '_zoom.pdf')))

def plot_limits_2d_summary(operators):
    colors_per_opcomb = {
        # ('cll2233', 'cle2233'): ROOT.kRed+1, 
        # ('cll2233', 'cee2233'): ROOT.kOrange, 
        # ('cle2233', 'cee2233'): ROOT.kOrange, 
        # ('cee2233', 'cle2233'): ROOT.kOrange, 
        # ('cle2332', 'cee2233'): ROOT.kAzure-2, 
        # ('cee2233', 'cle2332'): ROOT.kAzure-2, 
        # ('cll2233', 'cle2332'): ROOT.kGreen-2, 
        # ('cll2233', 'cll2332'): ROOT.kGreen-2
        ('cll2233', 'cle2233'): ROOT.TColor.GetColor('#d7191c'), 
        ('cee2233', 'cle2233'): ROOT.TColor.GetColor('#fdae61'), 
        ('cee2233', 'cle2332'): ROOT.TColor.GetColor('#abd9e9'), 
        ('cll2233', 'cle2332'): ROOT.TColor.GetColor('#2c7bb6'), 
        ('cll2233', 'cee2233'): ROOT.TColor.GetColor('#abd9e9'), 

    }
    
    limits_per_ops = OrderedDict()
    for (opx, opy) in operators:
        infilename = os.path.join(filefolder, 'limits_comb_%s_vs_%s.root' % (opx, opy))
        infile = ROOT.TFile(infilename, 'READ')

        limits_per_ops[(opx, opy)] = (infile.Get('g_obs_p'), infile.Get('g_obs_n'), infile.Get('g_exp_p'), infile.Get('g_exp_n'), infile.Get('g_68_p'), infile.Get('g_95_p'), infile.Get('g_68_n'), infile.Get('g_95_n'))

    maxdigits = 3
    lumitag='138 fb^{-1} (13 TeV)'
    # c, chist = tdrCanvas(canvName='c', x_min=-5.5E3, x_max=+5.5E3, y_min=-5.5E3, y_max=8.25E3, nameXaxis='C_{X}', nameYaxis='C_{Y}', square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), return_hist=True)
    c, chist = tdrCanvas(canvName='c', x_min=-6.0E3, x_max=+6.0E3, y_min=-6.0E3, y_max=8.0E3, nameXaxis='C_{X}', nameYaxis='C_{Y}', square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), return_hist=True)

    # # # STOP AT 2.5
    # yaxis = ROOT.TGaxis(ROOT.gPad.GetUxmin(),ROOT.gPad.GetUymin(),ROOT.gPad.GetUxmin(),(5.5E3/+8.25E3)*ROOT.gPad.GetUymax(),-5.5E3,5.5E3,510,'R')
    # yaxis.SetLabelSize(0)
    # yaxis.SetLabelFont(chist.GetYaxis().GetLabelFont())
    # yaxis.SetLabelColor(chist.GetYaxis().GetLabelColor())
    # yaxis.SetLabelSize(chist.GetYaxis().GetLabelSize())
    # yaxis.SetLabelOffset(chist.GetYaxis().GetLabelOffset())
    # yaxis.SetNdivisions(chist.GetYaxis().GetNdivisions())
    # yaxis.SetTitleOffset(chist.GetYaxis().GetTitleOffset())
    # yaxis.SetTextFont(42)
    # yaxis.SetTitleSize(chist.GetYaxis().GetTitleSize())
    # yaxis.SetTitle(chist.GetYaxis().GetTitle())
    # # chist.GetYaxis().SetLabelSize(0)
    # chist.GetYaxis().ChangeLabel(6, -1, 0, -1, -1, -1, '')
    # chist.GetYaxis().ChangeLabel(7, -1, 0, -1, -1, -1, '')
    # ROOT.TGaxis.SetExponentOffset(-0.095, -0.155, 'y')
    # yaxis.Draw()
    # chist.GetYaxis().SetTitle('')



    # # new x-axis on the bottom of the white box
    # xaxistop = ROOT.TGaxis(ROOT.gPad.GetUxmin(),(5.5E3/+8.25E3)*ROOT.gPad.GetUymax(),ROOT.gPad.GetUxmax(),(5.5E3/+8.25E3)*ROOT.gPad.GetUymax(),-5.5E3,+5.5E3,508,'-')
    # xaxistop.SetLabelFont(chist.GetYaxis().GetLabelFont())
    # xaxistop.SetLabelColor(chist.GetYaxis().GetLabelColor())
    # xaxistop.SetLabelSize(chist.GetYaxis().GetLabelSize())
    # xaxistop.SetLabelOffset(chist.GetXaxis().GetLabelOffset())
    # xaxistop.SetNdivisions(chist.GetYaxis().GetNdivisions())
    # xaxistop.SetLabelSize(0)
    # xaxistop.Draw()

    # # BOX COVERING TOP PART OF THE PLOT
    # box1 = ROOT.TBox(-5.5E3, +5.5E3, +5.5E3, 8.25E3) # fill
    # box2 = ROOT.TBox(-5.5E3, +5.5E3, +5.5E3, 8.25E3) # line
    # box1.SetFillStyle(1001)
    # box1.SetFillColor(0)
    # box2.SetFillStyle(0)
    # box2.SetLineStyle(1)
    # box2.SetLineWidth(1)
    # box2.SetLineColor(kBlack)

    # leg = tdrLeg(0.18,0.78,0.90,0.925, textSize=0.037)
    # leg = tdrLeg(0.185,0.67,0.90,0.85, textSize=0.040)
    leg = tdrLeg(0.31,0.73,0.89,0.92, textSize=0.036)
    leg.SetHeader('95% CL observed limits on C_{Y} vs. C_{X}')
    leg.SetNColumns(2)

    for (opx, opy) in limits_per_ops.keys():
        tdrDraw(limits_per_ops[(opx, opy)][0], "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=1, lwidth=303)
        tdrDraw(limits_per_ops[(opx, opy)][1], "L SAME", mcolor=colors_per_opcomb[(opx, opy)], lcolor=colors_per_opcomb[(opx, opy)], fcolor=colors_per_opcomb[(opx, opy)], fstyle=3013, lstyle=1, lwidth=-303)
    
        leg.AddEntry(limits_per_ops[(opx, opy)][0], '%s vs. %s' % (operatornames_pretty[opy], operatornames_pretty[opx]), 'L')


    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)
    smtext.DrawLatex(0., 4.5E2, 'SM')

    c.RedrawAxis()

    # box1.Draw()
    # box2.Draw()
    leg.Draw('SAME')
    # CMS_lumi( pad=c, iPosX=11, lumitag=lumitag)

    # cmstag = ROOT.TLatex()
    # cmstag.SetNDC()
    # cmstag.SetTextAlign(11)
    # cmstag.SetTextColor(kBlack)
    # cmstag.SetTextSize(0.045*0.95)
    # cmstag.SetTextFont(61)
    # cmstag.DrawLatex(0.16, 0.94, 'CMS')

    # supptag = ROOT.TLatex()
    # supptag.SetNDC()
    # supptag.SetTextAlign(11)
    # supptag.SetTextColor(kBlack)
    # supptag.SetTextSize(0.040*0.95)
    # supptag.SetTextFont(52)
    # supptag.DrawLatex(0.255, 0.94, 'Supplementary')

    outname = os.path.join(plotfolder, 'Limits2d_summaryall_comb.pdf')
    c.SaveAs(outname)


def plot_limits_2d_both(as_per_operator, bs_per_operator, opx, opy, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, abrs_per_operator=None, bbrs_per_operator=None, options=''):

    limits_r = get_limits_r_at_cl(cl=0.95)



    operatorname_sum = '%s+1p0x%s' % (opx, opy)
    if not operatorname_sum in bs_per_operator: raise ValueError('Could not find combined term c for operator-sum %s' % (operatorname_sum))

    ax = as_per_operator[opx]
    ay = as_per_operator[opy]
    if opx == 'cle2332': ax = 0.
    if opy == 'cle2332': ay = 0.
    bx = bs_per_operator[opx]
    by = bs_per_operator[opy]
    sumcoeff  = bs_per_operator[operatorname_sum]
    c = sumcoeff - bx - by

    acc_int_x = 1. if acceptances_int_per_operator == None else acceptances_int_per_operator[opx]
    acc_int_y = 1. if acceptances_int_per_operator == None else acceptances_int_per_operator[opy]
    if opx == 'cle2332': acc_int_x = 1.
    if opy == 'cle2332': acc_int_y = 1.
    acc_bsm_x = 1. if acceptances_bsm_per_operator == None else acceptances_bsm_per_operator[opx]
    acc_bsm_y = 1. if acceptances_bsm_per_operator == None else acceptances_bsm_per_operator[opy]
    acc_sm_x  = 1. if acceptances_sm_per_operator  == None else acceptances_sm_per_operator[opx]
    acc_sm_y  = 1. if acceptances_sm_per_operator  == None else acceptances_sm_per_operator[opy]
    if 1. in [acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_sm_x, acc_sm_y] and 'a' in options and not 'cle2332' in [opx, opy]: raise ValueError('Could not get all acceptances needed for using option \'f\'.')
    if opx == 'cle2332' and acc_int_y == 1: raise ValueError('Could not get all acceptances needed for using option \'f\'.')
    if opy == 'cle2332' and acc_int_x == 1: raise ValueError('Could not get all acceptances needed for using option \'f\'.')

    operatorname_sum_foracc = operatorname_sum if acceptances_bsm_per_operator != None and operatorname_sum in acceptances_bsm_per_operator else '%s+1p0x%s' % (opy, opx)
    if acceptances_bsm_per_operator == None or not operatorname_sum_foracc in acceptances_bsm_per_operator: raise ValueError('Could not find combined acceptance term for operator-sum %s' % (operatorname_sum_foracc))
    acc_bsm_withint = acceptances_bsm_per_operator[operatorname_sum_foracc]
    acc_bsm_intpart = get_acc_bsm_intpart(bx=bx, by=by, c=c, accx=acc_bsm_x, accy=acc_bsm_y, acctot=acc_bsm_withint)


    aincx = 0. if aincs_per_operator == None else aincs_per_operator[opx]
    aincy = 0. if aincs_per_operator == None else aincs_per_operator[opy]
    bincx = 0. if bincs_per_operator == None else bincs_per_operator[opx]
    bincy = 0. if bincs_per_operator == None else bincs_per_operator[opy]
    if 0. in [aincx, aincy, bincx, bincy] and 'f' in options and not 'cle2332' in [opx, opy]: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
    if opx == 'cle2332': 
        if aincy == 0.: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
        aincx = 0.
    if opy == 'cle2332': 
        if aincx == 0.: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
        aincy = 0.
    if not operatorname_sum in bincs_per_operator: raise ValueError('Could not find combined term c for inclusive operator-sum %s' % (operatorname_sum))
    sumcoeffinc  = bincs_per_operator[operatorname_sum]
    cinc = sumcoeffinc - bincx - bincy

    abrx = 0. if abrs_per_operator == None else abrs_per_operator[opx]
    abry = 0. if abrs_per_operator == None else abrs_per_operator[opy]
    bbrx = 0. if bbrs_per_operator == None else bbrs_per_operator[opx]
    bbry = 0. if bbrs_per_operator == None else bbrs_per_operator[opy]

    if not 'a' in options: 
        acc_int_x = acc_int_y = acc_bsm_x = acc_bsm_y = acc_sm_x = acc_sm_y = acc_bsm_intpart = 1.
    if not 'f' in options:
        aincx = ax
        aincy = ay
        bincx = bx
        bincy = by
        cinc  = c
    if not 'b' in options: 
        abrx = abry = bbrx = bbry = 0.


    # print 'opx:', opx
    # print 'opy:', opy
    # print 'ax:',  ax
    # print 'bx:',  bx
    # print 'ay:',  ay
    # print 'by:',  by
    # print 'c:',  c
    # print 'aincx:',  aincx
    # print 'bincx:',  bincx
    # print 'aincy:',  aincy
    # print 'bincy:',  bincy
    # print 'cinc:',  cinc
    # print 'acc_int_x:',  acc_int_x
    # print 'acc_int_y:',  acc_int_y
    # print 'acc_bsm_x:',  acc_bsm_x
    # print 'acc_bsm_y:',  acc_bsm_y
    # print 'acc_bsm_intpart:',  acc_bsm_intpart
    # print 'acc_sm_x:',  acc_sm_x
    # print 'acc_sm_y:',  acc_sm_y
    # print 'abrx:',  abrx
    # print 'abry:',  abry
    # print 'bbrx:',  bbrx
    # print 'bbry:',  bbry
    # print 'limits_r[0]:', limits_r[0]
    # print 'r_sm:', r_sm
    # raise ValueError('stop')

    # limit_x_1d_obs, limit_x_1d_exp, limit_x_1d_68low, limit_x_1d_68high, limit_x_1d_95low, limit_x_1d_95high = get_limits_at_cl(a=ax, b=bx, 'zttmm', 'comb', r_sm=r_sm, acc_int=acc_int_x, acc_sm=acc_sm_x, acc_bsm=acc_bsm_x, ainc=aincx, binc=bincx, options=options, cl=0.95)

    x_limit_obs = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[0]) 
    x_limit_obs_min, x_limit_obs_max = x_limit_obs[0], x_limit_obs[-1]
    x_limit_exp = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[1]) 
    x_limit_exp_min, x_limit_exp_max = x_limit_exp[0], x_limit_exp[-1]
    x_limit_low_68 = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[2]) 
    x_limit_low_68_min, x_limit_low_68_max = x_limit_low_68[0], x_limit_low_68[-1]
    x_limit_high_68 = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[3]) 
    x_limit_high_68_min, x_limit_high_68_max = x_limit_high_68[0], x_limit_high_68[-1]
    x_limit_low_95 = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[4]) 
    x_limit_low_95_min, x_limit_low_95_max = x_limit_low_95[0], x_limit_low_95[-1]
    x_limit_high_95 = get_limit_y_for_limit_x_comb_withacc_withf_withbr(x=0., ax=ay, bx=by, ay=ax, by=bx, c=c, aincx=aincy, bincx=bincy, aincy=aincx, bincy=bincx, cinc=cinc, acc_int_x=acc_int_y, acc_int_y=acc_int_x, acc_bsm_x=acc_bsm_y, acc_bsm_y=acc_bsm_x, acc_bsm_intpart=acc_bsm_intpart, acc_sm_x=acc_sm_y, acc_sm_y=acc_sm_y, abrx=abry, abry=abrx, bbrx=bbry, bbry=bbrx, limit=limits_r[5]) 
    x_limit_high_95_min, x_limit_high_95_max = x_limit_high_95[0], x_limit_high_95[-1]

    x_comb_obs_p     = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_obs_max*0.985, x_limit_obs_max*0.990, x_limit_obs_max*0.995, x_limit_obs_max*1.0, x_limit_obs_max*1.005, x_limit_obs_max*1.010, x_limit_obs_max*1.015] + [10, 25, 50, 75, 100, 175]
    x_comb_exp_p     = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_exp_max*0.985, x_limit_exp_max*0.990, x_limit_exp_max*0.995, x_limit_exp_max*1.0, x_limit_exp_max*1.005, x_limit_exp_max*1.010, x_limit_exp_max*1.015] + [10, 25, 50, 75, 100, 175]
    x_comb_low_68_p  = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_low_68_max*0.985, x_limit_low_68_max*0.990, x_limit_low_68_max*0.995, x_limit_low_68_max*1.0, x_limit_low_68_max*1.005, x_limit_low_68_max*1.010, x_limit_low_68_max*1.015] + [10, 25, 50, 75, 100, 175]
    x_comb_high_68_p = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_high_68_max*0.985, x_limit_high_68_max*0.990, x_limit_high_68_max*0.995, x_limit_high_68_max*1.0, x_limit_high_68_max*1.005, x_limit_high_68_max*1.010, x_limit_high_68_max*1.015] + [10, 25, 50, 75, 100, 175]
    x_comb_low_95_p  = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_low_95_max*0.985, x_limit_low_95_max*0.990, x_limit_low_95_max*0.995, x_limit_low_95_max*1.0, x_limit_low_95_max*1.005, x_limit_low_95_max*1.010, x_limit_low_95_max*1.015] + [10, 25, 50, 75, 100, 175]
    x_comb_high_95_p = slim_list(lst=[x for x in range(50000)], keep_every=250) + [x_limit_high_95_max*0.985, x_limit_high_95_max*0.990, x_limit_high_95_max*0.995, x_limit_high_95_max*1.0, x_limit_high_95_max*1.005, x_limit_high_95_max*1.010, x_limit_high_95_max*1.015] + [10, 25, 50, 75, 100, 175]

    x_comb_obs_p.sort()
    x_comb_exp_p.sort()
    x_comb_low_68_p.sort()
    x_comb_high_68_p.sort()
    x_comb_low_95_p.sort()
    x_comb_high_95_p.sort()

    x_comb_obs_n     = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_obs_min*0.985, x_limit_obs_min*0.990, x_limit_obs_min*0.995, x_limit_obs_min*1.0, x_limit_obs_min*1.005, x_limit_obs_min*1.010, x_limit_obs_min*1.015] + [-175, -100, -75, -50, -25, -10]
    x_comb_exp_n     = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_exp_min*0.985, x_limit_exp_min*0.990, x_limit_exp_min*0.995, x_limit_exp_min*1.0, x_limit_exp_min*1.005, x_limit_exp_min*1.010, x_limit_exp_min*1.015] + [-175, -100, -75, -50, -25, -10]
    x_comb_low_68_n  = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_low_68_min*0.985, x_limit_low_68_min*0.990, x_limit_low_68_min*0.995, x_limit_low_68_min*1.0, x_limit_low_68_min*1.005, x_limit_low_68_min*1.010, x_limit_low_68_min*1.015] + [-175, -100, -75, -50, -25, -10]
    x_comb_high_68_n = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_high_68_min*0.985, x_limit_high_68_min*0.990, x_limit_high_68_min*0.995, x_limit_high_68_min*1.0, x_limit_high_68_min*1.005, x_limit_high_68_min*1.010, x_limit_high_68_min*1.015] + [-175, -100, -75, -50, -25, -10]
    x_comb_low_95_n  = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_low_95_min*0.985, x_limit_low_95_min*0.990, x_limit_low_95_min*0.995, x_limit_low_95_min*1.0, x_limit_low_95_min*1.005, x_limit_low_95_min*1.010, x_limit_low_95_min*1.015] + [-175, -100, -75, -50, -25, -10]
    x_comb_high_95_n = slim_list(lst=[-x for x in range(50000)], keep_every=250) + [x_limit_high_95_min*0.985, x_limit_high_95_min*0.990, x_limit_high_95_min*0.995, x_limit_high_95_min*1.0, x_limit_high_95_min*1.005, x_limit_high_95_min*1.010, x_limit_high_95_min*1.015] + [-175, -100, -75, -50, -25, -10]

    x_comb_obs_n.sort(reverse=True)
    x_comb_exp_n.sort(reverse=True)
    x_comb_low_68_n.sort(reverse=True)
    x_comb_high_68_n.sort(reverse=True)
    x_comb_low_95_n.sort(reverse=True)
    x_comb_high_95_n.sort(reverse=True)

    # print 'ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]'
    # print ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]
    y_comb_obs_p =     [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[0]) for x in x_comb_obs_p ]
    y_comb_obs_n =     [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[0]) for x in x_comb_obs_n ]
    y_comb_exp_p =     [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[1]) for x in x_comb_exp_p ]
    y_comb_exp_n =     [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[1]) for x in x_comb_exp_n ]
    y_comb_low_68_p =  [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[2]) for x in x_comb_low_68_p ]
    y_comb_low_68_n =  [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[2]) for x in x_comb_low_68_n ]
    y_comb_high_68_p = [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[3]) for x in x_comb_high_68_p ]
    y_comb_high_68_n = [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[3]) for x in x_comb_high_68_n ]
    y_comb_low_95_p =  [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[4]) for x in x_comb_low_95_p ]
    y_comb_low_95_n =  [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[4]) for x in x_comb_low_95_n ]
    y_comb_high_95_p = [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[5]) for x in x_comb_high_95_p ]
    y_comb_high_95_n = [get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limits_r[5]) for x in x_comb_high_95_n ]

    x_int_all     = slim_list(lst=[x for x in range(-50000, 50000)], keep_every=10000)
    y_int_obs     = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[0]) for x in x_int_all]
    y_int_exp     = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[1]) for x in x_int_all]
    y_int_low_68  = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[2]) for x in x_int_all]
    y_int_high_68 = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[3]) for x in x_int_all]
    y_int_low_95  = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[4]) for x in x_int_all]
    y_int_high_95 = [get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limits_r[5]) for x in x_int_all]   


    x_comb_obs_pp = [x for (x, y)     in zip(x_comb_obs_p, y_comb_obs_p) if y is not None]
    y_comb_obs_pp = [y[-1] for (x, y) in zip(x_comb_obs_p, y_comb_obs_p) if y is not None]
    x_comb_obs_pn = [x for (x, y)     in zip(x_comb_obs_p, y_comb_obs_p) if y is not None]
    y_comb_obs_pn = [y[0] for (x, y)  in zip(x_comb_obs_p, y_comb_obs_p) if y is not None]
    x_comb_obs_np = [x for (x, y)     in zip(x_comb_obs_n, y_comb_obs_n) if y is not None]
    y_comb_obs_np = [y[-1] for (x, y) in zip(x_comb_obs_n, y_comb_obs_n) if y is not None]
    x_comb_obs_nn = [x for (x, y)     in zip(x_comb_obs_n, y_comb_obs_n) if y is not None]
    y_comb_obs_nn = [y[0] for (x, y)  in zip(x_comb_obs_n, y_comb_obs_n) if y is not None]

    x_comb_exp_pp = [x for (x, y)     in zip(x_comb_exp_p, y_comb_exp_p) if y is not None]
    y_comb_exp_pp = [y[-1] for (x, y) in zip(x_comb_exp_p, y_comb_exp_p) if y is not None]
    x_comb_exp_pn = [x for (x, y)     in zip(x_comb_exp_p, y_comb_exp_p) if y is not None]
    y_comb_exp_pn = [y[0] for (x, y)  in zip(x_comb_exp_p, y_comb_exp_p) if y is not None]
    x_comb_exp_np = [x for (x, y)     in zip(x_comb_exp_n, y_comb_exp_n) if y is not None]
    y_comb_exp_np = [y[-1] for (x, y) in zip(x_comb_exp_n, y_comb_exp_n) if y is not None]
    x_comb_exp_nn = [x for (x, y)     in zip(x_comb_exp_n, y_comb_exp_n) if y is not None]
    y_comb_exp_nn = [y[0] for (x, y)  in zip(x_comb_exp_n, y_comb_exp_n) if y is not None]

    x_comb_low_68_pp = [x for (x, y)     in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None]
    y_comb_low_68_pp = [y[-1] for (x, y) in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None]
    x_comb_low_68_pn = [x for (x, y)     in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None]
    y_comb_low_68_pn = [y[0] for (x, y)  in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None]
    x_comb_low_68_np = [x for (x, y)     in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None]
    y_comb_low_68_np = [y[-1] for (x, y) in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None]
    x_comb_low_68_nn = [x for (x, y)     in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None]
    y_comb_low_68_nn = [y[0] for (x, y)  in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None]

    x_comb_high_68_pp = [x for (x, y)     in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None]
    y_comb_high_68_pp = [y[-1] for (x, y) in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None]
    x_comb_high_68_pn = [x for (x, y)     in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None]
    y_comb_high_68_pn = [y[0] for (x, y)  in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None]
    x_comb_high_68_np = [x for (x, y)     in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None]
    y_comb_high_68_np = [y[-1] for (x, y) in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None]
    x_comb_high_68_nn = [x for (x, y)     in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None]
    y_comb_high_68_nn = [y[0] for (x, y)  in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None]

    x_comb_low_95_pp = [x for (x, y)     in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None]
    y_comb_low_95_pp = [y[-1] for (x, y) in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None]
    x_comb_low_95_pn = [x for (x, y)     in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None]
    y_comb_low_95_pn = [y[0] for (x, y)  in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None]
    x_comb_low_95_np = [x for (x, y)     in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None]
    y_comb_low_95_np = [y[-1] for (x, y) in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None]
    x_comb_low_95_nn = [x for (x, y)     in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None]
    y_comb_low_95_nn = [y[0] for (x, y)  in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None]

    x_comb_high_95_pp = [x for (x, y)     in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None]
    y_comb_high_95_pp = [y[-1] for (x, y) in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None]
    x_comb_high_95_pn = [x for (x, y)     in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None]
    y_comb_high_95_pn = [y[0] for (x, y)  in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None]
    x_comb_high_95_np = [x for (x, y)     in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None]
    y_comb_high_95_np = [y[-1] for (x, y) in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None]
    x_comb_high_95_nn = [x for (x, y)     in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None]
    y_comb_high_95_nn = [y[0] for (x, y)  in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None]

    g_comb_obs_p = ROOT.TGraph(len(y_comb_obs_pp+y_comb_obs_pn), array('d', x_comb_obs_pp+x_comb_obs_pn[::-1]), array('d', y_comb_obs_pp+y_comb_obs_pn[::-1])) 
    g_comb_obs_n = ROOT.TGraph(len(y_comb_obs_np+y_comb_obs_nn), array('d', x_comb_obs_np+x_comb_obs_nn[::-1]), array('d', y_comb_obs_np+y_comb_obs_nn[::-1])) 
    g_comb_exp_p = ROOT.TGraph(len(y_comb_exp_pp+y_comb_exp_pn), array('d', x_comb_exp_pp+x_comb_exp_pn[::-1]), array('d', y_comb_exp_pp+y_comb_exp_pn[::-1])) 
    g_comb_exp_n = ROOT.TGraph(len(y_comb_exp_np+y_comb_exp_nn), array('d', x_comb_exp_np+x_comb_exp_nn[::-1]), array('d', y_comb_exp_np+y_comb_exp_nn[::-1])) 

    if y_comb_obs_p[0] is not None and len(y_comb_obs_p[0]) >= 4:

        x_comb_obs_pp_inner = [x for (x, y)     in zip(x_comb_obs_p, y_comb_obs_p) if y is not None and len(y) >= 4]
        y_comb_obs_pp_inner = [y[-2] for (x, y) in zip(x_comb_obs_p, y_comb_obs_p) if y is not None and len(y) >= 4]
        x_comb_obs_pn_inner = [x for (x, y)     in zip(x_comb_obs_p, y_comb_obs_p) if y is not None and len(y) >= 4]
        y_comb_obs_pn_inner = [y[1] for (x, y)  in zip(x_comb_obs_p, y_comb_obs_p) if y is not None and len(y) >= 4]
        x_comb_obs_np_inner = [x for (x, y)     in zip(x_comb_obs_n, y_comb_obs_n) if y is not None and len(y) >= 4]
        y_comb_obs_np_inner = [y[-2] for (x, y) in zip(x_comb_obs_n, y_comb_obs_n) if y is not None and len(y) >= 4]
        x_comb_obs_nn_inner = [x for (x, y)     in zip(x_comb_obs_n, y_comb_obs_n) if y is not None and len(y) >= 4]
        y_comb_obs_nn_inner = [y[1] for (x, y)  in zip(x_comb_obs_n, y_comb_obs_n) if y is not None and len(y) >= 4]

        x_comb_exp_pp_inner = [x for (x, y)     in zip(x_comb_exp_p, y_comb_exp_p) if y is not None and len(y) >= 4]
        y_comb_exp_pp_inner = [y[-2] for (x, y) in zip(x_comb_exp_p, y_comb_exp_p) if y is not None and len(y) >= 4]
        x_comb_exp_pn_inner = [x for (x, y)     in zip(x_comb_exp_p, y_comb_exp_p) if y is not None and len(y) >= 4]
        y_comb_exp_pn_inner = [y[1] for (x, y)  in zip(x_comb_exp_p, y_comb_exp_p) if y is not None and len(y) >= 4]
        x_comb_exp_np_inner = [x for (x, y)     in zip(x_comb_exp_n, y_comb_exp_n) if y is not None and len(y) >= 4]
        y_comb_exp_np_inner = [y[-2] for (x, y) in zip(x_comb_exp_n, y_comb_exp_n) if y is not None and len(y) >= 4]
        x_comb_exp_nn_inner = [x for (x, y)     in zip(x_comb_exp_n, y_comb_exp_n) if y is not None and len(y) >= 4]
        y_comb_exp_nn_inner = [y[1] for (x, y)  in zip(x_comb_exp_n, y_comb_exp_n) if y is not None and len(y) >= 4]

        x_comb_low_68_pp_inner = [x for (x, y)     in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None and len(y) >= 4]
        y_comb_low_68_pp_inner = [y[-2] for (x, y) in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None and len(y) >= 4]
        x_comb_low_68_pn_inner = [x for (x, y)     in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None and len(y) >= 4]
        y_comb_low_68_pn_inner = [y[1] for (x, y)  in zip(x_comb_low_68_p, y_comb_low_68_p) if y is not None and len(y) >= 4]
        x_comb_low_68_np_inner = [x for (x, y)     in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None and len(y) >= 4]
        y_comb_low_68_np_inner = [y[-2] for (x, y) in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None and len(y) >= 4]
        x_comb_low_68_nn_inner = [x for (x, y)     in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None and len(y) >= 4]
        y_comb_low_68_nn_inner = [y[1] for (x, y)  in zip(x_comb_low_68_n, y_comb_low_68_n) if y is not None and len(y) >= 4]

        x_comb_high_68_pp_inner = [x for (x, y)     in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None and len(y) >= 4]
        y_comb_high_68_pp_inner = [y[-2] for (x, y) in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None and len(y) >= 4]
        x_comb_high_68_pn_inner = [x for (x, y)     in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None and len(y) >= 4]
        y_comb_high_68_pn_inner = [y[1] for (x, y)  in zip(x_comb_high_68_p, y_comb_high_68_p) if y is not None and len(y) >= 4]
        x_comb_high_68_np_inner = [x for (x, y)     in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None and len(y) >= 4]
        y_comb_high_68_np_inner = [y[-2] for (x, y) in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None and len(y) >= 4]
        x_comb_high_68_nn_inner = [x for (x, y)     in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None and len(y) >= 4]
        y_comb_high_68_nn_inner = [y[1] for (x, y)  in zip(x_comb_high_68_n, y_comb_high_68_n) if y is not None and len(y) >= 4]

        x_comb_low_95_pp_inner = [x for (x, y)     in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None and len(y) >= 4]
        y_comb_low_95_pp_inner = [y[-2] for (x, y) in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None and len(y) >= 4]
        x_comb_low_95_pn_inner = [x for (x, y)     in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None and len(y) >= 4]
        y_comb_low_95_pn_inner = [y[1] for (x, y)  in zip(x_comb_low_95_p, y_comb_low_95_p) if y is not None and len(y) >= 4]
        x_comb_low_95_np_inner = [x for (x, y)     in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None and len(y) >= 4]
        y_comb_low_95_np_inner = [y[-2] for (x, y) in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None and len(y) >= 4]
        x_comb_low_95_nn_inner = [x for (x, y)     in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None and len(y) >= 4]
        y_comb_low_95_nn_inner = [y[1] for (x, y)  in zip(x_comb_low_95_n, y_comb_low_95_n) if y is not None and len(y) >= 4]

        x_comb_high_95_pp_inner = [x for (x, y)     in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None and len(y) >= 4]
        y_comb_high_95_pp_inner = [y[-2] for (x, y) in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None and len(y) >= 4]
        x_comb_high_95_pn_inner = [x for (x, y)     in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None and len(y) >= 4]
        y_comb_high_95_pn_inner = [y[1] for (x, y)  in zip(x_comb_high_95_p, y_comb_high_95_p) if y is not None and len(y) >= 4]
        x_comb_high_95_np_inner = [x for (x, y)     in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None and len(y) >= 4]
        y_comb_high_95_np_inner = [y[-2] for (x, y) in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None and len(y) >= 4]
        x_comb_high_95_nn_inner = [x for (x, y)     in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None and len(y) >= 4]
        y_comb_high_95_nn_inner = [y[1] for (x, y)  in zip(x_comb_high_95_n, y_comb_high_95_n) if y is not None and len(y) >= 4]

        # redefine the graphs here, they need to be combined properly!
        # BREAK HERE! 

        g_comb_obs_p_inner = ROOT.TGraph(len(y_comb_obs_pp_inner+y_comb_obs_pn_inner), array('d', x_comb_obs_pp_inner+x_comb_obs_pn_inner[::-1]), array('d', y_comb_obs_pp_inner+y_comb_obs_pn_inner[::-1])) 
        g_comb_obs_n_inner = ROOT.TGraph(len(y_comb_obs_np_inner+y_comb_obs_nn_inner), array('d', x_comb_obs_np_inner+x_comb_obs_nn_inner[::-1]), array('d', y_comb_obs_np_inner+y_comb_obs_nn_inner[::-1])) 
        g_comb_exp_p_inner = ROOT.TGraph(len(y_comb_exp_pp_inner+y_comb_exp_pn_inner), array('d', x_comb_exp_pp_inner+x_comb_exp_pn_inner[::-1]), array('d', y_comb_exp_pp_inner+y_comb_exp_pn_inner[::-1])) 
        g_comb_exp_n_inner = ROOT.TGraph(len(y_comb_exp_np_inner+y_comb_exp_nn_inner), array('d', x_comb_exp_np_inner+x_comb_exp_nn_inner[::-1]), array('d', y_comb_exp_np_inner+y_comb_exp_nn_inner[::-1])) 

    x_comb_68_p  = x_comb_low_68_pp + x_comb_low_68_pn[::-1] + x_comb_high_68_pn + x_comb_high_68_pp[::-1]
    x_comb_68_n  = x_comb_low_68_np + x_comb_low_68_nn[::-1] + x_comb_high_68_nn + x_comb_high_68_np[::-1]
    x_comb_95_p  = x_comb_low_95_pp + x_comb_low_95_pn[::-1] + x_comb_high_95_pn + x_comb_high_95_pp[::-1]
    x_comb_95_n  = x_comb_low_95_np + x_comb_low_95_nn[::-1] + x_comb_high_95_nn + x_comb_high_95_np[::-1]

    y_comb_68_p = y_comb_low_68_pp + y_comb_low_68_pn[::-1] + y_comb_high_68_pn + y_comb_high_68_pp[::-1]
    y_comb_68_n = y_comb_low_68_np + y_comb_low_68_nn[::-1] + y_comb_high_68_nn + y_comb_high_68_np[::-1]
    y_comb_95_p = y_comb_low_95_pp + y_comb_low_95_pn[::-1] + y_comb_high_95_pn + y_comb_high_95_pp[::-1]
    y_comb_95_n = y_comb_low_95_np + y_comb_low_95_nn[::-1] + y_comb_high_95_nn + y_comb_high_95_np[::-1]

    g_comb_68_p = ROOT.TGraph(len(y_comb_68_p), array('d', x_comb_68_p), array('d', y_comb_68_p)) 
    g_comb_95_p = ROOT.TGraph(len(y_comb_95_p), array('d', x_comb_95_p), array('d', y_comb_95_p)) 
    g_comb_68_n = ROOT.TGraph(len(y_comb_68_n), array('d', x_comb_68_n), array('d', y_comb_68_n)) 
    g_comb_95_n = ROOT.TGraph(len(y_comb_95_n), array('d', x_comb_95_n), array('d', y_comb_95_n)) 

    if y_comb_obs_p[0] is not None and len(y_comb_obs_p[0]) >= 4:
        x_comb_68_p_inner  = x_comb_low_68_pp_inner + x_comb_low_68_pn_inner[::-1] + x_comb_high_68_pn_inner + x_comb_high_68_pp_inner[::-1]
        x_comb_68_n_inner  = x_comb_low_68_np_inner + x_comb_low_68_nn_inner[::-1] + x_comb_high_68_nn_inner + x_comb_high_68_np_inner[::-1]
        x_comb_95_p_inner  = x_comb_low_95_pp_inner + x_comb_low_95_pn_inner[::-1] + x_comb_high_95_pn_inner + x_comb_high_95_pp_inner[::-1]
        x_comb_95_n_inner  = x_comb_low_95_np_inner + x_comb_low_95_nn_inner[::-1] + x_comb_high_95_nn_inner + x_comb_high_95_np_inner[::-1]

        y_comb_68_p_inner = y_comb_low_68_pp_inner + y_comb_low_68_pn_inner[::-1] + y_comb_high_68_pn_inner + y_comb_high_68_pp_inner[::-1]
        y_comb_68_n_inner = y_comb_low_68_np_inner + y_comb_low_68_nn_inner[::-1] + y_comb_high_68_nn_inner + y_comb_high_68_np_inner[::-1]
        y_comb_95_p_inner = y_comb_low_95_pp_inner + y_comb_low_95_pn_inner[::-1] + y_comb_high_95_pn_inner + y_comb_high_95_pp_inner[::-1]
        y_comb_95_n_inner = y_comb_low_95_np_inner + y_comb_low_95_nn_inner[::-1] + y_comb_high_95_nn_inner + y_comb_high_95_np_inner[::-1]

        g_comb_68_p_inner = ROOT.TGraph(len(y_comb_68_p_inner), array('d', x_comb_68_p_inner), array('d', y_comb_68_p_inner)) 
        g_comb_95_p_inner = ROOT.TGraph(len(y_comb_95_p_inner), array('d', x_comb_95_p_inner), array('d', y_comb_95_p_inner)) 
        g_comb_68_n_inner = ROOT.TGraph(len(y_comb_68_n_inner), array('d', x_comb_68_n_inner), array('d', y_comb_68_n_inner)) 
        g_comb_95_n_inner = ROOT.TGraph(len(y_comb_95_n_inner), array('d', x_comb_95_n_inner), array('d', y_comb_95_n_inner)) 

    outfilename_comb = os.path.join(filefolder, 'limits_comb_%s_%s_vs_%s.root' % (options, opx, opy))
    outfile_comb = ROOT.TFile(outfilename_comb, 'RECREATE')
    g_comb_obs_p.SetName('g_obs_p')
    g_comb_obs_p.Write()
    g_comb_obs_n.SetName('g_obs_n')
    g_comb_obs_n.Write()
    g_comb_exp_p.SetName('g_exp_p')
    g_comb_exp_p.Write()
    g_comb_exp_n.SetName('g_exp_n')
    g_comb_exp_n.Write()
    g_comb_68_p.SetName('g_68_p')
    g_comb_68_p.Write()
    g_comb_95_p.SetName('g_95_p')
    g_comb_95_p.Write()
    g_comb_68_n.SetName('g_68_n')
    g_comb_68_n.Write()
    g_comb_95_n.SetName('g_95_n')
    g_comb_95_n.Write()
    if y_comb_obs_p[0] is not None and len(y_comb_obs_p[0]) >= 4:
        g_comb_68_p_inner.SetName('g_comb_68_p_inner')
        g_comb_68_p_inner.Write()
        g_comb_95_p_inner.SetName('g_comb_95_p_inner')
        g_comb_95_p_inner.Write()
        g_comb_68_n_inner.SetName('g_comb_68_n_inner')
        g_comb_68_n_inner.Write()
        g_comb_95_n_inner.SetName('g_comb_95_n_inner')
        g_comb_95_n_inner.Write()
    outfile_comb.Close()


    x_int_high_68 = copy.deepcopy(x_int_all)
    x_int_high_95 = copy.deepcopy(x_int_all)
    y_int_high_68.reverse()
    x_int_high_68.reverse()
    y_int_high_95.reverse()
    x_int_high_95.reverse()

    y_int_68 = y_int_low_68 + y_int_high_68
    y_int_95 = y_int_low_95 + y_int_high_95
    x_int_68 = x_int_all    + x_int_high_68
    x_int_95 = x_int_all    + x_int_high_95

    if not 'cle2332' in [opx, opy]:
        g_int_obs = ROOT.TGraph(len(x_int_all), array('d', x_int_all), array('d', y_int_obs)) 
        g_int_exp = ROOT.TGraph(len(x_int_all), array('d', x_int_all), array('d', y_int_exp)) 
        g_int_68  = ROOT.TGraph(len(y_int_68),  array('d', x_int_68),  array('d', y_int_68)) 
        g_int_95  = ROOT.TGraph(len(y_int_95),  array('d', x_int_95),  array('d', y_int_95)) 

        outfilename_interference = os.path.join(filefolder, 'limits_interference_%s_%s_vs_%s.root' % (options, opx, opy))
        outfile_interference = ROOT.TFile(outfilename_interference, 'RECREATE')
        g_int_obs.SetName('g_obs')
        g_int_obs.Write()
        g_int_exp.SetName('g_exp')
        g_int_exp.Write()
        g_int_68.SetName('g_68')
        g_int_68.Write()
        g_int_95.SetName('g_95')
        g_int_95.Write()
        outfile_interference.Close()


    maxdigits = 3
    xtitle = operatornames_pretty[opx]
    ytitle = operatornames_pretty[opy]
    
    marker_sm = ROOT.TMarker(0., 0., 20)

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)

    outname = os.path.join(plotfolder, 'Limits2d_summary_%s_vs_%s.pdf' % (opx, opy))
    if options == 'a': outname = outname.replace('_summary_', '_summary_withacc_')
    if options == 'f': outname = outname.replace('_summary_', '_summary_withf_')
    if options == 'af': outname = outname.replace('_summary_', '_summary_withacc_withf_')
    if options == 'afb': outname = outname.replace('_summary_', '_summary_withacc_withf_withbr_')

    # ccomb = tdrCanvas(canvName='ccomb', x_min=-1.0E4, x_max=+1.0E4, y_min=-1.0E4, y_max=+1.0E4, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    lumitag='138 fb^{-1} (13 TeV)'
    ccomb, chist = tdrCanvas(canvName='ccomb', x_min=-1.0E4, x_max=+1.0E4, y_min=-1.0E4, y_max=+1.0E4*1.5, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), return_hist=True, extraText='Supplementary', ndivisions=(505, 505))
    # ccomb, chist = tdrCanvas(canvName='ccomb', x_min=-4.0E4, x_max=+4.0E4, y_min=-2.0E3, y_max=+2.0E3*1.5, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, lumitag=lumitag, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits), return_hist=True, extraText='Supplementary', ndivisions=(505, 505))
    # chist.SetNdivisions(505)


    
    tdrDraw(g_comb_95_p, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_comb_68_p, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_comb_exp_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_comb_obs_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=303)

    tdrDraw(g_comb_95_n, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_comb_68_n, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_comb_exp_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_comb_obs_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-303)
    
    if y_comb_obs_p[0] is not None and len(y_comb_obs_p[0]) >= 4:
        tdrDraw(g_comb_95_p_inner, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        tdrDraw(g_comb_68_p_inner, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        tdrDraw(g_comb_exp_p_inner, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        tdrDraw(g_comb_obs_p_inner, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-303)

        tdrDraw(g_comb_95_n_inner, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        tdrDraw(g_comb_68_n_inner, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        tdrDraw(g_comb_exp_n_inner, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        tdrDraw(g_comb_obs_n_inner, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-303)
    
    # legcomb = tdrLeg(0.18,0.77,0.90,0.925, textSize=0.040)
    # legcomb = tdrLeg(0.185,0.70,0.89,0.85, textSize=0.040)
    # legcomb.SetNColumns(2)
    legcomb = tdrLeg(0.50,0.70,0.89,0.90, textSize=0.038)

    legcomb.SetHeader('95% CL limits')
    legcomb.AddEntry(g_comb_obs_p, 'Observed', 'L')
    legcomb.AddEntry(g_comb_exp_p, 'Median expected', 'L')
    legcomb.AddEntry(g_comb_68_p, '68% expected', 'LF')
    legcomb.AddEntry(g_comb_95_p, '95% expected', 'LF')

    # legcomb.Draw('SAME')
    marker_sm.Draw('SAME')
    smtext.DrawLatex(0., 8E2, 'SM')
    # ccomb.RedrawAxis()
    # box1.Draw()
    # box2.Draw()
    legcomb.Draw('SAME')

    # cmstag = ROOT.TLatex()
    # cmstag.SetNDC()
    # cmstag.SetTextAlign(11)
    # cmstag.SetTextColor(kBlack)
    # cmstag.SetTextSize(0.045*0.95)
    # cmstag.SetTextFont(61)
    # cmstag.DrawLatex(0.16, 0.94, 'CMS')

    # supptag = ROOT.TLatex()
    # supptag.SetNDC()
    # supptag.SetTextAlign(11)
    # supptag.SetTextColor(kBlack)
    # supptag.SetTextSize(0.040*0.95)
    # supptag.SetTextFont(52)
    # supptag.DrawLatex(0.255, 0.94, 'Supplementary')
    
    # CMS_lumi(pad=ccomb, iPosX=11, lumitag=lumitag)
    ccomb.SaveAs(outname.replace('_summary_', '_comb_'))
    ROOT.TGaxis.SetExponentOffset(0., 0., 'y')

    c = tdrCanvas(canvName='c', x_min=-1.0E4, x_max=+5.0E4, y_min=-1.0E4, y_max=+5.0E4, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))

    tdrDraw(g_comb_95_p, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_comb_68_p, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_comb_exp_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_comb_obs_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=303)

    tdrDraw(g_comb_95_n, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_comb_68_n, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_comb_exp_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_comb_obs_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-303)

    if not 'cle2332' in [opx, opy]:
        tdrDraw(g_int_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
        tdrDraw(g_int_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
        tdrDraw(g_int_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
        tdrDraw(g_int_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1, fstyle=3013, lwidth=303)

    leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    leg.SetHeader('95% CL limits')
    leg.AddEntry(g_comb_obs_p, 'Observed', 'L')
    leg.AddEntry(g_comb_exp_p, 'Median expected', 'L')
    leg.AddEntry(g_comb_68_p, '68% expected', 'LF')
    leg.AddEntry(g_comb_95_p, '95% expected', 'LF')

    leg.Draw('SAME')

    leg.Draw('SAME')
    marker_sm.Draw('SAME')
    smtext.DrawLatex(0., 8E2, 'SM')
    c.RedrawAxis()

    # outname = os.path.join(plotfolder, 'Limits2d_summary_%s_vs_%s.pdf' % (opx, opy))
    # if options == 'a': outname = outname.replace('_summary_', '_summary_withacc_')
    # if options == 'f': outname = outname.replace('_summary_', '_summary_withf_')
    # if options == 'af': outname = outname.replace('_summary_', '_summary_withacc_withf_')
    c.SaveAs(outname)

    if 'cle2332' in [opx, opy]: return
    cint = tdrCanvas(canvName='cint', x_min=-5.0E4, x_max=+5.0E4, y_min=-5.0E4, y_max=+5.0E4, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))

    tdrDraw(g_int_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_int_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_int_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_int_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1, fstyle=3013, lwidth=303)

    leg.Draw('SAME')
    marker_sm.Draw('SAME')
    smtext.DrawLatex(0., 8E2, 'SM')
    cint.RedrawAxis()
    cint.SaveAs(outname.replace('_summary_', '_interference_'))

def plot_limits_2d_comb(as_per_operator, bs_per_operator, proc_affected_per_operator, opx, opy, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, options=''):

    # # slightly less easy but still OK analytic relation for limits
    ax = as_per_operator[opx]
    ay = as_per_operator[opy]
    bx = bs_per_operator[opx]
    by = bs_per_operator[opy]
    operatorname_sum = '%s+1p0x%s' % (opx, opy)
    if not operatorname_sum in bs_per_operator:
        c = 0
    else:
        sumcoeff  = bs_per_operator[operatorname_sum]
        c = sumcoeff - bx - by
        # print 'x coefficient bx =', bx
        # print 'y coefficient by =', by
        # print 'interference coefficient c =', c
        
        if options == 'a' or options == 'af':
            # also get mixed acceptance in a similar way
            acc_int_x = acceptances_int_per_operator[opx]
            acc_int_y = acceptances_int_per_operator[opy]
            acc_bsm_x = acceptances_bsm_per_operator[opx]
            acc_bsm_y = acceptances_bsm_per_operator[opy]
            acc_sm_x = acceptances_sm_per_operator[opx]
            acc_sm_y = acceptances_sm_per_operator[opy]
            if operatorname_sum in acceptances_bsm_per_operator:
                acc_bsm_withint = acceptances_bsm_per_operator[operatorname_sum]
            elif '%s+1p0x%s' % (opy, opx) in acceptances_bsm_per_operator: # some operatornames are swapped, but for their dual insertion it doesn't matter
                acc_bsm_withint = acceptances_bsm_per_operator['%s+1p0x%s' % (opy, opx)]
            else: 
                raise ValueError('Trying to find combined operator %s (or swapped: %s) in list of acceptances, but it is not in. Need to skip this combination of operators?' % (operatorname_sum, '%s+1p0x%s' % (opy, opx)))
            acc_bsm_intpart = get_acc_bsm_intpart(bx=bx, by=by, c=c, accx=acc_bsm_x, accy=acc_bsm_y, acctot=acc_bsm_withint)

    if options == 'f' or options == 'af':
        if None in [aincs_per_operator, bincs_per_operator]: raise ValueError('When running plot_limits_2d_comb with option \'f\', the aincs and bincs must not be None.')
        aincx = aincs_per_operator[opx]
        aincy = aincs_per_operator[opy]
        bincx = bincs_per_operator[opx]
        bincy = bincs_per_operator[opy]
        if not operatorname_sum in bincs_per_operator:
            cinc = 0
        else:
            sumcoeffinc  = bincs_per_operator[operatorname_sum]
            cinc = sumcoeffinc - bincx - bincy

            print 'x coefficient bincx =', bincx
            print 'y coefficient bincy =', bincy
            print 'interference coefficient cinc =', cinc

    # c *= 1000

    # if there is an interference term, the solution to the master equation is (according to Wolfram):
    # y = -( (-math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) )

    def get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limit):
        result_plus  = -(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limit / r_sm - 1 - (ax * x + bx * x**2))/by)
        result_minus = -(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limit / r_sm - 1 - (ax * x + bx * x**2))/by)
        return (result_plus, result_minus)

    def get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limit):
        # for abc in [ax, bx, ay, by, c]: print abc
        result_plus  = -( (+math.sqrt(4*by*((limit / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) )
        result_minus = -( (-math.sqrt(4*by*((limit / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) )
        return (result_plus, result_minus)

    def get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limit):
        result_plus  = -(      (math.sqrt((4*by*(ax*x + bx*x**2 - limit/r_sm + 1) + ay**2*limit/r_sm)/(limit/r_sm)) + ay)/ (2*by) )
        result_minus = -( ay - (math.sqrt((4*by*(ax*x + bx*x**2 - limit/r_sm + 1) + ay**2*limit/r_sm)/(limit/r_sm)))     / (2*by) )
        return (result_plus, result_minus)

    # def get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limit):
    #     r = limit/r_sm
    #     acrx = acc_int_x / acc_sm_x
    #     acry = acc_int_y / acc_sm_y
    #     result = (math.sqrt(4*ax*acrx*r*x - 4*acry*r*(ax*x+1) + acry**2*r**2 + 4*r) - 2*ax*x + acry*r - 2)/(2*ay)
    #     return result

    def get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limit):
        if not acc_sm_x == acc_sm_y: raise ValueError('SM acceptances for the same process are different.')
        acc_sm = acc_sm_x

        r = limit/r_sm
        acrintx = acc_int_x / acc_sm
        acrinty = acc_int_y / acc_sm
        acrbsmx = acc_bsm_x / acc_sm
        acrbsmy = acc_bsm_y / acc_sm
        acrbsmintpart = acc_bsm_intpart / acc_sm

        equation = sp.Eq( (1+ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar) / ((ax*x)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrintx + (ay*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrinty + (bx*x**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmx + (by*yvar**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmy + (c*x*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmintpart), r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
        # if result != []: print 'solutions for x = %f:'%(x), result
        if len(result) == 0:
            # print 'No real solution, returning None'
            return None
        return result

    def get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limit):

        r = limit/r_sm

        equation = sp.Eq( (1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)**2/(1 + aincx*x + aincy*yvar + bincx*x**2 + bincy*yvar**2 + cinc*x*yvar), r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
        if len(result) == 0:
            # print 'No real solution, returning None'
            return None
        return result

    def get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limit):
        if not acc_sm_x == acc_sm_y: raise ValueError('SM acceptances for the same process are different.')
        acc_sm = acc_sm_x

        r = limit/r_sm
        acrintx = acc_int_x / acc_sm
        acrinty = acc_int_y / acc_sm
        acrbsmx = acc_bsm_x / acc_sm
        acrbsmy = acc_bsm_y / acc_sm
        acrbsmintpart = acc_bsm_intpart / acc_sm

        equation = sp.Eq( (1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)**2/(1 + aincx*x + aincy*yvar + bincx*x**2 + bincy*yvar**2 + cinc*x*yvar) / ((ax*x)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrintx + (ay*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrinty + (bx*x**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmx + (by*yvar**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmy + (c*x*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmintpart), r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
        if len(result) == 0:
            # print 'No real solution, returning None'
            return None
        return result
    
    limits_r = get_limits_r_at_cl(cl=0.95)

    if proc_affected_per_operator[opx] == 'zmmmm' and proc_affected_per_operator[opy] == 'zmmmm':
        raise ValueError('For 2d comb limits, both operators affect the 4mu process. We do not do this...')
    if proc_affected_per_operator[opx] == 'zmmmm':
        raise ValueError('For 2d comb limits, the x-operator affect the 4mu process. Use our convention and make it the y-operator')


    if options == '':
        if proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm': # easy case, same as before
    
            if c == 0: raise ValueError('Somehow there is no interference between the two operators in the quadratic case that affects the zttmm both times... there should...')
            x_obs_p     = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
            x_exp_p     = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
            x_low_68_p  = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
            x_high_68_p = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
            x_low_95_p  = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
            x_high_95_p = slim_list(lst=[x for x in range(50000) if 4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
    
            x_obs_n     = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[0] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
            x_exp_n     = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[1] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
            x_low_68_n  = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[2] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
            x_high_68_n = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[3] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
            x_low_95_n  = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[4] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
            x_high_95_n = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[5] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*x**2 > 0], keep_every=100)
    
            y_obs_pp     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[0])[0] for x in x_obs_p     ]
            y_exp_pp     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[1])[0] for x in x_exp_p     ]
            y_low_68_pp  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[2])[0] for x in x_low_68_p  ]
            y_high_68_pp = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[3])[0] for x in x_high_68_p ]
            y_low_95_pp  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[4])[0] for x in x_low_95_p  ]
            y_high_95_pp = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[5])[0] for x in x_high_95_p ]
    
            y_obs_pn     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[0])[1] for x in x_obs_p     ]
            y_exp_pn     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[1])[1] for x in x_exp_p     ]
            y_low_68_pn  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[2])[1] for x in x_low_68_p  ]
            y_high_68_pn = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[3])[1] for x in x_high_68_p ]
            y_low_95_pn  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[4])[1] for x in x_low_95_p  ]
            y_high_95_pn = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[5])[1] for x in x_high_95_p ]
    
            y_obs_np     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[0])[0] for x in x_obs_n     ]
            y_exp_np     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[1])[0] for x in x_exp_n     ]
            y_low_68_np  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[2])[0] for x in x_low_68_n  ]
            y_high_68_np = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[3])[0] for x in x_high_68_n ]
            y_low_95_np  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[4])[0] for x in x_low_95_n  ]
            y_high_95_np = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[5])[0] for x in x_high_95_n ]
    
            y_obs_nn     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[0])[1] for x in x_obs_n     ]
            y_exp_nn     = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[1])[1] for x in x_exp_n     ]
            y_low_68_nn  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[2])[1] for x in x_low_68_n  ]
            y_high_68_nn = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[3])[1] for x in x_high_68_n ]
            y_low_95_nn  = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[4])[1] for x in x_low_95_n  ]
            y_high_95_nn = [get_limit_y_for_limit_x_numnum_withint(x, ax, bx, ay, by, c, limits_r[5])[1] for x in x_high_95_n ] 

    elif options == 'a':
        if proc_affected_per_operator[opx] == 'zmmmm' or proc_affected_per_operator[opy] == 'zmmmm': raise ValueError('For 2d comb limits with acceptance, one of the two operators affects zmmmm, we do not do that...')
        if c == 0: raise ValueError('Somehow there is no interference between the two operators in the quadratic case that affects the zttmm both times... there should...')

        if None in [acceptances_int_per_operator, acceptances_bsm_per_operator, acceptances_sm_per_operator]: raise ValueError('When running plot_limits_2d_comb with option \'a\', the acceptances must none be None.')
         
        x_obs_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_exp_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_68_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_68_p = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_95_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_95_p = slim_list(lst=[x for x in range(20000)], keep_every=100)   

        x_obs_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_exp_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_68_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_68_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_95_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_95_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)    

        y_obs_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_obs_p ]
        y_obs_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_obs_n ]
        y_exp_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_exp_p ]
        y_exp_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_exp_n ]
        y_low_68_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_low_68_p ]
        y_low_68_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_low_68_n ]
        y_high_68_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_high_68_p ]
        y_high_68_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_high_68_n ]
        y_low_95_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_low_95_p ]
        y_low_95_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_low_95_n ]
        y_high_95_p = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_high_95_p ]
        y_high_95_n = [get_limit_y_for_limit_x_numnum_withint_withacc(x, ax, bx, ay, by, c, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_high_95_n ]


    elif options == 'f':
        if proc_affected_per_operator[opx] == 'zmmmm' or proc_affected_per_operator[opy] == 'zmmmm': raise ValueError('For 2d comb limits with f, one of the two operators affects zmmmm, we do not do that...')
        if c == 0 or cinc == 0: raise ValueError('Somehow there is no interference between the two operators in the quadratic case that affects the ttmm both times (either the zttmm nor the incttmm process) there should...')

        if None in [aincs_per_operator, bincs_per_operator]: raise ValueError('When running plot_limits_2d_comb with option \'f\', the aincs and bincs must none be None.')
         
        x_obs_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_exp_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_68_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_68_p = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_95_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_95_p = slim_list(lst=[x for x in range(20000)], keep_every=100)   

        x_obs_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_exp_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_68_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_68_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_95_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_95_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)    

        y_obs_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[0]) for x in x_obs_p ]
        y_obs_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[0]) for x in x_obs_n ]
        y_exp_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[1]) for x in x_exp_p ]
        y_exp_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[1]) for x in x_exp_n ]
        y_low_68_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[2]) for x in x_low_68_p ]
        y_low_68_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[2]) for x in x_low_68_n ]
        y_high_68_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[3]) for x in x_high_68_p ]
        y_high_68_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[3]) for x in x_high_68_n ]
        y_low_95_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[4]) for x in x_low_95_p ]
        y_low_95_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[4]) for x in x_low_95_n ]
        y_high_95_p = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[5]) for x in x_high_95_p ]
        y_high_95_n = [get_limit_y_for_limit_x_numnum_withint_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, limits_r[5]) for x in x_high_95_n ]


    elif options == 'af':
        if proc_affected_per_operator[opx] == 'zmmmm' or proc_affected_per_operator[opy] == 'zmmmm': raise ValueError('For 2d comb limits with f, one of the two operators affects zmmmm, we do not do that...')
        if c == 0 or cinc == 0: raise ValueError('Somehow there is no interference between the two operators in the quadratic case that affects the ttmm both times (either the zttmm nor the incttmm process) there should...')

        if None in [aincs_per_operator, bincs_per_operator, acceptances_int_per_operator, acceptances_bsm_per_operator, acceptances_sm_per_operator]: raise ValueError('When running plot_limits_2d_comb with option \'af\', the aincs and bincs and all acceptances must not be None.')
         
        x_obs_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_exp_p     = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_68_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_68_p = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_low_95_p  = slim_list(lst=[x for x in range(20000)], keep_every=100)
        x_high_95_p = slim_list(lst=[x for x in range(20000)], keep_every=100)   

        x_obs_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_exp_n     = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_68_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_68_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_low_95_n  = slim_list(lst=[-x for x in range(20000)], keep_every=100)
        x_high_95_n = slim_list(lst=[-x for x in range(20000)], keep_every=100)    

        y_obs_p =     [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_obs_p ]
        y_obs_n =     [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_obs_n ]
        y_exp_p =     [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_exp_p ]
        y_exp_n =     [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_exp_n ]
        y_low_68_p =  [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_low_68_p ]
        y_low_68_n =  [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_low_68_n ]
        y_high_68_p = [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_high_68_p ]
        y_high_68_n = [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_high_68_n ]
        y_low_95_p =  [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_low_95_p ]
        y_low_95_n =  [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_low_95_n ]
        y_high_95_p = [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_high_95_p ]
        y_high_95_n = [get_limit_y_for_limit_x_numnum_withint_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_high_95_n ]


 

    lwidth_base = 303
    if '2233' in opy or True:
        if options == '':
            lwidth = -lwidth_base

            x_high_68_p.reverse()
            x_high_95_p.reverse()   
            x_high_68_n.reverse()
            x_high_95_n.reverse() 

            y_high_68_pp.reverse()
            y_high_95_pp.reverse()
            y_high_68_pn.reverse()
            y_high_95_pn.reverse()
            y_high_68_np.reverse()
            y_high_95_np.reverse()
            y_high_68_nn.reverse()
            y_high_95_nn.reverse()

            x_68_p  = x_low_68_p + x_low_68_p[::-1] + x_high_68_p[::-1] + x_high_68_p
            x_95_p  = x_low_95_p + x_low_95_p[::-1] + x_high_95_p[::-1] + x_high_95_p
            x_68_n  = x_low_68_n + x_low_68_n[::-1] + x_high_68_n[::-1] + x_high_68_n
            x_95_n  = x_low_95_n + x_low_95_n[::-1] + x_high_95_n[::-1] + x_high_95_n

            y_68_p = y_low_68_pp + y_low_68_pn[::-1] + y_high_68_pn[::-1] + y_high_68_pp
            y_95_p = y_low_95_pp + y_low_95_pn[::-1] + y_high_95_pn[::-1] + y_high_95_pp
            y_68_n = y_low_68_np + y_low_68_nn[::-1] + y_high_68_nn[::-1] + y_high_68_np
            y_95_n = y_low_95_np + y_low_95_nn[::-1] + y_high_95_nn[::-1] + y_high_95_np

            g_obs_p = ROOT.TGraph(len(y_obs_pp+y_obs_pn), array('d', x_obs_p+x_obs_p[::-1]), array('d', y_obs_pp+y_obs_pn[::-1])) 
            g_exp_p = ROOT.TGraph(len(y_exp_pp+y_exp_pn), array('d', x_exp_p+x_exp_p[::-1]), array('d', y_exp_pp+y_exp_pn[::-1])) 
            g_obs_n = ROOT.TGraph(len(y_obs_np+y_obs_nn), array('d', x_obs_n+x_obs_n[::-1]), array('d', y_obs_np+y_obs_nn[::-1])) 
            g_exp_n = ROOT.TGraph(len(y_exp_np+y_exp_nn), array('d', x_exp_n+x_exp_n[::-1]), array('d', y_exp_np+y_exp_nn[::-1])) 
        
        elif options == 'a' or options == 'f' or options == 'af':
            lwidth = +lwidth_base

            x_obs_pp = [x for (x, y)     in zip(x_obs_p, y_obs_p) if y is not None]
            y_obs_pp = [y[-1] for (x, y) in zip(x_obs_p, y_obs_p) if y is not None]
            x_obs_pn = [x for (x, y)     in zip(x_obs_p, y_obs_p) if y is not None]
            y_obs_pn = [y[0] for (x, y)  in zip(x_obs_p, y_obs_p) if y is not None]
            x_obs_np = [x for (x, y)     in zip(x_obs_n, y_obs_n) if y is not None]
            y_obs_np = [y[-1] for (x, y) in zip(x_obs_n, y_obs_n) if y is not None]
            x_obs_nn = [x for (x, y)     in zip(x_obs_n, y_obs_n) if y is not None]
            y_obs_nn = [y[0] for (x, y)  in zip(x_obs_n, y_obs_n) if y is not None]

            x_exp_pp = [x for (x, y)     in zip(x_exp_p, y_exp_p) if y is not None]
            y_exp_pp = [y[-1] for (x, y) in zip(x_exp_p, y_exp_p) if y is not None]
            x_exp_pn = [x for (x, y)     in zip(x_exp_p, y_exp_p) if y is not None]
            y_exp_pn = [y[0] for (x, y)  in zip(x_exp_p, y_exp_p) if y is not None]
            x_exp_np = [x for (x, y)     in zip(x_exp_n, y_exp_n) if y is not None]
            y_exp_np = [y[-1] for (x, y) in zip(x_exp_n, y_exp_n) if y is not None]
            x_exp_nn = [x for (x, y)     in zip(x_exp_n, y_exp_n) if y is not None]
            y_exp_nn = [y[0] for (x, y)  in zip(x_exp_n, y_exp_n) if y is not None]

            x_low_68_pp = [x for (x, y)     in zip(x_low_68_p, y_low_68_p) if y is not None]
            y_low_68_pp = [y[-1] for (x, y) in zip(x_low_68_p, y_low_68_p) if y is not None]
            x_low_68_pn = [x for (x, y)     in zip(x_low_68_p, y_low_68_p) if y is not None]
            y_low_68_pn = [y[0] for (x, y)  in zip(x_low_68_p, y_low_68_p) if y is not None]
            x_low_68_np = [x for (x, y)     in zip(x_low_68_n, y_low_68_n) if y is not None]
            y_low_68_np = [y[-1] for (x, y) in zip(x_low_68_n, y_low_68_n) if y is not None]
            x_low_68_nn = [x for (x, y)     in zip(x_low_68_n, y_low_68_n) if y is not None]
            y_low_68_nn = [y[0] for (x, y)  in zip(x_low_68_n, y_low_68_n) if y is not None]

            x_high_68_pp = [x for (x, y)     in zip(x_high_68_p, y_high_68_p) if y is not None]
            y_high_68_pp = [y[-1] for (x, y) in zip(x_high_68_p, y_high_68_p) if y is not None]
            x_high_68_pn = [x for (x, y)     in zip(x_high_68_p, y_high_68_p) if y is not None]
            y_high_68_pn = [y[0] for (x, y)  in zip(x_high_68_p, y_high_68_p) if y is not None]
            x_high_68_np = [x for (x, y)     in zip(x_high_68_n, y_high_68_n) if y is not None]
            y_high_68_np = [y[-1] for (x, y) in zip(x_high_68_n, y_high_68_n) if y is not None]
            x_high_68_nn = [x for (x, y)     in zip(x_high_68_n, y_high_68_n) if y is not None]
            y_high_68_nn = [y[0] for (x, y)  in zip(x_high_68_n, y_high_68_n) if y is not None]

            x_low_95_pp = [x for (x, y)     in zip(x_low_95_p, y_low_95_p) if y is not None]
            y_low_95_pp = [y[-1] for (x, y) in zip(x_low_95_p, y_low_95_p) if y is not None]
            x_low_95_pn = [x for (x, y)     in zip(x_low_95_p, y_low_95_p) if y is not None]
            y_low_95_pn = [y[0] for (x, y)  in zip(x_low_95_p, y_low_95_p) if y is not None]
            x_low_95_np = [x for (x, y)     in zip(x_low_95_n, y_low_95_n) if y is not None]
            y_low_95_np = [y[-1] for (x, y) in zip(x_low_95_n, y_low_95_n) if y is not None]
            x_low_95_nn = [x for (x, y)     in zip(x_low_95_n, y_low_95_n) if y is not None]
            y_low_95_nn = [y[0] for (x, y)  in zip(x_low_95_n, y_low_95_n) if y is not None]

            x_high_95_pp = [x for (x, y)     in zip(x_high_95_p, y_high_95_p) if y is not None]
            y_high_95_pp = [y[-1] for (x, y) in zip(x_high_95_p, y_high_95_p) if y is not None]
            x_high_95_pn = [x for (x, y)     in zip(x_high_95_p, y_high_95_p) if y is not None]
            y_high_95_pn = [y[0] for (x, y)  in zip(x_high_95_p, y_high_95_p) if y is not None]
            x_high_95_np = [x for (x, y)     in zip(x_high_95_n, y_high_95_n) if y is not None]
            y_high_95_np = [y[-1] for (x, y) in zip(x_high_95_n, y_high_95_n) if y is not None]
            x_high_95_nn = [x for (x, y)     in zip(x_high_95_n, y_high_95_n) if y is not None]
            y_high_95_nn = [y[0] for (x, y)  in zip(x_high_95_n, y_high_95_n) if y is not None]

            g_obs_p = ROOT.TGraph(len(y_obs_pp+y_obs_pn), array('d', x_obs_pp+x_obs_pn[::-1]), array('d', y_obs_pp+y_obs_pn[::-1])) 
            g_obs_n = ROOT.TGraph(len(y_obs_np+y_obs_nn), array('d', x_obs_np+x_obs_nn[::-1]), array('d', y_obs_np+y_obs_nn[::-1])) 
            g_exp_p = ROOT.TGraph(len(y_exp_pp+y_exp_pn), array('d', x_exp_pp+x_exp_pn[::-1]), array('d', y_exp_pp+y_exp_pn[::-1])) 
            g_exp_n = ROOT.TGraph(len(y_exp_np+y_exp_nn), array('d', x_exp_np+x_exp_nn[::-1]), array('d', y_exp_np+y_exp_nn[::-1])) 

            x_68_p  = x_low_68_pp + x_low_68_pn[::-1] + x_high_68_pn + x_high_68_pp[::-1]
            x_68_n  = x_low_68_np + x_low_68_nn[::-1] + x_high_68_nn + x_high_68_np[::-1]
            x_95_p  = x_low_95_pp + x_low_95_pn[::-1] + x_high_95_pn + x_high_95_pp[::-1]
            x_95_n  = x_low_95_np + x_low_95_nn[::-1] + x_high_95_nn + x_high_95_np[::-1]

            y_68_p = y_low_68_pp + y_low_68_pn[::-1] + y_high_68_pn + y_high_68_pp[::-1]
            y_68_n = y_low_68_np + y_low_68_nn[::-1] + y_high_68_nn + y_high_68_np[::-1]
            y_95_p = y_low_95_pp + y_low_95_pn[::-1] + y_high_95_pn + y_high_95_pp[::-1]
            y_95_n = y_low_95_np + y_low_95_nn[::-1] + y_high_95_nn + y_high_95_np[::-1]
        else:
            raise ValueError('Illegal option(s).')

        

    # elif '2222' in opy:
    #     y_high_68_pp.reverse()
    #     y_high_95_pp.reverse()
    #     y_high_68_pn.reverse()
    #     y_high_95_pn.reverse()
    #     y_high_68_np.reverse()
    #     y_high_95_np.reverse()
    #     y_high_68_nn.reverse()
    #     y_high_95_nn.reverse()

    #     x_68_p  = x_low_68_p[::-1] + x_low_68_p + x_high_68_p + x_high_68_p[::-1]
    #     x_95_p  = x_low_95_p[::-1] + x_low_95_p + x_high_95_p + x_high_95_p[::-1]
    #     x_68_n  = x_low_68_n[::-1] + x_low_68_n + x_high_68_n + x_high_68_n[::-1]
    #     x_95_n  = x_low_95_n[::-1] + x_low_95_n + x_high_95_n + x_high_95_n[::-1]

    #     y_68_p = y_low_68_pp[::-1] + y_low_68_pn + y_high_68_pn + y_high_68_pp[::-1]
    #     y_95_p = y_low_95_pp[::-1] + y_low_95_pn + y_high_95_pn + y_high_95_pp[::-1]
    #     y_68_n = y_low_68_np[::-1] + y_low_68_nn + y_high_68_nn + y_high_68_np[::-1]
    #     y_95_n = y_low_95_np[::-1] + y_low_95_nn + y_high_95_nn + y_high_95_np[::-1]

    #     g_obs_p = ROOT.TGraph(len(y_obs_pp+y_obs_pn), array('d', x_obs_p[::-1]+x_obs_p), array('d', y_obs_pp[::-1]+y_obs_pn)) 
    #     g_exp_p = ROOT.TGraph(len(y_exp_pp+y_exp_pn), array('d', x_exp_p[::-1]+x_exp_p), array('d', y_exp_pp[::-1]+y_exp_pn)) 
    #     g_obs_n = ROOT.TGraph(len(y_obs_np+y_obs_nn), array('d', x_obs_n[::-1]+x_obs_n), array('d', y_obs_np[::-1]+y_obs_nn)) 
    #     g_exp_n = ROOT.TGraph(len(y_exp_np+y_exp_nn), array('d', x_exp_n[::-1]+x_exp_n), array('d', y_exp_np[::-1]+y_exp_nn)) 

    #     lwidth = -lwidth_base
    # else:
    #     raise ValueError('some operator other than 2222 or 2233 is supposed to be plotted on the y axis, not supported.')

    g_68_p = ROOT.TGraph(len(y_68_p), array('d', x_68_p), array('d', y_68_p)) 
    g_95_p = ROOT.TGraph(len(y_95_p), array('d', x_95_p), array('d', y_95_p)) 
    g_68_n = ROOT.TGraph(len(y_68_n), array('d', x_68_n), array('d', y_68_n)) 
    g_95_n = ROOT.TGraph(len(y_95_n), array('d', x_95_n), array('d', y_95_n)) 


    maxdigits = 3
    xtitle = operatornames_pretty[opx]
    ytitle = operatornames_pretty[opy]
    xmin = -2.0E4 
    xmax = +2.0E4
    ymin, ymax = xmin, xmax
    if not (proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm'):
        ymax /= 4.
        ymin /= 4.
    c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=ymin, y_max=ymax, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    
    if '2233' in opy:
        leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    else:
        leg = tdrLeg(0.37,0.68,0.72,0.9, textSize=0.038)
    leg.SetHeader('95% CL limits')
    tdrDraw(g_95_p, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_68_p, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_exp_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_obs_p, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=lwidth)

    tdrDraw(g_95_n, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_68_n, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_exp_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_obs_n, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, fstyle=3013, lstyle=1, lwidth=-lwidth)

    leg.AddEntry(g_obs_p, 'Observed', 'L')
    leg.AddEntry(g_exp_p, 'Median expected', 'L')
    leg.AddEntry(g_68_p, '68% expected', 'LF')
    leg.AddEntry(g_95_p, '95% expected', 'LF')

    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)
    smtext.DrawLatex(0., 8E2, 'SM')

    c.RedrawAxis()
    
    outname = os.path.join(plotfolder, 'Limits2d_%s_%s_vs_%s.pdf' % ('comb', opx, opy))
    if options == 'a': outname = outname.replace('comb', 'comb_withacc')
    if options == 'f': outname = outname.replace('comb', 'comb_withf')
    if options == 'af': outname = outname.replace('comb', 'comb_withacc_withf')
    c.SaveAs(outname)

    del c


def slim_list(lst, keep_every):
    if len(lst) < 2: raise ValueError('Given list is does not have 2 or more entries, just take it as is.')
    if len(lst) == 2:
        return lst

    result = [lst[0]]

    for idx in range(1, len(lst)-1):
        if idx % keep_every == 0: result.append(lst[idx])
    result.append(lst[-1])
    return result



def get_limits_r_at_cl(cl=0.95):
    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r(r_sm=r_sm)
    return (interpolate_x_from_y(value_list=r_obs, target_y=1.-cl), interpolate_x_from_y(value_list=r_exp, target_y=1.-cl), interpolate_x_from_y(value_list=r_68_low, target_y=1.-cl), interpolate_x_from_y(value_list=r_68_high, target_y=1.-cl), interpolate_x_from_y(value_list=r_95_low, target_y=1.-cl), interpolate_x_from_y(value_list=r_95_high, target_y=1.-cl))


def plot_limits_2d_interference(as_per_operator, proc_affected_per_operator, opx, opy, acceptances_int_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, options=''):
    limits_r = get_limits_r_at_cl(cl=0.95)
    

    def get_limit_y_for_limit_x_numnum(x, ax, ay, limit):
        result = (limit/r_sm - 1 - ax*x)/ay
        return result
    
    def get_limit_y_for_limit_x_numden(x, ax, ay, limit):
        result = ((1+ax*x)/(limit/r_sm)-1)/ay
        return result

    def get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limit):
        r = limit/r_sm
        acrx = acc_int_x / acc_sm_x
        acry = acc_int_y / acc_sm_y
        result = (math.sqrt(4*ax*acrx*r*x - 4*acry*r*(ax*x+1) + acry**2*r**2 + 4*r) - 2*ax*x + acry*r - 2)/(2*ay)
        return result

    def get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limit):
        r = limit/r_sm
        equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar), r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
        if len(result) == 0:
            return None
        return result[-1]

    def get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limit):
        r = limit/r_sm
        acrx = acc_int_x / acc_sm_x
        acry = acc_int_y / acc_sm_y
        equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar) / (ax*x/(1+ax*x+ay*yvar)*acrx + ay*yvar/(1+ax*x+ay*yvar)*acry + 1 - (ax*x+ay*yvar)/(1+ax*x+ay*yvar)), r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
        print x, ':  ', result
        if len(result) == 0:
            return None
        return result[-1]

    # easy analytic relation for limits, they will be a straight line
    ax = as_per_operator[opx]
    ay = as_per_operator[opy]
    if 0 in [ax, ay]:
        return 

    # make graphs for exp, obs, 68, 95 and fill the area between (like also done before)
    npoints = 5
    xmin = int(-5E4)
    xmax = int(5E4)

    # need to distinguish between cases that both affect the tau rate, and those that affect once the tau rate and once the mu rate (denominator)
    if proc_affected_per_operator[opx] == 'zmmmm' and proc_affected_per_operator[opy] == 'zmmmm':
        raise ValueError('For 2d interference limits, both operators affect the 4mu process. We do not do this...')
    if proc_affected_per_operator[opx] == 'zmmmm':
        raise ValueError('For 2d interference limits, the x-operator affect the 4mu process. Use our convention and make it the y-operator')
    if options == '':
        if proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm': # easy case, same as before
            x_all     = slim_list(lst=[x for x in range(xmin, xmax)], keep_every=100)
            y_obs     = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[0]) for x in x_all]
            y_exp     = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[1]) for x in x_all]
            y_low_68  = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[2]) for x in x_all]
            y_high_68 = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[3]) for x in x_all]
            y_low_95  = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[4]) for x in x_all]
            y_high_95 = [get_limit_y_for_limit_x_numnum(x, ax, ay, limits_r[5]) for x in x_all]    
        else: # opx affects the numerator and opy the denominator
            x_all     = slim_list(lst=[x for x in range(xmin, xmax)], keep_every=100)
            y_obs     = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[0]) for x in x_all]
            y_exp     = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[1]) for x in x_all]
            y_low_68  = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[2]) for x in x_all]
            y_high_68 = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[3]) for x in x_all]
            y_low_95  = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[4]) for x in x_all]
            y_high_95 = [get_limit_y_for_limit_x_numden(x, ax, ay, limits_r[5]) for x in x_all]    
    elif options == 'f':
        if aincs_per_operator == None: raise ValueError('When running plot_limits_2d_interference with option \'f\', the aincs must not be None.')
        aincx = aincs_per_operator[opx]
        aincy = aincs_per_operator[opy]

        if proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm': # easy case, same as before
            x_all     = slim_list(lst=[x for x in range(xmin, xmax)], keep_every=100)
            y_obs     = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[0]) for x in x_all]
            y_exp     = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[1]) for x in x_all]
            y_low_68  = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[2]) for x in x_all]
            y_high_68 = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[3]) for x in x_all]
            y_low_95  = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[4]) for x in x_all]
            y_high_95 = [get_limit_y_for_limit_x_numnum_withf(x, ax, ay, aincx, aincy, limits_r[5]) for x in x_all]    
        else:
            raise ValueError('In 2-d limits (interference) with acceptance (option \'f\'), we only consider operators that affect zttmm, but this one affects zmmmm. Nah.')
    elif options == 'a':
        if None in [acceptances_int_per_operator, acceptances_sm_per_operator]: raise ValueError('When running plot_limits_2d_interference with option \'a\', the acceptances must none be None.')
        acc_int_x = acceptances_int_per_operator[opx]
        acc_int_y = acceptances_int_per_operator[opy]
        acc_sm_x = acceptances_sm_per_operator[opx]
        acc_sm_y = acceptances_sm_per_operator[opy]

        if proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm': # easy case, same as before
            x_all     = slim_list(lst=[x for x in range(xmin, xmax)], keep_every=100)
            y_obs     = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_all]
            y_exp     = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_all]
            y_low_68  = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_all]
            y_high_68 = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_all]
            y_low_95  = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_all]
            y_high_95 = [get_limit_y_for_limit_x_numnum_withacc(x, ax, ay, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_all]    
        else:
            raise ValueError('In 2-d limits (interference) with acceptance (option \'a\'), we only consider operators that affect zttmm, but this one affects zmmmm. Nah.')
    elif options == 'af':
        if None in [acceptances_int_per_operator, acceptances_sm_per_operator, aincs_per_operator]: raise ValueError('When running plot_limits_2d_interference with option \'a\', the acceptances and aincs must none be None.')
        acc_int_x = acceptances_int_per_operator[opx]
        acc_int_y = acceptances_int_per_operator[opy]
        acc_sm_x = acceptances_sm_per_operator[opx]
        acc_sm_y = acceptances_sm_per_operator[opy]
        aincx = aincs_per_operator[opx]
        aincy = aincs_per_operator[opy]

        if proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm': # easy case, same as before
            x_all     = slim_list(lst=[x for x in range(xmin, xmax)], keep_every=10000)
            y_obs     = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[0]) for x in x_all]
            y_exp     = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[1]) for x in x_all]
            y_low_68  = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[2]) for x in x_all]
            y_high_68 = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[3]) for x in x_all]
            y_low_95  = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[4]) for x in x_all]
            y_high_95 = [get_limit_y_for_limit_x_numnum_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limits_r[5]) for x in x_all]    
        else:
            raise ValueError('In 2-d limits (interference) with acceptance (option \'f\'), we only consider operators that affect zttmm, but this one affects zmmmm. Nah.')
    else: raise ValueError('In 2-d limits (interference) computation, must either pass no option (\'\') or option \'a\'. Instead passed: %s' % (option))

    x_high_68 = copy.deepcopy(x_all)
    x_high_95 = copy.deepcopy(x_all)
    y_high_68.reverse()
    x_high_68.reverse()
    y_high_95.reverse()
    x_high_95.reverse()

    y_68 = y_low_68 + y_high_68
    y_95 = y_low_95 + y_high_95
    x_68 = x_all + x_high_68
    x_95 = x_all + x_high_95

    y_obs_larger_than_exp  = [y_obs[i] for i in range(len(y_obs)) if y_obs[i] >= y_exp[i]]
    y_obs_smaller_than_exp = [y_obs[i] for i in range(len(y_obs)) if y_obs[i] <= y_exp[i]]
    x_obs_larger_than_exp  = [x_all[i] for i in range(len(y_obs)) if y_obs[i] >= y_exp[i]]
    x_obs_smaller_than_exp = [x_all[i] for i in range(len(y_obs)) if y_obs[i] <= y_exp[i]]
    
    draw_above = len(y_obs_larger_than_exp) > 0
    draw_below = len(y_obs_smaller_than_exp) > 0
    print draw_above, draw_below
    if draw_above: g_obs_above = ROOT.TGraph(len(x_obs_larger_than_exp), array('d', x_obs_larger_than_exp), array('d', y_obs_larger_than_exp)) 
    if draw_below: g_obs_below = ROOT.TGraph(len(x_obs_smaller_than_exp), array('d', x_obs_smaller_than_exp), array('d', y_obs_smaller_than_exp)) 

    g_obs = ROOT.TGraph(len(x_all), array('d', x_all), array('d', y_obs)) 
    g_exp = ROOT.TGraph(len(x_all), array('d', x_all), array('d', y_exp)) 
    g_68 = ROOT.TGraph(len(y_68), array('d', x_68), array('d', y_68)) 
    g_95 = ROOT.TGraph(len(y_95), array('d', x_95), array('d', y_95)) 

    maxdigits = 3
    xtitle = operatornames_pretty[opx]
    ytitle = operatornames_pretty[opy]
    ymin, ymax = xmin, xmax
    if not (proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm'):
        ymax /= 10.
        ymin /= 10.
    c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=ymin, y_max=ymax, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))

    lwidth_base = 303
    
    leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    leg.SetHeader('95% CL limits')
    tdrDraw(g_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    # tdrDraw(g_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1, fstyle=3013, lwidth=lwidth_base)
    tdrDraw(g_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
    # if not (proc_affected_per_operator[opx] == 'zttmm' and proc_affected_per_operator[opy] == 'zttmm'):
    if draw_above: tdrDraw(g_obs_above, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1, fstyle=3013, lwidth=-lwidth_base)
    if draw_below: tdrDraw(g_obs_below, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1, fstyle=3013, lwidth=lwidth_base)
    leg.AddEntry(g_obs, 'Observed', 'L')
    leg.AddEntry(g_exp, 'Median expected', 'L')
    leg.AddEntry(g_68, '68% expected', 'LF')
    leg.AddEntry(g_95, '95% expected', 'LF')
    leg.Draw('SAME')
    
    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)
    smtext.DrawLatex(0., 0.075*ymax, 'SM')

    c.RedrawAxis()
    
    outname = os.path.join(plotfolder, 'Limits2d_%s_%s_vs_%s.pdf' % ('interference', opx, opy))
    if options == 'a': outname = outname.replace('interference', 'interference_withacc')
    if options == 'f': outname = outname.replace('interference', 'interference_withf')
    if options == 'af': outname = outname.replace('interference', 'interference_withacc_withf')
    c.SaveAs(outname)
    del c



def increase_distance_from_zero(number, amount):
    if number >= 0:
        return number + amount
    else:
        return number - amount

def round_away_from_zero(number):
    return math.ceil(number) if number > 0 else math.floor(number)

def ceil_power_of_10(n):
    is_neg = False
    if n < 0: is_neg = True
    if n == 0.: return 0.
    exp = math.log(abs(n), 10)
    exp = math.ceil(exp)
    return -10**exp if is_neg else 10**exp

def find_largest_common_y(*lists):
    ys = [set([y for (x,y) in l]) for l in lists]
    if len(ys) < 2:
        return max(ys[0])
    result = ys[0]
    for i, y in enumerate(ys):
        result = result & ys[i+1]
        if len(ys) == i+2: break
    return max(result)


def plot_limits_vs_cl(op, order, l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high, plottag=''):
    if 0 in [len(l_obs), len(l_exp), len(l_68_low), len(l_68_high), len(l_95_low), len(l_95_high)]: return

    l_obs.sort(key=lambda tup: tup[1])
    l_exp.sort(key=lambda tup: tup[1])
    l_68_low.sort(key=lambda tup: tup[1])
    l_68_high.sort(key=lambda tup: tup[1], reverse=True)
    l_95_low.sort(key=lambda tup: tup[1])
    l_95_high.sort(key=lambda tup: tup[1], reverse=True)

    # max_y = min(find_largest_common_y(l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high), 0.5)
    max_y = find_largest_common_y(l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high)

    l_obs = [p for p in l_obs if p[1] <= max_y]
    l_exp = [p for p in l_exp if p[1] <= max_y]
    l_68_low = [p for p in l_68_low if p[1] <= max_y]
    l_68_high = [p for p in l_68_high if p[1] <= max_y]
    l_95_low = [p for p in l_95_low if p[1] <= max_y]
    l_95_high = [p for p in l_95_high if p[1] <= max_y]
    if 0 in [len(l_obs), len(l_exp), len(l_68_low), len(l_68_high), len(l_95_low), len(l_95_high)]: return

    l_68 = l_68_low + l_68_high
    l_95 = l_95_low + l_95_high

    l_obs = [p for p in l_obs if p[1] <= max_y]
    l_exp = [p for p in l_exp if p[1] <= max_y]
    l_68_low = [p for p in l_68_low if p[1] <= max_y]
    l_68_high = [p for p in l_68_high if p[1] <= max_y]
    l_95_low = [p for p in l_95_low if p[1] <= max_y]
    l_95_high = [p for p in l_95_high if p[1] <= max_y]

    l_obs, was_obs_pos = ensure_limits_positive(limits=l_obs)
    l_exp, was_exp_pos = ensure_limits_positive(limits=l_exp)
    l_68, was_l68_pos = ensure_limits_positive(limits=l_68)
    l_95, was_l95_pos = ensure_limits_positive(limits=l_95)



    # make graphs
    g_exp = ROOT.TGraph(len(l_exp), array('d', zip(*l_exp)[0]), array('d', zip(*l_exp)[1]))
    g_obs = ROOT.TGraph(len(l_obs), array('d', zip(*l_obs)[0]), array('d', zip(*l_obs)[1]))
    g_68 = ROOT.TGraph(len(l_68), array('d', zip(*l_68)[0]), array('d', zip(*l_68)[1]))
    g_95 = ROOT.TGraph(len(l_95), array('d', zip(*l_95)[0]), array('d', zip(*l_95)[1]))

    x_max = 1E4 if abs(l_95_high[-1][0]) < 1E4 else 1E5
    xmax = max(0., x_max)
    xmin = min(0., x_max)
    maxdigits = 2 if abs(xmax) <= 1E4 else 3
    xtitle = operatornames_pretty[op]
    if not was_obs_pos:
        xtitle = '#minus %s' % (xtitle)
        xmax = xmax / 4.
    c = tdrCanvas(canvName='c', x_min=min(0., xmax), x_max=max(0., xmax), y_min=5E-3, y_max=10., nameXaxis=xtitle, nameYaxis='1 #minus CL', square=True, iPos=11, margins=(None, 0.10, 0.14, None), maxdigits=(maxdigits, None))
    c.SetLogy()
    
    leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    leg.SetHeader('Upper limits')
    tdrDraw(g_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
    leg.AddEntry(g_obs, 'Observed', 'L')
    leg.AddEntry(g_exp, 'Median expected', 'L')
    leg.AddEntry(g_68, '68% expected', 'LF')
    leg.AddEntry(g_95, '95% expected', 'LF')
    leg.Draw('SAME')

    line = ROOT.TLine(xmin, 0.05, xmax, 0.05)
    line.SetLineColor(ROOT.kGray+1)
    line.SetLineStyle(2)
    line.Draw()

    c.RedrawAxis()
    
    outname = os.path.join(plotfolder, 'Limits_%s_%s.pdf' % (op, order))
    if plottag != '':
        outname = outname.replace('.pdf', '_%s.pdf' % (plottag))
    c.SaveAs(outname)
    del c

def ensure_limits_positive(limits):
    result = []
    if len(limits) == 0 or limits[0][0] > 0:
        return limits, True
    for l, cl in limits:
        result.append((-l, cl))
    return result, False

def get_limits_at_cl(a, b, proc_affected, order, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=None, binc=None, options='', cl=0.95):
    limits_obs, limits_exp, limits_68_low, limits_68_high, limits_95_low, limits_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
    limit_obs_cl = interpolate_x_from_y(value_list=limits_obs, target_y=1.-cl) 
    return (interpolate_x_from_y(value_list=limits_obs, target_y=1.-cl), interpolate_x_from_y(value_list=limits_exp, target_y=1.-cl), interpolate_x_from_y(value_list=limits_68_low, target_y=1.-cl), interpolate_x_from_y(value_list=limits_68_high, target_y=1.-cl), interpolate_x_from_y(value_list=limits_95_low, target_y=1.-cl), interpolate_x_from_y(value_list=limits_95_high, target_y=1.-cl))



def interpolate_x_from_y(value_list, target_y):
    value_list.sort(key=lambda point: point[1])  # Sort the list of values based on y
    y_values = [v[1] for v in value_list]
    
    if len(value_list) == 0:
        return (None, None)
    if target_y <= y_values[0]:
        return (None, None)
    if target_y >= y_values[-1]:
        return (None, None)
    
    for i in range(len(y_values) - 1):
        if y_values[i] <= target_y <= y_values[i + 1]:
            x1, y1 = value_list[i]
            x2, y2 = value_list[i + 1]
            if y1 == y2:
                raise ValueError("y1 and y2 cannot be the same for interpolation.")
            return x1 + ((x2 - x1) / (y2 - y1)) * (target_y - y1)

    raise ValueError("Interpolation points not found.")


def get_limits_vs_cl(a, b, proc_affected, order, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=None, binc=None, abr=None, bbr=None, options=''):
    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r(r_sm=r_sm)

    limits_obs = []
    limits_exp = []
    limits_68_low  = []
    limits_68_high = []
    limits_95_low  = []
    limits_95_high = []

    for ipoint in range(len(r_obs)):
        if (1-r_obs[ipoint][1])*100. > 0. :
            # print 'for obs'
            limit_obs = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_obs[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            # print limit_obs
            if limit_obs is not None:
                limits_obs.append((limit_obs, r_obs[ipoint][1]))
        if (1-r_exp[ipoint][1])*100. > 0. :
            # print 'for exp'
            limit_exp =    get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_exp[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            if limit_exp is not None:
                limits_exp.append((limit_exp, r_exp[ipoint][1]))
        if (1-r_68_low[ipoint][1])*100. > 0. :
            limit_68_low = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_68_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            if limit_68_low is not None:
                limits_68_low.append((limit_68_low, r_68_low[ipoint][1]))
        if (1-r_68_high[ipoint][1])*100. > 0. :
            limit_68_high = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_68_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            if limit_68_high is not None:
                limits_68_high.append((limit_68_high, r_68_high[ipoint][1]))                    
        if (1-r_95_low[ipoint][1])*100. > 0. :
            limit_95_low = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_95_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            if limit_95_low is not None:
                limits_95_low.append((limit_95_low, r_95_low[ipoint][1]))
        if (1-r_95_high[ipoint][1])*100. > 0. :
            limit_95_high = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_95_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, abr=abr, bbr=bbr, options=options)
            if limit_95_high is not None:
                limits_95_high.append((limit_95_high, r_95_high[ipoint][1]))
    return (limits_obs, limits_exp, limits_68_low, limits_68_high, limits_95_low, limits_95_high)

def get_limit_on_c_for_target_r(a, b, proc_affected, order, r_target, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=None, binc=None, abr=None, bbr=None, options=''):
    # there is an easy expression for the effect on R: r_with_smeft = r_sm * (1 + a*C + b*C^2). Solving by C:
    r = r_target / r_sm
    if order == 'interference':
        if a == 0: return None
        if options == '':
            if proc_affected == 'zttmm': # formula: r = (1+aC) / (1+0) --> 0 = 1+aC - r
                return (r - 1)/a
            elif proc_affected == 'zmmmm': # formula: r = (1+0) / (1+aC)
                return (1 - r)/(a * r)
            else:
                raise ValueError('Affected process is neither \'zttmm\' nor \'zmmmm\', but instead illegal value: \'%s\'' % (proc_affected))
        elif options == 'f':
            if ainc == None: raise ValueError('Option \'f\' for interference 1-d requires ainc, currently it is None.')
            if proc_affected == 'zttmm': # formula: r/r = (1+aC) * (1+aC)/(1+ainc*C)
                if (4*a**2*r - 4*a*ainc*r + ainc**2*r**2) < 0: return None
                return ((math.sqrt(4*a**2*r - 4*a*ainc*r + ainc**2*r**2) - 2*a + ainc*r)/(2*a**2))
            else:
                raise ValueError('Only considering zttmm processes, please fix.')
        elif options == 'a':
            if None in [acc_int, acc_sm]: raise ValueError('Acceptances must not be None.')

            if proc_affected == 'zttmm': # formula: r = (1+aC) / (1+0) [N-correction in numerator] * 1 / ((aC)/(1+aC) * (ASMEFT/ASM) + 1 - (aC)/(1+aC)) / (1) [Ae-correction in denominator]
                if (acc_int/acc_sm*r)**2 - 4*acc_int/acc_sm*r + 4*r < 0: return None
                return (math.sqrt((acc_int/acc_sm*r)**2 - 4*acc_int/acc_sm*r + 4*r) + acc_int/acc_sm*r - 2) / (2*a)
            elif proc_affected == 'zmmmm': # formula: r = (1+0) / (1+aC) [N-correction in denominator] * ((aC)/(1+aC) * (ASMEFT/ASM) + 1 - (aC)/(1+aC)) [Ae-correction in numerator]
                raise ValueError('limits in the case of interference and with considering the acceptance effect on Z->4mu is not yet understood. Since the limits on C are generally negative, how does the acceptance change in this case? For positive interference it is easy, but for negative? Need to implement and think about, skip for now.')
            else:
                raise ValueError('Affected process is neither \'zttmm\' nor \'zmmmm\', but instead illegal value: \'%s\'' % (proc_affected))
        elif options == 'af':
            if None in [acc_int, acc_sm, ainc]: raise ValueError('Option \'af\' for interference 1-d requires ainc and acc_int and acc_sm, currently one or more is None.')
            if proc_affected == 'zttmm': # 
                equation = sp.Eq((1+a*xvar)**2/(1+ainc*xvar) / ((a*xvar)/(1+a*xvar) * acc_int/acc_sm + 1 - (a*xvar)/(1+a*xvar)), r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                if len(result) == 0:
                    return None
                return result[-1]
            else:
                raise ValueError('Only considering zttmm processes, please fix.')
        elif options == 'afb':
            if None in [acc_int, acc_sm, ainc, abr]: raise ValueError('Option \'afb\' for interference 1-d requires ainc and acc_int and acc_sm, currently one or more is None.')
            if proc_affected == 'zttmm': # 
                equation = sp.Eq((1+a*xvar)**2/(1+ainc*xvar) / ((a*xvar)/(1+a*xvar) * acc_int/acc_sm + 1 - (a*xvar)/(1+a*xvar)) / ((1+abr*xvar)/(1+br_sm*abr*xvar))**2, r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                # print result
                if len(result) == 0:
                    return None
                return result[-1]
            else:
                raise ValueError('Only considering zttmm processes, please fix.')

        else: raise ValueError('option string must either be empty or \'a\', need to implement the rest.')
    elif order == 'comb':
        if b == 0: return None
        if options == '':
            if proc_affected == 'zttmm': # formula: r = (1+aC+bC^2) / (1+0)
                if (a/(2*b))**2 < (1 - r)/b: return None
                return (-a/(2*b) + math.sqrt((a/(2*b))**2 - (1 - r)/b))
                
            elif proc_affected == 'zmmmm': # formula: r = (1+0) / (1+aC+bC^2)
                if (a*r)**2 - 4*b*(r)**2 + 4*b*r < 0: return None
                return ((math.sqrt((a*r)**2 - 4*b*(r)**2 + 4*b*r) - a*r)/(2*b*r))
            else: 
                raise ValueError('Affected process is neither \'zttmm\' nor \'zmmmm\', but instead illegal value: \'%s\'' % (proc_affected))
        elif options == 'f':
            if None in [ainc, binc]: raise ValueError('Option \'f\' for comb 1-d requires ainc and binc, currently at least one of them is None.')
            if proc_affected == 'zttmm': # formula: r/r = (1+aC+bC^2) * (1+aC+bC^2)/(1+ainc*C+binc*C^2)
                equation = sp.Eq((1+a*xvar+b*xvar**2)**2/(1+ainc*xvar+binc*xvar**2), r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                if len(result) == 0:
                    return None
                return result[-1]
            else:
                raise ValueError('Only considering zttmm processes, please fix.')
        elif options == 'a':
            if None in [acc_int, acc_bsm, acc_sm]: raise ValueError('Acceptances must not be None.')
            if proc_affected == 'zttmm': 
                equation = sp.Eq((1+a*xvar+b*xvar**2) * (1 / ((a*xvar)/(1+a*xvar+b*xvar**2) * (acc_int/acc_sm) + (b*xvar**2)/(1+a*xvar+b*xvar**2) * (acc_bsm/acc_sm) + 1 - (a*xvar+b*xvar**2)/(1+a*xvar+b*xvar**2))), r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                if len(result) == 0:
                    return None
                return result[-1]
            elif proc_affected == 'zmmmm': 
                raise ValueError('limits in the case of combined interference+pure and with considering the acceptance effect on Z->4mu is not yet understood. There won\'t be a limit on Z->4mu anyway, so better avoid this.')
            else:
                raise ValueError('Affected process is neither \'zttmm\' nor \'zmmmm\', but instead illegal value: \'%s\'' % (proc_affected))
        elif options == 'af':
            if None in [acc_int, acc_sm, acc_bsm, ainc, binc]: raise ValueError('Option \'af\' for comb 1-d requires ainc, binc, acc_int, acc_bsm, and acc_sm, currently one or more is None.')
            if proc_affected == 'zttmm': # 
                equation = sp.Eq((1+a*xvar+b*xvar**2)**2/(1+ainc*xvar+binc*xvar**2) / ((a*xvar)/(1+a*xvar+b*xvar**2)*acc_int/acc_sm + (b*xvar**2)/(1+a*xvar+b*xvar**2)*acc_bsm/acc_sm + 1 - (a*xvar+b*xvar**2)/(1+a*xvar+b*xvar**2)), r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                if len(result) == 0:
                    return None
                return result[-1]
            else:
                raise ValueError('Only considering zttmm processes, please fix.')
        elif options == 'afb':
            if None in [acc_int, acc_sm, acc_bsm, ainc, binc, abr, bbr]: raise ValueError('Option \'afb\' for comb 1-d requires ainc, binc, acc_int, acc_bsm, and acc_sm, abr, bbr, currently one or more is None.')
            if proc_affected == 'zttmm': # 
                equation = sp.Eq((1+a*xvar+b*xvar**2)**2/(1+ainc*xvar+binc*xvar**2) / ((a*xvar)/(1+a*xvar+b*xvar**2)*acc_int/acc_sm + (b*xvar**2)/(1+a*xvar+b*xvar**2)*acc_bsm/acc_sm + 1 - (a*xvar+b*xvar**2)/(1+a*xvar+b*xvar**2)) / ((1+abr*xvar+bbr*xvar**2)/(1+br_sm*(abr*xvar + bbr*xvar**2)))**2, r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
                # print result
                if len(result) == 0:
                    return None
                return result[-1]
            else:
                raise ValueError('Only considering zttmm processes, please fix.')
        else: raise ValueError('option string must either be empty (\'\') or \'a\', need to implement the rest.')
    raise ValueError('order %s not supported for limit calculation at the moment.' % (order))

def any_nonzero_value(dictionary):
    return any(value != 0 for value in dictionary.values())

def find_procs_affected(coefficients_per_proc):
    procs_affected = []
    for proc in coefficients_per_proc:
        if any_nonzero_value(dictionary=coefficients_per_proc[proc]) and proc != 'incttmm' and proc != 'tmvv':
            procs_affected.append(proc)
    return procs_affected

def affects_incttmm(coefficients_per_proc):
    if any_nonzero_value(dictionary=coefficients_per_proc['incttmm']):
        return True
    return False

def affects_tmvv(coefficients_per_proc):
    if any_nonzero_value(dictionary=coefficients_per_proc['tmvv']):
        return True
    return False

def find_orders_affected(coefficients_per_order):
    return [key for (key, value) in coefficients_per_order.items() if value != 0]

def get_acc_bsm_intpart(bx, by, c, accx, accy, acctot):
    return (bx*(acctot-accx) + by*(acctot-accy) + c*acctot) / (c)


def get_limit_y_for_limit_x_interference_withacc_withf(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, limit):
    r = limit/r_sm
    acrx = acc_int_x / acc_sm_x
    acry = acc_int_y / acc_sm_y
    equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar) / (ax*x/(1+ax*x+ay*yvar)*acrx + ay*yvar/(1+ax*x+ay*yvar)*acry + 1./(1+ax*x+ay*yvar)), r)
    result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
    if len(result) == 0:
        return None
    return result[-1]

def get_limit_y_for_limit_x_interference_withacc_withf_withbr(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limit):
    r = limit/r_sm
    acrx = acc_int_x / acc_sm_x
    acry = acc_int_y / acc_sm_y
    equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar) / (ax*x/(1+ax*x+ay*yvar)*acrx + ay*yvar/(1+ax*x+ay*yvar)*acry + 1./(1+ax*x+ay*yvar)) / ((1+abrx*x+abry*yvar)/(1+br_sm*(abrx*x+abry*yvar)))**2, r)
    result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
    # print result
    if len(result) == 0:
        return None
    return result[-1]

def get_limit_y_for_limit_x_interference_withacc_withf_withbr_allsolutions(x, ax, ay, aincx, aincy, acc_int_x, acc_int_y, acc_sm_x, acc_sm_y, abrx, abry, limit):
    r = limit/r_sm
    acrx = acc_int_x / acc_sm_x if acc_int_x != 0 else -1.
    acry = acc_int_y / acc_sm_y if acc_int_y != 0 else -1.
    # equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar) / (ax*x/(1+ax*x+ay*yvar)*acrx + ay*yvar/(1+ax*x+ay*yvar)*acry + 1./(1+ax*x+ay*yvar)) / ((1+abrx*x+abry*yvar)/(1+br_sm*(abrx*x+abry*yvar)))**2, r)
    equation = sp.Eq((1+ax*x+ay*yvar)**2/(1+aincx*x+aincy*yvar) * 
                     ((1+ax*x/acrx+ay*yvar/acry)/(1+ax*x+ay*yvar))
                     / ((1+abrx*x+abry*yvar)/(1+br_sm*(abrx*x+abry*yvar)))**2, r)
    result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
    # print result
    if len(result) == 0:
        return None
    return result

def get_limit_y_for_limit_x_comb_withacc_withf(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, limit):
    if not acc_sm_x == acc_sm_y: raise ValueError('SM acceptances for the same process are different.')
    acc_sm = acc_sm_x
    r = limit/r_sm
    acrintx = acc_int_x / acc_sm
    acrinty = acc_int_y / acc_sm
    acrbsmx = acc_bsm_x / acc_sm
    acrbsmy = acc_bsm_y / acc_sm
    acrbsmintpart = acc_bsm_intpart / acc_sm
    equation = sp.Eq( (1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)**2/(1 + aincx*x + aincy*yvar + bincx*x**2 + bincy*yvar**2 + cinc*x*yvar) / (1./(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar) + (ax*x)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrintx + (ay*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrinty + (bx*x**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmx + (by*yvar**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmy + (c*x*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmintpart), r)
    result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
    if len(result) == 0:
        return None
    return result

def get_limit_y_for_limit_x_comb_withacc_withf_withbr(x, ax, bx, ay, by, c, aincx, bincx, aincy, bincy, cinc, acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_bsm_intpart, acc_sm_x, acc_sm_y, abrx, abry, bbrx, bbry, limit):
    if not acc_sm_x == acc_sm_y: raise ValueError('SM acceptances for the same process are different.')
    acc_sm = acc_sm_x
    r = limit/r_sm
    acrintx = acc_int_x / acc_sm if acc_int_x != 0 else -1.
    acrinty = acc_int_y / acc_sm if acc_int_y != 0 else -1.
    acrbsmx = acc_bsm_x / acc_sm if acc_bsm_x != 0 else -1.
    acrbsmy = acc_bsm_y / acc_sm if acc_bsm_y != 0 else -1.
    acrbsmintpart = acc_bsm_intpart / acc_sm if acc_bsm_intpart != 0 else -1.

    # try:
    #     eq = lambda y: (1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)**2/(1 + aincx*x + aincy*y + bincx*x**2 + bincy*y**2 + cinc*x*y) / (1./(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y) + (ax*x)/(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)*acrintx + (ay*y)/(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)*acrinty + (bx*x**2)/(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)*acrbsmx + (by*y**2)/(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)*acrbsmy + (c*x*y)/(1 + ax*x + ay*y + bx*x**2 + by*y**2 + c*x*y)*acrbsmintpart) / ((1 + abrx*x + abry*y + bbrx*x**2 + bbry*y**2)/(1+br_sm*(abrx*x + abry*y + bbrx*x**2 + bbry*y**2)))**2 - r
    #     sol = root(eq, [-10000., 10000.])
    #     print sol
    # except:
    #     print 'something did not work for scipy.optimize.root, result empty'
    #     sol = []

    try:
        # equation = sp.Eq( (1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)**2/(1 + aincx*x + aincy*yvar + bincx*x**2 + bincy*yvar**2 + cinc*x*yvar) / (1./(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar) + (ax*x)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrintx + (ay*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrinty + (bx*x**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmx + (by*yvar**2)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmy + (c*x*yvar)/(1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)*acrbsmintpart) / ((1 + abrx*x + abry*yvar + bbrx*x**2 + bbry*yvar**2)/(1+br_sm*(abrx*x + abry*yvar + bbrx*x**2 + bbry*yvar**2)))**2, r)
        equation = sp.Eq( (1 + ax*x + ay*yvar + bx*x**2 + by*yvar**2 + c*x*yvar)/(1 + aincx*x + aincy*yvar + bincx*x**2 + bincy*yvar**2 + cinc*x*yvar) *
                         ((1 + ax*x/acrintx + ay*yvar/acrinty + bx*x**2/acrbsmx + by*yvar**2/acrbsmy + c*x*yvar/acrbsmintpart)) 
                         / ((1 + abrx*x + abry*yvar + bbrx*x**2 + bbry*yvar**2)/(1+br_sm*(abrx*x + abry*yvar + bbrx*x**2 + bbry*yvar**2)))**2, r)
        result = list(sp.solveset(equation, yvar, domain=sp.S.Reals))
    except:
        print '======== something did not work for sp.solveset, result empty'
        result = []
    # print abrx, abry, bbrx, bbry
    # print 'x:', x, result
    if len(result) == 0:
        return None
    return result


def get_coefficients_for_operators(opx, opy, options, as_per_operator, bs_per_operator, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, abrs_per_operator=None, bbrs_per_operator=None):

    operatorname_sum = '%s+1p0x%s' % (opx, opy)
    if not operatorname_sum in bs_per_operator: raise ValueError('Could not find combined term c for operator-sum %s' % (operatorname_sum))

    ax = as_per_operator[opx]
    ay = as_per_operator[opy]
    if opx == 'cle2332': ax = 0.
    if opy == 'cle2332': ay = 0.
    bx = bs_per_operator[opx]
    by = bs_per_operator[opy]
    sumcoeff  = bs_per_operator[operatorname_sum]
    c = sumcoeff - bx - by

    acc_int_x = 1. if acceptances_int_per_operator == None else acceptances_int_per_operator[opx]
    acc_int_y = 1. if acceptances_int_per_operator == None else acceptances_int_per_operator[opy]
    if opx == 'cle2332': acc_int_x = 1.
    if opy == 'cle2332': acc_int_y = 1.
    acc_bsm_x = 1. if acceptances_bsm_per_operator == None else acceptances_bsm_per_operator[opx]
    acc_bsm_y = 1. if acceptances_bsm_per_operator == None else acceptances_bsm_per_operator[opy]
    acc_sm_x  = 1. if acceptances_sm_per_operator  == None else acceptances_sm_per_operator[opx]
    acc_sm_y  = 1. if acceptances_sm_per_operator  == None else acceptances_sm_per_operator[opy]
    if 1. in [acc_int_x, acc_int_y, acc_bsm_x, acc_bsm_y, acc_sm_x, acc_sm_y] and 'a' in options and not 'cle2332' in [opx, opy]: raise ValueError('Could not get all acceptances needed for using option \'f\'.')
    if opx == 'cle2332' and acc_int_y == 1: raise ValueError('Could not get all acceptances needed for using option \'f\'.')
    if opy == 'cle2332' and acc_int_x == 1: raise ValueError('Could not get all acceptances needed for using option \'f\'.')

    operatorname_sum_foracc = operatorname_sum if acceptances_bsm_per_operator != None and operatorname_sum in acceptances_bsm_per_operator else '%s+1p0x%s' % (opy, opx)
    if acceptances_bsm_per_operator == None or not operatorname_sum_foracc in acceptances_bsm_per_operator: raise ValueError('Could not find combined acceptance term for operator-sum %s' % (operatorname_sum_foracc))
    acc_bsm_withint = acceptances_bsm_per_operator[operatorname_sum_foracc]
    acc_bsm_intpart = get_acc_bsm_intpart(bx=bx, by=by, c=c, accx=acc_bsm_x, accy=acc_bsm_y, acctot=acc_bsm_withint)


    aincx = 0. if aincs_per_operator == None else aincs_per_operator[opx]
    aincy = 0. if aincs_per_operator == None else aincs_per_operator[opy]
    bincx = 0. if bincs_per_operator == None else bincs_per_operator[opx]
    bincy = 0. if bincs_per_operator == None else bincs_per_operator[opy]
    if 0. in [aincx, aincy, bincx, bincy] and 'f' in options and not 'cle2332' in [opx, opy]: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
    if opx == 'cle2332': 
        if aincy == 0.: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
        aincx = 0.
    if opy == 'cle2332': 
        if aincx == 0.: raise ValueError('Could not get all inclusive coefficients needed for using option \'f\'.')
        aincy = 0.
    if not operatorname_sum in bincs_per_operator: raise ValueError('Could not find combined term c for inclusive operator-sum %s' % (operatorname_sum))
    sumcoeffinc  = bincs_per_operator[operatorname_sum]
    cinc = sumcoeffinc - bincx - bincy

    abrx = 0. if abrs_per_operator == None else abrs_per_operator[opx]
    abry = 0. if abrs_per_operator == None else abrs_per_operator[opy]
    bbrx = 0. if bbrs_per_operator == None else bbrs_per_operator[opx]
    bbry = 0. if bbrs_per_operator == None else bbrs_per_operator[opy]

    if not 'a' in options: 
        acc_int_x = acc_int_y = acc_bsm_x = acc_bsm_y = acc_sm_x = acc_sm_y = acc_bsm_intpart = 1.
    if not 'f' in options:
        aincx = ax
        aincy = ay
        bincx = bx
        bincy = by
        cinc  = c
    if not 'b' in options: 
        abrx = abry = bbrx = bbry = 0.


    if not acc_sm_x == acc_sm_y: raise ValueError('SM acceptances for the same process are different.')
    acc_sm = acc_sm_x
    acrintx = acc_int_x / acc_sm
    acrinty = acc_int_y / acc_sm
    acrbsmx = acc_bsm_x / acc_sm
    acrbsmy = acc_bsm_y / acc_sm
    acrbsmintpart = acc_bsm_intpart / acc_sm

    return (ax, ay, bx, by, c, aincx, aincy, bincx, bincy, cinc, acrintx, acrinty, acrbsmx, acrbsmy, acrbsmintpart, abrx, abry, bbrx, bbry)

if __name__ == '__main__':
    main()

