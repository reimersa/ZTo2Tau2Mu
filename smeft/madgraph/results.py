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

import os, sys, math, json
import subprocess
from yaml import safe_load
from limits_r_smp import get_limits_r
import copy


workarea = '/work/areimers/ZTo2Tau2Mu/smeft/madgraph'
plotfolder = os.path.join(workarea, 'plots/fixed_formula')
# plotfolder = os.path.join(workarea, 'plots/xcheck')
ensureDirectory(plotfolder)


r_sm = 0.9

# operators =  ['cll2222', 'cll2233', 'cll2332', 'cle2222', 'cle2233', 'cle3322', 'cle2332', 'cee2222', 'cee2233']
operators =  ['cll2233', 'cll2332', 'cle2233', 'cle3322', 'cle2332', 'cee2233']

operatornames_pretty = {
    'cll2233': 'C_{ll}^{2233}',
    'cll2222': 'C_{ll}^{2222}', 
    'cll2332': 'C_{ll}^{2332}', 
    'cee2222': 'C_{ee}^{2222}', 
    'cee2233': 'C_{ee}^{2233}', 
    'cle2222': 'C_{le}^{2222}', 
    'cle2233': 'C_{le}^{2233}', 
    'cle3322': 'C_{le}^{3322}', 
    'cle2332': 'C_{le}^{2332}',
}

def main():
    ROOT.gROOT.SetBatch(1)

    # always read in the json with coefficients
    with open(os.path.join(workarea, 'coefficients.json'), 'r') as j:
        coefficients = safe_load(j)
    with open(os.path.join(workarea, 'acceptances.json'), 'r') as j:
        acceptances = safe_load(j)


    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r()

        
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

        # limits_1d_per_operator[op]      = get_limits_at_cl(a=a, b=b, proc_affected=proc_affected, order='interference', r_sm=r_sm, cl=0.95)
        # limits_1d_per_operator_comb[op] = get_limits_at_cl(a=a, b=b, proc_affected=proc_affected, order='comb',         r_sm=r_sm, cl=0.95)
        as_per_operator[op] = a
        bs_per_operator[op] = b
        aincs_per_operator[op] = ainc
        bincs_per_operator[op] = binc
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
        proc_affected = proc_affected_per_operator[op]
        acc_int = acceptances_int_per_operator[op]
        acc_bsm = acceptances_bsm_per_operator[op]
        acc_sm = acceptances_sm_per_operator[op]
        for order in ['interference', 'comb']:
            # if order == 'comb': continue
            l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options='')
            plot_limits_vs_cl(op=op, order=order, l_obs=l_obs, l_exp=l_exp, l_68_low=l_68_low, l_68_high=l_68_high, l_95_low=l_95_low, l_95_high=l_95_high, plottag='')

            l_withacc_obs, l_withacc_exp, l_withacc_68_low, l_withacc_68_high, l_withacc_95_low, l_withacc_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options='a')
            plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_obs, l_exp=l_withacc_exp, l_68_low=l_withacc_68_low, l_68_high=l_withacc_68_high, l_95_low=l_withacc_95_low, l_95_high=l_withacc_95_high, plottag='withacc')

            l_withf_obs, l_withf_exp, l_withf_68_low, l_withf_68_high, l_withf_95_low, l_withf_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=ainc, binc=binc, options='f')
            plot_limits_vs_cl(op=op, order=order, l_obs=l_withf_obs, l_exp=l_withf_exp, l_68_low=l_withf_68_low, l_68_high=l_withf_68_high, l_95_low=l_withf_95_low, l_95_high=l_withf_95_high, plottag='withf')

            l_withacc_withf_obs, l_withacc_withf_exp, l_withacc_withf_68_low, l_withacc_withf_68_high, l_withacc_withf_95_low, l_withacc_withf_95_high = get_limits_vs_cl(a=a, b=b, proc_affected=proc_affected, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options='af')
            plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_withf_obs, l_exp=l_withacc_withf_exp, l_68_low=l_withacc_withf_68_low, l_68_high=l_withacc_withf_68_high, l_95_low=l_withacc_withf_95_low, l_95_high=l_withacc_withf_95_high, plottag='withacc_withf')

    # ['cll2233', 'cll2332', 'cle2233', 'cle3322', 'cle2332', 'cee2233']
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233')

    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')

    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', aincs_per_operator=aincs_per_operator, options='f')

    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')
    plot_limits_2d_interference(as_per_operator=as_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, options='af')


    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233')

    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, options='a')

    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='f')

    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cll2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cll2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle3322', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2233', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cle2332', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle3322', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')
    plot_limits_2d_comb(as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, proc_affected_per_operator=proc_affected_per_operator, opx='cle2332', opy='cee2233', acceptances_int_per_operator=acceptances_int_per_operator, acceptances_bsm_per_operator=acceptances_bsm_per_operator, acceptances_sm_per_operator=acceptances_sm_per_operator, aincs_per_operator=aincs_per_operator, bincs_per_operator=bincs_per_operator, options='af')






def plot_limits_2d_comb(as_per_operator, bs_per_operator, proc_affected_per_operator, opx, opy, acceptances_int_per_operator=None, acceptances_bsm_per_operator=None, acceptances_sm_per_operator=None, aincs_per_operator=None, bincs_per_operator=None, options=''):
    # if not '2233' in opx:
    #     raise ValueError('a non-2233 operator on the x-axis is not supported.')
    
    def get_acc_bsm_intpart(bx, by, c, accx, accy, acctot):
        return (bx*(acctot-accx) + by*(acctot-accy) + c*acctot) / (c)

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
            # x_obs_p     = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
            # x_exp_p     = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
            # x_low_68_p  = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
            # x_high_68_p = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
            # x_low_95_p  = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
            # x_high_95_p = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
        
            # x_obs_n     = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
            # x_exp_n     = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
            # x_low_68_n  = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
            # x_high_68_n = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
            # x_low_95_n  = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
            # x_high_95_n = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
        
            # y_obs_pp     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[0])[0] for x in x_obs_p     ]
            # y_exp_pp     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[1])[0] for x in x_exp_p     ]
            # y_low_68_pp  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[2])[0] for x in x_low_68_p  ]
            # y_high_68_pp = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[3])[0] for x in x_high_68_p ]
            # y_low_95_pp  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[4])[0] for x in x_low_95_p  ]
            # y_high_95_pp = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[5])[0] for x in x_high_95_p ]
        
            # y_obs_pn     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[0])[1] for x in x_obs_p     ]
            # y_exp_pn     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[1])[1] for x in x_exp_p     ]
            # y_low_68_pn  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[2])[1] for x in x_low_68_p  ]
            # y_high_68_pn = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[3])[1] for x in x_high_68_p ]
            # y_low_95_pn  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[4])[1] for x in x_low_95_p  ]
            # y_high_95_pn = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[5])[1] for x in x_high_95_p ]
        
            # y_obs_np     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[0])[0] for x in x_obs_n     ]
            # y_exp_np     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[1])[0] for x in x_exp_n     ]
            # y_low_68_np  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[2])[0] for x in x_low_68_n  ]
            # y_high_68_np = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[3])[0] for x in x_high_68_n ]
            # y_low_95_np  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[4])[0] for x in x_low_95_n  ]
            # y_high_95_np = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[5])[0] for x in x_high_95_n ]
        
            # y_obs_nn     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[0])[1] for x in x_obs_n     ]
            # y_exp_nn     = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[1])[1] for x in x_exp_n     ]
            # y_low_68_nn  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[2])[1] for x in x_low_68_n  ]
            # y_high_68_nn = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[3])[1] for x in x_high_68_n ]
            # y_low_95_nn  = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[4])[1] for x in x_low_95_n  ]
            # y_high_95_nn = [get_limit_y_for_limit_x_numnum(x, ax, bx, ay, by, limits_r[5])[1] for x in x_high_95_n ]
    
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
    
        # else: # num and denom
        #     x_obs_p     = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[0]/r_sm + 1) + ay**2*limits_r[0]/r_sm)/(limits_r[0]/r_sm) > 0], keep_every=100)
        #     x_exp_p     = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[1]/r_sm + 1) + ay**2*limits_r[1]/r_sm)/(limits_r[1]/r_sm) > 0], keep_every=100)
        #     x_low_68_p  = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[2]/r_sm + 1) + ay**2*limits_r[2]/r_sm)/(limits_r[2]/r_sm) > 0], keep_every=100)
        #     x_high_68_p = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[3]/r_sm + 1) + ay**2*limits_r[3]/r_sm)/(limits_r[3]/r_sm) > 0], keep_every=100)
        #     x_low_95_p  = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[4]/r_sm + 1) + ay**2*limits_r[4]/r_sm)/(limits_r[4]/r_sm) > 0], keep_every=100)
        #     x_high_95_p = slim_list(lst=[x for x in range(50000) if (4*by*(ax*x + bx*x**2 - limits_r[5]/r_sm + 1) + ay**2*limits_r[5]/r_sm)/(limits_r[5]/r_sm) > 0], keep_every=100)
        
        #     x_obs_n     = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[0]/r_sm + 1) + ay**2*limits_r[0]/r_sm)/(limits_r[0]/r_sm) > 0], keep_every=100)
        #     x_exp_n     = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[1]/r_sm + 1) + ay**2*limits_r[1]/r_sm)/(limits_r[1]/r_sm) > 0], keep_every=100)
        #     x_low_68_n  = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[2]/r_sm + 1) + ay**2*limits_r[2]/r_sm)/(limits_r[2]/r_sm) > 0], keep_every=100)
        #     x_high_68_n = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[3]/r_sm + 1) + ay**2*limits_r[3]/r_sm)/(limits_r[3]/r_sm) > 0], keep_every=100)
        #     x_low_95_n  = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[4]/r_sm + 1) + ay**2*limits_r[4]/r_sm)/(limits_r[4]/r_sm) > 0], keep_every=100)
        #     x_high_95_n = slim_list(lst=[-x for x in range(50000) if (4*by*(ax*(-x) + bx*x**2 - limits_r[5]/r_sm + 1) + ay**2*limits_r[5]/r_sm)/(limits_r[5]/r_sm) > 0], keep_every=100)
        
        #     y_obs_pp     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[0])[0] for x in x_obs_p     ]
        #     y_exp_pp     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[1])[0] for x in x_exp_p     ]
        #     y_low_68_pp  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[2])[0] for x in x_low_68_p  ]
        #     y_high_68_pp = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[3])[0] for x in x_high_68_p ]
        #     y_low_95_pp  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[4])[0] for x in x_low_95_p  ]
        #     y_high_95_pp = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[5])[0] for x in x_high_95_p ]
        
        #     y_obs_pn     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[0])[1] for x in x_obs_p     ]
        #     y_exp_pn     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[1])[1] for x in x_exp_p     ]
        #     y_low_68_pn  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[2])[1] for x in x_low_68_p  ]
        #     y_high_68_pn = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[3])[1] for x in x_high_68_p ]
        #     y_low_95_pn  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[4])[1] for x in x_low_95_p  ]
        #     y_high_95_pn = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[5])[1] for x in x_high_95_p ]
        
        #     y_obs_np     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[0])[0] for x in x_obs_n     ]
        #     y_exp_np     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[1])[0] for x in x_exp_n     ]
        #     y_low_68_np  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[2])[0] for x in x_low_68_n  ]
        #     y_high_68_np = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[3])[0] for x in x_high_68_n ]
        #     y_low_95_np  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[4])[0] for x in x_low_95_n  ]
        #     y_high_95_np = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[5])[0] for x in x_high_95_n ]
        
        #     y_obs_nn     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[0])[1] for x in x_obs_n     ]
        #     y_exp_nn     = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[1])[1] for x in x_exp_n     ]
        #     y_low_68_nn  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[2])[1] for x in x_low_68_n  ]
        #     y_high_68_nn = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[3])[1] for x in x_high_68_n ]
        #     y_low_95_nn  = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[4])[1] for x in x_low_95_n  ]
        #     y_high_95_nn = [get_limit_y_for_limit_x_numden(x, ax, bx, ay, by, limits_r[5])[1] for x in x_high_95_n ]  

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
    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r()
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


def get_limits_vs_cl(a, b, proc_affected, order, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=None, binc=None, options=''):
    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r()

    limits_obs = []
    limits_exp = []
    limits_68_low  = []
    limits_68_high = []
    limits_95_low  = []
    limits_95_high = []

    for ipoint in range(len(r_obs)):
        if (1-r_obs[ipoint][1])*100. > 0. :
            # print 'for obs'
            limit_obs = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_obs[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            # print limit_obs
            if limit_obs is not None:
                limits_obs.append((limit_obs, r_obs[ipoint][1]))
        if (1-r_exp[ipoint][1])*100. > 0. :
            # print 'for exp'
            limit_exp =    get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_exp[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            if limit_exp is not None:
                limits_exp.append((limit_exp, r_exp[ipoint][1]))
        if (1-r_68_low[ipoint][1])*100. > 0. :
            limit_68_low = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_68_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            if limit_68_low is not None:
                limits_68_low.append((limit_68_low, r_68_low[ipoint][1]))
        if (1-r_68_high[ipoint][1])*100. > 0. :
            limit_68_high = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_68_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            if limit_68_high is not None:
                limits_68_high.append((limit_68_high, r_68_high[ipoint][1]))                    
        if (1-r_95_low[ipoint][1])*100. > 0. :
            limit_95_low = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_95_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            if limit_95_low is not None:
                limits_95_low.append((limit_95_low, r_95_low[ipoint][1]))
        if (1-r_95_high[ipoint][1])*100. > 0. :
            limit_95_high = get_limit_on_c_for_target_r(a=a, b=b, proc_affected=proc_affected, order=order, r_target=r_95_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, ainc=ainc, binc=binc, options=options)
            if limit_95_high is not None:
                limits_95_high.append((limit_95_high, r_95_high[ipoint][1]))
    return (limits_obs, limits_exp, limits_68_low, limits_68_high, limits_95_low, limits_95_high)

def get_limit_on_c_for_target_r(a, b, proc_affected, order, r_target, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, ainc=None, binc=None, options=''):
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
            if None in [acc_int, acc_sm, acc_bsm, ainc, binc]: raise ValueError('Option \'af\' for interference 1-d requires ainc, binc, acc_int, acc_bsm, and acc_sm, currently one or more is None.')
            if proc_affected == 'zttmm': # 
                equation = sp.Eq((1+a*xvar+b*xvar**2)**2/(1+ainc*xvar+binc*xvar**2) / ((a*xvar)/(1+a*xvar+b*xvar**2)*acc_int/acc_sm + (b*xvar**2)/(1+a*xvar+b*xvar**2)*acc_bsm/acc_sm + 1 - (a*xvar+b*xvar**2)/(1+a*xvar+b*xvar**2)), r)
                result = list(sp.solveset(equation, xvar, domain=sp.S.Reals))
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
        if any_nonzero_value(dictionary=coefficients_per_proc[proc]) and proc != 'incttmm':
            procs_affected.append(proc)
    return procs_affected

def affects_incttmm(coefficients_per_proc):
    if any_nonzero_value(dictionary=coefficients_per_proc['incttmm']):
        return True
    return False

def find_orders_affected(coefficients_per_order):
    return [key for (key, value) in coefficients_per_order.items() if value != 0]


if __name__ == '__main__':
    main()

