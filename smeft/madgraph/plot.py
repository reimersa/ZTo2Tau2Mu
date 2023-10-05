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

import os, sys, math, json
import subprocess
from yaml import safe_load
from limits_r_smp import get_limits_r


description = """Running Madgraph locally and submitting."""
parser = ArgumentParser(prog="workspaces",description=description,epilog="Finished successfully!")
parser.add_argument('-c', "--cllimits",    dest="cllimits", default=False, action='store_true',
                                           help="plot limits on single operators for varying CL" )

args = parser.parse_args()

workarea = '/work/areimers/ZTo2Tau2Mu/smeft/madgraph'
plotfolder = os.path.join(workarea, 'plots')
ensureDirectory(plotfolder)


r_sm = 0.9

operators =  ['cll2222', 'cll2233', 'cll2332', 'cle2222', 'cle2233', 'cle3322', 'cle2332', 'cee2222', 'cee2233']
operators_interference = ['']

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

sign_per_proc = {
    'zttmm': 1.0,
    'zmmmm': -1.0
}

def main():
    ROOT.gROOT.SetBatch(1)

    # always read in the json with coefficients
    with open(os.path.join(workarea, 'coefficients.json'), 'r') as j:
        coefficients = safe_load(j)
    with open(os.path.join(workarea, 'acceptances.json'), 'r') as j:
        acceptances = safe_load(j)


    r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r()

    if args.cllimits:
        
        # get limits from known formula
        limits_1d_per_operator = {}
        limits_1d_per_operator_comb = {}
        as_per_operator = {}
        bs_per_operator = {}
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
                a = sign_per_proc[proc_affected] * coefficients_per_proc[proc_affected]['interference']
                b = sign_per_proc[proc_affected] * coefficients_per_proc[proc_affected]['purebsm']
                if op in operators:
                    acc_op_int = acceptances[op][proc_affected]['interference'] if 'interference' in acceptances[op][proc_affected] else 0.
                    acc_op_bsm = acceptances[op][proc_affected]['purebsm'] if 'purebsm' in acceptances[op][proc_affected] else 0.
                    acc_sm = acceptances['sm'][proc_affected]['sm']
            elif len(proc_affected) == 2:
                continue
            else:
                print '-- operator:', op
                print 'procs affected:', proc_affected
                raise ValueError('Not 1 or 2 process affected, huh?')

            limits_1d_per_operator[op] = get_limits_at_cl(a=a, b=b, order='interference', r_sm=r_sm, cl=0.95)
            limits_1d_per_operator_comb[op] = get_limits_at_cl(a=a, b=b, order='comb', r_sm=r_sm, cl=0.95)
            as_per_operator[op] = a
            bs_per_operator[op] = b
            if op in operators:
                acceptances_int_per_operator[op] = acc_op_int
                acceptances_bsm_per_operator[op] = acc_op_bsm
                acceptances_sm_per_operator[op] = acc_sm
        
        for op in operators:
            a = as_per_operator[op]
            b = bs_per_operator[op]
            acc_int = acceptances_int_per_operator[op]
            acc_bsm = acceptances_bsm_per_operator[op]
            acc_sm = acceptances_sm_per_operator[op]
            for order in ['interference', 'comb']:

                if order == 'comb': continue
                # print a, b
                # print acc_int, acc_bsm, acc_sm
                l_obs, l_exp, l_68_low, l_68_high, l_95_low, l_95_high = get_limits_vs_cl(a=a, b=b, order=order, r_sm=r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options='')
                plot_limits_vs_cl(op=op, order=order, l_obs=l_obs, l_exp=l_exp, l_68_low=l_68_low, l_68_high=l_68_high, l_95_low=l_95_low, l_95_high=l_95_high, plottag='')

                l_withacc_obs, l_withacc_exp, l_withacc_68_low, l_withacc_68_high, l_withacc_95_low, l_withacc_95_high = get_limits_vs_cl(a=a, b=b, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options='a')
                plot_limits_vs_cl(op=op, order=order, l_obs=l_withacc_obs, l_exp=l_withacc_exp, l_68_low=l_withacc_68_low, l_68_high=l_withacc_68_high, l_95_low=l_withacc_95_low, l_95_high=l_withacc_95_high, plottag='withacc')
                # print 'no acc:  ',get_limit_on_c_for_target_r(a=a, b=b, order='interference', r_target=6.2, r_sm=0.9, acc_int=None, acc_sm=None, acc_bsm=None, options='')
                # print 'with acc:',get_limit_on_c_for_target_r(a=a, b=b, order='interference', r_target=6.2, r_sm=0.9, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=None, options='a')


        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cll2233', opy='cle2233')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cll2233', opy='cee2233')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cle2233', opy='cee2233')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cll2233', opy='cll2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cll2233', opy='cle2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cll2233', opy='cee2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cle2233', opy='cll2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cle2233', opy='cle2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cle2233', opy='cee2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cee2233', opy='cll2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cee2233', opy='cle2222')
        plot_limits_2d_interference(limits_per_operator=limits_1d_per_operator, as_per_operator=as_per_operator, opx='cee2233', opy='cee2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2233')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2233')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2233')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cll2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cle2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cll2233', opy='cee2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cll2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cle2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cle2233', opy='cee2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cee2233', opy='cll2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cee2233', opy='cle2222')
        plot_limits_2d_comb(limits_per_operator=limits_1d_per_operator_comb, as_per_operator=as_per_operator, bs_per_operator=bs_per_operator, opx='cee2233', opy='cee2222')







def plot_limits_2d_comb(limits_per_operator, as_per_operator, bs_per_operator, opx, opy):
    if not '2233' in opx:
        raise ValueError('a non-2233 operator on the x-axis is not supported.')

    # # slightly less easy but still OK analytic relation for limits
    limits_x = limits_per_operator[opx]
    limits_y = limits_per_operator[opy]
    ax = as_per_operator[opx]
    ay = as_per_operator[opy]
    bx = bs_per_operator[opx]
    by = bs_per_operator[opy]
    operatorname_sum = '%s+1p0x%s' % (opx, opy)
    if not operatorname_sum in bs_per_operator:
        c = 0
    else:
        sumcoeff  = bs_per_operator['%s+1p0x%s' % (opx, opy)]
        c = sumcoeff - bx - by
        print '%s:' % (opx), bx
        print '%s:' % (opy), by
        print '%s:' % (operatorname_sum), sumcoeff
    print 'interference coefficient =', c

    # c *= 1000

    # if there is an interference term, the solution to the master equation is (according to Wolfram):
    # y = -( (-math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) )


    # y = -(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (R_smeft / R_sm - 1 - (ax * Cx + bx * Cx**2))/by)
    # R_smeft is the limit on R by SMP for (obs, exp, etc) at 95% CL
    limits_r = get_limits_r_at_cl(cl=0.95)
    x_obs_p     = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
    x_exp_p     = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
    x_low_68_p  = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
    x_high_68_p = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
    x_low_95_p  = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)
    x_high_95_p = slim_list(lst=[x for x in range(50000) if (ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by > 0], keep_every=100)

    x_obs_n     = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
    x_exp_n     = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
    x_low_68_n  = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
    x_high_68_n = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
    x_low_95_n  = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)
    x_high_95_n = slim_list(lst=[-x for x in range(50000) if (ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * (-x) + bx * (-x)**2))/by > 0], keep_every=100)

    y_obs_pp     = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_obs_p     ]
    y_exp_pp     = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_exp_p     ]
    y_low_68_pp  = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_68_p  ]
    y_high_68_pp = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_68_p ]
    y_low_95_pp  = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_95_p  ]
    y_high_95_pp = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_95_p ]

    y_obs_pn     = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_obs_p     ]
    y_exp_pn     = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_exp_p     ]
    y_low_68_pn  = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_68_p  ]
    y_high_68_pn = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_68_p ]
    y_low_95_pn  = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_95_p  ]
    y_high_95_pn = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_95_p ]

    y_obs_np     = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_obs_n]
    y_exp_np     = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_exp_n]
    y_low_68_np  = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_68_n]
    y_high_68_np = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_68_n]
    y_low_95_np  = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_95_n]
    y_high_95_np = [-(ay/(2*by)) + math.sqrt((ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_95_n]

    y_obs_nn     = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[0] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_obs_n]
    y_exp_nn     = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[1] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_exp_n]
    y_low_68_nn  = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[2] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_68_n]
    y_high_68_nn = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[3] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_68_n]
    y_low_95_nn  = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[4] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_low_95_n]
    y_high_95_nn = [-(ay/(2*by)) - math.sqrt((ay/(2*by))**2 + (limits_r[5] / r_sm - 1 - (ax * x + bx * x**2))/by) for x in x_high_95_n]

    if c != 0:
        x_obs_p     = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_exp_p     = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_low_68_p  = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_high_68_p = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_low_95_p  = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_high_95_p = slim_list(lst=[x for x in range(50000) if  4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2 > 0], keep_every=100)
        x_obs_n     = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[0] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)
        x_exp_n     = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[1] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)
        x_low_68_n  = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[2] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)
        x_high_68_n = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[3] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)
        x_low_95_n  = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[4] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)
        x_high_95_n = slim_list(lst=[-x for x in range(50000) if 4*by*((limits_r[5] / r_sm - 1) - (-x)*(ax + bx*(-x))) + ay**2 + 2*ay*c*(-x) + c**2*(-x)**2 > 0], keep_every=100)

        y_obs_pp     = [-( (-math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_obs_p     ]
        y_exp_pp     = [-( (-math.sqrt(4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_exp_p     ]
        y_low_68_pp  = [-( (-math.sqrt(4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_68_p  ]
        y_high_68_pp = [-( (-math.sqrt(4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_68_p ]
        y_low_95_pp  = [-( (-math.sqrt(4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_95_p  ]
        y_high_95_pp = [-( (-math.sqrt(4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_95_p ]
        y_obs_pn     = [-( (+math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_obs_p     ]
        y_exp_pn     = [-( (+math.sqrt(4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_exp_p     ]
        y_low_68_pn  = [-( (+math.sqrt(4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_68_p  ]
        y_high_68_pn = [-( (+math.sqrt(4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_68_p ]
        y_low_95_pn  = [-( (+math.sqrt(4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_95_p  ]
        y_high_95_pn = [-( (+math.sqrt(4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_95_p ]
        y_obs_np     = [-( (-math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_obs_n]
        y_exp_np     = [-( (-math.sqrt(4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_exp_n]
        y_low_68_np  = [-( (-math.sqrt(4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_68_n]
        y_high_68_np = [-( (-math.sqrt(4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_68_n]
        y_low_95_np  = [-( (-math.sqrt(4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_95_n]
        y_high_95_np = [-( (-math.sqrt(4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_95_n]
        y_obs_nn     = [-( (+math.sqrt(4*by*((limits_r[0] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_obs_n]
        y_exp_nn     = [-( (+math.sqrt(4*by*((limits_r[1] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_exp_n]
        y_low_68_nn  = [-( (+math.sqrt(4*by*((limits_r[2] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_68_n]
        y_high_68_nn = [-( (+math.sqrt(4*by*((limits_r[3] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_68_n]
        y_low_95_nn  = [-( (+math.sqrt(4*by*((limits_r[4] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_low_95_n]
        y_high_95_nn = [-( (+math.sqrt(4*by*((limits_r[5] / r_sm - 1) - x*(ax + bx*x)) + ay**2 + 2*ay*c*x + c**2*x**2) + ay + c*x)/(2*by) ) for x in x_high_95_n]

    x_high_68_p.reverse()
    x_high_95_p.reverse()   
    x_high_68_n.reverse()
    x_high_95_n.reverse()    

    lwidth_base = 303
    if '2233' in opy:
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

        lwidth = +lwidth_base

    elif '2222' in opy:
        y_high_68_pp.reverse()
        y_high_95_pp.reverse()

        y_high_68_pn.reverse()
        y_high_95_pn.reverse()
        
        y_high_68_np.reverse()
        y_high_95_np.reverse()

        y_high_68_nn.reverse()
        y_high_95_nn.reverse()

        x_68_p  = x_low_68_p[::-1] + x_low_68_p + x_high_68_p + x_high_68_p[::-1]
        x_95_p  = x_low_95_p[::-1] + x_low_95_p + x_high_95_p + x_high_95_p[::-1]
        x_68_n  = x_low_68_n[::-1] + x_low_68_n + x_high_68_n + x_high_68_n[::-1]
        x_95_n  = x_low_95_n[::-1] + x_low_95_n + x_high_95_n + x_high_95_n[::-1]

        y_68_p = y_low_68_pp[::-1] + y_low_68_pn + y_high_68_pn + y_high_68_pp[::-1]
        y_95_p = y_low_95_pp[::-1] + y_low_95_pn + y_high_95_pn + y_high_95_pp[::-1]
        y_68_n = y_low_68_np[::-1] + y_low_68_nn + y_high_68_nn + y_high_68_np[::-1]
        y_95_n = y_low_95_np[::-1] + y_low_95_nn + y_high_95_nn + y_high_95_np[::-1]

        g_obs_p = ROOT.TGraph(len(y_obs_pp+y_obs_pn), array('d', x_obs_p[::-1]+x_obs_p), array('d', y_obs_pp[::-1]+y_obs_pn)) 
        g_exp_p = ROOT.TGraph(len(y_exp_pp+y_exp_pn), array('d', x_exp_p[::-1]+x_exp_p), array('d', y_exp_pp[::-1]+y_exp_pn)) 
        g_obs_n = ROOT.TGraph(len(y_obs_np+y_obs_nn), array('d', x_obs_n[::-1]+x_obs_n), array('d', y_obs_np[::-1]+y_obs_nn)) 
        g_exp_n = ROOT.TGraph(len(y_exp_np+y_exp_nn), array('d', x_exp_n[::-1]+x_exp_n), array('d', y_exp_np[::-1]+y_exp_nn)) 

        lwidth = +lwidth_base
    else:
        raise ValueError('some operator other than 2222 or 2233 is supposed to be plotted on the y axis, not supported.')

    g_68_p = ROOT.TGraph(len(y_68_p), array('d', x_68_p), array('d', y_68_p)) 
    g_95_p = ROOT.TGraph(len(y_95_p), array('d', x_95_p), array('d', y_95_p)) 

    g_68_n = ROOT.TGraph(len(y_68_n), array('d', x_68_n), array('d', y_68_n)) 
    g_95_n = ROOT.TGraph(len(y_95_n), array('d', x_95_n), array('d', y_95_n)) 


    maxdigits = 3
    xtitle = operatornames_pretty[opx]
    ytitle = operatornames_pretty[opy]
    xmin = -2.0E4 
    xmax = +2.0E4
    c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=xmin, y_max=xmax, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    
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
    leg.Draw('SAME')

    # line = ROOT.TLine(xmin, 0, xmax, 0)
    # line.SetLineColor(ROOT.kBlack)
    # line.SetLineStyle(2)
    # line.Draw()
    # line2 = ROOT.TLine(0, xmin, 0, xmax)
    # line2.SetLineColor(ROOT.kBlack)
    # line2.SetLineStyle(2)
    # line2.Draw()

    marker_sm = ROOT.TMarker(0., 0., 20)
    marker_sm.Draw('SAME')

    smtext = ROOT.TLatex()
    smtext.SetTextFont(42)
    smtext.SetTextSize(0.035)
    smtext.SetTextAlign(21)
    smtext.DrawLatex(0., 8E2, 'SM')

    c.RedrawAxis()
    
    outname = os.path.join(plotfolder, 'Limits2d_%s_%s_vs_%s.pdf' % ('comb', opx, opy))
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


def plot_limits_2d_interference(limits_per_operator, as_per_operator, opx, opy):
    # easy analytic relation for limits, they will be a straight line
    limits_x = limits_per_operator[opx]
    limits_y = limits_per_operator[opy]
    ax = as_per_operator[opx]
    ay = as_per_operator[opy]

    # y = lim_y  - ax/ay * lim_x

    # make graphs for exp, obs, 68, 95 and fill the area between (like also done before)
    npoints = 5
    # eval_min = round_away_from_zero(xmin/limits_x[0])
    # eval_max = math.ceil(xmax/limits_x[0])
    xmin = int(-5E4)
    xmax = int(5E4)
    y_obs     = [limits_y[0] - ax/ay * x*limits_x[0] for x in np.arange(round_away_from_zero(xmin/limits_x[0]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[0]), abs((round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints)), (round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints)]
    y_exp     = [limits_y[1] - ax/ay * x*limits_x[1] for x in np.arange(round_away_from_zero(xmin/limits_x[1]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[1]), abs((round_away_from_zero(xmax/limits_x[1])-round_away_from_zero(xmin/limits_x[1]))/npoints)), (round_away_from_zero(xmax/limits_x[1])-round_away_from_zero(xmin/limits_x[1]))/npoints)]
    y_low_68  = [limits_y[2] - ax/ay * x*limits_x[2] for x in np.arange(round_away_from_zero(xmin/limits_x[2]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[2]), abs((round_away_from_zero(xmax/limits_x[2])-round_away_from_zero(xmin/limits_x[2]))/npoints)), (round_away_from_zero(xmax/limits_x[2])-round_away_from_zero(xmin/limits_x[2]))/npoints)]
    y_high_68 = [limits_y[3] - ax/ay * x*limits_x[3] for x in np.arange(round_away_from_zero(xmin/limits_x[3]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[3]), abs((round_away_from_zero(xmax/limits_x[3])-round_away_from_zero(xmin/limits_x[3]))/npoints)), (round_away_from_zero(xmax/limits_x[3])-round_away_from_zero(xmin/limits_x[3]))/npoints)]
    y_low_95  = [limits_y[4] - ax/ay * x*limits_x[4] for x in np.arange(round_away_from_zero(xmin/limits_x[4]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[4]), abs((round_away_from_zero(xmax/limits_x[4])-round_away_from_zero(xmin/limits_x[4]))/npoints)), (round_away_from_zero(xmax/limits_x[4])-round_away_from_zero(xmin/limits_x[4]))/npoints)]
    y_high_95 = [limits_y[5] - ax/ay * x*limits_x[5] for x in np.arange(round_away_from_zero(xmin/limits_x[5]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[5]), abs((round_away_from_zero(xmax/limits_x[5])-round_away_from_zero(xmin/limits_x[5]))/npoints)), (round_away_from_zero(xmax/limits_x[5])-round_away_from_zero(xmin/limits_x[5]))/npoints)]
    x_obs     = [x*limits_x[0] for x in np.arange(round_away_from_zero(xmin/limits_x[0]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[0]), abs((round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints)), (round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints)]
    x_exp     = [x*limits_x[1] for x in np.arange(round_away_from_zero(xmin/limits_x[1]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[1]), abs((round_away_from_zero(xmax/limits_x[1])-round_away_from_zero(xmin/limits_x[1]))/npoints)), (round_away_from_zero(xmax/limits_x[1])-round_away_from_zero(xmin/limits_x[1]))/npoints)]
    x_low_68  = [x*limits_x[2] for x in np.arange(round_away_from_zero(xmin/limits_x[2]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[2]), abs((round_away_from_zero(xmax/limits_x[2])-round_away_from_zero(xmin/limits_x[2]))/npoints)), (round_away_from_zero(xmax/limits_x[2])-round_away_from_zero(xmin/limits_x[2]))/npoints)]
    x_high_68 = [x*limits_x[3] for x in np.arange(round_away_from_zero(xmin/limits_x[3]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[3]), abs((round_away_from_zero(xmax/limits_x[3])-round_away_from_zero(xmin/limits_x[3]))/npoints)), (round_away_from_zero(xmax/limits_x[3])-round_away_from_zero(xmin/limits_x[3]))/npoints)]
    x_low_95  = [x*limits_x[4] for x in np.arange(round_away_from_zero(xmin/limits_x[4]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[4]), abs((round_away_from_zero(xmax/limits_x[4])-round_away_from_zero(xmin/limits_x[4]))/npoints)), (round_away_from_zero(xmax/limits_x[4])-round_away_from_zero(xmin/limits_x[4]))/npoints)]
    x_high_95 = [x*limits_x[5] for x in np.arange(round_away_from_zero(xmin/limits_x[5]), increase_distance_from_zero(round_away_from_zero(xmax/limits_x[5]), abs((round_away_from_zero(xmax/limits_x[5])-round_away_from_zero(xmin/limits_x[5]))/npoints)), (round_away_from_zero(xmax/limits_x[5])-round_away_from_zero(xmin/limits_x[5]))/npoints)]
    y_high_68.reverse()
    x_high_68.reverse()
    y_high_95.reverse()
    x_high_95.reverse()
    # print limits_x[0]
    # print round_away_from_zero(xmin/limits_x[0]), round_away_from_zero(xmax/limits_x[0]), (round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints
    # print np.arange(round_away_from_zero(xmin/limits_x[0]), round_away_from_zero(xmax/limits_x[0]), (round_away_from_zero(xmax/limits_x[0])-round_away_from_zero(xmin/limits_x[0]))/npoints)

    y_68 = y_low_68 + y_high_68
    y_95 = y_low_95 + y_high_95
    x_68 = x_low_68 + x_high_68
    x_95 = x_low_95 + x_high_95

    g_obs = ROOT.TGraph(len(x_obs), array('d', x_obs), array('d', y_obs)) 
    g_exp = ROOT.TGraph(len(x_exp), array('d', x_exp), array('d', y_exp)) 
    g_68 = ROOT.TGraph(len(y_68), array('d', x_68), array('d', y_68)) 
    g_95 = ROOT.TGraph(len(y_95), array('d', x_95), array('d', y_95)) 

    # xmin = min(x_obs+x_exp+x_68+x_95)
    # xmax = max(x_obs+x_exp+x_68+x_95)
    # ymin = min(y_obs+y_exp+y_68+y_95)
    # ymax = max(y_obs+y_exp+y_68+y_95)
    # xmax = ceil_power_of_10(xmax)
    # ymax = ceil_power_of_10(ymax)
    # xmin = ceil_power_of_10(xmin)
    # ymin = ceil_power_of_10(ymin)
    # ymax_abs = max(abs(ymin), abs(ymax))
    # xmax_abs = max(abs(xmin), abs(xmax))

    maxdigits = 3
    xtitle = operatornames_pretty[opx]
    ytitle = operatornames_pretty[opy]
    # c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=ymin, y_max=ymax, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    # c = tdrCanvas(canvName='c', x_min=-xmax_abs, x_max=xmax_abs, y_min=-ymax_abs, y_max=ymax_abs, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    c = tdrCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=xmin, y_max=xmax, nameXaxis=xtitle, nameYaxis=ytitle, square=True, iPos=11, margins=(None, 0.10, 0.14, 0.16), maxdigits=(maxdigits, maxdigits))
    
    leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    leg.SetHeader('95% CL limits')
    tdrDraw(g_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    tdrDraw(g_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    tdrDraw(g_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    tdrDraw(g_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
    leg.AddEntry(g_obs, 'Observed', 'L')
    leg.AddEntry(g_exp, 'Median expected', 'L')
    leg.AddEntry(g_68, '68% expected', 'LF')
    leg.AddEntry(g_95, '95% expected', 'LF')
    leg.Draw('SAME')

    line = ROOT.TLine(xmin, 0, xmax, 0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.Draw()
    line2 = ROOT.TLine(0, xmin, 0, xmax)
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineStyle(2)
    line2.Draw()

    c.RedrawAxis()
    
    outname = os.path.join(plotfolder, 'Limits2d_%s_%s_vs_%s.pdf' % ('interference', opx, opy))
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
    c = tdrCanvas(canvName='c', x_min=min(0., xmax), x_max=max(0., xmax), y_min=5E-3, y_max=5., nameXaxis=xtitle, nameYaxis='1 #minus CL', square=True, iPos=11, margins=(None, 0.10, 0.14, None), maxdigits=(maxdigits, None))
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

def get_limits_at_cl(a, b, order, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options='', cl=0.95):
    limits_obs, limits_exp, limits_68_low, limits_68_high, limits_95_low, limits_95_high = get_limits_vs_cl(a=a, b=b, order=order, r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
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


def get_limits_vs_cl(a, b, order, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options=''):
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
            limit_obs = get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_obs[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_obs is not None:
                limits_obs.append((limit_obs, r_obs[ipoint][1]))
        if (1-r_exp[ipoint][1])*100. > 0. :
            # print 'for exp'
            limit_exp =    get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_exp[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_exp is not None:
                limits_exp.append((limit_exp, r_exp[ipoint][1]))
        if (1-r_68_low[ipoint][1])*100. > 0. :
            limit_68_low = get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_68_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_68_low is not None:
                limits_68_low.append((limit_68_low, r_68_low[ipoint][1]))
        if (1-r_68_high[ipoint][1])*100. > 0. :
            limit_68_high = get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_68_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_68_high is not None:
                limits_68_high.append((limit_68_high, r_68_high[ipoint][1]))                    
        if (1-r_95_low[ipoint][1])*100. > 0. :
            limit_95_low = get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_95_low[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_95_low is not None:
                limits_95_low.append((limit_95_low, r_95_low[ipoint][1]))
        if (1-r_95_high[ipoint][1])*100. > 0. :
            limit_95_high = get_limit_on_c_for_target_r(a=a, b=b, order=order, r_target=r_95_high[ipoint][0], r_sm=r_sm, acc_int=acc_int, acc_sm=acc_sm, acc_bsm=acc_bsm, options=options)
            if limit_95_high is not None:
                limits_95_high.append((limit_95_high, r_95_high[ipoint][1]))
    return (limits_obs, limits_exp, limits_68_low, limits_68_high, limits_95_low, limits_95_high)

def get_limit_on_c_for_target_r(a, b, order, r_target, r_sm, acc_int=None, acc_sm=None, acc_bsm=None, options=''):
    # there is an easy expression for the effect on R: r_with_smeft = r_sm * (1 + a*C + b*C^2). Solving by C:
    if order == 'interference':
        if a == 0: return None
        if options == '':
            return (r_target/r_sm - 1)/a
        elif options == 'a':
            if None in [acc_int, acc_sm]: raise ValueError('Acceptances must not be None.')

            if ((acc_int/acc_sm)**2 + 2*(acc_int/acc_sm)*(r_target/r_sm - 2) + (r_target/r_sm)**2 + 4) < 0: return None
            return ( +math.sqrt((acc_int/acc_sm)**2 + 2*(acc_int/acc_sm)*(r_target/r_sm - 2) + (r_target/r_sm)**2 + 4) + acc_int/acc_sm + r_target/r_sm - 2)/(2*a)


        else: raise ValueError('option string must either be empty or \'a\', need to implement the rest.')
    if order == 'comb':
        if b == 0: return None
        if (a/(2*b))**2 < (1 - r_target/r_sm)/b: return None

        if options == '':
            result = -a/(2*b) + math.sqrt((a/(2*b))**2 - (1 - r_target/r_sm)/b)
        elif options == 'a':
            if None in [acc_int, acc_bsm, acc_sm]: raise ValueError('Acceptances must not be None.')
            pass
        else: raise ValueError('option string must either be empty (\'\') or \'a\', need to implement the rest.')
        return result
    raise ValueError('order %s not supported for limit calculation at the moment.' % (order))

def any_nonzero_value(dictionary):
    return any(value != 0 for value in dictionary.values())

def find_procs_affected(coefficients_per_proc):
    procs_affected = []
    for proc in coefficients_per_proc:
        if any_nonzero_value(dictionary=coefficients_per_proc[proc]):
            procs_affected.append(proc)
    return procs_affected

def find_orders_affected(coefficients_per_order):
    return [key for (key, value) in coefficients_per_order.items() if value != 0]


if __name__ == '__main__':
    main()
