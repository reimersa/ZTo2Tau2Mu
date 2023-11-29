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

import os, sys, math, json
import subprocess


description = """Running Madgraph locally and submitting."""
parser = ArgumentParser(prog="workspaces",description=description,epilog="Finished successfully!")
# parser.add_argument('-l', "--local",       dest="tag", type=str, default='_newffsuman2302_bounduncor_unnorm_fixedttunc_tauidinflatedpog_tauiddmuncor_finalresults', action='store',
#                                            help="tag that gets appended to folders and filenames" )
parser.add_argument('-s', "--submit",      dest="submit", default=False, action='store_true',
                                           help="submit jobs to the cluster" )
parser.add_argument('-l', "--local",       dest="local", default=False, action='store_true',
                                           help="run jobs locally" )
parser.add_argument('-r', "--rerun",       dest="rerun", default=False, action='store_true',
                                           help="rerun crashed jobs (locally or on cluster), currently only works for generation" )
parser.add_argument('-d', "--diagrams",    dest="diagrams", default=False, action='store_true',
                                           help="recreate folders and find diagrams" )
parser.add_argument('-g', "--generate",    dest="generate", default=False, action='store_true',
                                           help="generate events in existing folders" )
# parser.add_argument('-x', "--xsec",        dest="xsec", default=False, action='store_true',
#                                            help="extract the cross sections from the generation logs" )
parser.add_argument('-e', "--extract",     dest="extract", default=False, action='store_true',
                                           help="extract the cross sections from the generation logs and store them in a json." )
# parser.add_argument('-n', "--ntuples",     dest="ntuples", default=False, action='store_true',
#                                            help="produce n-tuples from LHE files" )
parser.add_argument('-c', "--convert",     dest="convert", default=False, action='store_true',
                                           help="convert GENSIM to n-tuples." )
parser.add_argument('-p', "--plot",        dest="plot", default=False, action='store_true',
                                           help="plot distributions from n-tuples" )
args = parser.parse_args()
if args.local and args.submit:
    raise ValueError('Both \'local\' and \'submit\' arguments used, please decide for one.')
if not any([args.local, args.submit]) and not args.extract:
    raise ValueError('Neither \'local\' nor \'submit\' arguments used, this is only allowed when extracting the x-sec. Please decide for one.')

mgfolder = '/work/areimers/MG5_aMC_v2_7_2'
workarea = '/work/areimers/ZTo2Tau2Mu/smeft/madgraph'
plotfolder = os.path.join(workarea, 'plots')
ensureDirectory(plotfolder)




operators = []
# operators +=  ['cll2222', 'cle2222', 'cee2222']
operators += ['cll2233', 'cle2233', 'cee2233']
operators += ['cll2332', 'cle3322']
operators += [           'cle2332']

operators += [                       'cll2233+1p0xcle2233', 'cll2233+1p0xcee2233', 'cll2233+1p0xcll2332', 'cll2233+1p0xcle3322', 'cll2233+1p0xcle2332']
operators += ['cle2233+1p0xcll2233',                        'cle2233+1p0xcee2233', 'cle2233+1p0xcll2332', 'cle2233+1p0xcle3322', 'cle2233+1p0xcle2332']
operators += ['cee2233+1p0xcll2233', 'cee2233+1p0xcle2233',                        'cee2233+1p0xcll2332', 'cee2233+1p0xcle3322', 'cee2233+1p0xcle2332']
operators += ['cll2332+1p0xcll2233', 'cll2332+1p0xcle2233', 'cll2332+1p0xcee2233',                        'cll2332+1p0xcle3322', 'cll2332+1p0xcle2332']
operators += ['cle3322+1p0xcll2233', 'cle3322+1p0xcle2233', 'cle3322+1p0xcee2233', 'cle3322+1p0xcll2332',                        'cle3322+1p0xcle2332']
operators += ['cle2332+1p0xcll2233', 'cle2332+1p0xcle2233', 'cle2332+1p0xcee2233', 'cle2332+1p0xcll2332', 'cle2332+1p0xcle3322'                       ]

# operators += ['cll2233+1p0xcle2233', 'cll2233+1p0xcee2233', 'cll2233+1p0xcll2332', 'cll2233+1p0xcle3322', 'cll2233+1p0xcle2332']
# operators += [                       'cle2233+1p0xcee2233', 'cle2233+1p0xcll2332', 'cle2233+1p0xcle3322', 'cle2233+1p0xcle2332']
# operators += [                                              'cee2233+1p0xcll2332', 'cee2233+1p0xcle3322', 'cee2233+1p0xcle2332']
# operators += [                                                                     'cll2332+1p0xcle3322', 'cll2332+1p0xcle2332']
# operators += [                                                                                            'cle3322+1p0xcle2332']
# operators += ['cee2233+1p0xcle2233']


# baseprocs = ['pp_zttmm', 'pp_zmmmm']
# baseprocs = ['pp_incttmm']

# baseprocs = ['pp_zttmm'] # for convert/plot acceptances: only zttmm, not incttmm!
baseprocs = ['pp_incttmm', 'pp_zttmm', 'tmvv']
# baseprocs = ['tmvv']
orders = ['interference', 'purebsm']
# orders = ['interference']

translation_to_gensim = {
    'pp_zttmm_smeftsim': 'ZTo2Mu2Tau_SMEFT/ZTo2Mu2Tau_2TauTo2Mu_SMEFT',
    'pp_zmmmm_smeftsim': 'ZTo4Mu_SMEFT/ZTo4Mu_SMEFT',
    'pp_zttmm_oldsample_sm': 'ZTo4L/ZTo2Mu2Tau_2TauTo2Mu',
}

procnames = []

# procnames += ['pp_zttmm_smeftsim_sm', 'pp_zttmm_sm_sm']
# procnames += ['pp_zmmmm_smeftsim_sm', 'pp_zmmmm_sm_sm']
# procnames += ['pp_zttmm_smeftsim_sm']
# procnames += ['pp_zttmm_oldsample_sm']
# procnames += ['pp_zmmmm_smeftsim_sm']
# procnames += ['pp_incttmm_smeftsim_sm']
procnames += ['pp_incttmm_smeftsim_sm', 'pp_zttmm_smeftsim_sm', 'tmvv_smeftsim_sm']
# procnames += ['tmvv_smeftsim_sm']
procnames += ['%s_smeftsim_%s-%s' % (p, order, op) for p in baseprocs for order in orders for op in operators]
print procnames

r_upper_limit_abs = 6.2
r_sm = 0.9

n_gensimfiles = 100


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

    commandfiles_per_step = defaultdict(list)

    procnames_to_use_per_step = defaultdict(list)
    for procname in procnames:
        gensetup = GenerationSetup(procname=procname, mgfolder=mgfolder, workarea=workarea)
        gensetup.create_environment()
        if args.diagrams or args.generate:
            gensetup.make_cards()

        n_gensimfiles_thisproc = n_gensimfiles
        if procname in ['pp_zttmm_smeftsim_sm', 'pp_zttmm_oldsample_sm']:
            n_gensimfiles_thisproc = 1000

        # generating the diagrams and setting up the MG output folder for each process
        if args.diagrams:
            commandfiles_per_step['diagrams'].append(gensetup.make_command_diagrams())
            procnames_to_use_per_step['diagrams'].append(procname)

        # actually generating events for each process
        if args.rerun: # only append crashed jobs
            if args.generate:
                if not check_generation_successful(procname=procname, workarea=workarea):
                    commandfiles_per_step['generate'].append(gensetup.make_command_generate())
                    procnames_to_use_per_step['generate'].append(procname)
            
            missing_files, missing_indices = missing_conversion_files(procname=procname, workarea=workarea, nfiles_to_check=n_gensimfiles_thisproc)
            commands_convert_all = gensetup.make_command_convert(nfiles=n_gensimfiles_thisproc)
            commands_convert_resubmit = [x.replace('/ntuples/', '/commands/').replace('events_gensim', 'convert').replace('.root', '.sh') for x in missing_files]
            procnames_convert_resubmit = ['%s_%i' % (procname, idx) for idx in missing_indices]
            commandfiles_per_step['convert'] += commands_convert_resubmit
            procnames_to_use_per_step['convert'] += procnames_convert_resubmit

        elif not args.extract:
            if args.generate:
                commandfiles_per_step['generate'].append(gensetup.make_command_generate())
                procnames_to_use_per_step['generate'].append(procname)

            # command_ntuplize = gensetup.make_command_ntuplize()
            # if command_ntuplize is not None:
            #     commandfiles_per_step['ntuplize'].append(command_ntuplize)
            #     procnames_to_use_per_step['ntuplize'].append(procname)
            
            commands_convert = gensetup.make_command_convert(nfiles=n_gensimfiles_thisproc)
            if commands_convert is not None:
                commandfiles_per_step['convert'] += commands_convert
                procnames_to_use_per_step['convert'] += ['%s_%i' % (procname, idx+1) for idx in range(len(commands_convert))]

            command_plot = gensetup.make_command_plot()
            if command_plot is not None:
                commandfiles_per_step['plot'].append(command_plot)
                procnames_to_use_per_step['plot'].append(procname)



    if args.local:
        if args.diagrams:
            run(scriptname=commandfiles_per_step['diagrams'])
        if args.generate:
            run(scriptname=commandfiles_per_step['generate'])
    elif args.submit:
        if args.generate:
            submit(scriptname=commandfiles_per_step['generate'], procnames=procnames_to_use_per_step['generate'], submissionscript='submit_generic.sh', arguments='%s %s'%(os.environ['MGDIR'], os.environ['CMSSW_BASE']), runtime=(1,00,00), ncores=4)
        # elif args.ntuples:
        #     submit(scriptname=commandfiles_per_step['ntuplize'], procnames=procnames_to_use_per_step['ntuplize'], submissionscript='submit_generic.sh', arguments='%s %s'%(os.environ['MGDIR'], os.environ['CMSSW_BASE']), runtime=(0,10,00), ncores=1)
        elif args.convert:
            # print commandfiles_per_step['convert']
            # print procnames_to_use_per_step['convert']
            submit(scriptname=commandfiles_per_step['convert'], procnames=procnames_to_use_per_step['convert'], submissionscript='submit_generic.sh', arguments='%s %s'%(os.environ['MGDIR'], os.environ['CMSSW_BASE']), runtime=(1,00,00), ncores=2)
        elif args.plot:
            # print commandfiles_per_step['plot']
            submit(scriptname=commandfiles_per_step['plot'], procnames=procnames_to_use_per_step['plot'], submissionscript='submit_generic.sh', arguments='%s %s'%(os.environ['MGDIR'], os.environ['CMSSW_BASE']), runtime=(1,00,00), ncores=1)
        else:
            raise ValueError('Trying to submit something other than the \'generate\' step, which is not implemented yet.')




    if args.extract:
        # with the generation done, now extract the corrections to the SM operator-by-operator
        results = defaultdict(dict)
        coefficients = {}
        for o in operators:
            results[o] = defaultdict(dict)
        
        results['sm'] = defaultdict(dict)
        print '  --> Now reading out cross sections for all requested processes.'
        pbar = tqdm(procnames, desc="Cross sections read out")
        for p in pbar:

            # get values
            if p.startswith('tmvv'):
                (value, value_unc) = get_width_from_procname(procname=p, workarea=workarea)
                syst_up = syst_down = 0.
            else:
                (value, value_unc, syst_up, syst_down) = get_xsec_from_procname(procname=p, workarea=workarea, with_systs=False)

            # inject them into dict
            (proc, operator) = (p.split('_')[1], p.split('_')[-1])
            if p.startswith('tmvv'):
                proc = p.split('_')[0]
            if operator == 'sm':
                results[operator][proc]['sm'] = (value, value_unc, syst_up, syst_down)
            else:
                (order, operator) = (operator.split('-')[0], operator.split('-')[1])
                results[operator][proc][order] = (value, value_unc, syst_up, syst_down)
                

        
        for operator in results:
            if operator == 'sm': continue
            if not operator in coefficients: coefficients[operator] = {}
            for proc in results[operator]:
                if not proc in coefficients[operator]: coefficients[operator][proc] = {}
                for order in results[operator][proc]:
                    if order in coefficients[operator][proc]: raise ValueError('trying to overwrite a value, abort. Please fix, this should never happen.')
                    coefficients[operator][proc][order] = results[operator][proc][order][0] / results['sm'][proc]['sm'][0]

        with open(os.path.join(workarea, 'coefficients.json'), 'w') as j:
            json.dump(obj=coefficients, fp=j, indent=2, sort_keys=True)


    # if args.xsec:
    #     # with the generation done, now extract the corrections to the SM operator-by-operator
    #     results = defaultdict(dict)
    #     for o in operators:
    #         results[o] = defaultdict(dict)
        
    #     results['sm'] = defaultdict(dict)
    #     print '  --> Now reading out cross sections for all requested processes.'
    #     pbar = tqdm(procnames, desc="Cross sections read out")
    #     for p in pbar:

    #         # get values
    #         (xsec, xsec_unc, syst_up, syst_down) = get_xsec_from_procname(procname=p, workarea=workarea)

    #         # inject them into dict
    #         (proc, operator) = (p.split('_')[1], p.split('_')[-1])
    #         if operator == 'sm':
    #             results[operator][proc]['sm'] = (xsec, xsec_unc, syst_up, syst_down)
    #         else:
    #             (order, operator) = (operator.split('-')[0], operator.split('-')[1])
    #             results[operator][proc][order] = (xsec, xsec_unc, syst_up, syst_down)

    #     for op in operators:
    #         ratios_mmmm = defaultdict(dict)
    #         ratios_ttmm = defaultdict(dict)
    #         total_correction = defaultdict(dict)
    #         operator_values_around_limit = {}
    #         r_with_operator = OrderedDict()

    #         # operator_stepsize = 0.25 if '2222' not in op else -0.25
    #         operator_stepsize = 0.25 
    #         operator_values_fine = [float('%.2f' % (x*operator_stepsize)) for x in range(int(100000./abs(operator_stepsize)))]
    #         # if '2222' in op:
    #         pbar = tqdm(operator_values_fine, desc="Operator values tested")
    #         for operator_value in pbar:
    #             # if operator_value % 1000. == 0: print 'at value: %f' % operator_value

    #             if results['sm']['zmmmm']['sm'][0] != 0:
    #                 ratios_mmmm[operator_value]['interference'] = (operator_value) * results[op]['zmmmm']['interference'][0] / results['sm']['zmmmm']['sm'][0]
    #                 ratios_mmmm[operator_value]['purebsm'] = (operator_value)**2 * results[op]['zmmmm']['purebsm'][0] / results['sm']['zmmmm']['sm'][0]
    #             else: 
    #                 ratios_mmmm[operator_value]['interference'] = 0.
    #                 ratios_mmmm[operator_value]['purebsm'] = 0.

    #             if results['sm']['zttmm']['sm'][0] != 0:
    #                 ratios_ttmm[operator_value]['interference'] = operator_value * results[op]['zttmm']['interference'][0] / results['sm']['zttmm']['sm'][0]
    #                 ratios_ttmm[operator_value]['purebsm'] = operator_value**2 * results[op]['zttmm']['purebsm'][0] / results['sm']['zttmm']['sm'][0]
    #             else: 
    #                 ratios_ttmm[operator_value]['interference'] = 0.
    #                 ratios_ttmm[operator_value]['purebsm'] = 0.

    #             total_correction[operator_value]['interference'] = ratios_ttmm[operator_value]['interference'] - ratios_mmmm[operator_value]['interference']
    #             total_correction[operator_value]['purebsm'] = ratios_ttmm[operator_value]['purebsm'] - ratios_mmmm[operator_value]['purebsm']
    #             # total_correction[operator_value]['comb'] = total_correction[operator_value]['interference'] + total_correction[operator_value]['purebsm']

    #             # also print the limit on this operator by comparing its contribution with the exclusion of CMS-SMP-22-016
    #             # r_with_operator[operator_value] = {order: r_sm * (1+total_correction[operator_value][order]) for order in total_correction[operator_value]}
    #             if total_correction[operator_value][order] >= 0:
    #                 r_with_operator[operator_value] = {}
    #             else:
    #                 r_with_operator[-operator_value] = {}
    #             for order in total_correction[operator_value]:
    #                 if total_correction[operator_value][order] >= 0:
    #                     r_with_operator[operator_value][order] = r_sm * (1+total_correction[operator_value][order])
    #                 else:
    #                     r_with_operator[-operator_value][order] = r_sm * (1-total_correction[operator_value][order])


    #         # if '2222' not in op: val_to_print = 1.0
    #         # else: val_to_print = -1.0
    #         val_to_print = 1.0
    #         print ''
    #         print '  --> Corrections (SMEFT/SM) for operator %s = %f:' % (op, val_to_print)
    #         print '    --> Interference-only:', total_correction[abs(val_to_print)]['interference']
    #         # print '    --> Quadratic-only:   ', total_correction[abs(val_to_print)]['purebsm']
    #         # print '    --> Combined:         ', total_correction[abs(val_to_print)]['comb']
            
    #         # now find intersection with R_limit from CMS-SMP-22-016
    #         r_obs, r_exp, r_68_low, r_68_high, r_95_low, r_95_high = get_limits_r()
    #         # orders = ['interference', 'comb'] if '2222' not in op else ['interference']
    #         orders = ['interference', 'comb']
    #         if '2222' in op and 'x' not in op:
    #             orders = ['interference']
    #         #  if all(('2222' not in op, '+' not in op, 'x' not in op)) else ['interference']
    #         for order in orders:
    #         # for order in ['interference']:
    #             if total_correction[val_to_print][order] == 0.: continue
    #             limits_obs = []
    #             limits_exp = []
    #             limits_68_low  = []
    #             limits_68_high = []
    #             limits_95_low  = []
    #             limits_95_high = []
    #             pbar = tqdm(range(len(r_obs)), desc="CLs checked")
    #             for ipoint in pbar:
    #                 # print '  --> at CL-point no. %i' % (ipoint)
    #                 # print '  --> checking limit for CL = %f, corresponding limit on R = %f' % (cl, r_obs[ipoint][0])
    #                 # if (1-r_obs[ipoint][1])*100. > 50:
    #                 if (1-r_obs[ipoint][1])*100. > 0.:
    #                     limit_obs = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_obs[ipoint][0])
    #                     if limit_obs is not None:
    #                         limits_obs.append((limit_obs, r_obs[ipoint][1]))

    #                 # if (1-r_exp[ipoint][1])*100. > 50:
    #                 if (1-r_exp[ipoint][1])*100. > 0.:
    #                     limit_exp = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_exp[ipoint][0])
    #                     if limit_exp is not None:
    #                         limits_exp.append((limit_exp, r_exp[ipoint][1]))

    #                 # if (1-r_68_low[ipoint][1])*100. > 50:
    #                 if (1-r_68_low[ipoint][1])*100. > 0.:
    #                     limit_68_low = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_68_low[ipoint][0])
    #                     if limit_68_low is not None:
    #                         limits_68_low.append((limit_68_low, r_68_low[ipoint][1]))

    #                 # if (1-r_68_high[ipoint][1])*100. > 50:
    #                 if (1-r_68_high[ipoint][1])*100. > 0.:
    #                     limit_68_high = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_68_high[ipoint][0])
    #                     if limit_68_high is not None:
    #                         limits_68_high.append((limit_68_high, r_68_high[ipoint][1]))
                    
    #                 # if (1-r_95_low[ipoint][1])*100. > 50:
    #                 if (1-r_95_low[ipoint][1])*100. > 0.:
    #                     limit_95_low = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_95_low[ipoint][0])
    #                     if limit_95_low is not None:
    #                         limits_95_low.append((limit_95_low, r_95_low[ipoint][1]))

    #                 # if (1-r_95_high[ipoint][1])*100. > 50:
    #                 if (1-r_95_high[ipoint][1])*100. > 0.:
    #                     limit_95_high = get_upper_limit_on_operator(r_with_operator=r_with_operator, order=order, r_upper_limit=r_95_high[ipoint][0])
    #                     if limit_95_high is not None:
    #                         limits_95_high.append((limit_95_high, r_95_high[ipoint][1]))
    #                 # if limit_exp is not None:
    #                 #     print '  --> Upper limit from %s at %f%% CL on %s: %f + %f' % (order, (1-r_exp[ipoint][1])*100., op, limit_exp, limit_95_high)
    #                 # else:
    #                 #     print '  --> No upper limit found from %s at %f%% CL on %s' % (order, (1-r_obs[ipoint][1])*100., op)

    #             limits_obs.sort(key=lambda tup: tup[1])
    #             limits_exp.sort(key=lambda tup: tup[1])
    #             limits_68_low.sort(key=lambda tup: tup[1])
    #             limits_68_high.sort(key=lambda tup: tup[1], reverse=True)
    #             limits_95_low.sort(key=lambda tup: tup[1])
    #             limits_95_high.sort(key=lambda tup: tup[1], reverse=True)
    #             limits_68 = limits_68_low + limits_68_high
    #             limits_95 = limits_95_low + limits_95_high
    #             limits_obs_cleaned = []
    #             limits_exp_cleaned = []
    #             limits_68_cleaned = []
    #             limits_95_cleaned = []

    #             # is_limit_negative = False
    #             for il, l in enumerate(limits_obs):
    #                 if il == 0: 
    #                     if l[0] < 0: 
    #                         limits_obs_cleaned.append((-l[0], l[1]))
    #                         # is_limit_negative = True
    #                     else: limits_obs_cleaned.append(l)
    #                 if l[0] == operator_values_fine[1] and limits_obs[il-1][0] == l[0]: continue
    #                 if l[0] < 0: 
    #                     limits_obs_cleaned.append((-l[0], l[1]))
    #                     # is_limit_negative = True
    #                 else: limits_obs_cleaned.append(l)
    #             for il, l in enumerate(limits_exp):
    #                 if il == 0: 
    #                     if l[0] < 0: limits_exp_cleaned.append((-l[0], l[1]))
    #                     else:        limits_exp_cleaned.append(l)
    #                 if l[0] == operator_values_fine[1] and limits_exp[il-1][0] == l[0]: continue
    #                 if l[0] < 0: limits_exp_cleaned.append((-l[0], l[1]))
    #                 else:        limits_exp_cleaned.append(l)
    #             for l in limits_68:
    #                 if l[0] < 0: limits_68_cleaned.append((-l[0], l[1]))
    #                 else:        limits_68_cleaned.append(l)
    #             for l in limits_95:
    #                 if l[0] < 0: limits_95_cleaned.append((-l[0], l[1]))
    #                 else:        limits_95_cleaned.append(l)

    #             # scan x-points. At each x-value, find upper and lower limit on "1-CL"
    #             g_exp = ROOT.TGraph(len(limits_exp_cleaned), array('d', zip(*limits_exp_cleaned)[0]), array('d', zip(*limits_exp_cleaned)[1]))
    #             g_obs = ROOT.TGraph(len(limits_obs_cleaned), array('d', zip(*limits_obs_cleaned)[0]), array('d', zip(*limits_obs_cleaned)[1]))
    #             g_68 = ROOT.TGraph(len(limits_68_cleaned), array('d', zip(*limits_68_cleaned)[0]), array('d', zip(*limits_68_cleaned)[1]))
    #             g_95 = ROOT.TGraph(len(limits_95_cleaned), array('d', zip(*limits_95_cleaned)[0]), array('d', zip(*limits_95_cleaned)[1]))

    #             xmax = 1E4 if abs(limits_95_high[-1][0]) < 1E4 else 1E5
    #             # if limits_95_high[-1][0] < 0:
    #             #     xmax = -6E3 if limits_95_high[-1][0] > -4E3 else -6E4
    #             maxdigits = 2 if abs(xmax) <= 1E4 else 3
    #             # xtitle = get_pretty_operatorname_with_factors(input_string=op) if limits_95_high[-1][0] >= 0 else '#minus %s' % (get_pretty_operatorname_with_factors(input_string=op))
    #             xtitle = get_pretty_operatorname_with_factors(input_string=op)
    #             if op in ['cll2222', 'cle2222', 'cee2222']:
    #                 xtitle = '#minus %s' % (xtitle)
    #             c = tdrCanvas(canvName='c', x_min=min(0., xmax), x_max=max(0., xmax), y_min=5E-3, y_max=5., nameXaxis=xtitle, nameYaxis='1 #minus CL', square=True, iPos=11, margins=(None, 0.10, 0.14, None), maxdigits=(maxdigits, None))
    #             c.SetLogy()
    #             # c.SetLogx()

    #             if xmax > 0: 
    #                 leg = tdrLeg(0.52,0.68,0.87,0.9, textSize=0.038)
    #             else:
    #                 leg = tdrLeg(0.18,0.53,0.53,0.75, textSize=0.038)
    #             tdrHeader(leg, 'Upper limits')
    #             tdrDraw(g_95, "F SAME", mcolor=ROOT.kOrange, lcolor=ROOT.kOrange, fcolor=ROOT.kOrange)
    #             tdrDraw(g_68, "F SAME", mcolor=ROOT.kGreen+1, lcolor=ROOT.kGreen+1, fcolor=ROOT.kGreen+1)
    #             tdrDraw(g_exp, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=2)
    #             tdrDraw(g_obs, "L SAME", mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, fcolor=ROOT.kBlack, lstyle=1)
    #             leg.AddEntry(g_obs, 'Observed', 'L')
    #             leg.AddEntry(g_exp, 'Median expected', 'L')
    #             leg.AddEntry(g_68, '68% expected', 'LF')
    #             leg.AddEntry(g_95, '95% expected', 'LF')
    #             leg.Draw('SAME')
    #             c.RedrawAxis()
                
    #             outname = os.path.join(plotfolder, 'Limits_%s_%s.pdf' % (op, order))
    #             c.SaveAs(outname)

        


def get_upper_limit_on_operator(r_with_operator, order, r_upper_limit):
    for iop, operator_value in enumerate(r_with_operator):
        if iop == 0: continue
        # if r_with_operator[operator_value][order] > r_upper_limit:
        #     return 0
        # if iop % 1000 == 0 or iop == 1: 
        #     print '  --> Trying operator value: %f' % (operator_value)
        #     print r_with_operator[operator_value][order], r_with_operator[r_with_operator.keys()[iop-1]][order], r_upper_limit
        if r_with_operator[operator_value][order] > r_upper_limit and (r_with_operator[r_with_operator.keys()[iop-1]][order] <= r_upper_limit or iop == 1):
            # if iop == 1: return 0.
            # print '  --> Found limit: %f (R here: %f, point before: %f, SMP: %f)' % (operator_value, r_with_operator[operator_value][order], r_with_operator[r_with_operator.keys()[iop-1]][order], r_upper_limit)
            return operator_value
    return None


def get_xsec_from_procname(procname, workarea, with_systs=False):
    file_to_parse = os.path.join(workarea, 'logs', procname, 'generate.log')
    if not check_generation_successful(procname=procname, workarea=workarea): return (None, None, None, None)
    with open(file_to_parse, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'Survey return zero cross section' in line: return (0., 0., 0., 0.)
        if 'Cross-section :' in line: line_crosssection = line
        elif 'scale variation:' in line: line_scale = line
        elif 'PDF variation:' in line: line_pdf = line

    # find correct lines to parse and set up the patterns
    pattern_crosssection = '{}:   {} +- {} pb\n'
    pattern_scale        = '{}: +{}% -{}%\n'
    pattern_pdf          = '{}: +{}% -{}%\n'

    # do the parsing
    parser_xsec  = parse.compile(pattern_crosssection)
    parser_scale = parse.compile(pattern_scale)
    parser_pdf   = parse.compile(pattern_pdf)

    pre, xsec, xsec_unc = parser_xsec.parse(line_crosssection)
    xsec, xsec_unc      = float(xsec), float(xsec_unc)
    if with_systs:
        pre, scale_up, scale_down = parser_scale.parse(line_scale)
        pre, pdf_up, pdf_down = parser_pdf.parse(line_pdf)
        scale_up, scale_down, pdf_up, pdf_down = float(scale_up), float(scale_down), float(pdf_up), float(pdf_down)
        syst_up = math.sqrt((scale_up/100.*xsec)**2 + (pdf_up/100.*xsec)**2)
        syst_down = math.sqrt((scale_down/100.*xsec)**2 + (pdf_down/100.*xsec)**2)
    else: 
        syst_up = syst_down = None

    # xsec, xsec_unc, scale_up, scale_down, pdf_up, pdf_down = float(xsec), float(xsec_unc), float(scale_up), float(scale_down), float(pdf_up), float(pdf_down)

    # syst_up = math.sqrt((scale_up/100.*xsec)**2 + (pdf_up/100.*xsec)**2)
    # syst_down = math.sqrt((scale_down/100.*xsec)**2 + (pdf_down/100.*xsec)**2)
    return (xsec, xsec_unc, syst_up, syst_down)


def get_width_from_procname(procname, workarea):
    file_to_parse = os.path.join(workarea, 'logs', procname, 'generate.log')
    if not check_generation_successful_width(procname=procname, workarea=workarea): return (None, None, None, None)
    with open(file_to_parse, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'Survey return zero cross section' in line: return (0., 0.)
        if 'Width :' in line: line_width = line

    # find correct lines to parse and set up the patterns
    pattern_width = '{}:   {} +- {} GeV\n'

    # do the parsing
    parser_width  = parse.compile(pattern_width)

    pre, width, width_unc = parser_width.parse(line_width)
    width, width_unc      = float(width), float(width_unc)

    print (width, width_unc)
    return (width, width_unc)


def check_generation_successful_width(procname, workarea):
    # checks generate.log for unexpected crashes, i.e. being unable to extract the xsec. However, crashing because of 0 xsec is OK.

    file_to_parse = os.path.join(workarea, 'logs', procname, 'generate.log')
    with open(file_to_parse, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'Survey return zero cross section' in line: return True
        if 'Width :' in line: return True
    return False


def check_generation_successful(procname, workarea, need_systs=False):
    # checks generate.log for unexpected crashes, i.e. being unable to extract the xsec. However, crashing because of 0 xsec is OK.

    file_to_parse = os.path.join(workarea, 'logs', procname, 'generate.log')
    with open(file_to_parse, 'r') as f:
        lines = f.readlines()

    found_xsec  = False
    found_scale = False
    found_pdf   = False
    for line in lines:
        if 'Survey return zero cross section' in line: return True
        if 'Cross-section :' in line: found_xsec = True
        elif 'scale variation:' in line: found_scale = True
        elif 'PDF variation:' in line: found_pdf = True

    if not need_systs: return found_xsec
    else:              return all([found_xsec, found_scale, found_pdf])

def missing_conversion_files(procname, workarea, nfiles_to_check=100):
    # checks generate.log for unexpected crashes, i.e. being unable to extract the xsec. However, crashing because of 0 xsec is OK.
    missing_files = []
    missing_indices = []
    for idx in range(nfiles_to_check):
        file_to_check = os.path.join(workarea, 'ntuples', procname, 'events_gensim_%i.root' % (idx+1))
        if not is_file_large_enough(filename=file_to_check, min_size=51200):
            missing_files.append(file_to_check)
            missing_indices.append(idx+1)
    return missing_files, missing_indices

def is_file_large_enough(filename, min_size):
    if os.path.exists(filename):
        file_size = os.path.getsize(filename)
        if file_size >= min_size:
            return True
    return False

def get_pretty_operatorname_with_factors(input_string, translation_dict=operatornames_pretty):
    # Delimiters for separating operators and factors within each operator
    operator_delimiter = '+'
    factor_delimiter = 'x'

    # Split input string into individual operators
    operator_parts = input_string.split(operator_delimiter)
    translated_and_converted_operators = []

    for operator in operator_parts:
        # Split operator into individual factors
        factor_parts = operator.split(factor_delimiter)
        translated_and_converted_factors = []

        for factor in factor_parts:
            if factor in translation_dict:
                # Translate factor if found in the dictionary
                translated_and_converted_factors.append(translation_dict[factor])
            elif 'p' in factor:
                # Convert factor with 'p' to float format (replace 'p' with '.')
                float_value = factor.replace('p', '.')
                translated_and_converted_factors.append(float_value)
            else:
                # Append unchanged factor
                translated_and_converted_factors.append(factor)

        # Join translated and converted factors with the factor delimiter
        translated_and_converted_operators.append(factor_delimiter.join(translated_and_converted_factors))

    # Join translated and converted operators with the operator delimiter
    translated_and_converted_string = operator_delimiter.join(translated_and_converted_operators)
    return translated_and_converted_string


class GenerationSetup:
    def __init__(self, procname, mgfolder, workarea):
        self.procname = procname
        self.mgfolder = mgfolder
        self.mgexec   = os.path.join(self.mgfolder, 'bin', 'mg5_aMC')
        self.workarea = workarea
        self.cardfolder_baseline = os.path.join(self.workarea, 'cards', 'baseline') # contains pre-made cards that can be copied and then modified
        self.cardfolder = os.path.join(self.workarea, 'cards', self.procname)
        self.outputfolder = os.path.join(self.workarea, 'output', self.procname)
        self.commandfolder = os.path.join(self.workarea, 'commands', self.procname)
        self.logfolder = os.path.join(self.workarea, 'logs', self.procname)
        self.lhefolder = os.path.join(self.workarea, 'lheevents', self.procname)
        self.ntuplefolder = os.path.join(self.workarea, 'ntuples', self.procname)
        # self.gensimfolder = os.path.join('root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/areimers/GENSIM', self.procname.replace('pp_zmmmm_smeftsim', translation_to_gensim['pp_zmmmm_smeftsim']).replace('pp_zttmm_smeftsim', translation_to_gensim['pp_zttmm_smeftsim']).replace('pp_zttmm_oldsample_sm', translation_to_gensim['pp_zttmm_oldsample_sm']).replace('-c', '_c'))
        self.gensimfolder = os.path.join('root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/areimers/GENSIM', self.procname.replace('pp_zmmmm_smeftsim', translation_to_gensim['pp_zmmmm_smeftsim']).replace('pp_zttmm_smeftsim', translation_to_gensim['pp_zttmm_smeftsim']).replace('pp_zttmm_oldsample_sm', translation_to_gensim['pp_zttmm_oldsample_sm'])).split('-')[0]
        operatorpart = self.procname.split('-')[-1]
        
        
        # if not self.gensimfolder.endswith('_sm') and not 'ZTo4L' in self.gensimfolder:
        #     self.gensimfolder += '-1p0'
        if not self.gensimfolder.endswith('_sm') and not 'ZTo4L' in self.gensimfolder:
            operators_in_operatorpart = operatorpart.split('+')
            operatorstring = ''
            for o in operators_in_operatorpart:
                if 'x' in o and not '1p0x' in o:
                    raise ValueError('there is a factor in one of the operator strings (like 0p5xcll2233), but it is not 1p0x. This is what we assume at the moment, please check.')
                o = o.replace('1p0x', '')
                operatorstring += '_%s-1p0' % (o)
            self.gensimfolder += operatorstring
        self.plotfolder = os.path.join(self.workarea, 'plots', self.procname)
        self.genexec = os.path.join(self.outputfolder, 'bin', 'generate_events')
        self.operatorsettings = self.get_operatorsettings()


        self.genstrings_per_proc = {
            'pp_zttmm_smeftsim_sm':   'generate p p > z > mu+ mu- ta+ ta- / h NP=0 SMHLOOP=0 NPprop=0', 
            'pp_incttmm_smeftsim_sm': 'generate p p > mu+ mu- ta+ ta- / h NP=0 SMHLOOP=0 NPprop=0', 
            'pp_zttmm_sm_sm':         'generate p p > z > mu+ mu- ta+ ta- / h ', 
            'pp_zmmmm_smeftsim_sm':   'generate p p > z > mu+ mu- mu+ mu- / h NP=0 SMHLOOP=0 NPprop=0', 
            'pp_zmmmm_sm_sm':         'generate p p > z > mu+ mu- mu+ mu- / h ', 
            'tmvv_smeftsim_sm':       'generate ta- > mu- vt vm~ / h NP=0 SMHLOOP=0 NPprop=0', 
            'tmvv_sm_sm':             'generate ta- > mu- vt vm~ / h ', 

            'pp_zttmm_smeftsim_interference':   'generate p p > z > mu+ mu- ta+ ta- / h NP<=1 NP^2==1 SMHLOOP=0 NProp=0',
            'pp_zttmm_smeftsim_purebsm':        'generate p p > z > mu+ mu- ta+ ta- / h NP==1 SMHLOOP=0 NProp=0',
            'pp_incttmm_smeftsim_interference': 'generate p p > mu+ mu- ta+ ta- / h NP<=1 NP^2==1 SMHLOOP=0 NProp=0',
            'pp_incttmm_smeftsim_purebsm':      'generate p p > mu+ mu- ta+ ta- / h NP==1 SMHLOOP=0 NProp=0',
            'pp_zmmmm_smeftsim_interference':   'generate p p > z > mu+ mu- mu+ mu- / h NP<=1 NP^2==1 SMHLOOP=0 NProp=0',
            'pp_zmmmm_smeftsim_purebsm':        'generate p p > z > mu+ mu- mu+ mu- / h NP==1 SMHLOOP=0 NProp=0',
            'tmvv_smeftsim_interference':       'generate ta- > mu- vt vm~ / h NP<=1 NP^2==1 SMHLOOP=0 NProp=0',
            'tmvv_smeftsim_purebsm':            'generate ta- > mu- vt vm~ / h NP==1 SMHLOOP=0 NProp=0',

        }
        self.modelnames_per_proc = {
            'pp_zttmm_smeftsim_sm':   'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_incttmm_smeftsim_sm': 'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_zttmm_sm_sm':         'sm', 
            'pp_zmmmm_smeftsim_sm':   'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_zmmmm_sm_sm':         'sm', 
            'tmvv_smeftsim_sm':       'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'tmvv_sm_sm':             'sm', 

            'pp_zttmm_smeftsim_interference':   'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_zttmm_smeftsim_purebsm':        'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_incttmm_smeftsim_interference': 'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_incttmm_smeftsim_purebsm':      'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_zmmmm_smeftsim_interference':   'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'pp_zmmmm_smeftsim_purebsm':        'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'tmvv_smeftsim_interference':       'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
            'tmvv_smeftsim_purebsm':            'SMEFTsim_top_alphaScheme_UFO-Zmmtt', 
        }

        self.param_settings_common = OrderedDict() # All MG_2.7.2 default from model 'sm'
        self.param_settings_common['mass 5'] =    4.700000e+00 # MB
        self.param_settings_common['mass 6'] =    1.730000e+02 # MT 
        self.param_settings_common['mass 15'] =   1.777000e+00 # MTA
        self.param_settings_common['mass 23'] =   9.118800e+01 # MZ
        self.param_settings_common['mass 25'] =   1.250000e+02 # MH  

        self.param_settings_common['yukawa 5'] =  4.700000e+00 # ymb 
        self.param_settings_common['yukawa 6'] =  1.730000e+02 # ymt 
        self.param_settings_common['yukawa 15'] = 1.777000e+00 # ymtau 

        self.param_settings_common['DECAY 6'] =   1.491500e+00 # WT 
        self.param_settings_common['DECAY 23'] =  2.441404e+00 # WZ 
        self.param_settings_common['DECAY 24'] =  2.047600e+00 # WW 
        self.param_settings_common['DECAY 25'] =  6.382339e-03 # WH 

        # SMEFT Operators, 0 by default. Will be overridden by specific settings, because those will be applied after these common ones.
        self.param_settings_common['smeft 5'] =  0   # chdd
        self.param_settings_common['smeft 9'] =  0   # chwb
        self.param_settings_common['smeft 110'] =  0 # chl122
        self.param_settings_common['smeft 111'] =  0 # chl133
        self.param_settings_common['smeft 112'] =  0 # chl311
        self.param_settings_common['smeft 113'] =  0 # chl322
        self.param_settings_common['smeft 114'] =  0 # chl333
        self.param_settings_common['smeft 116'] =  0 # che22
        self.param_settings_common['smeft 117'] =  0 # che33
        self.param_settings_common['smeft 119'] =  0 # cll2222
        self.param_settings_common['smeft 123'] =  0 # cll2233
        self.param_settings_common['smeft 124'] =  0 # cll1221
        self.param_settings_common['smeft 126'] =  0 # cll2332
        self.param_settings_common['smeft 140'] =  0 # cee2222
        self.param_settings_common['smeft 144'] =  0 # cee2233
        self.param_settings_common['smeft 176'] =  0 # cle2222
        self.param_settings_common['smeft 180'] =  0 # cle2233
        self.param_settings_common['smeft 183'] =  0 # cle3322
        self.param_settings_common['smeft 186'] =  0 # cle2332

        self.operator_to_param_translation = {
            'cll2222': 'smeft 119',
            'cll2233': 'smeft 123',
            'cll2332': 'smeft 126',
            'cee2222': 'smeft 140',
            'cee2233': 'smeft 144',
            'cle2222': 'smeft 176',
            'cle2233': 'smeft 180',
            'cle3322': 'smeft 183',
            'cle2332': 'smeft 186'
        }

        self.param_settings_specific = OrderedDict()
        if self.operatorsettings != 'sm':
            for operator in self.operatorsettings:
                self.param_settings_specific[self.operator_to_param_translation[operator]] = self.operatorsettings[operator]



    def create_environment(self):
        ensureDirectory(self.cardfolder)
        ensureDirectory(self.outputfolder)
        ensureDirectory(self.commandfolder)
        ensureDirectory(self.logfolder)
        ensureDirectory(self.lhefolder)
        ensureDirectory(self.ntuplefolder)     
        ensureDirectory(self.plotfolder)    

    def get_genstring(self):
        count, matching_keys = count_and_get_matching_keys(dictionary=self.genstrings_per_proc, target_string=self.procname)
        if count == 1:
            return self.genstrings_per_proc[matching_keys[0]]
        elif count == 0:
            raise ValueError('For process \'%s\', no genstring (MG command for proc card) has been defined yet, please add.' % (self.procname))
        else:
            print 'more than one matching key: ', matching_keys
            raise ValueError('For process \'%s\', more than one genstring key matches, please resolve.' % (self.procname))
        

    def get_modelname(self):
        count, matching_keys = count_and_get_matching_keys(dictionary=self.modelnames_per_proc, target_string=self.procname)
        if count == 1:
            return self.modelnames_per_proc[matching_keys[0]]
        elif count == 0:
            raise ValueError('For process \'%s\', no modelname (UFO to be imported) has been defined yet, please add.' % (self.procname))
        else:
            print 'more than one matching key: ', matching_keys
            raise ValueError('For process \'%s\', more than one modelname key matches, please resolve.' % (self.procname))
            
        # if not self.procname in self.modelnames_per_proc:
        #     raise ValueError('For process \'%s\', no modelname (to be used by MG) has been defined yet, please add.' % (self.procname))
        # return self.modelnames_per_proc[self.procname]

    def get_operatorsettings(self):
        operatorstring = self.procname.split('_')[-1]
        if operatorstring == 'sm':
            return 'sm'
        operatorstring = operatorstring.split('-')[-1]
        
        operatorsettings = {}
        terms = operatorstring.split('+')
        for t in terms:
            if len(t.split('x')) > 1:
                factor = convert_string_to_float(input_string = t.split('x')[0])
                op = t.split('x')[1]
            else:
                factor = 1.0
                op = t
            if op in operatorsettings:
                raise ValueError('trying to set operator %s to value %f but it already exists with value %f. Complete string was: %s' % (op, factor, operatorsettings[op], operatorstring))
            operatorsettings[op] = factor
        return operatorsettings

    def make_cards(self):
        self.make_proccard()
        self.make_runcard()
        self.make_gencard()

    def make_proccard(self):
        cardname = os.path.join(self.cardfolder, 'proc_card.dat')

        proccard_text = []
        proccard_text.append('import %s\n' % (self.get_modelname()))
        proccard_text.append('define p = 21 2 4 1 3 -2 -4 -1 -3 5 -5 # pass to 5 flavors\n')
        proccard_text.append('define j = p\n')
        proccard_text.append('%s\n' % (self.get_genstring()))
        proccard_text.append('output %s\n' % (self.outputfolder))
        proccard_text.append('y\n')
        proccard_text.append('exit')
        with open(cardname, 'w') as f:
            f.writelines(proccard_text)
        self.proccard = cardname

    def make_gencard(self):
        cardname = os.path.join(self.cardfolder, 'gen_card.dat')

        gencard_text = []
        gencard_text.append('done\n')

        # now "customizecard" content modifying run_card and param_card: common settings
        gencard_text += ['set param_card %s %f\n' % (key, self.param_settings_common[key]) for key in self.param_settings_common.keys()]

        # now "customizecard" content modifying run_card and param_card: proc-specific settings
        # gencard_text += ['set param_card %s %f\n' % (key, self.param_settings_specific[self.procname][key]) for key in self.param_settings_specific[self.procname].keys()]
        gencard_text += ['set param_card %s %f\n' % (key, self.param_settings_specific[key]) for key in self.param_settings_specific.keys()]

        # for the tau decay, remove the PDF (reset to MG default) and switch off cuts
        if self.procname.startswith('tmvv'):
            gencard_text.append('set run_card pdlabel nn23lo1\n')
            gencard_text.append('set run_card lhaid 230000\n')
            gencard_text.append('set run_card ptl 0.0\n')
            gencard_text.append('set run_card mmnl 0.0\n')
            gencard_text.append('set run_card mmnlmax -1.0\n')
            gencard_text.append('set run_card mmll 0.0\n')


        # print self.operatorsettings
        # print gencard_text


        gencard_text.append('done')
        with open(cardname, 'w') as f:
            f.writelines(gencard_text)
        self.gencard = cardname
    
    def make_runcard(self):
        runcard_to_use = 'run_card.dat'
        if 'forkinematics' in self.procname:
            runcard_to_use = 'run_card_forkinematics.dat'
        cardname = os.path.join(self.cardfolder, 'run_card.dat')
        cardname_baseline = os.path.join(self.cardfolder_baseline, runcard_to_use)
        execute_command_silent(command=get_copy_command(from_filename=cardname_baseline, to_filename=cardname))

    def make_command_diagrams(self):
        command = 'cat %s | %s' % (self.proccard, self.mgexec)
        commandfilename = os.path.join(self.commandfolder, 'diagrams.sh')
        with open(commandfilename, 'w') as f:
            f.write(command)
        return commandfilename

    def make_command_generate(self):
        commands = []

        # copy run_card into MG output folder for generating
        from_filename = os.path.join(self.cardfolder, 'run_card.dat') # here all the custom-made cards are stored
        to_filename   = os.path.join(self.outputfolder, 'Cards', 'run_card.dat') # this is the 'Cards' folder in the MG output directory
        commands.append('%s\n' % (get_copy_command(from_filename=from_filename, to_filename=to_filename)))
        commands.append('cd %s\n' % (self.outputfolder))
        
        # run the generation
        commands.append('cat %s | %s\n' % (self.gencard, self.genexec))

        # delete the LHE events
        commands.append('rm -rf %s\n' % (os.path.join(self.outputfolder, 'Events')))
        commands.append('rm -rf %s\n' % (os.path.join(self.outputfolder, 'SubProcesses')))
        commands.append('rm -rf %s\n' % (os.path.join(self.outputfolder, 'lib')))
        commands.append('rm -rf %s\n' % (os.path.join(self.outputfolder, 'bin')))
        commandfilename = os.path.join(self.commandfolder, 'generate.sh')
        with open(commandfilename, 'w') as f:
            f.writelines(commands)
        return commandfilename

    # def make_command_ntuplize(self):
    #     commands = []
    #     packed_events = os.path.join(self.outputfolder, 'Events', 'run_01', 'unweighted_events.lhe.gz') 
    #     unpacked_events = os.path.join(self.lhefolder, 'events.lhe')
    #     root_events = os.path.join(self.ntuplefolder, 'events_lhe.root')

    #     if not os.path.isfile(packed_events): return None
    #     if os.path.isfile(unpacked_events): commands.append('rm %s\n' % unpacked_events)
    #     if os.path.isfile(root_events):     commands.append('rm %s\n' % root_events)

    #     command_unpack = 'gunzip -c %s > %s\n' % (packed_events, unpacked_events)
    #     commands.append(command_unpack)

    #     command_convert = 'python %s -i %s -o %s\n' % (os.path.join(self.workarea, 'convert_lhe_root.py'), unpacked_events, root_events)
    #     commands.append(command_convert)

    #     if os.path.isfile(unpacked_events): commands.append('rm %s\n' % unpacked_events)
        
    #     commandfilename = os.path.join(self.commandfolder, 'ntuplize.sh')
    #     with open(commandfilename, 'w') as f:
    #         f.writelines(commands)
    #     return commandfilename

    # def make_command_convert(self, nfiles=100):
    #     commands = []
    #     #  if os.path.isfile(os.path.join(self.gensimfolder+'-1p0', 'GENSIM_%i.root' % (ifile+1)))
    #     # filenames_gensim = [os.path.join(self.gensimfolder+'-1p0', 'GENSIM_%i.root' % (ifile+1)) for ifile in range(nfiles)]
    #     filenames_gensim = [os.path.join(self.gensimfolder, 'GENSIM_%i.root' % (ifile+1)) for ifile in range(nfiles)]
    #     gensim_events = ' '.join(filenames_gensim)
    #     root_events = os.path.join(self.ntuplefolder, 'events_gensim.root')
    #     # print self.ntuplefolder
    #     # print filenames_gensim
    #     # print gensim_events

    #     if os.path.isfile(root_events):     commands.append('rm %s\n' % root_events)

    #     command_convert = 'python %s -i %s -o %s\n' % (os.path.join(self.workarea, 'convert_gensim_root.py'), gensim_events, root_events)

    #     # print command_convert
    #     commands.append(command_convert)

    #     commandfilename = os.path.join(self.commandfolder, 'convert.sh')
    #     with open(commandfilename, 'w') as f:
    #         f.writelines(commands)
    #     return commandfilename

    def make_command_convert(self, nfiles=100):
        commandfilenames = []
        for i in range(nfiles):
            commands = []
            gensim_events = os.path.join(self.gensimfolder, 'GENSIM_%i.root' % (i+1))
            root_events = os.path.join(self.ntuplefolder, 'events_gensim_%i.root' % (i+1))
            if os.path.isfile(root_events):     
                commands.append('rm %s\n' % root_events)
            command_convert = 'python %s -i %s -o %s\n' % (os.path.join(self.workarea, 'convert_gensim_root.py'), gensim_events, root_events)

            # print command_convert
            commands.append(command_convert)

            commandfilename = os.path.join(self.commandfolder, 'convert_%i.sh' % (i+1))
            with open(commandfilename, 'w') as f:
                f.writelines(commands)
            commandfilenames.append(commandfilename)
        return commandfilenames

    def make_command_plot(self, filetype='gensim', nfiles=1000):
        commands = []
        if filetype == 'gensim':
            filename_sm = []
            filename_smeft = []
            for idx in range(nfiles):
                fn_sm = os.path.join('_'.join(self.ntuplefolder.split('_')[:-1] + ['sm']), 'events_%s_%i.root' % (filetype, idx+1)).replace('_oldsample', '_smeftsim')
                fn_smeft = os.path.join(self.ntuplefolder, 'events_%s_%i.root' % (filetype, idx+1))
                if is_file_large_enough(filename=fn_sm, min_size=51200):    filename_sm.append(fn_sm)
                if is_file_large_enough(filename=fn_smeft, min_size=51200): filename_smeft.append(fn_smeft)

        else:
            filename_sm    = os.path.join('_'.join(self.ntuplefolder.split('_')[:-1] + ['sm']), 'events_%s.root' % (filetype))
            filename_smeft = os.path.join(self.ntuplefolder, 'events_%s.root' % (filetype))
            if not os.path.isfile(filename_sm) or not os.path.isfile(filename_smeft): return None

        if isinstance(filename_sm, list):
            filename_sm = ' '.join(filename_sm)
        if isinstance(filename_smeft, list):
            filename_smeft = ' '.join(filename_smeft)
        command_plot = 'python %s -r %s -a %s -o %s -b %s\n' % (os.path.join(self.workarea, 'plot_ntuples.py'), filename_sm, filename_smeft, self.plotfolder, self.workarea)
        commands.append(command_plot)

        commandfilename = os.path.join(self.commandfolder, 'plot.sh')
        with open(commandfilename, 'w') as f:
            f.writelines(commands)
        return commandfilename

def run(scriptname):
    if not isinstance(scriptname, list):
        scriptname = [scriptname]
    commands = ['source %s' % (sn) for sn in scriptname]
    logfiles = [os.path.join(os.path.dirname(sn).replace('commands', 'logs'), sn.split('/')[-1].replace('.sh', '.log')) for sn in scriptname]
    print(blue('--> Executing %i commands\n\n' % (len(commands))))
    parallelize(commands=commands, logfiles=logfiles, niceness=None, remove_temp_files=False)
    print(green('\n\n--> Done executing %i commands' % (len(commands))))

def submit(scriptname, procnames, submissionscript, arguments, runtime=(1,00,00), ncores=1):
    runtime_str, queue = format_runtime(runtime)
    if not isinstance(scriptname, list):
        scriptname = [scriptname]
    if len(scriptname) == 0:
        print(yellow('  --> Gave a list of length 0, no jobs to be submitted.'))
        return
    logfiles = [os.path.join(os.path.dirname(sn).replace('commands', 'logs'), sn.split('/')[-1].replace('.sh', '.log')) for sn in scriptname]
    print(blue('--> Submitting %i jobs\n\n' % (len(scriptname))))
    submitcommands = [ 'sbatch --parsable -J %s -p %s -t %s --mem-per-cpu 8000 --cpus-per-task %i --ntasks-per-core 1 -o %s -e %s %s %s %s' % ('generate_%s' % (pn), queue, runtime_str, ncores, l, l, submissionscript, arguments, sn) for sn,l,pn in zip(scriptname, logfiles, procnames)]
    for sc, pn in zip(submitcommands, procnames):
        jobid = int(subprocess.check_output(sc.split(' ')))
        print(green('  --> Submitted generation job \'generate_%s\' with ID %i' % (pn, jobid)))
    print(green('\n\n--> Done submitting %i jobs' % (len(scriptname))))

def get_copy_command(from_filename, to_filename):
    return 'cp %s %s' % (from_filename, to_filename)

def count_and_get_matching_keys(dictionary, target_string):
    keys_list = list(dictionary.keys())
    matching_keys = [key for key in keys_list if key in target_string]
    count = len(matching_keys)
    return count, matching_keys

def convert_string_to_float(input_string):
    if 'p' in input_string:
        integer_part, decimal_part = input_string.split('p')
        result = float(integer_part) + float(decimal_part) / 10
    else:
        result = float(input_string)
    return result





    

def get_limits_r():
    # Numbers copied from "../files/UpperLimit_ratio_log.C" made by Fanqiang. To be changed whenever the limits change (if ever).
    y_95 = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498, 0.9900498, 0.9390532, 0.8906833, 0.8448049, 0.8012897, 0.7600159, 0.7208681, 0.6837368, 0.6485181, 0.6151134, 0.5834295, 0.5533775, 0.5248735, 0.4978377, 0.4721944, 0.4478721, 0.4248026, 0.4029213, 0.3821672, 0.3624821, 0.3438109, 0.3261015, 0.3093043, 0.2933723, 0.2782609, 0.2639279, 0.2503332, 0.2374387, 0.2252084, 0.2136081, 0.2026054, 0.1921693, 0.1822708, 0.1728822, 0.1639772, 0.1555309, 0.1475196, 0.139921, 0.1327138, 0.1258778, 0.1193939, 0.113244, 0.1074109, 0.1018783, 0.09663062, 0.09165325, 0.08693227, 0.08245446, 0.0782073, 0.0741789, 0.07035801, 0.06673392, 0.06329651, 0.06003616, 0.05694375, 0.05401062, 0.05122858, 0.04858984, 0.04608701, 0.04371311, 0.04146148, 0.03932584, 0.03730019, 0.03537889, 0.03355655, 0.03182808, 0.03018864, 0.02863365, 0.02715876, 0.02575983, 0.02443296, 0.02317444, 0.02198074, 0.02084853, 0.01977464, 0.01875607, 0.01778996, 0.01687361, 0.01600446, 0.01518009, 0.01439817, 0.01365654, 0.0129531, 0.0122859, 0.01165306, 0.01105282, 0.0104835, 0.009943501, 0.009431319, 0.00894552, 0.008484743, 0.008047701, 0.007633171, 0.007239992, 0.006867066, 0.006513349, 0.006177852, 0.005859636, 0.005557811, 0.005271532, 0.005]
    x_95 = [26.7274, 26.5588, 26.485, 26.377, 26.234, 26.0911, 26.0543, 25.912, 25.746, 25.5951, 25.4538, 25.3822, 25.2182, 25.0778, 24.9844, 24.8449, 24.7057, 24.5667, 24.4884, 24.328, 24.1903, 24.123, 23.955, 23.8187, 23.6827, 23.5472, 23.412, 23.2772, 23.1627, 22.9585, 22.8257, 22.6933, 22.5615, 22.3597, 22.2549, 22.151, 21.951, 21.8631, 21.6644, 21.4665, 21.3808, 21.1844, 21.0603, 20.9059, 20.7838, 20.6314, 20.4395, 20.2612, 20.1149, 19.8949, 19.7476, 19.6016, 19.4889, 19.2736, 19.073, 18.9328, 18.7217, 18.5241, 18.3888, 18.2237, 18.0195, 17.8275, 17.5951, 17.4062, 17.2517, 17.066, 16.8435, 16.6295, 16.4125, 16.2035, 15.9963, 15.715, 15.5148, 15.2836, 15.0891, 14.8273, 14.5694, 14.3157, 14.0668, 13.7896, 13.535, 13.2452, 12.9728, 12.6636, 12.3319, 12.0066, 11.6825, 11.335, 10.9481, 10.5529, 10.1315, 9.6958, 9.1984, 8.6632, 8.1326, 7.4833, 6.7259, 5.9077, 4.773, 3.2971, 1.1258, 0.25, 0.25, 0.5, 0.5, 0.75, 0.875, 1, 1.125, 0.6875, 1.125, 1.2188, 1.3125, 1.4062, 1.4531, 1.8047, 1.9141, 2.0234, 2.0781, 2.1875, 2.2969, 2.3516, 2.4609, 2.7246, 2.8125, 2.9004, 2.9883, 3.0762, 3.1641, 3.3301, 3.3398, 3.5117, 3.7188, 3.8125, 3.875, 3.9688, 4.0312, 4.125, 4.1875, 4.2812, 4.4795, 4.5762, 4.6406, 4.7051, 4.7695, 4.8662, 5.0801, 5.1465, 5.2461, 5.3125, 5.3789, 5.4453, 5.5928, 5.6602, 5.7444, 5.8286, 5.9814, 6.0498, 6.1182, 6.1865, 6.3442, 6.4136, 6.4829, 6.5522, 6.6216, 6.7852, 6.8555, 6.9082, 6.9609, 7.1289, 7.2002, 7.2715, 7.3428, 7.4141, 7.5007, 7.624, 7.6963, 7.7686, 7.8408, 7.9131, 8.0029, 8.0757, 8.2031, 8.2397, 8.313, 8.3862, 8.4595, 8.5327, 8.6836, 8.7578, 8.7949, 8.8691, 8.9434, 9.0176, 9.1143, 9.189, 9.2866, 9.3618, 9.437, 9.4746, 9.5498, 9.6497]
    (y_95_high, y_95_low) = y_95[:len(y_95)/2], y_95[len(y_95)/2:]
    (x_95_high, x_95_low) = x_95[:len(x_95)/2], x_95[len(x_95)/2:]
    y_68 = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498, 0.9900498, 0.9390532, 0.8906833, 0.8448049, 0.8012897, 0.7600159, 0.7208681, 0.6837368, 0.6485181, 0.6151134, 0.5834295, 0.5533775, 0.5248735, 0.4978377, 0.4721944, 0.4478721, 0.4248026, 0.4029213, 0.3821672, 0.3624821, 0.3438109, 0.3261015, 0.3093043, 0.2933723, 0.2782609, 0.2639279, 0.2503332, 0.2374387, 0.2252084, 0.2136081, 0.2026054, 0.1921693, 0.1822708, 0.1728822, 0.1639772, 0.1555309, 0.1475196, 0.139921, 0.1327138, 0.1258778, 0.1193939, 0.113244, 0.1074109, 0.1018783, 0.09663062, 0.09165325, 0.08693227, 0.08245446, 0.0782073, 0.0741789, 0.07035801, 0.06673392, 0.06329651, 0.06003616, 0.05694375, 0.05401062, 0.05122858, 0.04858984, 0.04608701, 0.04371311, 0.04146148, 0.03932584, 0.03730019, 0.03537889, 0.03355655, 0.03182808, 0.03018864, 0.02863365, 0.02715876, 0.02575983, 0.02443296, 0.02317444, 0.02198074, 0.02084853, 0.01977464, 0.01875607, 0.01778996, 0.01687361, 0.01600446, 0.01518009, 0.01439817, 0.01365654, 0.0129531, 0.0122859, 0.01165306, 0.01105282, 0.0104835, 0.009943501, 0.009431319, 0.00894552, 0.008484743, 0.008047701, 0.007633171, 0.007239992, 0.006867066, 0.006513349, 0.006177852, 0.005859636, 0.005557811, 0.005271532, 0.005]
    x_68 = [20.9055, 20.7649, 20.6314, 20.58, 20.4468, 20.3136, 20.2628, 20.1299, 19.9907, 19.947, 19.8145, 19.682, 19.544, 19.412, 19.2743, 19.1427, 19.0113, 18.8799, 18.8371, 18.7011, 18.5703, 18.4396, 18.3931, 18.2629, 18.1329, 18.0031, 17.8734, 17.7439, 17.5667, 17.4379, 17.3094, 17.1811, 17.053, 16.9252, 16.7962, 16.6675, 16.5408, 16.3701, 16.244, 16.1183, 15.9487, 15.8237, 15.6991, 15.5306, 15.4069, 15.2393, 15.1166, 14.9505, 14.8288, 14.6632, 14.4982, 14.3337, 14.2141, 14.0509, 13.8877, 13.7258, 13.5647, 13.4027, 13.2431, 13.04, 12.882, 12.7228, 12.5224, 12.3644, 12.1661, 12.0096, 11.8137, 11.6147, 11.4213, 11.2245, 11.0285, 10.8398, 10.6394, 10.4102, 10.2197, 9.9861, 9.7542, 9.5243, 9.2964, 9.0261, 8.8349, 8.583, 8.3077, 8.0332, 7.7488, 7.4373, 7.1732, 6.8224, 6.5145, 6.1758, 5.8429, 5.424, 5.0124, 4.6664, 4.1865, 3.7194, 3.1762, 2.5964, 1.9095, 1.3174, 0.8132, 0.375, 0.375, 0.75, 0.75, 1.125, 1.0938, 1.25, 1.4062, 1.7188, 1.5938, 1.7266, 2.1328, 2.2852, 2.3613, 2.6748, 2.8369, 2.8364, 3.0801, 3.2422, 3.4043, 3.4854, 3.6475, 3.7861, 4.0078, 4.1331, 4.2583, 4.3835, 4.5088, 4.6594, 4.7593, 4.9136, 4.9971, 5.123, 5.3281, 5.457, 5.543, 5.6719, 5.7578, 5.8867, 6.0575, 6.1882, 6.2754, 6.3625, 6.5197, 6.6519, 6.761, 6.8494, 6.9819, 7.0703, 7.2328, 7.3221, 7.4608, 7.5507, 7.6631, 7.7754, 7.9174, 8.0079, 8.0984, 8.1889, 8.3338, 8.4249, 8.516, 8.607, 8.7396, 8.8466, 8.9383, 9.049, 9.1604, 9.227, 9.3193, 9.4543, 9.547, 9.6397, 9.7602, 9.8407, 9.934, 10.0273, 10.1653, 10.2589, 10.3377, 10.4317, 10.5128, 10.6052, 10.6995, 10.7938, 10.888, 10.9823, 11.0973, 11.1921, 11.2396, 11.3344, 11.4293, 11.5241, 11.6069, 11.702, 11.7854, 11.8808, 11.9763, 12.073, 12.1688, 12.2532]
    (y_68_high, y_68_low) = y_68[:len(y_68)/2], y_68[len(y_68)/2:]
    (x_68_high, x_68_low) = x_68[:len(x_68)/2], x_68[len(x_68)/2:]
    y_exp = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498]
    x_exp = [15.9375, 15.875, 15.75, 15.6875, 15.5625, 15.4375, 15.375, 15.25, 15.1875, 15.0625, 14.9375, 14.8125, 14.75, 14.625, 14.5625, 14.4375, 14.3125, 14.1875, 14.0625, 14, 13.875, 13.75, 13.6875, 13.5625, 13.4375, 13.3125, 13.1875, 13.0625, 13, 12.875, 12.75, 12.625, 12.5, 12.375, 12.2812, 12.1875, 12.0625, 11.9375, 11.8125, 11.6875, 11.5625, 11.4375, 11.3125, 11.1875, 11.0625, 10.9375, 10.8125, 10.6562, 10.5, 10.375, 10.25, 10.125, 10, 9.875, 9.6875, 9.5625, 9.4375, 9.25, 9.125, 9, 8.875, 8.6875, 8.5625, 8.375, 8.25, 8.0625, 7.9375, 7.75, 7.625, 7.4375, 7.25, 7.125, 6.875, 6.75, 6.5625, 6.375, 6.1875, 6, 5.8125, 5.625, 5.375, 5.25, 5, 4.75, 4.625, 4.375, 4.125, 3.875, 3.75, 3.5, 3.25, 3, 2.75, 2.25, 2, 1.75, 1.5, 1, 1, 0.5, 0.5]
    x_obs = [10.9389, 10.8679, 10.7812, 10.6943, 10.6071, 10.5196, 10.4318, 10.3438, 10.2305, 10.1447, 10.078, 9.9677, 9.8954, 9.7986, 9.7113, 9.6238, 9.506, 9.4148, 9.3234, 9.2316, 9.1395, 9.0471, 8.9548, 8.8772, 8.7838, 8.6901, 8.596, 8.5016, 8.4068, 8.3117, 8.1924, 8.0996, 7.9963, 7.9182, 7.8187, 7.7242, 7.6136, 7.5205, 7.4217, 7.3225, 7.2229, 7.1228, 7.0222, 6.9207, 6.804, 6.7021, 6.5997, 6.4969, 6.3935, 6.2897, 6.1831, 6.0776, 5.972, 5.8655, 5.7584, 5.6566, 5.5483, 5.4402, 5.3315, 5.2223, 5.1124, 5.002, 4.8751, 4.7634, 4.6511, 4.5382, 4.4245, 4.3103, 4.2, 4.0757, 3.9665, 3.8498, 3.7237, 3.6044, 3.4904, 3.3634, 3.2417, 3.1192, 2.9959, 2.8749, 2.7499, 2.6239, 2.4971, 2.3693, 2.2406, 2.111, 1.9798, 1.8481, 1.7154, 1.5793, 1.4443, 1.3082, 1.171, 1.0333, 0.8938, 0.7518, 0.6098, 0.4665, 0.3219, 0.1755, 0.0278]
    y_obs = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498]

    return (zip(x_obs, y_obs), zip(x_exp, y_exp), zip(x_68_low, y_68_low), zip(x_68_high, y_68_high), zip(x_95_low, y_95_low), zip(x_95_high, y_95_high))


if __name__ == '__main__':
    main()
