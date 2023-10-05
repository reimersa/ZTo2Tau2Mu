#! /usr/bin/env python

from argparse import ArgumentParser

import ROOT
from tqdm import tqdm
from array import array
from tdrstyle_all import *
import subprocess
import os, math, json
from yaml import safe_load


ROOT.gROOT.SetBatch(1)



description = """Plotting variables from ntuples."""
parser = ArgumentParser(prog="plotter",description=description,epilog="Finished successfully!")
#                                            help="tag that gets appended to folders and filenames" )
parser.add_argument('-r', "--reference",   dest="reference", nargs='+', default=None, action='store', required=True,
                                           help="Name of the SM file(s) for a given process" )
parser.add_argument('-a', "--alternative", dest="alternative", nargs='+', default=None, action='store', required=True,
                                           help="Name of the alternative file(s) to compare to the reference." )
parser.add_argument('-o', "--outfolder",   dest="outfolder", action='store', required=True,
                                           help="Name of the existing folder to store plots in.")
parser.add_argument('-b', "--basefolder",  dest="basefolder", action='store', required=True,
                                           help="Name of the existing folder to store JSON file in, should be the base directory of this script.")
args = parser.parse_args()



def main():
    
    print '--> Starting to plot variables from ntuples.'

    make_plots(infilenames_ref=args.reference, infilenames_alt=args.alternative)

    print '--> Done with plotting variables from ntuples.'


def make_plots(infilenames_ref, infilenames_alt, normalize_to_binwidth=True, nfiles=None, maxevents=None):
    print '--> Plotting ...'

    histholder_ref = HistHolder()    
    histholder_alt = HistHolder()
    histholder_ref.book_default_hists()
    histholder_alt.book_default_hists()

    nsel_ref, ntot_ref = fill_histograms(histholder=histholder_ref, infilenames=infilenames_ref, nfiles=nfiles, maxevents=maxevents)
    nsel_alt, ntot_alt = fill_histograms(histholder=histholder_alt, infilenames=infilenames_alt, nfiles=nfiles, maxevents=maxevents)

    infopart = args.outfolder.split('_')[-1]
    process = args.outfolder.split('_')[1]
    if infopart == 'sm':
        operator = 'sm'
        order = 'sm'
    else:
        if len(infopart.split('-')) != 2: raise ValueError('Unexpected format of \'infopart\' being \'%s\', should have ==2 entries.' % (infopart))
        operator = infopart.split('-')[-1]
        order = infopart.split('-')[0]
    
    eff_alt = float(nsel_alt) / float(ntot_alt)

    # Read dict from JSON
    if os.path.isfile(os.path.join(args.basefolder, 'acceptances.json')):
        with open(os.path.join(args.basefolder, 'acceptances.json'), 'r') as j:
            acceptances = safe_load(j)
    if acceptances == None: acceptances = {}
    
    if operator in acceptances:
        if process in acceptances[operator]:
            acceptances[operator][process][order] = eff_alt
        else:
            acceptances[operator][process] = {order: eff_alt}
    else:
        acceptances[operator] = {process: {order: eff_alt}}

    # write back to json - not thread-safe though... good enough for now anyway
    with open(os.path.join(args.basefolder, 'acceptances.json'), 'w') as j:
        json.dump(obj=acceptances, fp=j, indent=2, sort_keys=True)
    print '  --> Updated acceptance-dictionary.'

    print '  --> Selected %i events from the reference sample and %i from the alternative.' % (nsel_ref, nsel_alt)


    # make plots, one for each histogram in the histfolder
    histnames = histholder_ref.histdict.keys()
    for histname in histnames:

        # get corresponding histograms from all samples
        firsthist = histholder_ref.histdict[histname]
        xmin = firsthist.GetXaxis().GetXmin()
        xmax = firsthist.GetXaxis().GetXmax()
        nameXaxis = firsthist.GetXaxis().GetTitle()
        nameYaxis = firsthist.GetYaxis().GetTitle()
        if normalize_to_binwidth: nameYaxis = 'Events / GeV'
        ymax = 5

        leg = tdrLeg(0.45,0.75,0.90,0.85, textSize=0.040)
        c = tdrDiCanvas(canvName='c', x_min=xmin, x_max=xmax, y_min=5E-4, y_max=ymax, y_min2=0.3, y_max2=1.7, nameXaxis=nameXaxis, nameYaxis=nameYaxis, nameYaxis2='Alt. / Ref.', square=False, iPos=11)

        # upper part (distributions)
        c.cd(1)
        hist_ref = histholder_ref.histdict[histname]
        if normalize_to_binwidth: normalize_content_to_bin_width(histogram=hist_ref)
        if hist_ref.Integral() > 0: hist_ref.Scale(1/hist_ref.Integral())
        tdrDraw(hist_ref, 'E HIST', mcolor=ROOT.kBlack, lcolor=ROOT.kBlack, marker=1, fstyle=0, lstyle=1)
        hist_ref.SetLineWidth(2)
        eff_ref = 1.
        if not ntot_ref == 0:
            eff_ref = float(nsel_ref)/ntot_ref
        eff_ref_err = math.sqrt((1./math.sqrt(nsel_ref))**2 + (1./math.sqrt(ntot_ref))**2) * eff_ref
        leg.AddEntry(hist_ref, 'Reference, #epsilon = %.4f +- %.4f' % (eff_ref, eff_ref_err), 'L')

        hist_alt = histholder_alt.histdict[histname]
        if normalize_to_binwidth: normalize_content_to_bin_width(histogram=hist_alt)
        if hist_alt.Integral() > 0: hist_alt.Scale(1/hist_alt.Integral())
        tdrDraw(hist_alt, 'E HIST', mcolor=ROOT.kRed+1, lcolor=ROOT.kRed+1, marker=1, fstyle=0, lstyle=1)
        hist_alt.SetLineWidth(2)
        eff_alt = 1.
        if not ntot_alt == 0:
            eff_alt = float(nsel_alt)/ntot_alt
        
        eff_alt_err = math.sqrt((1./math.sqrt(nsel_alt))**2 + (1./math.sqrt(ntot_alt))**2) * eff_alt
        leg.AddEntry(hist_alt, 'Alternative, #epsilon = %.4f +- %.4f' % (eff_alt, eff_alt_err), 'L')

        # ratio
        c.cd(2)
        ratio = ROOT.TGraphAsymmErrors(hist_alt, hist_ref, 'pois')
        tdrDraw(ratio, 'E0Z', mcolor=ROOT.kRed+1, lcolor=ROOT.kRed+1, marker=1, fstyle=0, lstyle=1)
        ratio.SetLineWidth(2)

        l_unity = ROOT.TLine(xmin, 1, xmax, 1)
        l_unity.SetLineColor(ROOT.kBlack)
        l_unity.SetLineWidth(2)
        l_unity.SetLineStyle(2)
        l_unity.Draw('SAME')

        c.cd(1)
        leg.Draw('SAME')
        ROOT.gPad.SetLogy(1)

        c.SaveAs(os.path.join(args.outfolder, histname+'.pdf'))
        del c



def fill_histograms(histholder, infilenames, nfiles=None, maxevents=None):

    ievent = 0
    nbaseline = 0
    nselected = 0
    chain = ROOT.TChain('Events')
    nfiles_loaded = 0
    for infilename in infilenames:
        if nfiles is not None and nfiles_loaded >= nfiles: break
        chain.Add(infilename)
        nfiles_loaded += 1

    print '  --> Loaded %i events' % chain.GetEntries()
    for event in chain:
        if ievent%10000 == 0: print '    --> Filling event no. %i' % ievent
        if maxevents is not None and ievent >= maxevents:
            break

        ievent += 1
        # Selection
        keep_event = True

        if event.m_mupmum_1 > 4 and event.m_mupmum_2 > 4 and event.m_mupmum_3 > 4 and event.m_mupmum_4 > 4 and 12. < event.m_mupmum_max < 75. and 40 < event.m_mumumumu < 100 and event.mu1_charge+event.mu2_charge+event.mu3_charge+event.mu4_charge == 0:
            nbaseline += 1
        else:
            keep_event = keep_event and False

        keep_event = keep_event and all((event.mu1_pt > 26, event.mu2_pt > 3.5, event.mu3_pt > 3.5, abs(event.mu1_eta) < 2.4, abs(event.mu2_eta) < 2.4, abs(event.mu3_eta) < 2.4, abs(event.mu4_eta) < 2.4))
        if abs(event.mu4_eta) < 1.2:
            keep_event = keep_event and event.mu4_pt > 3.5
        else:
            keep_event = keep_event and event.mu4_pt > 2.5
        # keep_event = keep_event and event.mu1_charge+event.mu2_charge+event.mu3_charge+event.mu4_charge == 0


        if not keep_event: continue

        histholder.fill('tau1pt', event.tau1_pt, 1)
        histholder.fill('tau2pt', event.tau2_pt, 1)
        histholder.fill('tau1charge', event.tau1_charge, 1)
        histholder.fill('tau2charge', event.tau2_charge, 1)
        histholder.fill('mu1pt', event.mu1_pt, 1)
        histholder.fill('mu2pt', event.mu2_pt, 1)
        histholder.fill('mu3pt', event.mu3_pt, 1)
        histholder.fill('mu4pt', event.mu4_pt, 1)
        histholder.fill('mu1eta', event.mu1_eta, 1)
        histholder.fill('mu2eta', event.mu2_eta, 1)
        histholder.fill('mu3eta', event.mu3_eta, 1)
        histholder.fill('mu4eta', event.mu4_eta, 1)
        histholder.fill('mu1charge', event.mu1_charge, 1)
        histholder.fill('mu2charge', event.mu2_charge, 1)
        histholder.fill('mu3charge', event.mu3_charge, 1)
        histholder.fill('mu4charge', event.mu4_charge, 1)

        histholder.fill('m_tautaumumu', event.m_tautaumumu, 1)
        histholder.fill('m_mumumumu', event.m_mumumumu, 1)
        histholder.fill('m_mumumumu_40to100', event.m_mumumumu, 1)
        histholder.fill('m_mumumumu_sensitive', event.m_mumumumu, 1)
        histholder.fill('m_mupmum_1', event.m_mupmum_1, 1)
        histholder.fill('m_mupmum_2', event.m_mupmum_2, 1)
        histholder.fill('m_mupmum_3', event.m_mupmum_3, 1)
        histholder.fill('m_mupmum_4', event.m_mupmum_4, 1)
        histholder.fill('m_mupmum_max', event.m_mupmum_max, 1)

        histholder.fill('n_tau', event.n_tau, 1)
        histholder.fill('n_mu', event.n_mu, 1)
        nselected += 1

    return nselected, nbaseline


class HistHolder():
    def __init__(self):
        self.histdict = {}

    def book_hist(self, name, *args):
        hist = ROOT.TH1F(name, *args)
        hist.SetDirectory(0)
        self.histdict[name] = hist

    def fill(self, name, *args):
        self.histdict[name].Fill(*args)

    def book_default_hists(self):
        self.book_hist('tau1pt', ';p_{T}^{gen. #tau 1} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('tau2pt', ';p_{T}^{gen. #tau 2} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('tau1charge', ';charge (gen. #tau 1);Events / bin', 3, -1.5, 1.5)
        self.book_hist('tau2charge', ';charge (gen. #tau 2);Events / bin', 3, -1.5, 1.5)
        self.book_hist('mu1pt', ';p_{T}^{gen. #mu 1} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('mu2pt', ';p_{T}^{gen. #mu 2} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('mu3pt', ';p_{T}^{gen. #mu 3} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('mu4pt', ';p_{T}^{gen. #mu 4} [GeV];Events / bin', 20, 0, 100)
        self.book_hist('mu1eta', ';#eta^{gen. #mu 1} ;Events / bin', 100, -5, 5)
        self.book_hist('mu2eta', ';#eta^{gen. #mu 2} ;Events / bin', 100, -5, 5)
        self.book_hist('mu3eta', ';#eta^{gen. #mu 3} ;Events / bin', 100, -5, 5)
        self.book_hist('mu4eta', ';#eta^{gen. #mu 4} ;Events / bin', 100, -5, 5)
        self.book_hist('mu1charge', ';charge (gen. #mu 1) ;Events / bin', 3, -1.5, 1.5)
        self.book_hist('mu2charge', ';charge (gen. #mu 2) ;Events / bin', 3, -1.5, 1.5)
        self.book_hist('mu3charge', ';charge (gen. #mu 3) ;Events / bin', 3, -1.5, 1.5)
        self.book_hist('mu4charge', ';charge (gen. #mu 4) ;Events / bin', 3, -1.5, 1.5)

        self.book_hist('m_tautaumumu', ';m_{#tau#tau#mu#mu} [GeV];Events / bin', 40, 0, 200)
        self.book_hist('m_mumumumu', ';m_{#mu#mu#mu#mu} [GeV];Events / bin', 40, 0, 200)
        m_mumumumu_40to100_bins   = [ 40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100 ]
        self.book_hist('m_mumumumu_40to100', ';m_{#mu#mu#mu#mu} [GeV];Events / bin', len(m_mumumumu_40to100_bins)-1, array('f', m_mumumumu_40to100_bins))
        m_mumumumu_sensitive   = [ 50, 60, 65, 70, 75, 80, 85, 90, 95, 100 ]
        self.book_hist('m_mumumumu_sensitive', ';m_{#mu#mu#mu#mu} (sensitive bins) [GeV];Events / bin', 6, 50, 80)
        self.book_hist('m_mupmum_1', ';m_{#mu^{+}#mu^{-}} pair 1 [GeV];Events / bin', 40, 0, 200)
        self.book_hist('m_mupmum_2', ';m_{#mu^{+}#mu^{-}} pair 2 [GeV];Events / bin', 40, 0, 200)
        self.book_hist('m_mupmum_3', ';m_{#mu^{+}#mu^{-}} pair 3 [GeV];Events / bin', 40, 0, 200)
        self.book_hist('m_mupmum_4', ';m_{#mu^{+}#mu^{-}} pair 4 [GeV];Events / bin', 40, 0, 200)
        self.book_hist('m_mupmum_max', ';m_{#mu^{+}#mu^{-}} maximum [GeV];Events / bin', 40, 0, 200)

        self.book_hist('n_tau', ';N_{#tau};Events / bin', 11, -0.5, 10.5)
        self.book_hist('n_mu', ';N_{#mu};Events / bin', 11, -0.5, 10.5)
        




def normalize_content_to_bin_width(histogram):

    # Normalize the content of each bin to the bin width
    for bin in range(1, histogram.GetNbinsX() + 1):
        content = histogram.GetBinContent(bin)
        error = histogram.GetBinError(bin)
        bin_width = histogram.GetBinWidth(bin)
        normalized_content = content / bin_width
        normalized_error = error / bin_width
        histogram.SetBinContent(bin, normalized_content)
        histogram.SetBinError(bin, normalized_error)

if __name__ == '__main__':
    main()
