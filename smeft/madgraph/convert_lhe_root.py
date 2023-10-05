#! /usr/bin/env python

from argparse import ArgumentParser

import ROOT
from tqdm import tqdm
from array import array
ROOT.gROOT.SetBatch(1)



description = """Converting LHE files to ROOT files for analysis."""
parser = ArgumentParser(prog="converter",description=description,epilog="Finished successfully!")
# parser.add_argument('-l', "--local",       dest="tag", type=str, default='_newffsuman2302_bounduncor_unnorm_fixedttunc_tauidinflatedpog_tauiddmuncor_finalresults', action='store',
parser.add_argument('-i', "--infilename",  dest="infilename", default=None, action='store',
                                           help="Name of the LHE file" )
parser.add_argument('-o', "--outfilename", dest="outfilename", default=None, action='store',
                                           help="Name of the output ROOT file" )
args = parser.parse_args()



def main():
    
    print '--> Starting LHE -> ROOT conversion.'
    print '--> Converting \'%s\' to \'%s\'.' % (args.infilename, args.outfilename)

    file_lhe = open(args.infilename, 'r')
    
    file_root = ROOT.TFile(args.outfilename, 'RECREATE')
    outtree = ROOT.TTree('Events', 'Some variables converted from LHE to ROOT format')
    
    tau1_pt = array('f', [ 0. ])
    tau1_eta = array('f', [ 0. ])
    tau1_phi = array('f', [ 0. ])
    tau1_e = array('f', [ 0. ])
    
    tau2_pt  = array('f', [ 0. ])
    tau2_eta = array('f', [ 0. ])
    tau2_phi = array('f', [ 0. ])
    tau2_e   = array('f', [ 0. ])
    
    mu1_pt   = array('f', [ 0. ])
    mu1_eta  = array('f', [ 0. ])
    mu1_phi  = array('f', [ 0. ])
    mu1_e    = array('f', [ 0. ])
    
    mu2_pt   = array('f', [ 0. ])
    mu2_eta  = array('f', [ 0. ])
    mu2_phi  = array('f', [ 0. ])
    mu2_e    = array('f', [ 0. ])
    
    mu3_pt   = array('f', [ 0. ])
    mu3_eta  = array('f', [ 0. ])
    mu3_phi  = array('f', [ 0. ])
    mu3_e    = array('f', [ 0. ])
    
    mu4_pt   = array('f', [ 0. ])
    mu4_eta  = array('f', [ 0. ])
    mu4_phi  = array('f', [ 0. ])
    mu4_e    = array('f', [ 0. ])
    
    m_tautaumumu = array('f', [ 0. ])
    m_mumumumu   = array('f', [ 0. ])

    outtree.Branch('tau1_pt',  tau1_pt,  'tau1_pt/F')
    outtree.Branch('tau1_eta', tau1_eta, 'tau1_eta/F')
    outtree.Branch('tau1_phi', tau1_phi, 'tau1_phi/F')
    outtree.Branch('tau1_e',   tau1_e,   'tau1_e/F')
    
    outtree.Branch('tau2_pt',  tau2_pt,  'tau2_pt/F')
    outtree.Branch('tau2_eta', tau2_eta, 'tau2_eta/F')
    outtree.Branch('tau2_phi', tau2_phi, 'tau2_phi/F')
    outtree.Branch('tau2_e',   tau2_e,   'tau2_e/F')
    
    outtree.Branch('mu1_pt',  mu1_pt,  'mu1_pt/F')
    outtree.Branch('mu1_eta', mu1_eta, 'mu1_eta/F')
    outtree.Branch('mu1_phi', mu1_phi, 'mu1_phi/F')
    outtree.Branch('mu1_e',   mu1_e,   'mu1_e/F')
    
    outtree.Branch('mu2_pt',  mu2_pt,  'mu2_pt/F')
    outtree.Branch('mu2_eta', mu2_eta, 'mu2_eta/F')
    outtree.Branch('mu2_phi', mu2_phi, 'mu2_phi/F')
    outtree.Branch('mu2_e',   mu2_e,   'mu2_e/F')
    
    outtree.Branch('mu3_pt',  mu3_pt,  'mu3_pt/F')
    outtree.Branch('mu3_eta', mu3_eta, 'mu3_eta/F')
    outtree.Branch('mu3_phi', mu3_phi, 'mu3_phi/F')
    outtree.Branch('mu3_e',   mu3_e,   'mu3_e/F')
    
    outtree.Branch('mu4_pt',  mu4_pt,  'mu4_pt/F')
    outtree.Branch('mu4_eta', mu4_eta, 'mu4_eta/F')
    outtree.Branch('mu4_phi', mu4_phi, 'mu4_phi/F')
    outtree.Branch('mu4_e',   mu4_e,   'mu4_e/F')
    
    outtree.Branch('m_tautaumumu', m_tautaumumu, 'm_tautaumumu/F')
    outtree.Branch('m_mumumumu', m_mumumumu, 'm_mumumumu/F')
    
    lines = file_lhe.readlines()
    
    # Strips the newline character
    tau1 = ROOT.TLorentzVector()
    tau2 = ROOT.TLorentzVector()
    mu1  = ROOT.TLorentzVector()
    mu2  = ROOT.TLorentzVector()
    mu3  = ROOT.TLorentzVector()
    mu4  = ROOT.TLorentzVector()
    
    nevt = 0
    found_event = False
    filled_tau1 = False
    filled_tau2 = False
    filled_mu1  = False
    filled_mu2  = False
    filled_mu3  = False
    filled_mu4  = False
    
    for line in lines:
        word_counter = 0
        at_tau1 = False
        at_tau2 = False
        at_mu1  = False
        at_mu2  = False
        at_mu3  = False
        at_mu4  = False
        px = 0
        py = 0
        pz = 0
        enegry = 0
        particle_id = 0
        for word in line.split():
            if ('<event>' in word):
                found_event = True
            if not found_event:
                continue
        
            if (word == '15' or word == '-15' or at_tau1) and (not filled_tau1):
                at_tau1 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10): 
                    energy = float(word)
                    tau1.SetPxPyPzE(px, py, pz, energy)
                    filled_tau1 = True
                    at_tau1 = False
      
            if (word == '15' or word == '-15' or at_tau2) and not filled_tau2 and filled_tau1:
                at_tau2 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10):
                    energy = float(word)
                    tau2.SetPxPyPzE(px, py, pz, energy)
                    filled_tau2 = True
                    at_tau2 = False
      
            if (word == '13' or word == '-13' or at_mu1) and not filled_mu1:
                at_mu1 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10):
                    energy = float(word)
                    mu1.SetPxPyPzE(px, py, pz, energy)
                    filled_mu1 = True
                    at_mu1 = False
      
            if (word == '13' or word == '-13' or at_mu2) and not filled_mu2 and filled_mu1:
                at_mu2 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10):
                    energy = float(word)
                    mu2.SetPxPyPzE(px, py, pz, energy)
                    filled_mu2 = True
                    at_mu2 = False
      
            if (word == '13' or word == '-13' or at_mu3) and not filled_mu3 and all((filled_mu1, filled_mu2)):
                at_mu3 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10):
                    energy = float(word)
                    mu3.SetPxPyPzE(px, py, pz, energy)
                    filled_mu3 = True
                    at_mu3 = False
      
            if (word == '13' or word == '-13' or at_mu4) and not filled_mu4 and all((filled_mu1, filled_mu2, filled_mu3)):
                at_mu4 = True
                word_counter += 1
                if (word_counter==1): particle_id = word
                if (word_counter==7): px = float(word)
                if (word_counter==8): py = float(word)
                if (word_counter==9): pz = float(word)
                if (word_counter==10):
                    energy = float(word)
                    mu4.SetPxPyPzE(px, py, pz, energy)
                    filled_mu4 = True
                    at_mu4 = False
      
      
      
      
            if (word == '</event>'):
                nevt += 1
                if nevt % 1000 == 0: print '  --> At event %i' % (nevt)
                
                tau1_pt[0]   = tau1.Pt()
                tau1_eta[0]  = tau1.Eta()
                tau1_phi[0]  = tau1.Phi()
                tau1_e[0]    = tau1.E()
                
                tau2_pt[0]   = tau2.Pt()
                tau2_eta[0]  = tau2.Eta()
                tau2_phi[0]  = tau2.Phi()
                tau2_e[0]    = tau2.E()
                
                mu1_pt[0]   = mu1.Pt()
                mu1_eta[0]  = mu1.Eta()
                mu1_phi[0]  = mu1.Phi()
                mu1_e[0]    = mu1.E()
                
                mu2_pt[0]   = mu2.Pt()
                mu2_eta[0]  = mu2.Eta()
                mu2_phi[0]  = mu2.Phi()
                mu2_e[0]    = mu2.E()
                
                mu3_pt[0]   = mu3.Pt()
                mu3_eta[0]  = mu3.Eta()
                mu3_phi[0]  = mu3.Phi()
                mu3_e[0]    = mu3.E()
                
                mu4_pt[0]   = mu4.Pt()
                mu4_eta[0]  = mu4.Eta()
                mu4_phi[0]  = mu4.Phi()
                mu4_e[0]    = mu4.E()
          
                if all((filled_tau1, filled_tau2, filled_mu1, filled_mu2)):
                    m_tautaumumu[0] = (tau1 + tau2 + mu1 + mu2).M()      
                if all((filled_mu1, filled_mu2, filled_mu3, filled_mu4)):
                    m_mumumumu[0] = (mu1 + mu2 + mu3 + mu4).M()
                outtree.Fill()
          
                found_event = False
                filled_tau1 = False
                filled_tau2 = False
                filled_mu1  = False
                filled_mu2  = False
                filled_mu3  = False
                filled_mu4  = False
    
    outtree.Write()
    file_root.Close()
    file_lhe.close()

    print '--> Converted \'%s\' to \'%s\'.' % (args.infilename, args.outfilename)
    print '--> Done with LHE -> ROOT conversion.'







if __name__ == '__main__':
    main()
