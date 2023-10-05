#! /usr/bin/env python

from argparse import ArgumentParser

import ROOT
from DataFormats.FWLite import Events, Handle

from tqdm import tqdm
from array import array
ROOT.gROOT.SetBatch(1)



description = """Converting GENSIM files to flat ROOT files for analysis."""
parser = ArgumentParser(prog="converter",description=description,epilog="Finished successfully!")
parser.add_argument('-i', "--infilenames", dest="infilenames", default=None, action='store', nargs='+',
                                           help="Name of the GENSIM file(s)" )
parser.add_argument('-o', "--outfilename", dest="outfilename", default=None, action='store',
                                           help="Name of the output ROOT file" )
args = parser.parse_args()



def main():
    
    print '--> Starting GENSIM -> ROOT conversion.'
    
    # events = Events(args.infilenames)
    existing_files = []
    # print args.infilenames
    for idx in range(len(args.infilenames)):
        try:
            f = ROOT.TFile.Open(args.infilenames[idx], 'READ')
            f.Close()
            existing_files.append(args.infilenames[idx])
        except:
            print '  --> File %s does not exist.' % (args.infilenames[idx])
    events = Events(existing_files)
    print '  --> Loaded %i files.' % (len(existing_files))

    file_root = ROOT.TFile(args.outfilename, 'RECREATE')
    outtree = ROOT.TTree('Events', 'Some variables converted from GENSIM to flat ROOT format')
    
    tau1_pt = array('f',  [ 0. ])
    tau1_eta = array('f', [ 0. ])
    tau1_phi = array('f', [ 0. ])
    tau1_e = array('f',   [ 0. ])
    tau1_charge = array('f',   [ 0. ])
    
    tau2_pt  = array('f', [ 0. ])
    tau2_eta = array('f', [ 0. ])
    tau2_phi = array('f', [ 0. ])
    tau2_e   = array('f', [ 0. ])
    tau2_charge   = array('f', [ 0. ])
    
    mu1_pt   = array('f', [ 0. ])
    mu1_eta  = array('f', [ 0. ])
    mu1_phi  = array('f', [ 0. ])
    mu1_e    = array('f', [ 0. ])
    mu1_charge    = array('f', [ 0. ])
    
    mu2_pt   = array('f', [ 0. ])
    mu2_eta  = array('f', [ 0. ])
    mu2_phi  = array('f', [ 0. ])
    mu2_e    = array('f', [ 0. ])
    mu2_charge    = array('f', [ 0. ])
    
    mu3_pt   = array('f', [ 0. ])
    mu3_eta  = array('f', [ 0. ])
    mu3_phi  = array('f', [ 0. ])
    mu3_e    = array('f', [ 0. ])
    mu3_charge    = array('f', [ 0. ])
    
    mu4_pt   = array('f', [ 0. ])
    mu4_eta  = array('f', [ 0. ])
    mu4_phi  = array('f', [ 0. ])
    mu4_e    = array('f', [ 0. ])
    mu4_charge    = array('f', [ 0. ])
    
    m_tautaumumu = array('f', [ 0. ])
    m_mumumumu   = array('f', [ 0. ])
    m_mupmum_1   = array('f', [ 0. ])
    m_mupmum_2   = array('f', [ 0. ])
    m_mupmum_3   = array('f', [ 0. ])
    m_mupmum_4   = array('f', [ 0. ])
    m_mupmum_max = array('f', [ 0. ])

    n_tau        = array('f', [ 0. ])
    n_mu         = array('f', [ 0. ])

    outtree.Branch('tau1_pt',  tau1_pt,  'tau1_pt/F')
    outtree.Branch('tau1_eta', tau1_eta, 'tau1_eta/F')
    outtree.Branch('tau1_phi', tau1_phi, 'tau1_phi/F')
    outtree.Branch('tau1_e',   tau1_e,   'tau1_e/F')
    outtree.Branch('tau1_charge',   tau1_charge,   'tau1_charge/F')
    
    outtree.Branch('tau2_pt',  tau2_pt,  'tau2_pt/F')
    outtree.Branch('tau2_eta', tau2_eta, 'tau2_eta/F')
    outtree.Branch('tau2_phi', tau2_phi, 'tau2_phi/F')
    outtree.Branch('tau2_e',   tau2_e,   'tau2_e/F')
    outtree.Branch('tau2_charge',   tau2_charge,   'tau2_charge/F')
    
    outtree.Branch('mu1_pt',  mu1_pt,  'mu1_pt/F')
    outtree.Branch('mu1_eta', mu1_eta, 'mu1_eta/F')
    outtree.Branch('mu1_phi', mu1_phi, 'mu1_phi/F')
    outtree.Branch('mu1_e',   mu1_e,   'mu1_e/F')
    outtree.Branch('mu1_charge',   mu1_charge,   'mu1_charge/F')
    
    outtree.Branch('mu2_pt',  mu2_pt,  'mu2_pt/F')
    outtree.Branch('mu2_eta', mu2_eta, 'mu2_eta/F')
    outtree.Branch('mu2_phi', mu2_phi, 'mu2_phi/F')
    outtree.Branch('mu2_e',   mu2_e,   'mu2_e/F')
    outtree.Branch('mu2_charge',   mu2_charge,   'mu2_charge/F')
    
    outtree.Branch('mu3_pt',  mu3_pt,  'mu3_pt/F')
    outtree.Branch('mu3_eta', mu3_eta, 'mu3_eta/F')
    outtree.Branch('mu3_phi', mu3_phi, 'mu3_phi/F')
    outtree.Branch('mu3_e',   mu3_e,   'mu3_e/F')
    outtree.Branch('mu3_charge',   mu3_charge,   'mu3_charge/F')
    
    outtree.Branch('mu4_pt',  mu4_pt,  'mu4_pt/F')
    outtree.Branch('mu4_eta', mu4_eta, 'mu4_eta/F')
    outtree.Branch('mu4_phi', mu4_phi, 'mu4_phi/F')
    outtree.Branch('mu4_e',   mu4_e,   'mu4_e/F')
    outtree.Branch('mu4_charge',   mu4_charge,   'mu4_charge/F')
    
    outtree.Branch('m_tautaumumu', m_tautaumumu, 'm_tautaumumu/F')
    outtree.Branch('m_mumumumu', m_mumumumu, 'm_mumumumu/F')
    outtree.Branch('m_mupmum_1', m_mupmum_1, 'm_mupmum_1/F')
    outtree.Branch('m_mupmum_2', m_mupmum_2, 'm_mupmum_2/F')
    outtree.Branch('m_mupmum_3', m_mupmum_3, 'm_mupmum_3/F')
    outtree.Branch('m_mupmum_4', m_mupmum_4, 'm_mupmum_4/F')
    outtree.Branch('m_mupmum_max', m_mupmum_max, 'm_mupmum_max/F')

    outtree.Branch('n_tau', n_tau, 'n_tau/F')
    outtree.Branch('n_mu',  n_mu,  'n_mu/F' )

    handle_gps, label_gps = Handle('std::vector<reco::GenParticle>'), 'genParticles'

    ie = 0
    for e in events:
        # print '\n\n\n===== NEW EVENT ====='
        if ie%500 == 0: print 'New event no. %i' % (ie)
        e.getByLabel(label_gps,handle_gps)
        gps = handle_gps.product()

        mu_final = [p for p in gps if isFinal(p) and abs(p.pdgId()) == 13 and p.status() == 1]
        mu_final.sort(key=lambda p: p.pt(), reverse=True)
        tau_final = [p for p in gps if isFinal(p) and abs(p.pdgId()) == 15]
        tau_final.sort(key=lambda p: p.pt(), reverse=True)
        gps_hard = [p for p in gps if p.isHardProcess()]
        # gps_hard_lep = [p for p in gps_hard if abs(p.pdgId()) in [13, 15]]
        gps_hard_mu = [p for p in gps_hard if abs(p.pdgId()) == 13]
        gps_hard_tau = [p for p in gps_hard if abs(p.pdgId()) == 15]
        # if len(gps_hard_mu) != 2 or len(gps_hard_tau) != 2:
        #     raise ValueError('Not exactly two muons and two taus in hard process')
        if gps_hard[0].numberOfDaughters() != 1:
            continue
        # if len(zbosons) != 1:
        #     for i in range(len(gps_hard)):
        #         printParticle(gps_hard[i])
        #     for i in range(len(mu_final)):
        #         print 'muon number %i, mother pdgId: %i, pt: %f, status: %i' % (i, mu_final[i].mother(0).pdgId(), mu_final[i].pt(), mu_final[i].status())
        #     raise ValueError('Not exactly one Z boson in hard process')

        # for i in range(len(gps_hard)):
        #     printParticle(gps_hard[i])
        # for i in range(len(mu_final)):
        #     print 'muon number %i, mother pdgId: %i, pt: %f, status: %i' % (i, mu_final[i].mother(0).pdgId(), mu_final[i].pt(), mu_final[i].status())
        
        lepv4 = ROOT.TLorentzVector()
        for l in gps_hard_mu + gps_hard_tau:
            thisv4 = ROOT.TLorentzVector()
            thisv4.SetPxPyPzE(l.p4().Px(), l.p4().Py(), l.p4().Pz(), l.p4().E())
            lepv4 += thisv4
        minv_hard = lepv4.M()
        # print 'invariant mass of 4 leptons: %f' % (minv_hard)

        # printDecayChain(gps_hard)

        ie += 1
    
        tau1 = ROOT.TLorentzVector()
        tau2 = ROOT.TLorentzVector()
        mu1  = ROOT.TLorentzVector()
        mu2  = ROOT.TLorentzVector()
        mu3  = ROOT.TLorentzVector()
        mu4  = ROOT.TLorentzVector()

        mus_pos = []
        mus_neg = []

        if len(tau_final) > 0: tau1.SetPxPyPzE(tau_final[0].p4().Px(), tau_final[0].p4().Py(), tau_final[0].p4().Pz(), tau_final[0].p4().E())
        if len(tau_final) > 1: tau2.SetPxPyPzE(tau_final[1].p4().Px(), tau_final[1].p4().Py(), tau_final[1].p4().Pz(), tau_final[1].p4().E())
        
        mu1.SetPxPyPzE(mu_final[0].p4().Px(), mu_final[0].p4().Py(), mu_final[0].p4().Pz(), mu_final[0].p4().E())
        mu2.SetPxPyPzE(mu_final[1].p4().Px(), mu_final[1].p4().Py(), mu_final[1].p4().Pz(), mu_final[1].p4().E())
        mu3.SetPxPyPzE(mu_final[2].p4().Px(), mu_final[2].p4().Py(), mu_final[2].p4().Pz(), mu_final[2].p4().E())
        mu4.SetPxPyPzE(mu_final[3].p4().Px(), mu_final[3].p4().Py(), mu_final[3].p4().Pz(), mu_final[3].p4().E())

        if len(tau_final) > 0: 
          tau1_pt[0]   = tau1.Pt()
          tau1_eta[0]  = tau1.Eta()
          tau1_phi[0]  = tau1.Phi()
          tau1_e[0]    = tau1.E()
          tau1_charge[0]= -1. if tau_final[0].pdgId() > 0 else +1.

        if len(tau_final) > 1: 
          tau2_pt[0]   = tau2.Pt()
          tau2_eta[0]  = tau2.Eta()
          tau2_phi[0]  = tau2.Phi()
          tau2_e[0]    = tau2.E()
          tau2_charge[0]= -1. if tau_final[1].pdgId() > 0 else +1.
        
        mu1_pt[0]   = mu1.Pt()
        mu1_eta[0]  = mu1.Eta()
        mu1_phi[0]  = mu1.Phi()
        mu1_e[0]    = mu1.E()
        mu1_charge[0]= -1. if mu_final[0].pdgId() > 0 else +1.
        if mu1_charge[0] < 0: 
            mus_neg.append(mu1)
        else:                 
            mus_pos.append(mu1)
        
        mu2_pt[0]   = mu2.Pt()
        mu2_eta[0]  = mu2.Eta()
        mu2_phi[0]  = mu2.Phi()
        mu2_e[0]    = mu2.E()
        mu2_charge[0]= -1. if mu_final[1].pdgId() > 0 else +1.
        if mu2_charge[0] < 0: 
            mus_neg.append(mu2)
        else:                 
            mus_pos.append(mu2)
        
        mu3_pt[0]   = mu3.Pt()
        mu3_eta[0]  = mu3.Eta()
        mu3_phi[0]  = mu3.Phi()
        mu3_e[0]    = mu3.E()
        mu3_charge[0]= -1. if mu_final[2].pdgId() > 0 else +1.
        if mu3_charge[0] < 0: 
            mus_neg.append(mu3)
        else:                 
            mus_pos.append(mu3)
        
        mu4_pt[0]   = mu4.Pt()
        mu4_eta[0]  = mu4.Eta()
        mu4_phi[0]  = mu4.Phi()
        mu4_e[0]    = mu4.E()
        mu4_charge[0]= -1. if mu_final[3].pdgId() > 0 else +1.
        if mu4_charge[0] < 0: 
            mus_neg.append(mu4)
        else:                 
            mus_pos.append(mu4)

        m_mupmums = []
        for ip in range(len(mus_pos)):
            for im in range(len(mus_neg)):
                m_mupmums.append((mus_pos[ip] + mus_neg[im]).M())
                # m_mupmum_max = max(m_mupmum_max, m_mupmum)
        if len(m_mupmums) > 0:
            m_mupmum_max[0] = max(m_mupmums)
        else: 
            m_mupmum_max[0] = -1.
        if len(m_mupmums) > 0: m_mupmum_1[0] = m_mupmums[0]
        else: m_mupmum_1[0] = -1.
        if len(m_mupmums) > 1: m_mupmum_2[0] = m_mupmums[1]
        else: m_mupmum_2[0] = -1.
        if len(m_mupmums) > 2: m_mupmum_3[0] = m_mupmums[2]
        else: m_mupmum_3[0] = -1.
        if len(m_mupmums) > 3: m_mupmum_4[0] = m_mupmums[3]
        else: m_mupmum_4[0] = -1.
    
        m_tautaumumu[0] = minv_hard      
        m_mumumumu[0]   = (mu1 + mu2 + mu3 + mu4).M()
        n_tau[0] = len(tau_final)
        n_mu[0]  = len(mu_final)
        outtree.Fill()
    

    file_root.cd()
    outtree.Write()
    file_root.Close()

    print '--> Output written to: %s' % (args.outfilename)
    print '--> Done with GENSIM -> ROOT conversion.'







def isFinal(p):
  # check if one daughter is final and has same PID, then it's not final
  return not (p.numberOfDaughters()==1 and p.daughter(0).pdgId()==p.pdgId())

def printParticle(p):
  string = "Particle with pdgId %9d: status=%2d, pt=%7.2f, eta=%5.2f, phi=%5.2f, final=%5s"%(p.pdgId(),p.status(),p.pt(),p.eta(),p.phi(),isFinal(p))
  if p.numberOfMothers()>=2:
    string += ", mothers %s, %s"%(p.mother(0).pdgId(),p.mother(1).pdgId())
  elif p.numberOfMothers()==1:
    string += ", mother %s"%(p.mother(0).pdgId())
  if p.numberOfDaughters()>=2:
    string += ", daughters %s, %s"%(p.daughter(0).pdgId(),p.daughter(1).pdgId())
  elif p.numberOfDaughters()==1:
    string += ", daughter %s"%(p.daughter(0).pdgId())
  print string

def printDecayChain(mothers):
  for m in mothers:
    daus   = []
    gdaus  = []
    ggdaus = []
    gggdaus = []
    print '='*10 + ' New mother:'
    printParticle(m)
    for d_idx in range(m.numberOfDaughters()):
      d = m.daughter(d_idx)
      daus.append(d)
      # printParticle(d)
    print '-'*10 + ' Daughters:'
    for d in daus:
      printParticle(d)
      for g_idx in range(d.numberOfDaughters()):
        g = d.daughter(g_idx)
        gdaus.append(g)
    print '-'*10 + ' Grand daughters:'
    for g in gdaus:
      printParticle(g)
      for gg_idx in range(g.numberOfDaughters()):
        gg = g.daughter(gg_idx)
        ggdaus.append(gg)
    print '-'*10 + ' Great grand daughters:'
    for gg in ggdaus:
      printParticle(gg)

  for g_idx in range(d.numberOfDaughters()):
    g = d.daughter(g_idx)
    printParticle(g)
    for gg_idx in range(d.numberOfDaughters()):
      gg = d.daughter(gg_idx)
      printParticle(gg)

def getFinalDaughter(p):
    thisp = p
    while not isFinal(thisp):
        thisp = thisp.daughter(0)
    return thisp


if __name__ == '__main__':
    main()
