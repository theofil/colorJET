#!/usr/bin/env python
#python test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
#ipython -- test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Handle, Events
#from PhysicsTools.HeppyCore.utils.deltar import *
#from PhysicsTools.HeppyCore.statistics.tree import *

from array import array

import os
from sys import argv

from math import *


from itertools import combinations 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib

import pdb

### define here filename and max number of events to process
#fileName = "/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/VBF_HToInvisible_PU0/inputs104X_VBF_HToInvisible_PU0_job1.root"
#fileName = '/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/150720.done/VBF_HToInvisible_PU200/inputs110X_1.root'
fileName = '/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/150720.done/QCD_Pt15to3000_PU200/inputs110X_1.root'
#fileName = '/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/WToLNu_PU200/inputs104X_WToLNu_PU200_job14.root'
#fileName = '/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/VBFHToBB_PU200/inputs104X_VBFHToBB_PU200_job1.root'
#fileName = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/685D035E-FBE3-E811-83B4-48FD8E282493.root'
#fileName = 'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_1_0/RelValQCD_FlatPt15to3000_13UP18_RD/GEN-SIM/111X_upgrade2018_realistic_RunDep_v1-v1/10000/005600A2-F19A-2B44-92A8-CF227D764F73.root'
#fileName = 'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/VBF_HToInvisible_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-DIGI-RAW/PU140_110X_mcRun4_realistic_v3_ext1-v2/270000/011556F3-DF8D-244F-80EB-B7DAAF20AD63.root'
#fileName = '/afs/cern.ch/work/t/theofil/CMSSW/CMSSW_10_1_7/src/FastPUPPI/NtupleProducer/python/scripts/COlOR/ROOTs/ZZqqnunu_small.root'
#fileName = '/afs/cern.ch/work/t/theofil/CMSSW/CMSSW_10_1_7/src/FastPUPPI/NtupleProducer/python/scripts/COlOR/ROOTs/ZZqqnunu.root'

# configuration
#edmStyle = 'MiniAOD'
edmStyle  = 'AOD'
makePlots = False
folders   = ['ZZ', 'VBFHinv', 'Minbias','QCD', 'VBFHbb']
folder    = folders[3]
maxEvents =  50e6
debug = False



# start analysis
events = Events(fileName)
genmet      = Handle("vector<reco::GenMET>") # vector<reco::GenMET>                  "genMetTrue"                ""                "SIM"
lhe         = Handle("LHEEventProduct") # https://github.com/amarini/ChargedHiggs/blob/nano/test/vbs_pol.py
evtInfo     = Handle("GenEventInfoProduct")     # https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
genj        = Handle("std::vector<reco::GenJet>")

if edmStyle == 'AOD':
    genp = Handle("std::vector<reco::GenParticle>")  

if edmStyle == 'MiniAOD':
    genp    = Handle("std::vector<pat::PackedGenParticle>") 
    genpacked   = Handle("std::vector<pat::PackedGenParticle>")
    genpruned   = Handle("std::vector<reco::GenParticle>") # genprunedGenParticles

# make a simple counter
class counter:
    """counting events"""
    pass
count = counter()
count.maxeve   = maxEvents
count.alleve   = 0





### ntpl vars
filename = 'tree.root' if folder =='' else 'tree_'+folder+'.root'
treefile = ROOT.TFile(filename, 'recreate')

mcWeights = ROOT.TH1I('mcWeights','mcWeights', 4, -2,2);


tree = ROOT.TTree("events", "events")

nJetsMax = 5
tvars = []
t_nJets       = array('i', [0]); tvars += [t_nJets]
t_jetPt       = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPt]
t_jetEta      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetEta]
t_jetPhi      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPhi]
t_jetM        = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetM] 
t_jetPV1      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPV1]
t_jetPV2      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPV2]     
t_jetPVA      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVA]     
t_jetPVM      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVM]     
t_jetFlag     = array('i', [0   for i in range(nJetsMax)]); tvars += [t_jetFlag]    
t_jetMuIndex  = array('i', [-1 for i in range(nJetsMax)]) ; tvars += [t_jetMuIndex] 
t_jetBtag     = array('B', [False  for i in range(nJetsMax)]); tvars += [t_jetBtag]
t_mjj         = array('f', [0.]) ; tvars += [t_mjj]      
t_ptjj        = array('f', [0.]) ; tvars += [t_ptjj]     
t_dYjj        = array('f', [0.]) ; tvars += [t_dYjj]    
t_dPhijj      = array('f', [0.]) ; tvars += [t_dPhijj]   
t_c21         = array('f', [-3.]); tvars += [t_c21]      
t_c12         = array('f', [-3.]); tvars += [t_c12]      
t_met         = array('f', [0.]) ; tvars += [t_met]      
t_metphi      = array('f', [0.]) ; tvars += [t_metphi]   
t_weight      = array('f', [1.]) ; tvars += [t_weight]   
t_LHEweight   = array('f', [1.]) ; tvars += [t_LHEweight]


def reset():
    global tvars
    for var in tvars:
        typecode = var.typecode
        for i in range(len(var)):
            if typecode == 'f': var[i] = -99999.9
            if typecode == 'i': var[i] = int(-99999)
            if typecode == 'B': var[i] = False

tree.Branch("nJets",      t_nJets,      "nJets/I")
tree.Branch("jetPt",      t_jetPt,      "jetPt[nJets]/F")
tree.Branch("jetEta",     t_jetEta,     "jetEta[nJets]/F")
tree.Branch("jetPhi",     t_jetPhi,     "jetPhi[nJets]/F")
tree.Branch("jetM",       t_jetM,       "jetM[nJets]/F")
tree.Branch("jetPV1",     t_jetPV1,     "jetPV1[nJets]/F")
tree.Branch("jetPV2",     t_jetPV2,     "jetPV2[nJets]/F")
tree.Branch("jetPVA",     t_jetPVA,     "jetPVA[nJets]/F")
tree.Branch("jetPVM",     t_jetPVM,     "jetPVM[nJets]/F")
tree.Branch("jetFlag",    t_jetFlag,    "jetFlag[nJets]/I")
tree.Branch("jetBtag",    t_jetBtag,    "jetBtag[nJets]/O")
tree.Branch("jetMuIndex", t_jetMuIndex, "jetMuIndex[nJets]/I")
tree.Branch("mjj",        t_mjj,        "mjj/F")
tree.Branch("c21",        t_c21,        "c21/F")
tree.Branch("c12",        t_c12,        "c12/F")
tree.Branch("ptjj",       t_ptjj,       "ptjj/F")
tree.Branch("dYjj",       t_dYjj,       "dYjj/F")
tree.Branch("dPhijj",     t_dPhijj,     "dPhijj/F")
tree.Branch("met",        t_met,        "met/F")
tree.Branch("metphi",     t_metphi,     "metphi/F")
tree.Branch("weight",     t_weight,     "weight/F")
tree.Branch("LHEweight",  t_LHEweight,   "LHEweight/F")

def fillJets(collection):
    nmax = min(len(collection), nJetsMax)
    t_nJets[0] = nmax 
    for iobj, obj in enumerate(collection[0:nmax]):
        t_jetPt[iobj]      = round(obj.pt() ,  1)
        t_jetEta[iobj]     = round(obj.eta(),  2)
        t_jetPhi[iobj]     = round(obj.phi(),  2)
        t_jetM[iobj]       = round(obj.mass(), 1)
        t_jetPV1[iobj]     = round(obj.pv1,     7)
        t_jetPV2[iobj]     = round(obj.pv2,     7)
        t_jetPVA[iobj]     = round(obj.pva,     5)
        t_jetPVM[iobj]     = round(obj.pvm,     7)
        t_jetFlag[iobj]    = obj.flag
        t_jetBtag[iobj]    = obj.btag
        t_jetMuIndex[iobj] = obj.muIndex

nMuonsMax = 10
t_nMuons        = array('i', [0])
t_muonPt        = array('f', [0 for i in range(nMuonsMax)])
t_muonEta       = array('f', [0 for i in range(nMuonsMax)])
t_muonPhi       = array('f', [0 for i in range(nMuonsMax)])
t_muonM         = array('f', [0 for i in range(nMuonsMax)])
t_muonPtRel     = array('f', [0 for i in range(nMuonsMax)])
t_muonPt03      = array('f', [0 for i in range(nMuonsMax)])
tree.Branch("nMuons",     t_nMuons,   "nMuons/I")
tree.Branch("muonPt",     t_muonPt,     "muonPt[nMuons]/F")
tree.Branch("muonEta",    t_muonEta,    "muonEta[nMuons]/F")
tree.Branch("muonPhi",    t_muonPhi,    "muonPhi[nMuons]/F")
tree.Branch("muonM",      t_muonM,      "muonM[nMuons]/F")
tree.Branch("muonPtRel",  t_muonPtRel,   "muonPtRel[nMuons]/F")
tree.Branch("muonPt03",   t_muonPt03,   "muonPt03[nMuons]/F")
def fillMuons(collection):
    nmax = min(nMuonsMax, len(collection))
    t_nMuons[0] = nmax 
    for iobj, obj in enumerate(collection[0:nmax]):
        t_muonPt[iobj]    = round(obj.pt() ,  1)
        t_muonEta[iobj]   = round(obj.eta(),  2)
        t_muonPhi[iobj]   = round(obj.phi(),  2)
        t_muonM[iobj]     = round(obj.mass(), 2)
        t_muonPtRel[iobj] = obj.PtRel
        t_muonPt03[iobj]  = obj.Pt03
    

# check if "a" is mother of "p" 
def isAncestor(a,p) :
    if a == p :
            return True
    for i in xrange(0,p.numberOfMothers()) :
            if isAncestor(a,p.mother(i)) :
                     return True
    return False

def Pt03(p, particles):
    myP4 = toTLV(p)
    sumPt = 0
    for aParticle in particles:
        aParticleP4 = toTLV(aParticle)
        DR = myP4.DeltaR(aParticleP4)
        if 0 < DR < 0.3:
            sumPt +=  aParticleP4.Pt()

    return sumPt
         
def PtRel(p, jets):
    ptrel = 0
    minDR = 1.e9
    myP4 = toTLV(p)
        
    for aJet in jets:
        aJetP4 = toTLV(aJet)
        DR = myP4.DeltaR(aJetP4)
        if 0 < DR < minDR:
            DR = minDR
            ptrel = myP4.Perp(aJetP4.Vect())
            
    return ptrel            


#def mujets(mus, jets):
    
 ##   ptrel = 0
 ##   minDR = 1.e9
 ##   myP4 = toTLV(p)
 ##       
 ##   for aJet in jets:
 ##       aJetP4 = toTLV(aJet)
 ##       DR = myP4.DeltaR(aJetP4)
 ##       if 0 < DR < minDR:
 ##           DR = minDR
 ##           ptrel = myP4.Perp(aJetP4.Vect())
 ##           
 ##   return ptrel            
    


# translate a LorentzVector to TLorentzVector
toTLV = lambda x: ROOT.TLorentzVector(x.px(), x.py(), x.pz(), x.energy())  

# vector sum of four vectors
def sumP4(fvecs):
    megaJ = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D('double'))(0,0,0,0)
    for jet in fvecs:
        megaJ += jet.p4()
    return megaJ

class miniEvent:
    """mini event class"""
    
    pass


# from https://gitlab.cern.ch/cms-sw/cmssw/blob/e303d9f2c3d4f25397db5feb7ad59d2f20c842f2/PhysicsTools/HeppyCore/python/utils/deltar.py
def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > np.pi:
        res -= 2*np.pi
    while res < -np.pi:
        res += 2*np.pi
    return res



# start the bloody event loop
for iev,event in enumerate(events):
    count.alleve += 1
    reset()
    if iev > maxEvents: break
    idev = "%d:%d:%d" % ( event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())
    if debug: print "Event %s" % idev

    # get pruned collection

    jetPtMin   = 30
    jetEtaMax  = 4.7 

    try:
        event.getByLabel("generator", evtInfo)
        t_weight[0] = evtInfo.product().weights()[0]         
        if debug: print('evtInfo %2.1f'%t_weight[0])
    except:
        t_weight[0] = 1.0         
        if debug: print('eventInfo product does not exist, seeting weight %2.1f'%t_weight[0])

    if t_weight[0] < 0: mcWeights.Fill(-0.999)
    else: mcWeights.Fill(0.999) 

    try:
        res = event.getByLabel("genMetTrue", genmet)
        t_met[0] =  genmet.product()[0].p4().pt()
        t_metphi[0] =  genmet.product()[0].p4().phi()
    except:
        t_met[0] = -1
        t_metphi[0] = -10

  
    try:
        event.getByLabel("externalLHEProducer", lhe)
        t_LHEweight[0] = lhe.product().weights()[0].wgt
        if debug :print('LHE weight %2.1f'%t_LHEweight[0])
    except:
        t_LHEweight[0] = 1.0
        if debug :print('LHE product does not exist setting LHE weight %2.1f'%t_LHEweight[0])

    if edmStyle == 'AOD':
        event.getByLabel("genParticles", genp)           
        event.getByLabel("ak4GenJetsNoNu", genj)         

    if edmStyle == 'MiniAOD':
        event.getByLabel("packedGenParticles", genp)            
        event.getByLabel("slimmedGenJets", genj)   
        event.getByLabel("prunedGenParticles", genpruned)
        bhadrons = [p for p in genpruned.product() if abs(p.pdgId()) > 500 and abs(p.pdgId()) < 600]


    # get genJets
    genjs   = [ g for g in genj.product() if abs(g.eta())<jetEtaMax  and g.pt() > jetPtMin] ### study all jets within pseudorapidity 4.7 and pt > 30

    # check if the pointer to the daughter exists
    dauPtrExists = lambda j, i: j.daughterPtr(i).isAvailable() and j.daughterPtr(i).isNonnull()              

    genps            = [p for p in genp.product() if p.status()==1]
    #gennus           = [p for p in genp.product() if p.status()==1 and abs(p.pdgId()) in (12,14,15)] 
    genmus           = [p for p in genp.product() if p.status()==1 and abs(p.pdgId()) == (13)] 
    # genels           = [p for p in genp.product() if p.status()==1 and abs(p.pdgId()) == (11)] 
    #genpromptps      = [p for p in genp.product() if p.status()==1 and p.isPromptFinalState()] 
    if edmStyle != 'MiniAOD': bhadrons         = [p for p in genp.product() if abs(p.pdgId()) > 500 and abs(p.pdgId()) < 600]

    


    if debug: print('---------')
    for genmu in genmus:
         if debug: print('pt %2.1f   eta %2.1f   phi %2.1f   pt03 %2.1f   ptrel %2.1f'%(genmu.pt(), genmu.eta(), genmu.phi(), Pt03(genmu, genps), PtRel(genmu, genjs)))
         genmu.PtRel = PtRel(genmu, genjs)
         genmu.Pt03  = Pt03(genmu, genps)

    # make a set of bjets {}
    bjets = {j for j in genjs for i in range(j.numberOfDaughters()) for bhad in bhadrons if dauPtrExists(j,i) if isAncestor(bhad, j.daughter(i)) }

    jcs = [] # will store here all jet constitutents
     
    # flag jets
    for j in genjs: # initialize added members 
        j.flag    =  0 
        j.btag    =  False 
        j.muIndex = -1
        j.pv1     = 0.0
        j.pv2     = 0.0
    for j in genjs: 
        if j in bjets:
            j.flag = 1 << 0
            j.btag = True

    for j in genjs: 
        jetconst = [j.daughter(d) for d in range(j.numberOfDaughters()) if dauPtrExists(j, d)]
        for jc in jetconst:
            jcs += [jc]
            if abs(jc.pdgId()) == 13:
                j.flag = 1 << 1
                matchedMuon = jc
                for iobj, obj in enumerate(genmus):
                    if obj == matchedMuon: j.muIndex = iobj

        # calculate pull vector
        pullV = np.array([0., 0.])
        for jc in jetconst:
            if jc.pt() < 1.e-3: continue
            dY    = jc.rapidity() - j.rapidity()
            dPhi  =  deltaPhi(jc.phi(), j.phi())
            r     = (dY**2 + dPhi**2)**0.5
            pullV +=  (jc.pt()/j.pt())*r*np.array([dY, dPhi])
        j.pv1 = pullV[0]
        j.pv2 = pullV[1]
        r = (j.pv1**2 + j.pv2**2)**0.5
        theta  = np.arctan2(j.pv2/r, j.pv1/r)
        j.pvm = r
        j.pva = theta
          

    ### calculate Relative Pull Angle
    def RPA(j1, j2):
        '''cos21 needs j1 and j2 in that order, do not invert'''
        p  = np.array([[j1.pv1], [j1.pv2]])
        v2 = np.array([[j2.rapidity()] , [j2.phi()] ])
        v1 = np.array([[j1.rapidity()] , [j1.phi()] ])
        r = v2 - v1
        mag_p = (p.T).dot(p)[0][0]**0.5
        mag_r = (r.T).dot(r)[0][0]**0.5
        cos21 = (r.T).dot(p)[0][0]/(mag_p*mag_r)
        if debug:
            print('p = ', p)
            print('v2 = ', v2)
            print('v1 = ', v1)
            print('r = ', r)
            print('mag_r = ', mag_r)
            print('mag_rp = ', mag_p)
            print('(r.T).dot(p)[0][0] = ', (r.T).dot(p)[0][0])
        return cos21
    



    ### fill tree
    fillJets  (genjs)
    fillMuons (genmus)
    if len(genjs) >= 2:
        t_mjj[0]  = round(sumP4(genjs[0:2]).mass() , 1)
        t_ptjj[0] = round(sumP4(genjs[0:2]).pt() , 1)
        t_dYjj[0]   = 0. if len(genjs) < 2 else round(genjs[0].rapidity() - genjs[1].rapidity()  , 1)
        t_dPhijj[0] = 0. if len(genjs) < 2 else round(deltaPhi(genjs[0].phi(), genjs[1].phi())   , 1)
        j1 = genjs[0]
        j2 = genjs[1]
        cos21 = RPA(j1, j2)
        cos12 = RPA(j2, j1)
        t_c21[0]  = round(cos21 , 2)
        t_c12[0]  = round(cos12 , 2)
    tree.Fill()

    
    if makePlots:
        plt.rc('font', size=16)
        ### make first subplot
        f, (ax1) = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(11,7))
        cmap = matplotlib.cm.get_cmap('viridis')
    
        collection = [p for p in genps if p.pt()>0.0]
        #collection = [p for p in jcs]
        etas     = map(lambda p:p.y()  , collection) ### it's rapidity that is used in reality
        phis     = map(lambda p:p.phi(), collection)
        pts      = map(lambda p:p.pt() , collection)
        energies = map(lambda p:p.energy() , collection)
        print('jcs lenght = %d'%len(jcs))
        #normalize = matplotlib.colors.Normalize(vmin=min(size), vmax=max(size))
        normalize = matplotlib.colors.LogNorm(vmin=2, vmax=7000)
        colors = [cmap(normalize(value)) for value in pts]
        #cax, _ = matplotlib.colorbar.make_axes(ax1)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
       
        ax1.scatter(etas, phis, s=pts, marker='o', color = colors)
        plt.xlabel('Rapidity')   
        plt.ylabel('Phi')
        
        theta = np.arange(0, 2*np.pi , 0.01)
    
        if len(genjs)>0:
            myjet = genjs[0]
            ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
    
        if len(genjs)>1:
            myjet = genjs[1]
            ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
        
        for myjet in genjs[2:]:
            ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'g')
    
        for myjet in bjets:
            ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'r', linestyle = 'dotted')

        def draw(x1, y1, x2, y2, color = 'y'):
            k    = (y2 - y1)/(x2 - x1)
    
            xs   = np.linspace( min(x1, x2), max(x1, x2), 101)
            ys =  k*(xs - x1) + y1
            ax1.plot( xs , ys , color)
            #ax1.plot( x1 + 0.4 * np.cos(theta), y1 + 0.4* np.sin(theta), color = 'y', linestyle = 'dashed')
    
            r = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
            cosT  = (x2 - x1)/r
            sinT  = (y2 - y1)/r
            txs   = np.linspace(min(x1, x1 + cosT*0.4), max(x1, x1 + cosT*0.4), 101)
            tys = k*(txs - x1) + y1
            ax1.plot( txs , tys , 'r')
            if debug: print('iev %d  yj %2.3f phij %2.3f dY %2.3f dPhi %2.3f r %2.5f'%(iev, x1,  y1, x2, y2, r))
    
        
    
        for myjet in genjs:
            Yj   = myjet.rapidity()
            Phij = myjet.phi()
            Yc   = Yj + myjet.pv1
            Phic = Phij + myjet.pv2
            draw(Yj, Phij, Yc, Phic)
    
    
        # test block
        if False:
            Yj   = 0
            Phij = 0
            Yc   = 0.5*0.4*2**0.5
            Phic = 0.5*0.4*2**0.5
            draw(Yj, Phij, Yc, Phic)
    
            Yj   = -1.5
            Phij = -1.5
            Yc   = -1.5 -0.5*0.4*2**0.5
            Phic = -1.5 -0.5*0.4*2**0.5
            draw(Yj, Phij, Yc, Phic)
    
        if False:
            dY   = 0.5*0.4*2**0.5
            dPhi = 0.5*0.4*2**0.5
            Ys   = np.arange(0 , dY, 0.05)
            Phis = (dPhi/dY)*(Ys - 0) + 0
            ax1.plot( Ys , Phis , color = 'y')
            ax1.plot( 0 + 0.4 * np.cos(theta), 0 + 0.4* np.sin(theta), color = 'y', linestyle = 'dashed')
    
        ax1.set_xlim(-4.7,4.7)
        ax1.set_ylim(-3.2, 3.2)
         
        rad2degrees = 180./np.pi
        if len(genjs)>1: f.suptitle('mjj %2.1f  ptjj %2.1f  dYjj %2.1f dPjj = %2.1f  \n cosTheta21 = %2.3f cosTheta12 = %2.3f  \n j1  = (%2.1f, %2.1f, %2.1f, %2.2f) j2  = (%2.1f, %2.1f, %2.1f, %2.2f)'%(t_mjj[0], t_ptjj[0], 
        t_dYjj[0], t_dPhijj[0],
        t_c21[0], t_c12[0],
        genjs[0].pt(), genjs[0].rapidity(), genjs[0].phi(), t_jetPVA[0]*rad2degrees, 
        genjs[1].pt(), genjs[1].rapidity(), genjs[1].phi(), t_jetPVA[1]*rad2degrees, 
        ), size = 12)
        plt.savefig('/afs/cern.ch/user/t/theofil/www/images/'+folder+'/neve%d.png'%iev)
        plt.close('all')


print("events.size    = %d"%events.size())
print("count.alleve   = %d"%count.alleve)


treefile.cd()
tree.Write()
mcWeights.Write()
treefile.Close()
