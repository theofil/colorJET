#!/usr/bin/env python
#python test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
#ipython -- test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

from array import array

import os
from sys import argv

from math import *

from DataFormats.FWLite import Handle, Events
from PhysicsTools.HeppyCore.utils.deltar import *

from itertools import combinations 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib

import pdb
### define here filename and max number of events to process
#fileName = "/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/VBF_HToInvisible_PU0/inputs104X_VBF_HToInvisible_PU0_job1.root"
fileName = '/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/150720.done/VBF_HToInvisible_PU200/inputs110X_1.root'
#fileName = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/685D035E-FBE3-E811-83B4-48FD8E282493.root"
maxEvents =  7
debug = False



# start analysis
events = Events(fileName)
genj        = Handle("std::vector<reco::GenJet>")
genmet      = Handle("vector<reco::GenMET>")
genp        = Handle("std::vector<reco::GenParticle>") 

# make a simple counter
class counter:
    """counting events"""
    pass
count = counter()
count.maxeve   = maxEvents
count.alleve   = 0


### ntpl vars
treefile = ROOT.TFile("tree.root", 'recreate')
tree = ROOT.TTree("events", "events")

t_nbjets    = array('i', [0])
t_njets     = array('i', [0])
t_b1pt      = array('f', [0])
t_b2pt      = array('f', [0])
t_b1eta     = array('f', [0])
t_b2eta     = array('f', [0])
t_mbb       = array('f', [0])
t_dphibb    = array('f', [0])
t_detabb    = array('f', [0])
t_drbb      = array('f', [0])
t_j1pt      = array('f', [0])
t_j2pt      = array('f', [0])
t_j1eta     = array('f', [0])
t_j2eta     = array('f', [0])

tree.Branch("njets",  t_njets,   "njets/I")
tree.Branch("nbjets", t_nbjets,  "nbjets/I")
tree.Branch("b1pt",   t_b1pt,    "b1pt/F")
tree.Branch("b2pt",   t_b2pt,    "b2pt/F")
tree.Branch("b1eta",  t_b1eta,   "b1eta/F")
tree.Branch("b2eta",  t_b2eta,   "b2eta/F")
tree.Branch("mbb",    t_mbb,     "mbb/F")
tree.Branch("dphibb", t_dphibb,  "dphibb/F")
tree.Branch("detabb", t_detabb,  "detabb/F")
tree.Branch("drbb",   t_drbb,    "drbb/F")
tree.Branch("j1pt",   t_j1pt,    "j1pt/F")
tree.Branch("j2pt",   t_j2pt,    "j2pt/F")
tree.Branch("j1eta",  t_j1eta,   "j1eta/F")
tree.Branch("j2eta",  t_j2eta,   "j2eta/F")


# check if "a" is mother of "p" 
def isAncestor(a,p) :
        if a == p :
                return True
        for i in xrange(0,p.numberOfMothers()) :
                if isAncestor(a,p.mother(i)) :
                         return True
        return False

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

# start the bloody event loop
for iev,event in enumerate(events):
    count.alleve += 1
    if iev > maxEvents: break
    idev = "%d:%d:%d" % ( event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())
    if debug: print "Event %s" % idev

    # get pruned collection
#    event.getByLabel("prunedGenParticles", genpruned)
#    bhadrons = [p for p in genpruned.product() if abs(p.pdgId()) > 500 and abs(p.pdgId()) < 600]

    # get genjets
    event.getByLabel("ak4GenJetsNoNu", genj)
    genjs   = [ g for g in genj.product() if abs(g.eta())<5] ### study all jets within pseudorapidity 3

    # check if the pointer to the daughter exists
    dauPtrExists = lambda j, i: j.daughterPtr(i).isAvailable() and j.daughterPtr(i).isNonnull()              

    # make a set of bjets {}
#    bjets = {j for j in genjs for i in range(j.numberOfDaughters()) for bhad in bhadrons if dauPtrExists(j,i) if isAncestor(bhad, j.daughter(i)) }


    event.getByLabel("genParticles", genp)
    genps            = [p for p in genp.product() if p.status()==1]
    gennus           = [p for p in genp.product() if p.status()==1 and abs(p.pdgId()) in (12,14,15)] 
    genpromptps      = [p for p in genp.product() if p.status()==1 and p.isPromptFinalState()] 
    genpromptnus     = [p for p in genp.product() if p.status()==1 and p.isPromptFinalState() and abs(p.pdgId()) in (12,14,15)] 
   
    genV             = sumP4(genpromptnus)                         # the prompt V
    genpsnopromptnus = [p for p in genps if p not in genpromptnus] # all but prompt nu



    collection = [p for p in genpsnopromptnus if p.pt()>0.5]
    etas     = map(lambda p:p.y()  , collection) ### it's rapidity that is used in reality
    phis     = map(lambda p:p.phi(), collection)
    size     = map(lambda p:p.pt() , collection)
    energies = map(lambda p:p.energy() , collection)

    f, (ax1) = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(11,7))
    
    ### make first subplot
    cmap = matplotlib.cm.get_cmap('viridis')
    #normalize = matplotlib.colors.Normalize(vmin=min(size), vmax=max(size))
    normalize = matplotlib.colors.LogNorm(vmin=2, vmax=7000)
    colors = [cmap(normalize(value)) for value in size]
    #cax, _ = matplotlib.colorbar.make_axes(ax1)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
    ax1.scatter(etas, phis, s=size, marker='o', color = colors)
    plt.xlabel('Rapidity')   
    plt.ylabel('Phi')
    
    theta = np.arange(0, 2*np.pi , 0.01)
    

    myjet = genjs[0]
    ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
    myjet = genjs[1]
    ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
    
    for myjet in genjs[2:]:
        ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'r')

    for myjet in genjs[2:]:
        ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'r')
    

    ax1.plot( genV.y() + 0.3 * np.cos(theta), genV.phi() + 0.3* np.sin(theta), color = 'g')
    
    ax1.set_xlim(-5,5)
    ax1.set_ylim(-3.2, 3.2)
    dijetmass  = sumP4(genjs[0:2]).mass()
     
    bjetLead = 0
    if genjs[0] in genjs[0:2] and  genjs[1] in genjs[0:2] : bjetLead = 1

    f.suptitle('mjj = %2.1f ptjj = %2.1f j_pt = (%2.1f, %2.1f) mV = %2.1f ptV = %2.1f yV = %2.1f  phiV = %2.1f'%(dijetmass, sumP4(genjs[0:2]).Pt(), 
    genjs[0].pt(), genjs[1].pt(), genV.mass(), genV.pt(), genV.Y(), genV.Phi()))
    plt.savefig('/afs/cern.ch/user/t/theofil/www/images/figs/neve%d.png'%iev)
    #plt.savefig('../figs/eve%d.png'%iev)
    plt.close('all')
    #plt.show()

print("events.size    = %d"%events.size())
print("count.alleve   = %d"%count.alleve)

plt.close('all')
