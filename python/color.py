#!/usr/bin/env python
#python test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
#ipython -- test.py -n 1 -g inputs_GluGluToHHTo4B_PU0_job1.root
# http://theofil.web.cern.ch/theofil/images/figs/
# https://cpp.hotexamples.com/examples/-/GenParticle/-/cpp-genparticle-class-examples.html
# https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/PackedGenParticle.h
# https://cmssdt.cern.ch/lxr/source/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
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
#fileName = "inputs_GluGluToHHTo4B_PU0_job1.root"
#fileName = "inputs_GluGluToHHTo4B_PU0_job1.root"
fileName = "../ROOTs/ZZqqnunu_small.root"
#fileName = "../ROOTs/ZZqqnunu.root"
#fileName = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/685D035E-FBE3-E811-83B4-48FD8E282493.root"
maxEvents =  5
debug = False



# start analysis
events = Events(fileName)
genj        = Handle("std::vector<reco::GenJet>")
genpacked   = Handle("std::vector<pat::PackedGenParticle>")
genpruned   = Handle("std::vector<reco::GenParticle>") # genprunedGenParticles

# make a simple counter
class counter:
    """counting events"""
    pass
count = counter()
count.maxeve   = maxEvents
count.alleve   = 0
count.fpeve    = 0
count.twobhad  = 0
count.manybhad = 0
count.goodeve  = 0 
count.njl2     = 0


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
    if iev >= maxEvents: break
    idev = "%d:%d:%d" % ( event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())
    if debug: print "Event %s" % idev

    # get pruned collection
    event.getByLabel("prunedGenParticles", genpruned)
    bhadrons = [p for p in genpruned.product() if abs(p.pdgId()) > 500 and abs(p.pdgId()) < 600]

    # count the number of bhadrons
    if len(bhadrons) <2  : continue
    if len(bhadrons) == 2: count.twobhad +=1
    if len(bhadrons) > 2 : 
        if debug: print ("len(bhadrons) = %d"%len(bhadrons))
        count.manybhad +=1

    # get genjets
    event.getByLabel("slimmedGenJets", genj)
    allGenJ   = [ g for g in genj.product() if abs(g.eta())<3] ### study all jets within pseudorapidity 3

    # count the number of events with less than two jets
    if len(allGenJ) <2: 
        count.njl2 +=0
        continue     

    # check if the pointer to the daughter exists
    dauPtrExists = lambda j, i: j.daughterPtr(i).isAvailable() and j.daughterPtr(i).isNonnull()              

    # make a set of bjets {}
    bjets = {j for j in allGenJ for i in range(j.numberOfDaughters()) for bhad in bhadrons if dauPtrExists(j,i) if isAncestor(bhad, j.daughter(i)) }

    # if there are no b-jets, skip the event
    if len(bjets) < 2: continue

    # alljets but the b-jets
    allGenJnoB = [j for j in allGenJ if j not in bjets]
    
    count.goodeve += 1

    #event.getByLabel("genParticles", genpacked)
    event.getByLabel("packedGenParticles", genpacked)
    allparticles   = [p for p in genpacked.product() if p.status()==1]
    allprompt      = [p for p in genpacked.product() if p.status()==1 and p.isPromptFinalState()] 
    allNu          = [p for p in genpacked.product() if p.status()==1 and p.isPromptFinalState() and abs(p.pdgId()) in (12,14,15)] 
    nonprompt      = [p for p in genpacked.product() if p.status()==1 and not p.isPromptFinalState()]
    allpromptNoNu  = [p for p in allprompt if p not in allNu]
    allNoNu        = [p for p in allparticles if p not in allNu]


#    collection = allNoNu
    collection = [p for p in allNoNu if p.pt()>0.5]
    #collection = allpromptNoNu
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
    
    bjlist = [b for b in bjets]

    myjet = bjlist[0]
    ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
    myjet = bjlist[1]
    ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'b')
    
    for myjet in allGenJnoB[2:]:
        ax1.plot( myjet.y() + 0.4 * np.cos(theta), myjet.phi() + 0.4* np.sin(theta), color = 'r')
    
    
    ax1.set_xlim(-5,5)
    ax1.set_ylim(-3.2, 3.2)
    dijetmass  = sumP4(bjlist[0:2]).mass()
     
    bjetLead = 0
    if bjlist[0] in allGenJ[0:2] and  bjlist[1] in allGenJ[0:2] : bjetLead = 1

    f.suptitle('M = %2.1f Pt = %2.1f bj_pt = (%2.1f, %2.1f) N = %d Nbj = %d bjL = %d maxPt = %2.1f'%(dijetmass, sumP4(bjlist[0:2]).Pt() ,bjlist[0].pt(), bjlist[1].pt(), len(allparticles), len(bjlist), bjetLead, max(size)))
    #plt.savefig('/afs/cern.ch/user/t/theofil/www/images/figs/eve%d.png'%iev)
    #plt.savefig('../figs/eve%d.png'%iev)
    plt.close('all')
    #plt.show()

print("events.size    = %d"%events.size())
print("count.alleve   = %d"%count.alleve)
print("count.fpeve    = %d"%count.fpeve)
print("count.goodeve  = %d"%count.goodeve)
print("count.twobhad  = %d"%count.twobhad)
print("count.manybhad = %d"%count.manybhad)
print("count.njl2     = %d"%count.njl2)

plt.close('all')
