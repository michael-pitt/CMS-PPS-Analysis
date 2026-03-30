#!/usr/bin/env python3
import os, sys, math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class AsymmetryModule(Module):
    def __init__(self, channel="mu", year=2026):
        self.channel = channel
        self.year = year
        self.rp_ids = {"45": [3, 23], "56": [103, 123]}

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        # Event / Pileup Summaries
        self.out.branch("nano_nTrkPV05", "I") # nTracks at PV with pT > 0.5
        self.out.branch("nano_nTrkPV09", "I") # nTracks at PV with pT > 0.9     
        
        # Proton Summaries
        self.out.branch("nano_nProtons", "I")
        self.out.branch("nano_pps_arm", "I", lenVar="nProtons")
        self.out.branch("nano_pps_rpid", "I", lenVar="nProtons")
        self.out.branch("nano_pps_x", "F", lenVar="nProtons")
        self.out.branch("nano_pps_y", "F", lenVar="nProtons")
        
        # Jet Summaries
        self.out.branch("nano_nJets", "I")
        self.out.branch("nano_jet_pt", "F", lenVar="nJets")
        self.out.branch("nano_jet_eta", "F", lenVar="nJets")
        self.out.branch("nano_jet_phi", "F", lenVar="nJets")
        self.out.branch("nano_mJets", "F")
        self.out.branch("nano_yJets", "F")
        
        # Lepton Summaries
        self.out.branch("nano_nLeptons", "I")
        self.out.branch("nano_lep_pt", "F", lenVar="nLeptons")
        self.out.branch("nano_lep_eta", "F", lenVar="nLeptons")
        self.out.branch("nano_lep_phi", "F", lenVar="nLeptons")
        self.out.branch("nano_mll", "F")
        self.out.branch("nano_yll", "F")
        self.out.branch("nano_ptll", "F")

        self.out.branch("nano_w_mT", "F")
        self.out.branch("nano_w_pt", "F")
        self.out.branch("nano_w_phi", "F")

    def analyze(self, event):
           
        # 1. Load Collections
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        jets = Collection(event, "Jet")
        protons = Collection(event, "PPSLocalTrack") # Access tracks for x, y for now
        
        # 2. Basic Selections
        
        sel_mu = [m for m in muons if m.pt > 15 and abs(m.eta) < 2.4 and m.looseId]
        
        sel_el = [e for e in electrons if e.pt > 15 and abs(e.eta) < 2.5]

        sel_jets = []
        for j in jets:
            # Standard Tight ID criteria for 13TeV/Run3
            pass_id = (j.neHEF < 0.9 and j.neEmEF < 0.9 and j.nConstituents > 1)
            if abs(j.eta) < 2.4:
                pass_id = pass_id and (j.chHEF > 0 and j.chMultiplicity > 0)
            
            if j.pt > 25 and abs(j.eta) < 4.7 and pass_id:
                sel_jets.append(j)
        jet_sum = ROOT.TLorentzVector()
        for j in sel_jets: jet_sum += j.p4()

        # 3. Event filtering
        if self.channel == "mu":
            if len(sel_mu) == 0: return False
            leptons = sel_mu
        elif self.channel == "el":
            if len(sel_el) == 0: return False
            leptons = sel_el
        elif self.channel == "mj":
            if len(sel_jets) < 2: return False
            leptons = sel_mu
        elif self.channel == "zb":
            leptons = sel_mu + sel_el # Save everything
        else:
            return False

        # 4. Calculations: Dileptons & W variables
        mll = yll = ptll = -1.0
        w_mT = w_pt = w_phi = -1.0
        
        if len(leptons) >= 2:
            dilep = leptons[0].p4() + leptons[1].p4()
            mll, yll, ptll = dilep.M(), dilep.Rapidity(), dilep.Pt()
            
        if len(leptons) >= 1:
            # Transverse Mass calculation using MET
            met_pt = event.PuppiMET_pt
            met_phi = event.PuppiMET_phi
            dphi = ROOT.TVector2.Phi_mpi_pi(leptons[0].phi - met_phi)
            w_mT = math.sqrt(2 * leptons[0].pt * met_pt * (1 - math.cos(dphi)))
            
            # W Vector
            w_vec = ROOT.TVector2(leptons[0].pt * math.cos(leptons[0].phi) + met_pt * math.cos(met_phi),
                                  leptons[0].pt * math.sin(leptons[0].phi) + met_pt * math.sin(met_phi))
            w_pt, w_phi = w_vec.Mod(), w_vec.Phi()

        # 5. Protons Logic (LocalTrack mapping)
        pps_arm = []
        pps_rpid = []
        pps_x = []
        pps_y = []
        for p in protons:
            # Map decRPId to arm (0 for 45, 1 for 56)
            arm = 0 if p.decRPId < 100 else 1
            pps_arm.append(arm)
            pps_rpid.append(p.decRPId)
            pps_x.append(p.x)
            pps_y.append(p.y)

        # 6. Fill Branches
        self.out.fillBranch("nano_nTrkPV05", getattr(event, "PV_ntrk0p5", 0))
        self.out.fillBranch("nano_nTrkPV09", getattr(event, "PV_ntrk0p9", 0))

        self.out.fillBranch("nano_nProtons", len(pps_arm))
        self.out.fillBranch("nano_pps_arm", pps_arm)
        self.out.fillBranch("nano_pps_rpid", pps_rpid)
        self.out.fillBranch("nano_pps_x", pps_x)
        self.out.fillBranch("nano_pps_y", pps_y)

        self.out.fillBranch("nano_nJets", len(sel_jets))
        self.out.fillBranch("nano_jet_pt", [j.pt for j in sel_jets])
        self.out.fillBranch("nano_jet_eta", [j.eta for j in sel_jets])
        self.out.fillBranch("nano_jet_phi", [j.phi for j in sel_jets])
        self.out.fillBranch("nano_mJets", jet_sum.M())
        self.out.fillBranch("nano_yJets", jet_sum.Rapidity())

        self.out.fillBranch("nano_nLeptons", len(leptons))
        self.out.fillBranch("nano_lep_pt", [l.pt for l in leptons])
        self.out.fillBranch("nano_lep_eta", [l.eta for l in leptons])
        self.out.fillBranch("nano_lep_phi", [l.phi for l in leptons])
        self.out.fillBranch("nano_mll", mll)
        self.out.fillBranch("nano_yll", yll)
        self.out.fillBranch("nano_ptll", ptll)
        self.out.fillBranch("nano_w_mT", w_mT)
        self.out.fillBranch("nano_w_pt", w_pt)
        self.out.fillBranch("nano_w_phi", w_phi)

        return True

asymmetry_mu  = lambda : AsymmetryModule(channel="mu")
asymmetry_el  = lambda : AsymmetryModule(channel="el")
asymmetry_mj  = lambda : AsymmetryModule(channel="mj")
asymmetry_zb  = lambda : AsymmetryModule(channel="zb")

