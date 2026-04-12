#!/usr/bin/env python3
import os, sys, math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# --- Define here the helpers ---
def deltaR(obj1, obj2):
  deta = obj1.eta - obj2.eta
  dphi = ROOT.TVector2.Phi_mpi_pi(obj1.phi - obj2.phi)
  return math.hypot(deta, dphi)

class JetCEPModule(Module):
    def __init__(self, channel="mj", year=2026):
        self.channel = channel
        self.year = year
        self.rp_ids = {"45": [3, 23], "56": [103, 123]}
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        # Event / Pileup Summaries
        self.out.branch("nano_nPV", "I")
        self.out.branch("nano_PV0_ntrk05", "I") # nTracks at PV with pT > 0.5
        self.out.branch("nano_PV0_ntrk09", "I") # nTracks at PV with pT > 0.9
        self.out.branch("nano_PV0_z", "F")
        self.out.branch("nano_NMPI05", "I")
        self.out.branch("nano_NMPI09", "I")
        
        # Proton Summaries
        self.out.branch("nano_nProtons", "I")
        self.out.branch("nano_pps_arm", "I", lenVar="nProtons")
        self.out.branch("nano_pps_rpid", "I", lenVar="nProtons")
        self.out.branch("nano_pps_x", "F", lenVar="nProtons")
        self.out.branch("nano_pps_y", "F", lenVar="nProtons")
        self.out.branch("nano_pps_t", "F", lenVar="nProtons")
        
        # Jet Summaries
        self.out.branch("nano_nJets", "I")
        self.out.branch("nano_jet_pt", "F", lenVar="nJets")
        self.out.branch("nano_jet_eta", "F", lenVar="nJets")
        self.out.branch("nano_jet_phi", "F", lenVar="nJets")
        self.out.branch("nano_Jet_ntrk05", "I", lenVar="nJets")
        self.out.branch("nano_Jet_ntrk09", "I", lenVar="nJets")
        self.out.branch("nano_mJets", "F")
        self.out.branch("nano_yJets", "F")
        self.out.branch("nano_ptJets", "F")
        self.out.branch("nano_dphiJets", "F")
        
        
    def analyze(self, event):
       
        # Load Collections
        jets = Collection(event, "Jet")
        protons = Collection(event, "PPSLocalTrack")
        
        # Object selection
        sel_jets = [j for j in jets if j.pt > 25 and abs(j.eta) < 4.7]
        sel_jets = sorted(sel_jets,key=lambda j: j.pt,reverse=True)
                      
        # event selection
        if len(sel_jets) < 2: return False

        jet_sum = ROOT.TLorentzVector()
        for j in sel_jets: jet_sum += j.p4()
        
        # Protons Logic (LocalTrack mapping)
        pps_arm = []
        pps_rpid = []
        pps_x = []
        pps_y = []
        pps_t = []
        for p in protons:
            # Map decRPId to arm (0 for 45, 1 for 56)
            arm = 0 if p.decRPId < 100 else 1
            pps_arm.append(arm)
            pps_rpid.append(p.decRPId)
            pps_x.append(p.x)
            pps_y.append(p.y)
            pps_t.append(p.time)

        
       
        # Read PV0 tracks
        pv0_ntrk05 = event.PV_ntrk0p5 if hasattr(event, "PV_ntrk0p5") else 0
        pv0_ntrk09 = event.PV_ntrk0p9 if hasattr(event, "PV_ntrk0p9") else 0
        
        # Extract object footprints (Default to 0 if branch missing)
        jet_trk05 = [getattr(j, "ntrk0p5", 0) for j in sel_jets]
        jet_trk09 = [getattr(j, "ntrk0p9", 0) for j in sel_jets]

        
        # N_trk^MPI = N_trk^PV - Sum(N_trk^objects)
        nmpi05 = pv0_ntrk05 - sum(jet_trk05)
        nmpi09 = pv0_ntrk09 - sum(jet_trk09)
        
        # Safety catch: prevent negative MPI in case of overlapping object cones
        if nmpi05 < 0: nmpi05 = 0
        if nmpi09 < 0: nmpi09 = 0
        
        # 7. Fill Branches
        self.out.fillBranch("nano_nPV", getattr(event, "PV_npvsGood", 0))
        self.out.fillBranch("nano_PV0_ntrk05", pv0_ntrk05)
        self.out.fillBranch("nano_PV0_ntrk09", pv0_ntrk09)
        self.out.fillBranch("nano_PV0_z", getattr(event, "PV_z", 0))
        self.out.fillBranch("nano_NMPI05", nmpi05)
        self.out.fillBranch("nano_NMPI09", nmpi09)

        self.out.fillBranch("nano_nProtons", len(pps_arm))
        self.out.fillBranch("nano_pps_arm", pps_arm)
        self.out.fillBranch("nano_pps_rpid", pps_rpid)
        self.out.fillBranch("nano_pps_x", pps_x)
        self.out.fillBranch("nano_pps_y", pps_y)
        self.out.fillBranch("nano_pps_t", pps_t)

        self.out.fillBranch("nano_nJets", len(sel_jets))
        self.out.fillBranch("nano_jet_pt", [j.pt for j in sel_jets])
        self.out.fillBranch("nano_jet_eta", [j.eta for j in sel_jets])
        self.out.fillBranch("nano_jet_phi", [j.phi for j in sel_jets])
        self.out.fillBranch("nano_Jet_ntrk05", jet_trk05)
        self.out.fillBranch("nano_Jet_ntrk09", jet_trk09)
        self.out.fillBranch("nano_mJets", jet_sum.M())
        self.out.fillBranch("nano_yJets", jet_sum.Rapidity())
        self.out.fillBranch("nano_ptJets", jet_sum.Pt())
        self.out.fillBranch("nano_dphiJets", abs(ROOT.TVector2.Phi_mpi_pi(sel_jets[0].phi - sel_jets[1].phi)))

        return True


search_cep_mj  = lambda : JetCEPModule(channel="mj")


