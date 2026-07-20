#!/usr/bin/env python3
import os, sys, math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# --- Define here the helpers ---
def safe_get(obj, attr_name, default=0):
    """Safely retrieves attributes from NanoAOD events or collections,
       catching the NanoAODTools RuntimeError if the branch is missing."""
    try:
        return getattr(obj, attr_name)
    except (RuntimeError, AttributeError):
        return default


def deltaR(obj1, obj2):
    deta = obj1.eta - obj2.eta
    dphi = ROOT.TVector2.Phi_mpi_pi(obj1.phi - obj2.phi)
    return math.hypot(deta, dphi)


JES_VARIATION_SIGNS = {
    "nominal": 0,
    "jestotalup": 1,
    "jestotaldown": -1,
}


def normalize_jes_variation(value):
    """Return the canonical name used by the three JES production passes."""
    key = str(value or "nominal").strip().replace("_", "").lower()
    aliases = {
        "nom": "nominal",
        "central": "nominal",
        "up": "jestotalup",
        "jesup": "jestotalup",
        "totalup": "jestotalup",
        "down": "jestotaldown",
        "jesdown": "jestotaldown",
        "totaldown": "jestotaldown",
    }
    key = aliases.get(key, key)
    if key not in JES_VARIATION_SIGNS:
        allowed = "nominal, jesTotalUp, jesTotalDown"
        raise ValueError(f"Unknown JES variation '{value}'. Expected one of: {allowed}.")
    return {
        "nominal": "nominal",
        "jestotalup": "jesTotalUp",
        "jestotaldown": "jesTotalDown",
    }[key]


def shift_jet_pt_mass(pt, mass, uncertainty, variation_sign):
    """Apply a fractional JES uncertainty to a jet four-vector scale."""
    if not math.isfinite(uncertainty) or uncertainty < 0.0 or uncertainty >= 1.0:
        raise RuntimeError(f"Invalid fractional JES uncertainty: {uncertainty}")
    scale = 1.0 + variation_sign * uncertainty
    return pt * scale, mass * scale


class AnalysisJet:
    def __init__(self, source, pt, eta, phi, mass, jes_uncertainty=0.0):
        self.source = source
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass
        self.jes_uncertainty = jes_uncertainty

        self._p4 = ROOT.TLorentzVector()
        self._p4.SetPtEtaPhiM(pt, eta, phi, mass)

    def p4(self):
        return self._p4

    def __getattr__(self, name):
        return getattr(self.source, name)
  
def get_nu_p4(lep_vec, met_pt, met_phi):
    """Reconstructs the neutrino 4-vector using the W mass constraint."""
    MW = 80.379
    px_nu = met_pt * math.cos(met_phi)
    py_nu = met_pt * math.sin(met_phi)

    Lambda = (MW**2) / 2.0 + lep_vec.Px() * px_nu + lep_vec.Py() * py_nu
    A = lep_vec.Pt()**2
    B = -2.0 * Lambda * lep_vec.Pz()
    C = (lep_vec.E()**2) * (met_pt**2) - Lambda**2

    delta = B**2 - 4 * A * C

    # Solve the quadratic equation for pz
    if delta >= 0:
        pz1 = (-B + math.sqrt(delta)) / (2.0 * A)
        pz2 = (-B - math.sqrt(delta)) / (2.0 * A)
        # Take the solution with the smallest absolute value
        pz_nu = pz1 if abs(pz1) < abs(pz2) else pz2
    else:
        # If complex, take the real part
        pz_nu = -B / (2.0 * A)

    nu_vec = ROOT.TLorentzVector()
    e_nu = math.sqrt(met_pt**2 + pz_nu**2)
    nu_vec.SetPxPyPzE(px_nu, py_nu, pz_nu, e_nu)
    return nu_vec  

class AsymmetryModule(Module):
    def __init__(self, channel="mu", year=2026, jes_variation=None,
                 jes_json=None, jes_uncertainty_name=None, jes_jec_name=None):
        self.channel = channel
        self.year = year
        self.rp_ids = {"45": [3, 23], "56": [103, 123]}

        requested_variation = (
            os.environ.get("PROTON_ASYM_JES_VARIATION", "nominal")
            if jes_variation is None else jes_variation
        )
        self.jes_variation = normalize_jes_variation(requested_variation)
        self.jes_sign = JES_VARIATION_SIGNS[
            self.jes_variation.replace("_", "").lower()
        ]
        self.jes_json = jes_json or os.environ.get("PROTON_ASYM_JES_JSON", "")
        self.jes_uncertainty_name = jes_uncertainty_name or os.environ.get(
            "PROTON_ASYM_JES_UNCERTAINTY_NAME", ""
        )
        self.jes_jec_name = jes_jec_name or os.environ.get(
            "PROTON_ASYM_JES_JEC_NAME", ""
        )
        self.jes_uncertainty_evaluator = None
        self.jes_jec_evaluator = None
        
        # --- HARDCODED KINEMATIC PARAMETERS ---
        self.min_muon_pt = 15.0      # For W/Z Control Regions
        self.min_ele_pt = 15.0       # For W/Z Control Regions
        self.min_soft_muon_pt = 3.0  # For Inclusive Dimuon Region
        self.min_jet_pt = 25.0       # For Jet/MJ selections

    @property
    def has_jes_shift(self):
        return self.jes_sign != 0

    def _load_jes(self):
        if not self.has_jes_shift:
            return
        if not self.jes_json or not self.jes_uncertainty_name or not self.jes_jec_name:
            raise RuntimeError(
                "A non-nominal JES pass requires PROTON_ASYM_JES_JSON, "
                "PROTON_ASYM_JES_UNCERTAINTY_NAME, and "
                "PROTON_ASYM_JES_JEC_NAME. The uncertainty correction supplies "
                "delta_JES; the nominal compound JEC is needed for low-pT "
                "CorrT1METJet objects in the Type-1 PUPPI MET propagation."
            )

        if not os.path.isfile(self.jes_json):
            raise RuntimeError(f"JES payload does not exist: {self.jes_json}")
        if self.jes_json.endswith(".gz"):
            with open(self.jes_json, "rb") as payload:
                if payload.read(2) != b"\x1f\x8b":
                    raise RuntimeError(
                        f"JES payload is not a gzip file: {self.jes_json}. "
                        "Check that an HTML login page was not downloaded instead."
                    )

        try:
            import correctionlib
        except ImportError as exc:
            raise RuntimeError(
                "JES requested but correctionlib is not available in this "
                "CMSSW environment."
            ) from exc

        cset = correctionlib.CorrectionSet.from_file(self.jes_json)
        try:
            self.jes_uncertainty_evaluator = cset[self.jes_uncertainty_name]
        except KeyError as exc:
            raise RuntimeError(
                f"JES uncertainty correction '{self.jes_uncertainty_name}' "
                f"was not found in {self.jes_json}."
            ) from exc

        try:
            self.jes_jec_evaluator = cset[self.jes_jec_name]
        except KeyError:
            compound = getattr(cset, "compound", {})
            try:
                self.jes_jec_evaluator = compound[self.jes_jec_name]
            except KeyError as exc:
                raise RuntimeError(
                    f"Nominal JEC correction '{self.jes_jec_name}' was not "
                    f"found as a correction or compound correction in {self.jes_json}."
                ) from exc

        print("[ProtonAsymModule] JES systematic enabled")
        print(f"[ProtonAsymModule]   variation: {self.jes_variation}")
        print(f"[ProtonAsymModule]   JSON: {self.jes_json}")
        print(
            "[ProtonAsymModule]   uncertainty: "
            f"{self.jes_uncertainty_name}"
        )
        print(f"[ProtonAsymModule]   nominal JEC: {self.jes_jec_name}")

    def _correction_input_value(self, input_info, obj, event, pt, mass=0.0):
        name = input_info.name
        key = name.lower().replace("_", "").replace(" ", "")
        input_type = getattr(input_info, "type", "")

        if input_type == "string":
            raise RuntimeError(
                f"Correction input '{name}' is a string. JES uncertainty and "
                "nominal JEC evaluators must be selected by their correction "
                "names, not by passing a variation string."
            )
        if key in ("jetpt", "pt"):
            return pt
        if key in ("jetmass", "mass"):
            return mass
        if key in ("jeteta", "eta"):
            return obj.eta
        if key in ("jetphi", "phi"):
            return obj.phi
        if key in ("jetarea", "area", "jeta"):
            area = safe_get(obj, "area", None)
            if area is None:
                raise RuntimeError(f"Correction input '{name}' requires jet area.")
            return area
        if key == "rho":
            rho = safe_get(event, "fixedGridRhoFastjetAll", None)
            if rho is None:
                raise RuntimeError(
                    "Correction requires fixedGridRhoFastjetAll, but the branch "
                    "is missing from the input selection."
                )
            return rho
        if key == "run":
            run = safe_get(event, "run", None)
            if run is None:
                raise RuntimeError(
                    "Correction requires the run number, but the branch is missing."
                )
            return int(run)

        raise RuntimeError(
            f"Do not know how to provide correction input '{name}'."
        )

    def _evaluate(self, evaluator, obj, event, pt, mass=0.0):
        inputs = [
            self._correction_input_value(inp, obj, event, pt, mass)
            for inp in evaluator.inputs
        ]
        return evaluator.evaluate(*inputs)

    def _jes_uncertainty(self, obj, event, pt):
        if not self.has_jes_shift:
            return 0.0
        uncertainty = float(
            self._evaluate(self.jes_uncertainty_evaluator, obj, event, pt)
        )
        if not math.isfinite(uncertainty) or uncertainty < 0.0 or uncertainty >= 1.0:
            raise RuntimeError(
                f"Invalid JES uncertainty {uncertainty} for jet "
                f"(pt={pt}, eta={obj.eta})."
            )
        return uncertainty

    def _analysis_jet(self, jet, event):
        uncertainty = self._jes_uncertainty(jet, event, jet.pt)
        pt, mass = shift_jet_pt_mass(
            jet.pt, jet.mass, uncertainty, self.jes_sign
        )

        return AnalysisJet(
            jet,
            pt,
            jet.eta,
            jet.phi,
            mass,
            jes_uncertainty=uncertainty,
        )

    def _stored_jet_type1_pt(self, jet):
        """Return the nominal Type-1 jet pT and its MET eligibility.

        The 15 GeV threshold is applied after removing the raw muon component,
        matching the NanoAODTools Type-1 MET recipe.
        """
        raw_factor = safe_get(jet, "rawFactor", None)
        muon_subtr = safe_get(jet, "muonSubtrFactor", None)
        ne_em = safe_get(jet, "neEmEF", None)
        ch_em = safe_get(jet, "chEmEF", None)
        if None in (raw_factor, muon_subtr, ne_em, ch_em):
            raise RuntimeError(
                "Proper JES propagation to PUPPI MET requires Jet_rawFactor, "
                "Jet_muonSubtrFactor, Jet_neEmEF, and Jet_chEmEF."
            )

        raw_pt = jet.pt * (1.0 - raw_factor)
        corrected_no_mu_pt = jet.pt * (1.0 - muon_subtr)
        muon_pt = raw_pt * muon_subtr
        type1_pt = corrected_no_mu_pt + muon_pt
        eligible = corrected_no_mu_pt > 15.0 and (ne_em + ch_em) < 0.9
        return type1_pt, eligible

    def _low_pt_jet_type1_pt(self, jet, event):
        """Reconstruct a low-pT Type-1 jet with the configured nominal JEC."""
        raw_pt = safe_get(jet, "rawPt", None)
        muon_subtr = safe_get(jet, "muonSubtrFactor", None)
        if raw_pt is None or muon_subtr is None:
            raise RuntimeError(
                "Proper JES propagation to PUPPI MET requires "
                "CorrT1METJet_rawPt and CorrT1METJet_muonSubtrFactor."
            )

        jec = float(self._evaluate(self.jes_jec_evaluator, jet, event, raw_pt))
        if not math.isfinite(jec) or jec <= 0.0:
            raise RuntimeError(
                f"Invalid nominal JEC factor {jec} for CorrT1METJet "
                f"(rawPt={raw_pt}, eta={jet.eta})."
            )

        corrected_no_mu_pt = raw_pt * (1.0 - muon_subtr) * jec
        muon_pt = raw_pt * muon_subtr
        type1_pt = corrected_no_mu_pt + muon_pt
        return type1_pt, corrected_no_mu_pt > 15.0

    def _propagate_jes_to_puppi_met(self, event, jets):
        """Propagate the JES delta from all Type-1 jets to nominal PUPPI MET."""
        met_pt = safe_get(event, "PuppiMET_pt", None)
        met_phi = safe_get(event, "PuppiMET_phi", None)
        if met_pt is None or met_phi is None:
            raise RuntimeError("PuppiMET_pt/PuppiMET_phi are missing.")
        if not self.has_jes_shift:
            return met_pt, met_phi

        met_px = met_pt * math.cos(met_phi)
        met_py = met_pt * math.sin(met_phi)

        type1_jets = []
        for jet in jets:
            type1_pt, eligible = self._stored_jet_type1_pt(jet)
            if eligible:
                type1_jets.append((jet, type1_pt))

        n_low_pt = safe_get(event, "nCorrT1METJet", None)
        if n_low_pt is None:
            raise RuntimeError(
                "nCorrT1METJet is missing. Keep nCorrT1METJet and "
                "CorrT1METJet_* for non-nominal JES production."
            )
        for jet in Collection(event, "CorrT1METJet"):
            type1_pt, eligible = self._low_pt_jet_type1_pt(jet, event)
            if eligible:
                type1_jets.append((jet, type1_pt))

        for jet, type1_pt in type1_jets:
            uncertainty = self._jes_uncertainty(jet, event, type1_pt)
            delta_pt = self.jes_sign * uncertainty * type1_pt
            met_px -= delta_pt * math.cos(jet.phi)
            met_py -= delta_pt * math.sin(jet.phi)

        return math.hypot(met_px, met_py), math.atan2(met_py, met_px)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self._load_jes()
        
        # MPI Summaries 
        self.out.branch("nano_NMPI05", "I") # nMPI at PV with pT > 0.5
        self.out.branch("nano_NMPI09", "I") # nMPI at PV with pT > 0.9   
        
        # Proton Summaries
        self.out.branch("nano_nProtons", "I")
        self.out.branch("nano_pps_arm", "I", lenVar="nano_nProtons")
        self.out.branch("nano_pps_rpid", "I", lenVar="nano_nProtons")
        self.out.branch("nano_pps_x", "F", lenVar="nano_nProtons")
        self.out.branch("nano_pps_y", "F", lenVar="nano_nProtons")
        
        # Jet Summaries
        self.out.branch("nano_nJets", "I")
        self.out.branch("nano_jet_pt", "F", lenVar="nano_nJets")
        self.out.branch("nano_jet_eta", "F", lenVar="nano_nJets")
        self.out.branch("nano_jet_phi", "F", lenVar="nano_nJets")
        self.out.branch("nano_jet_jesUncertainty", "F", lenVar="nano_nJets")
        self.out.branch("nano_Jet_ntrk05", "I", lenVar="nano_nJets")
        self.out.branch("nano_Jet_ntrk09", "I", lenVar="nano_nJets")
        self.out.branch("nano_mJets", "F")
        self.out.branch("nano_yJets", "F")
        
        # Lepton Summaries
        self.out.branch("nano_nLeptons", "I")
        self.out.branch("nano_lep_pt", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_eta", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_phi", "F", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_charge", "I", lenVar="nano_nLeptons")
        self.out.branch("nano_lep_ntrk05", "I", lenVar="nano_nLeptons") 
        self.out.branch("nano_lep_ntrk09", "I", lenVar="nano_nLeptons")
        self.out.branch("nano_w_mT", "F")
        self.out.branch("nano_w_pt", "F")
        self.out.branch("nano_w_phi", "F")
        self.out.branch("nano_w_y", "F")
        self.out.branch("nano_w_m", "F")
        self.out.branch("nano_mll", "F")
        self.out.branch("nano_yll", "F")
        self.out.branch("nano_ptll", "F")
        
        # event branches
        self.out.branch("nano_jesVariation", "I")
        self.out.branch("nano_puppiMET_pt", "F")
        self.out.branch("nano_puppiMET_phi", "F")
        self.out.branch("nano_Mall", "F")
        self.out.branch("nano_Yall", "F")
        
        # cutflow histogram
        self.h_cutflow = ROOT.TH1D("Cutflow", "Event Cutflow", 3, 0, 3)
        self.h_cutflow.GetXaxis().SetBinLabel(1, "All Initial Events")
        self.h_cutflow.GetXaxis().SetBinLabel(2, "Pass HLT")
        self.h_cutflow.GetXaxis().SetBinLabel(3, "Pass Module Filter")
        
        self.h_cutflow.SetBinContent(1, inputTree.GetEntries())
        
        self.events_seen = 0
        self.events_passed = 0

    def analyze(self, event):
        
        # count events that survived the "-c CUT" in nano_postproc.py
        self.events_seen += 1
        
        # Load Collections
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        jets = Collection(event, "Jet")
        analysis_jets = [self._analysis_jet(j, event) for j in jets]
        analysis_met_pt, analysis_met_phi = self._propagate_jes_to_puppi_met(
            event, jets
        )
        protons = Collection(event, "PPSLocalTrack")
                
        # Object Selections & IDs
        # ----------------------------------------------------------------------
        # Muons: Loose (for ZCR/veto) and Tight (for WCR)
        loose_mu = [m for m in muons if m.pt > self.min_muon_pt and abs(m.eta) < 2.5 and m.looseId and m.pfRelIso04_all < 0.25]
        tight_mu = [m for m in loose_mu if m.tightId and m.pfRelIso04_all < 0.15]
        
        # Ensure soft muons are sorted by pT to correctly grab the leading pair
        soft_mu = sorted([m for m in muons if m.pt > self.min_soft_muon_pt and abs(m.eta) < 2.4 and m.looseId], key=lambda x: x.pt, reverse=True)
        
        # Electrons: Loose (cutBased=2) and Tight (cutBased=4)
        loose_el = [e for e in electrons if e.pt > self.min_ele_pt and abs(e.eta) < 2.5 and not (1.4442 < abs(e.eta) < 1.566) and e.cutBased >= 2]
        tight_el = [e for e in loose_el if abs(e.eta) < 2.5 and e.cutBased >= 4]

        # Combine and sort by pT
        loose_leps = sorted(loose_mu + loose_el, key=lambda x: x.pt, reverse=True)
        tight_leps = sorted(tight_mu + tight_el, key=lambda x: x.pt, reverse=True)

        # Object Overlap Removal (Jets vs Leptons)
        # ----------------------------------------------------------------------
        # Remove jets that fall within dR < 0.4 of any loose lepton
        raw_jets = sorted(
            [
                j for j in analysis_jets
                if j.pt > self.min_jet_pt and abs(j.eta) < 4.7
            ],
            key=lambda jet: jet.pt,
            reverse=True,
        )
        sel_jets = []
        for j in raw_jets:
            has_overlap = False
            for l in loose_leps:
                if deltaR(j, l) < 0.4:
                    has_overlap = True
                    break
            if not has_overlap:
                sel_jets.append(j)

        jet_sum = ROOT.TLorentzVector()
        for j in sel_jets: jet_sum += j.p4()

        # Event Overlap Removal & Signal Region Definitions
        # ----------------------------------------------------------------------
        # Booleans calculated independently to prevent ELIF cross-talk bugs
        
        # DY: Exactly 2 loose leptons, opposite sign, same flavor
        is_ZCR = (len(loose_leps) == 2 and loose_leps[0].charge != loose_leps[1].charge and loose_leps[0].pdgId == -loose_leps[1].pdgId)
        
        # W+jets: Exactly 1 tight lepton, AND exactly 1 loose lepton (vetoes events with a 2nd loose lepton)
        is_WCR = (len(tight_leps) == 1 and len(loose_leps) == 1)
        
        # MJ Control Region: >= 2 isolated jets, strictly 0 loose leptons
        is_mj = (len(sel_jets) >= 2 and sel_jets[0].pt > 140.0)
        
        # Soft Dimuon: >= 2 soft muons, leading pair is opposite sign
        is_dimuon_inclusive = (len(soft_mu) >= 2 and soft_mu[0].charge != soft_mu[1].charge)

        # 5. Event Filtering based on Channel
        # ----------------------------------------------------------------------
        leptons_to_save = []
        
        if self.channel == "mu":
            if is_WCR and abs(tight_leps[0].pdgId) == 13: leptons_to_save = tight_leps
            else: return False

        elif self.channel == "dimuon":
            if is_ZCR and abs(loose_leps[0].pdgId) == 13: leptons_to_save = loose_leps
            else: return False

        elif self.channel == "softmm":
            if is_dimuon_inclusive: leptons_to_save = soft_mu
            else: return False

        elif self.channel == "el":
            if is_ZCR and abs(loose_leps[0].pdgId) == 11: leptons_to_save = loose_leps
            elif is_WCR and abs(tight_leps[0].pdgId) == 11: leptons_to_save = tight_leps
            else: return False
            
        elif self.channel == "mj":
            if not is_mj: return False
            
        elif self.channel == "zb":
            leptons_to_save = loose_leps
            
        else:
            return False

        # Calculations: Dileptons & W variables
        # ----------------------------------------------------------------------
        mll = ptll = -999.0
        yll = -999.0
        w_mT = w_pt = w_phi = w_m = -999.0
        w_y = -999.0
        Mall = Yall = -999.0
        
        v_all = ROOT.TLorentzVector()
        
        # 1. Add all jets to the global system
        for j in sel_jets:
            v_all += j.p4()
        
        if len(leptons_to_save) >= 2:
            dilep = leptons_to_save[0].p4() + leptons_to_save[1].p4()
            mll, yll, ptll = dilep.M(), dilep.Rapidity(), dilep.Pt()
            v_all += leptons_to_save[0].p4()
            v_all += leptons_to_save[1].p4()
            
        if len(leptons_to_save) >= 1 and (self.channel in ["mu", "el"] and is_WCR):
            met_pt = analysis_met_pt
            met_phi = analysis_met_phi
            dphi = ROOT.TVector2.Phi_mpi_pi(leptons_to_save[0].phi - met_phi)
            w_mT = math.sqrt(2 * leptons_to_save[0].pt * met_pt * (1 - math.cos(dphi)))
            
            w_vec = ROOT.TVector2(leptons_to_save[0].pt * math.cos(leptons_to_save[0].phi) + met_pt * math.cos(met_phi),
                                  leptons_to_save[0].pt * math.sin(leptons_to_save[0].phi) + met_pt * math.sin(met_phi))
            w_pt, w_phi = w_vec.Mod(), w_vec.Phi()
            
            # Full reconstruction using W mass constraint
            lep_p4 = ROOT.TLorentzVector()
            lep_p4.SetPtEtaPhiM(leptons_to_save[0].pt, leptons_to_save[0].eta, leptons_to_save[0].phi, leptons_to_save[0].mass)            
            nu_p4 = get_nu_p4(lep_p4, met_pt, met_phi)
            w_p4 = lep_p4 + nu_p4
            w_pt, w_phi, w_m = w_p4.Pt(), w_p4.Phi(), w_p4.M()
            w_y = w_p4.Rapidity() if w_p4.E() > abs(w_p4.Pz()) else -999.0
            
            v_all += lep_p4
            v_all += nu_p4
            
        elif is_mj or self.channel == "zb":
            for l in leptons_to_save:
                v_all += l.p4()
                
        # Calculate Global Variables
        if v_all.E() > 0:
            Mall = v_all.M()
            if v_all.E() > abs(v_all.Pz()): # Protection for rapidity calculation
                Yall = v_all.Rapidity()

        
        # Protons Logic (LocalTrack mapping)
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

        
        # TRACK MULTIPLICITIES & MPI LOGIC
        
        # Read PV0 tracks
        pv0_ntrk05 = safe_get(event, "PV_ntrk0p5", 0)
        pv0_ntrk09 = safe_get(event, "PV_ntrk0p9", 0)
        
        # Extract object tracks (Default to 0 if branch missing)
        jet_trk05 = [safe_get(j, "ntrk0p5", 0) for j in sel_jets]
        jet_trk09 = [safe_get(j, "ntrk0p9", 0) for j in sel_jets]
        
        lep_trk05 = [safe_get(l, "ntrk0p5", 0) for l in leptons_to_save]
        lep_trk09 = [safe_get(l, "ntrk0p9", 0) for l in leptons_to_save]
        
        # N_trk^MPI = N_trk^PV - Sum(N_trk^objects)
        nmpi05 = pv0_ntrk05 - sum(jet_trk05) - sum(lep_trk05)
        nmpi09 = pv0_ntrk09 - sum(jet_trk09) - sum(lep_trk09)
        
        # Safety catch: prevent negative MPI in case of overlapping object cones
        if nmpi05 < 0: nmpi05 = 0
        if nmpi09 < 0: nmpi09 = 0
        
        # 7. Fill Branches
        self.out.fillBranch("nano_NMPI05", nmpi05)
        self.out.fillBranch("nano_NMPI09", nmpi09)

        self.out.fillBranch("nano_nProtons", len(pps_arm))
        self.out.fillBranch("nano_pps_arm", pps_arm)
        self.out.fillBranch("nano_pps_rpid", pps_rpid)
        self.out.fillBranch("nano_pps_x", pps_x)
        self.out.fillBranch("nano_pps_y", pps_y)

        self.out.fillBranch("nano_nJets", len(sel_jets))
        self.out.fillBranch("nano_jet_pt", [j.pt for j in sel_jets])
        self.out.fillBranch("nano_jet_eta", [j.eta for j in sel_jets])
        self.out.fillBranch("nano_jet_phi", [j.phi for j in sel_jets])
        self.out.fillBranch(
            "nano_jet_jesUncertainty",
            [j.jes_uncertainty for j in sel_jets],
        )
        self.out.fillBranch("nano_Jet_ntrk05", jet_trk05)
        self.out.fillBranch("nano_Jet_ntrk09", jet_trk09)
        self.out.fillBranch("nano_mJets", jet_sum.M() if len(sel_jets) > 0 else -999.0)
        self.out.fillBranch("nano_yJets", jet_sum.Rapidity() if jet_sum.E() > abs(jet_sum.Pz()) else -999.0)

        self.out.fillBranch("nano_nLeptons", len(leptons_to_save))
        self.out.fillBranch("nano_lep_pt", [l.pt for l in leptons_to_save])
        self.out.fillBranch("nano_lep_eta", [l.eta for l in leptons_to_save])
        self.out.fillBranch("nano_lep_phi", [l.phi for l in leptons_to_save])
        self.out.fillBranch("nano_lep_charge", [l.charge for l in leptons_to_save])
        self.out.fillBranch("nano_lep_ntrk05", lep_trk05)
        self.out.fillBranch("nano_lep_ntrk09", lep_trk09)
        self.out.fillBranch("nano_mll", mll)
        self.out.fillBranch("nano_yll", yll)
        self.out.fillBranch("nano_ptll", ptll)
        self.out.fillBranch("nano_w_mT", w_mT)
        self.out.fillBranch("nano_w_pt", w_pt)
        self.out.fillBranch("nano_w_phi", w_phi)
        self.out.fillBranch("nano_w_y", w_y)
        self.out.fillBranch("nano_w_m", w_m)

        self.out.fillBranch("nano_jesVariation", self.jes_sign)
        self.out.fillBranch("nano_puppiMET_pt", analysis_met_pt)
        self.out.fillBranch("nano_puppiMET_phi", analysis_met_phi)
        self.out.fillBranch("nano_Mall", Mall)
        self.out.fillBranch("nano_Yall", Yall)
        
        # total accepted events
        self.events_passed += 1

        return True

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        # Fill the final bins and write to the file
        self.h_cutflow.SetBinContent(2, self.events_seen)
        self.h_cutflow.SetBinContent(3, self.events_passed)
        
        outputFile.cd()
        self.h_cutflow.Write()

def _asymmetry_factory(channel, variation):
    return lambda: AsymmetryModule(channel=channel, jes_variation=variation)


# Nominal factories preserve the existing nano_postproc.py interface.
asymmetry_mu     = _asymmetry_factory("mu", "nominal")
asymmetry_dimuon = _asymmetry_factory("dimuon", "nominal")
asymmetry_softmm = _asymmetry_factory("softmm", "nominal")
asymmetry_el     = _asymmetry_factory("el", "nominal")
asymmetry_mj     = _asymmetry_factory("mj", "nominal")
asymmetry_zb     = _asymmetry_factory("zb", "nominal")

# Explicit factories make the three production passes reproducible without
# relying on a mutable variation environment variable.
asymmetry_mu_jesTotalUp         = _asymmetry_factory("mu", "jesTotalUp")
asymmetry_mu_jesTotalDown       = _asymmetry_factory("mu", "jesTotalDown")
asymmetry_dimuon_jesTotalUp     = _asymmetry_factory("dimuon", "jesTotalUp")
asymmetry_dimuon_jesTotalDown   = _asymmetry_factory("dimuon", "jesTotalDown")
asymmetry_softmm_jesTotalUp     = _asymmetry_factory("softmm", "jesTotalUp")
asymmetry_softmm_jesTotalDown   = _asymmetry_factory("softmm", "jesTotalDown")
asymmetry_el_jesTotalUp         = _asymmetry_factory("el", "jesTotalUp")
asymmetry_el_jesTotalDown       = _asymmetry_factory("el", "jesTotalDown")
asymmetry_mj_jesTotalUp         = _asymmetry_factory("mj", "jesTotalUp")
asymmetry_mj_jesTotalDown       = _asymmetry_factory("mj", "jesTotalDown")
asymmetry_zb_jesTotalUp         = _asymmetry_factory("zb", "jesTotalUp")
asymmetry_zb_jesTotalDown       = _asymmetry_factory("zb", "jesTotalDown")
