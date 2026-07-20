import importlib.util
import math
import pathlib
import sys
import tempfile
import types
import unittest
from unittest import mock


class FakeLorentzVector:
    def SetPtEtaPhiM(self, pt, eta, phi, mass):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass


def load_module():
    root = types.ModuleType("ROOT")
    root.PyConfig = types.SimpleNamespace(IgnoreCommandLineOptions=False)
    root.TLorentzVector = FakeLorentzVector
    root.TVector2 = types.SimpleNamespace(
        Phi_mpi_pi=lambda value: math.atan2(math.sin(value), math.cos(value))
    )

    eventloop = types.ModuleType(
        "PhysicsTools.NanoAODTools.postprocessing.framework.eventloop"
    )
    eventloop.Module = object
    datamodel = types.ModuleType(
        "PhysicsTools.NanoAODTools.postprocessing.framework.datamodel"
    )
    datamodel.Collection = lambda event, name: []

    mocked = {
        "ROOT": root,
        "PhysicsTools": types.ModuleType("PhysicsTools"),
        "PhysicsTools.NanoAODTools": types.ModuleType("PhysicsTools.NanoAODTools"),
        "PhysicsTools.NanoAODTools.postprocessing": types.ModuleType(
            "PhysicsTools.NanoAODTools.postprocessing"
        ),
        "PhysicsTools.NanoAODTools.postprocessing.framework": types.ModuleType(
            "PhysicsTools.NanoAODTools.postprocessing.framework"
        ),
        "PhysicsTools.NanoAODTools.postprocessing.framework.eventloop": eventloop,
        "PhysicsTools.NanoAODTools.postprocessing.framework.datamodel": datamodel,
    }
    sys.modules.update(mocked)

    path = pathlib.Path(__file__).parents[1] / "python" / "ProtonAsymModule.py"
    spec = importlib.util.spec_from_file_location("proton_asym_module_for_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MODULE = load_module()


class Input:
    def __init__(self, name, input_type="real"):
        self.name = name
        self.type = input_type


class ConstantEvaluator:
    def __init__(self, value, inputs):
        self.value = value
        self.inputs = [Input(name) for name in inputs]

    def evaluate(self, *values):
        return self.value


class Obj:
    def __init__(self, **values):
        self.__dict__.update(values)


class JESVariationTest(unittest.TestCase):
    def make_module(self, variation):
        result = MODULE.AsymmetryModule(channel="mu", jes_variation=variation)
        result.jes_uncertainty_evaluator = ConstantEvaluator(
            0.10, ["JetEta", "JetPt"]
        )
        result.jes_jec_evaluator = ConstantEvaluator(
            2.0, ["JetA", "JetEta", "JetPt", "Rho"]
        )
        result.jes_l1_evaluator = ConstantEvaluator(
            1.0, ["JetA", "JetEta", "JetPt", "Rho"]
        )
        return result

    def test_variation_aliases_and_invalid_value(self):
        self.assertEqual(MODULE.normalize_jes_variation("up"), "jesTotalUp")
        self.assertEqual(MODULE.normalize_jes_variation("total_down"), "jesTotalDown")
        self.assertEqual(MODULE.normalize_jes_variation("10up"), "jetPt10Up")
        with self.assertRaises(ValueError):
            MODULE.normalize_jes_variation("flavorQCDUp")

    def test_analysis_jet_scales_pt_and_mass_only(self):
        module = self.make_module("jesTotalUp")
        event = Obj()
        jet = Obj(pt=100.0, mass=20.0, eta=1.2, phi=-0.4)

        shifted = module._analysis_jet(jet, event)

        self.assertAlmostEqual(shifted.pt, 110.0)
        self.assertAlmostEqual(shifted.mass, 22.0)
        self.assertEqual(shifted.eta, jet.eta)
        self.assertEqual(shifted.phi, jet.phi)
        self.assertAlmostEqual(shifted.jes_uncertainty, 0.10)

    def test_fixed_ten_percent_stress_shift_needs_no_payload(self):
        module = MODULE.AsymmetryModule(channel="mj", jes_variation="jetPt10Down")
        event = Obj(PuppiMET_pt=37.0, PuppiMET_phi=0.7)
        jet = Obj(pt=100.0, mass=20.0, eta=1.2, phi=-0.4)

        shifted = module._analysis_jet(jet, event)

        self.assertFalse(module.has_jes_shift)
        self.assertTrue(module.has_jet_pt10_shift)
        self.assertAlmostEqual(shifted.pt, 90.0)
        self.assertAlmostEqual(shifted.mass, 18.0)
        self.assertEqual(shifted.jes_uncertainty, 0.0)
        self.assertEqual(
            module._propagate_jes_to_puppi_met(event, []), (37.0, 0.7)
        )

    def test_mj_fixed_shift_factories_are_exported(self):
        self.assertEqual(MODULE.asymmetry_mj_jetPt10Up().jet_pt10_sign, 1)
        self.assertEqual(MODULE.asymmetry_mj_jetPt10Down().jet_pt10_sign, -1)
        self.assertTrue(MODULE.asymmetry_mj().apply_central_jec)
        self.assertTrue(MODULE.asymmetry_mj_jetPt10Up().apply_central_jec)
        self.assertFalse(MODULE.asymmetry_mj_promptReco().apply_central_jec)

    def test_central_jec_undoes_stored_correction_before_recorrecting(self):
        module = MODULE.AsymmetryModule(
            channel="mj", jes_variation="jetPt10Up", apply_central_jec=True
        )
        module.central_jec_evaluator = ConstantEvaluator(
            1.5, ["JetA", "JetEta", "JetPt", "Rho"]
        )
        event = Obj(Rho_fixedGridRhoFastjetAll=8.0)
        jet = Obj(
            pt=100.0,
            mass=20.0,
            eta=0.3,
            phi=0.0,
            area=0.5,
            rawFactor=0.2,
        )

        shifted = module._analysis_jet(jet, event)

        # Raw jet: (80, 16); new central JEC: x1.5; stress shift: x1.1.
        self.assertAlmostEqual(shifted.pt, 132.0)
        self.assertAlmostEqual(shifted.mass, 26.4)
        self.assertAlmostEqual(shifted.central_jec_factor, 1.5)

    def test_promptreco_factory_keeps_the_stored_nanoaod_jec(self):
        module = MODULE.asymmetry_mj_promptReco()
        event = Obj()
        jet = Obj(
            pt=100.0,
            mass=20.0,
            eta=0.3,
            phi=0.0,
            rawFactor=0.2,
        )

        shifted = module._analysis_jet(jet, event)

        self.assertAlmostEqual(shifted.pt, 100.0)
        self.assertAlmostEqual(shifted.mass, 20.0)
        self.assertFalse(module.apply_central_jec)

    def test_type1_met_is_recomputed_from_raw_met_and_l1_reference(self):
        jet = Obj(
            pt=100.0,
            mass=20.0,
            eta=0.3,
            phi=0.0,
            area=0.5,
            rawFactor=0.2,
            muonSubtrFactor=0.1,
            neEmEF=0.1,
            chEmEF=0.1,
        )
        event = Obj(
            PuppiMET_pt=50.0,
            PuppiMET_phi=0.0,
            RawPuppiMET_pt=200.0,
            RawPuppiMET_phi=0.0,
            Rho_fixedGridRhoFastjetAll=8.0,
            nCorrT1METJet=0,
        )
        MODULE.Collection = lambda current_event, name: []

        up_pt, up_phi = self.make_module("jesTotalUp")._propagate_jes_to_puppi_met(
            event, [jet]
        )
        down_pt, down_phi = self.make_module(
            "jesTotalDown"
        )._propagate_jes_to_puppi_met(event, [jet])

        # raw-no-muon pT = 100*(1-.2)*(1-.1) = 72 GeV.
        # L1 pT = 72; nominal pT = 144; varied pT = 158.4/129.6.
        # MET_x = RawMET_x - (pT_varied - pT_L1).
        self.assertAlmostEqual(up_pt, 113.6)
        self.assertAlmostEqual(down_pt, 142.4)
        self.assertAlmostEqual(up_phi, 0.0)
        self.assertAlmostEqual(down_phi, 0.0)

    def test_low_pt_type1_jet_uses_nominal_jec_before_met_shift(self):
        low_pt_jet = Obj(
            rawPt=10.0,
            eta=0.2,
            phi=0.0,
            area=0.5,
            muonSubtrFactor=0.0,
        )
        event = Obj(
            PuppiMET_pt=50.0,
            PuppiMET_phi=0.0,
            RawPuppiMET_pt=50.0,
            RawPuppiMET_phi=0.0,
            nCorrT1METJet=1,
            Rho_fixedGridRhoFastjetAll=8.0,
        )
        MODULE.Collection = lambda current_event, name: [low_pt_jet]

        shifted_pt, shifted_phi = self.make_module(
            "jesTotalUp"
        )._propagate_jes_to_puppi_met(event, [])

        # pT_L1=10, pT_varied=22, so Type-1 removes 12 GeV from raw MET.
        self.assertAlmostEqual(shifted_pt, 38.0)
        self.assertAlmostEqual(shifted_phi, 0.0)

    def test_type1_threshold_is_applied_after_the_jes_shift(self):
        low_pt_jet = Obj(
            rawPt=7.4,
            eta=0.2,
            phi=0.0,
            area=0.5,
            muonSubtrFactor=0.0,
        )
        event = Obj(
            PuppiMET_pt=50.0,
            PuppiMET_phi=0.0,
            RawPuppiMET_pt=50.0,
            RawPuppiMET_phi=0.0,
            nCorrT1METJet=1,
            Rho_fixedGridRhoFastjetAll=8.0,
        )
        MODULE.Collection = lambda current_event, name: [low_pt_jet]

        up_pt, _ = self.make_module("jesTotalUp")._propagate_jes_to_puppi_met(
            event, []
        )
        down_pt, _ = self.make_module("jesTotalDown")._propagate_jes_to_puppi_met(
            event, []
        )

        # Nominal pT=14.8: JES-up migrates above 15 GeV, JES-down stays below.
        self.assertAlmostEqual(up_pt, 41.12)
        self.assertAlmostEqual(down_pt, 50.0)

    def test_nominal_pass_needs_no_payload_or_corrt1_collection(self):
        module = MODULE.AsymmetryModule(channel="mu", jes_variation="nominal")
        event = Obj(PuppiMET_pt=37.0, PuppiMET_phi=0.7)
        self.assertEqual(module._propagate_jes_to_puppi_met(event, []), (37.0, 0.7))

    def test_2026_defaults_are_the_official_cat_payload_and_tags(self):
        module = MODULE.AsymmetryModule(channel="mu", jes_variation="jesTotalUp")
        self.assertIn(
            "Run3-26Prompt-Summer24-NanoAODv15/2026-07-15/jet_jerc.json.gz",
            module.jes_json,
        )
        self.assertEqual(
            module.jes_uncertainty_name,
            "Summer24Prompt26_V1_MC_Total_AK4PFPuppi",
        )
        self.assertEqual(
            module.jes_jec_name,
            "Summer24Prompt26_V1_MC_L1L2L3Res_AK4PFPuppi",
        )
        self.assertEqual(
            module.jes_l1_name,
            "Summer24Prompt26_V1_MC_L1FastJet_AK4PFPuppi",
        )

    def test_payload_loader_resolves_compound_nominal_jec(self):
        uncertainty = ConstantEvaluator(0.10, ["JetEta", "JetPt"])
        nominal_jec = ConstantEvaluator(1.0, ["JetA", "JetEta", "JetPt", "Rho"])
        l1_jec = ConstantEvaluator(0.9, ["JetA", "JetEta", "JetPt", "Rho"])

        class FakeCorrectionSet:
            compound = {"nominal_jec": nominal_jec}

            def __getitem__(self, name):
                if name == "total_uncertainty":
                    return uncertainty
                if name == "l1_jec":
                    return l1_jec
                # correctionlib 2.7 raises IndexError("map::at") for a name
                # absent from the ordinary correction map. The requested
                # nominal JEC is intentionally in the compound map instead.
                raise IndexError("map::at")

        correctionlib = types.ModuleType("correctionlib")
        correctionlib.CorrectionSet = types.SimpleNamespace(
            from_file=lambda path: FakeCorrectionSet()
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            payload = pathlib.Path(tmpdir) / "payload.json.gz"
            payload.write_bytes(b"\x1f\x8b")
            module = MODULE.AsymmetryModule(
                channel="mu",
                jes_variation="jesTotalUp",
                jes_json=str(payload),
                jes_uncertainty_name="total_uncertainty",
                jes_jec_name="nominal_jec",
                jes_l1_name="l1_jec",
            )
            with mock.patch.dict(sys.modules, {"correctionlib": correctionlib}):
                module._load_jes()

        self.assertIs(module.jes_uncertainty_evaluator, uncertainty)
        self.assertIs(module.jes_jec_evaluator, nominal_jec)
        self.assertIs(module.jes_l1_evaluator, l1_jec)

    def test_payload_loader_rejects_html_saved_as_gzip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = pathlib.Path(tmpdir) / "payload.json.gz"
            payload.write_bytes(b"<!DOCTYPE html>")
            module = MODULE.AsymmetryModule(
                channel="mu",
                jes_variation="jesTotalUp",
                jes_json=str(payload),
                jes_uncertainty_name="total_uncertainty",
                jes_jec_name="nominal_jec",
                jes_l1_name="l1_jec",
            )
            with self.assertRaisesRegex(RuntimeError, "not a gzip file"):
                module._load_jes()


if __name__ == "__main__":
    unittest.main()
