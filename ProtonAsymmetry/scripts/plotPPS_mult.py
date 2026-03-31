import ROOT
import os
import math

# 1. Performance & Batch Setup
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT()

# --- SMART CONFIGURATION ---
# Options: 'MUON', 'JET', 'ZEROBIAS'
MODE = 'MUON' 

config = {
    'MUON': {'label': 'Muon', 'path_key': 'mu', 'hlt': 'HLT_Mu15', 'filter': 'nano_nLeptons >= 1'},
    'JET':  {'label': 'JetMET', 'path_key': 'mj', 'hlt': 'HLT_PFJet80_L1Jet60', 'filter': 'nano_nJets >= 2'},
    'ZEROBIAS': {'label': 'ZeroBias', 'path_key': 'zb', 'hlt': 'HLT_ZeroBias', 'filter': '1'}
}

cfg = config[MODE]
BASE_DIR = "/eos/cms/store/group/phys_diffraction/CMSLowPU2026/PilotRun/ntuples/"

# 2. Local File Discovery (Crawling EOS structure)
def get_local_files(cfg):
    input_path = os.path.join(BASE_DIR, cfg['path_key'])
    files = []
    if not os.path.exists(input_path):
        print(f"X Warning: Path {input_path} does not exist.")
        return files
    
    for root, dirs, filenames in os.walk(input_path):
        for f in filenames:
            if f.endswith(".root"):
                files.append(os.path.join(root, f))
    return files

all_files = get_local_files(cfg)
if not all_files: exit(f"X ERROR: No files found in {BASE_DIR}/{cfg['path_key']}")
print(f"--> Found {len(all_files)} files for {cfg['label']}")

# 3. RDataFrame Setup
df = ROOT.RDataFrame("Events", all_files)

# Use custom branches produced by Asymmetry.py
# Note: nProtons logic adjusted to use the 'nano_pps_rpid' branch
rp_ids = [3, 23, 103, 123]
for rid in rp_ids:
    # Use the custom branch 'nano_pps_rpid' created in our Asymmetry module
    df = df.Define(f"n_{rid}", f"Sum(nano_pps_rpid == {rid})")

df_filtered = df.Define("has_45", "n_3 > 0 || n_23 > 0") \
                .Define("has_56", "n_103 > 0 || n_123 > 0") \
                .Define("has_either", "has_45 || has_56") \
                .Define("has_both", "has_45 && has_56")

# 4. Book Histograms
h_probes = [df_filtered.Histo1D(ROOT.RDF.TH1DModel(f"h_{rid}", "", 7, -0.5, 6.5), f"n_{rid}") for rid in rp_ids]
h_45 = df_filtered.Histo1D(ROOT.RDF.TH1DModel("h45", "", 2, -0.5, 1.5), "has_45")
h_56 = df_filtered.Histo1D(ROOT.RDF.TH1DModel("h56", "", 2, -0.5, 1.5), "has_56")
h_either = df_filtered.Histo1D(ROOT.RDF.TH1DModel("heither", "", 2, -0.5, 1.5), "has_either")
h_both = df_filtered.Histo1D(ROOT.RDF.TH1DModel("hboth", "", 2, -0.5, 1.5), "has_both")

# 5. Execute Plotting (2x2)
c1 = ROOT.TCanvas("c1", "Multiplicities", 1200, 1000); c1.Divide(2, 2)
fractions = []

for i, prob in enumerate(h_probes):
    pad = c1.cd(i + 1); pad.SetLogy(1)
    hist = prob.GetValue()
    total = hist.GetEntries()
    # Frac >= 1 is total minus the zero-bin
    frac = (total - hist.GetBinContent(1)) / total if total > 0 else 0
    fractions.append(frac)
    
    hist.SetStats(0)
    hist.SetMarkerStyle(20)
    hist.GetXaxis().SetTitle(f"Number of tracks in RP {rp_ids[i]}")
    hist.GetYaxis().SetTitle("Events")
    hist.Draw("PE")
    
    # Poisson Fit
    f_p = ROOT.TF1(f"f_{i}", "[0]*TMath::Poisson(x, [1])", 0, 6)
    f_p.SetParameters(total, hist.GetMean())
    f_p.SetLineColor(ROOT.kRed); f_p.SetLineStyle(2)
    f_p.Draw("SAME")
    ROOT.SetOwnership(f_p, False)
    
    # Info Labels
    ltx = ROOT.TLatex(); ltx.SetNDC(); ltx.SetTextSize(0.04)
    ltx.DrawLatex(0.45, 0.85, f"#bf{{RP {rp_ids[i]}}}")
    ltx.SetTextColor(ROOT.kRed+1); ltx.DrawLatex(0.45, 0.80, f"#mu = {hist.GetMean():.4f}")
    ltx.SetTextColor(ROOT.kBlack); ltx.DrawLatex(0.45, 0.75, f"Frac #geq 1: {frac:.4f}")
    ltx.SetTextSize(0.03); ltx.DrawLatex(0.45, 0.70, f"Stream: {cfg['label']}")

# 6. Execute Summary Bar Plot
c2 = ROOT.TCanvas("c2", "Summary", 1000, 600)
h_sum = ROOT.TH1F("h_sum", f"PPS Occupancy Summary ({cfg['label']})", 8, 0, 8)
h_sum.SetStats(0)
h_sum.GetYaxis().SetTitle("Fraction of Events with #geq 1 Track")

labels = ["RP 3", "RP 23", "RP 103", "RP 123", "Arm 45", "Arm 56", "Either", "Both"]
total_ev = h_45.GetEntries()
sum_vals = fractions + [
    h_45.GetBinContent(2)/total_ev if total_ev > 0 else 0, 
    h_56.GetBinContent(2)/total_ev if total_ev > 0 else 0,
    h_either.GetBinContent(2)/total_ev if total_ev > 0 else 0, 
    h_both.GetBinContent(2)/total_ev if total_ev > 0 else 0
]

for i, (val, lab) in enumerate(zip(sum_vals, labels)):
    h_sum.SetBinContent(i+1, val)
    h_sum.GetXaxis().SetBinLabel(i+1, lab)

h_sum.SetFillColor(ROOT.kAzure+7); h_sum.SetLineColor(ROOT.kBlack)
h_sum.SetMinimum(0); h_sum.SetMaximum(max(sum_vals)*1.4)
h_sum.Draw("BAR TEXT0")

# 7. Save
output_suffix = cfg['label'].lower()
if not os.path.exists("plots"): os.makedirs("plots")
c1.SaveAs(f"plots/pps_multiplicity_{output_suffix}.png")
c2.SaveAs(f"plots/pps_summary_{output_suffix}.png")

print(f"--> Success. Plots saved in 'plots/' using ntuples from EOS.")