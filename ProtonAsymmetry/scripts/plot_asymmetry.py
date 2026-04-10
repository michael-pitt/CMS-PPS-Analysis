import ROOT
import os
import re

# 1. Performance & Batch Setup
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT()

# --- SMART CONFIGURATION ---
MODE = "DIMUON"       # Options: "MUON", "DIMUON", or "JET"
STREAM = None         # e.g. "Muon", "JetMET"; if None, taken from MODE
RUNS = None           # e.g. [402575, 402576]; if None, auto-discover all runs for the stream

BASE_DIR = "/eos/cms/store/group/phys_diffraction/CMSLowPU2026/ntuples/data_new"
OUT_DIR = "./plots"

config = {
    "MUON": {
        "label": "1 tight/iso lepton",
        "stream": "Muon",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 1",
        "var_name": "lep_eta",
        "var_expr": "nano_lep_eta[0]",
        "x_range": (-2.5, 2.5),
        "x_title": "Leading Lepton #eta",
        "nbins": 40,
    },
    "MUONP": {
        "label": "W+ selection",
        "stream": "Muon",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 1 && nano_lep_charge[0]>0 && nano_w_mT>50",
        "var_name": "lepp_eta",
        "var_expr": "nano_lep_eta[0]",
        "x_range": (-2.5, 2.5),
        "x_title": "Leading Lepton #eta",
        "nbins": 40,
    },
    "MUONM": {
        "label": "W- selection",
        "stream": "Muon",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 1 && nano_lep_charge[0]<0 && nano_w_mT>50",
        "var_name": "lepm_eta",
        "var_expr": "nano_lep_eta[0]",
        "x_range": (-2.5, 2.5),
        "x_title": "Leading Lepton #eta",
        "nbins": 40,
    },
    "DIMUON": {
        "label": "Z #rightarrow #mu#mu",
        "stream": "Muon",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 2 && abs(nano_mll-91)<20",
        "var_name": "y_ll",
        "var_expr": "nano_yll",
        "x_range": (-2.5, 2.5),
        "x_title": "Dimuon Rapidity y_{ll}",
        "nbins": 40,
    },
    "JET": {
        "label": "MultiJet",
        "stream": "JetMET",
        "hlt": "HLT_PFJet80_L1Jet60",
        "filter": "nano_nJets >= 2",
        "var_name": "y_jj",
        "var_expr": "nano_yJets",
        "x_range": (-4.7, 4.7),
        "x_title": "Dijet Rapidity y_{jj}",
        "nbins": 40,
    },
}

cfg = dict(config[MODE])
if STREAM is None:
    STREAM = cfg["stream"]

# =========================================
# 3. Discover runs and files
# =========================================
def discover_runs(base_dir, stream):
    """
    Return a sorted list of run numbers for which base_dir/RunXXXXXX/<stream> exists
    and contains at least one ROOT file.
    """
    runs = []
    if not os.path.exists(base_dir):
        print(f"X ERROR: base dir does not exist: {base_dir}")
        return runs

    for entry in sorted(os.listdir(base_dir)):
        full_run_dir = os.path.join(base_dir, entry)
        if not os.path.isdir(full_run_dir):
            continue

        m = re.match(r"Run(\d+)$", entry)
        if not m:
            continue

        run_number = int(m.group(1))
        stream_dir = os.path.join(full_run_dir, stream)
        if not os.path.isdir(stream_dir):
            continue

        has_root = any(f.endswith(".root") for f in os.listdir(stream_dir))
        if has_root:
            runs.append(run_number)

    return runs


def get_local_files(base_dir, runs, stream):
    """
    Return all ROOT files under:
      base_dir/RunXXXXXX/<stream>/*.root
    """
    files = []
    for run in runs:
        stream_dir = os.path.join(base_dir, f"Run{run}", stream)
        if not os.path.isdir(stream_dir):
            print(f"X Warning: missing directory {stream_dir}")
            continue

        for fname in sorted(os.listdir(stream_dir)):
            if fname.endswith(".root"):
                files.append(os.path.join(stream_dir, fname))

    return files


if RUNS is None:
    RUNS = discover_runs(BASE_DIR, STREAM)

if not RUNS:
    raise RuntimeError(f"No runs found in {BASE_DIR} for stream {STREAM}")

all_files = get_local_files(BASE_DIR, RUNS, STREAM)
if not all_files:
    raise RuntimeError(f"No ROOT files found for stream {STREAM} in runs {RUNS}")

print(f"--> Found {len(all_files)} files for stream {STREAM}")
print(f"--> Runs: {RUNS}")

if len(RUNS) == 1:
    run_str = f"Run{RUNS[0]}"
else:
    run_str = f"Run{RUNS[0]}to{RUNS[-1]}"


# =========================================
# 3. RDataFrame & Physics Selection
# =========================================
df = ROOT.RDataFrame("Events", all_files)

total_events = df.Count().GetValue()
print(f"--> Total events to process: {total_events:,}")

# Enable ROOT's native RDataFrame progress bar
if hasattr(ROOT.RDF.Experimental, "AddProgressBar"):
    ROOT.RDF.Experimental.AddProgressBar(df)
else:
    print("--> (Live progress bar requires ROOT 6.28+)")
    
df_phys = (
    df.Filter(f"{cfg['hlt']} == 1", "Trigger")
      .Filter(cfg["filter"], "Physics selection")
      .Define("plot_var", cfg["var_expr"])
      .Define("is_Zplus",  "Sum(nano_pps_rpid == 3 || nano_pps_rpid == 23) > 0")
      .Define("is_Zminus", "Sum(nano_pps_rpid == 103 || nano_pps_rpid == 123) > 0")
)

# 4. Book Histograms
h_model = ROOT.RDF.TH1DModel("h",f";{cfg['x_title']};Events",cfg["nbins"],cfg["x_range"][0],cfg["x_range"][1])
h_inc = df_phys.Histo1D(h_model, "plot_var")
h_zp  = df_phys.Filter("is_Zplus && !is_Zminus","z+ only").Histo1D(h_model, "plot_var")
h_zm  = df_phys.Filter("is_Zminus && !is_Zplus","z- only").Histo1D(h_model, "plot_var")

# Force execution
hi = h_inc.GetValue()
hp = h_zp.GetValue()
hm = h_zm.GetValue()

def style_hist(h, color, width=2):
    h.SetStats(0)
    h.SetTitle("")
    h.SetLineColor(color)
    h.SetLineWidth(width)

for h, c in [(hi, ROOT.kBlack), (hp, ROOT.kRed), (hm, ROOT.kBlue)]:
    style_hist(h, c)

def normalized_ratio(h_num, h_den, name):
    """
    Build ratio of normalized shapes:
      (h_num / integral_num) / (h_den / integral_den)
    """
    r_num = h_num.Clone(f"{name}_num")
    r_den = h_den.Clone(f"{name}_den")

    if r_num.Integral() > 0:
        r_num.Scale(1.0 / r_num.Integral())
    if r_den.Integral() > 0:
        r_den.Scale(1.0 / r_den.Integral())

    ratio = r_num.Clone(name)
    ratio.Divide(r_den)
    ratio.SetStats(0)
    return ratio

r_p = normalized_ratio(hp, hi, "ratio_zplus")
r_m = normalized_ratio(hm, hi, "ratio_zminus")
    
# 6. Plotting Logic
c1 = ROOT.TCanvas("c1", "", 900, 1000)
c1.Divide(1, 2)

pad1 = c1.cd(1)
pad2 = c1.cd(2)

pad1.SetPad(0.0, 0.35, 1.0, 1.0)
pad2.SetPad(0.0, 0.00, 1.0, 0.35)

pad1.SetLogy(True)
pad2.SetGridy(True)

pad1.SetBottomMargin(0.02)
pad1.SetTopMargin(0.08)
pad1.SetLeftMargin(0.14)
pad1.SetRightMargin(0.04)

pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.30)
pad2.SetLeftMargin(0.14)
pad2.SetRightMargin(0.04)

pad1.cd()
max_y = max(hi.GetMaximum(), hp.GetMaximum(), hm.GetMaximum())
hi.SetMaximum(max(10.0, max_y * 500.0))
hi.SetMinimum(0.5)
hi.GetXaxis().SetLabelSize(0); hi.GetXaxis().SetTitleSize(0)
hi.GetYaxis().SetTitle("Events")
hi.GetYaxis().SetTitleSize(0.06); hi.GetYaxis().SetTitleOffset(1.1); hi.GetYaxis().SetLabelSize(0.05)

hi.SetLineColor(ROOT.kBlack); hi.SetLineWidth(2); hi.Draw("HIST")
hp.SetLineColor(ROOT.kRed);   hp.SetLineWidth(2); hp.Draw("HIST SAME")
hm.SetLineColor(ROOT.kBlue);  hm.SetLineWidth(2); hm.Draw("HIST SAME")

# Legend & Info
leg = ROOT.TLegend(0.18, 0.65, 0.45, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04)
leg.AddEntry(hi, "Inclusive", "l"); leg.AddEntry(hp, "PPS Tag z+", "l"); leg.AddEntry(hm, "PPS Tag z-", "l")
leg.Draw()

ltx = ROOT.TLatex(); ltx.SetNDC(); ltx.SetTextSize(0.035)
ltx.DrawLatex(0.55, 0.86, f"{cfg['label']}")
ltx.DrawLatex(0.55, 0.81, f"#bf{{Stream:}} {cfg['stream']}")
ltx.DrawLatex(0.55, 0.76, f"#bf{{Trigger:}} {cfg['hlt']}")
ltx.SetTextColor(ROOT.kBlack); ltx.DrawLatex(0.55, 0.70, f"Incl: {int(hi.GetEntries()):,}")
ltx.SetTextColor(ROOT.kRed);   ltx.DrawLatex(0.55, 0.65, f"z+: {int(hp.GetEntries()):,}")
ltx.SetTextColor(ROOT.kBlue);  ltx.DrawLatex(0.55, 0.60, f"z-: {int(hm.GetEntries()):,}")

# Bottom Pad (Ratio)
pad2.cd()
def get_norm_ratio(h_num, h_den):
    h_n = h_num.Clone(); h_d = h_den.Clone()
    if h_n.Integral() > 0: h_n.Scale(1.0/h_n.Integral())
    if h_d.Integral() > 0: h_d.Scale(1.0/h_d.Integral())
    h_n.Divide(h_d)
    return h_n

r_p = get_norm_ratio(hp, hi); r_m = get_norm_ratio(hm, hi)
r_p.GetYaxis().SetTitle("PPS / Incl (Norm.)")
r_p.GetYaxis().CenterTitle(True); r_p.GetYaxis().SetRangeUser(0.75, 1.25)
r_p.GetXaxis().SetTitle(cfg['x_title']); r_p.GetXaxis().CenterTitle(True)
r_p.GetXaxis().SetTitleSize(0.12); r_p.GetXaxis().SetLabelSize(0.1); r_p.GetXaxis().SetTitleOffset(0.9)
r_p.GetYaxis().SetTitleSize(0.08); r_p.GetYaxis().SetLabelSize(0.08); r_p.GetYaxis().SetTitleOffset(0.8)

r_p.SetMarkerStyle(20); r_p.SetMarkerColor(ROOT.kRed); r_p.SetLineColor(ROOT.kRed); r_p.Draw("P")
r_m.SetMarkerStyle(24); r_m.SetMarkerColor(ROOT.kBlue); r_m.SetLineColor(ROOT.kBlue); r_m.Draw("P SAME")

# Save
os.makedirs(OUT_DIR, exist_ok=True)
out_path = os.path.join(
    OUT_DIR,
    f"correlation_{STREAM.lower()}_{cfg['var_name']}_{run_str}.png"
)

c1.SaveAs(out_path)
print(f"--> Correlation plot saved to {out_path}")