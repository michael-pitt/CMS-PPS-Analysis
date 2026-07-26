import argparse
import math
import os
import re

import ROOT
import cmsstyle as CMS

# 1. Performance & Batch Setup
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT()
ROOT.TH1.SetDefaultSumw2(True)

# --- SMART CONFIGURATION ---
DEFAULT_MODE = "MUON"
STREAM = None         # e.g. "Muon", "JetMET"; if None, taken from MODE
RUNS = None           # e.g. [402575, 402576]; if None, auto-discover all runs for the stream
EXCLUDED_RUNS = {403194}
TOTAL_LUMI = 2.04     # fb^-1
SQRT_S_TEV = 13.6

BASE_DIR = "/eos/cms/store/group/phys_diffraction/CMSLowPU2026/ntuples/data"
OUT_DIR = "./plots"

# Project convention for internal AN/PAS/paper-draft figures.
CMS.SetExtraText("Internal")
CMS.SetLumi(TOTAL_LUMI, unit="fb", run="", round_lumi=2)
CMS.SetEnergy(SQRT_S_TEV)
CMS.setCMSStyle()
CMS.cmsGrid(False)

COLOR_ZPLUS = ROOT.kBlue + 2
COLOR_ZMINUS = ROOT.kRed + 2

config = {
    "MUON": {
        "label": "W boson selection",
        "stream": "Muon_mu",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 1 && nano_w_mT > 50",
        "var_name": "lep_eta",
        "var_expr": "nano_lep_eta[0]",
        "x_range": (-2.5, 2.5),
        "x_title": "Leading Lepton #eta",
        # "var_name": "Yall",
        # "var_expr": "nano_Yall",
        # "x_range": (-5.0, 5.0),
        # "x_title": "Central System Rapidity Y_{all}",
        "nbins": 40,
    },
    "MUONP": {
        "label": "W+ selection",
        "stream": "Muon_mu",
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
        "stream": "Muon_mu",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 1 && nano_lep_charge[0]<0 && nano_w_mT>50",
        "var_name": "lepm_eta",
        "var_expr": "nano_lep_eta[0]",
        "x_range": (-2.5, 2.5),
        "x_title": "Leading Lepton #eta",
        "nbins": 40,
    },
    "DIMUON": {
        "label": "Z boson selection",
        "stream": "Muon_dimuon",
        "hlt": "HLT_Mu15",
        "filter": "nano_nLeptons == 2 && abs(nano_mll-91)<20",
        "var_name": "y_ll",
        "var_expr": "nano_yll",
        "x_range": (-2.5, 2.5),
        "x_title": "Dimuon Rapidity y_{ll}",
        # "var_name": "Yall",
        # "var_expr": "nano_Yall",
        # "x_range": (-3.0, 3.0),
        # "x_title": "Central System Rapidity Y_{all}",
        "nbins": 40,
    },
    "JET": {
        "label": "Multijet selection",
        "stream": "JetMET",
        "hlt": "HLT_PFJet110",
        "filter": "nano_nJets >= 2 && nano_jet_pt[0] > 140.0",
        "var_name": "y_jets",
        "var_expr": "nano_yJets",
        "x_range": (-5.0, 5.0),
        "x_title": "Multijet rapidity y_{jets}",
        "nbins": 40,
    },
}

parser = argparse.ArgumentParser(
    description="Plot normalized proton-tagged and inclusive shapes by signal region."
)
parser.add_argument(
    "--mode",
    choices=tuple(config),
    default=DEFAULT_MODE,
    help=f"signal region to plot (default: {DEFAULT_MODE})",
)
parser.add_argument(
    "--double-arm",
    choices=("keep", "exclude"),
    default="keep",
    help=(
        "treatment of events with MultiRP tags in both PPS arms: "
        "'keep' includes them in both arm histograms (default), while "
        "'exclude' removes them from both arm histograms"
    ),
)
args = parser.parse_args()
MODE = args.mode
DOUBLE_ARM_POLICY = args.double_arm

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


def get_local_files(base_dir, run, stream):
    """
    Return all ROOT files under base_dir/RunXXXXXX/<stream>/ for one run.
    """
    stream_dir = os.path.join(base_dir, f"Run{run}", stream)
    if not os.path.isdir(stream_dir):
        print(f"X Warning: missing directory {stream_dir}", flush=True)
        return []

    files = []
    for fname in sorted(os.listdir(stream_dir)):
        if fname.endswith(".root"):
            files.append(os.path.join(stream_dir, fname))

    return files


if RUNS is None:
    RUNS = discover_runs(BASE_DIR, STREAM)

if EXCLUDED_RUNS:
    RUNS = [run for run in RUNS if run not in EXCLUDED_RUNS]
    print(f"--> Excluded runs: {sorted(EXCLUDED_RUNS)}")
    
if not RUNS:
    raise RuntimeError(f"No runs found in {BASE_DIR} for stream {STREAM}")

print(f"--> Processing {len(RUNS)} runs for stream {STREAM}", flush=True)
print(f"--> Runs: {RUNS}", flush=True)
print(f"--> Double-arm MultiRP policy: {DOUBLE_ARM_POLICY}", flush=True)


# =========================================
# 3. Process and accumulate one run at a time
# =========================================
def add_histogram(total, current, name):
    if total is None:
        total = current.Clone(name)
        total.SetDirectory(0)
    else:
        total.Add(current)
    return total


hi = None
hp = None
hm = None
total_events = 0
total_selected = 0
total_double_arm = 0
total_files = 0
processed_runs = []
has_progress_bar = hasattr(ROOT.RDF.Experimental, "AddProgressBar")

if not has_progress_bar:
    print("--> Live event progress requires ROOT 6.28+; run-level progress remains enabled")

for run_index, run in enumerate(RUNS, start=1):
    run_files = get_local_files(BASE_DIR, run, STREAM)
    if not run_files:
        print(
            f"--> [{run_index}/{len(RUNS)}] Run {run}: no ROOT files, skipping",
            flush=True,
        )
        continue

    total_files += len(run_files)
    print(
        f"--> [{run_index}/{len(RUNS)}] Run {run}: "
        f"reading {len(run_files)} ROOT files",
        flush=True,
    )

    df = ROOT.RDataFrame("Events", run_files)
    if has_progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(df)

    df_phys = (
        df.Filter(f"{cfg['hlt']} == 1", "Trigger")
          .Filter(cfg["filter"], "Physics selection")
          .Define("plot_var", cfg["var_expr"])
          .Define(
              "is_Zplus",
              "Sum(nano_pps_rpid == 3) > 0 && Sum(nano_pps_rpid == 23) > 0",
          )
          .Define(
              "is_Zminus",
              "Sum(nano_pps_rpid == 103) > 0 && Sum(nano_pps_rpid == 123) > 0",
          )
          .Define("is_double_arm_multirp", "is_Zplus && is_Zminus")
    )
    if DOUBLE_ARM_POLICY == "keep":
        df_zp = df_phys.Filter("is_Zplus", "z+ MultiRP (double-arm kept)")
        df_zm = df_phys.Filter("is_Zminus", "z- MultiRP (double-arm kept)")
    else:
        df_zp = df_phys.Filter("is_Zplus && !is_Zminus", "z+ MultiRP only")
        df_zm = df_phys.Filter("is_Zminus && !is_Zplus", "z- MultiRP only")

    total_count = df.Count()
    selected_count = df_phys.Count()
    zp_count = df_zp.Count()
    zm_count = df_zm.Count()
    double_arm_count = df_phys.Filter(
        "is_double_arm_multirp", "double-arm MultiRP"
    ).Count()

    model_args = (
        cfg["nbins"],
        cfg["x_range"][0],
        cfg["x_range"][1],
    )
    h_inc = df_phys.Histo1D(
        ROOT.RDF.TH1DModel(
            f"h_inc_{run}", f";{cfg['x_title']};Events", *model_args
        ),
        "plot_var",
    )
    h_zp = df_zp.Histo1D(
        ROOT.RDF.TH1DModel(
            f"h_zp_{run}", f";{cfg['x_title']};Events", *model_args
        ),
        "plot_var",
    )
    h_zm = df_zm.Histo1D(
        ROOT.RDF.TH1DModel(
            f"h_zm_{run}", f";{cfg['x_title']};Events", *model_args
        ),
        "plot_var",
    )

    # All actions above are booked before the first GetValue(), so ROOT executes
    # one event loop for this run and then returns the remaining results.
    run_total = int(total_count.GetValue())
    run_selected = int(selected_count.GetValue())
    run_zp = int(zp_count.GetValue())
    run_zm = int(zm_count.GetValue())
    run_double_arm = int(double_arm_count.GetValue())

    hi = add_histogram(hi, h_inc.GetValue(), "h_inclusive_all_runs")
    hp = add_histogram(hp, h_zp.GetValue(), "h_zplus_all_runs")
    hm = add_histogram(hm, h_zm.GetValue(), "h_zminus_all_runs")

    total_events += run_total
    total_selected += run_selected
    total_double_arm += run_double_arm
    processed_runs.append(run)
    print(
        f"--> [{run_index}/{len(RUNS)}] Run {run}: done; "
        f"events={run_total:,}, selected={run_selected:,}, "
        f"PPS +z={run_zp:,}, PPS -z={run_zm:,}, "
        f"double-arm MultiRP={run_double_arm:,}",
        flush=True,
    )

if not processed_runs:
    raise RuntimeError(f"No ROOT files found for stream {STREAM} in runs {RUNS}")

print(
    f"--> Completed {len(processed_runs)} runs and {total_files} files: "
    f"events={total_events:,}, selected={total_selected:,}, "
    f"double-arm MultiRP={total_double_arm:,}",
    flush=True,
)

if len(processed_runs) == 1:
    run_str = f"Run{processed_runs[0]}"
else:
    run_str = f"Run{processed_runs[0]}to{processed_runs[-1]}"

def style_hist(h, color, width=2):
    h.SetStats(0)
    h.SetTitle("")
    h.SetLineColor(color)
    h.SetLineWidth(width)


for h, c in [(hi, ROOT.kBlack), (hp, COLOR_ZPLUS), (hm, COLOR_ZMINUS)]:
    style_hist(h, c)


def in_range_events(histogram):
    """Return the unweighted number of entries drawn inside the x-axis range."""
    return histogram.Integral(1, histogram.GetNbinsX())


print("--> Plot event summary (in-range entries):", flush=True)
print(f"    inclusive : {in_range_events(hi):,.0f}", flush=True)
print(f"    PPS +z    : {in_range_events(hp):,.0f}", flush=True)
print(f"    PPS -z    : {in_range_events(hm):,.0f}", flush=True)


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


def ratio_axis_range(*ratios):
    """Choose 0.1-spaced bounds, using 0.8--1.2 unless data require more room."""
    lower_values = []
    upper_values = []
    for ratio in ratios:
        for bin_index in range(1, ratio.GetNbinsX() + 1):
            value = ratio.GetBinContent(bin_index)
            error = ratio.GetBinError(bin_index)
            if not math.isfinite(value) or value <= 0:
                continue
            lower_values.append(value - error)
            upper_values.append(value + error)

    if not upper_values:
        return 0.8, 1.2

    lower = min(0.8, math.floor(10.0 * min(lower_values)) / 10.0)
    upper = max(1.2, math.ceil(10.0 * max(upper_values)) / 10.0)
    return lower, upper


ratio_y_min, ratio_y_max = ratio_axis_range(r_p, r_m)
ratio_major_intervals = max(1, int(round(10.0 * (ratio_y_max - ratio_y_min))))
print(
    f"--> Ratio y-axis range: {ratio_y_min:.1f} to {ratio_y_max:.1f} "
    "with 0.1 major spacing",
    flush=True,
)

# 6. Plotting Logic
c1 = ROOT.TCanvas("c1", "", 900, 900)
c1.Divide(1, 2)

pad1 = c1.cd(1)
pad2 = c1.cd(2)

pad1.SetPad(0.0, 0.35, 1.0, 1.0)
pad2.SetPad(0.0, 0.00, 1.0, 0.35)

pad1.SetLogy(True)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad2.SetTickx(1)
pad2.SetTicky(1)

pad1.SetBottomMargin(0.02)
pad1.SetTopMargin(0.12)
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.04)

pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.32)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.04)

pad1.cd()
max_y = max(hi.GetMaximum(), hp.GetMaximum(), hm.GetMaximum())
hi.SetMaximum(max(10.0, max_y * 500.0))
hi.SetMinimum(0.5)
hi.GetXaxis().SetLabelSize(0); hi.GetXaxis().SetTitleSize(0)
hi.GetYaxis().SetTitle("Events")
hi.GetYaxis().SetTitleSize(0.070); hi.GetYaxis().SetTitleOffset(1.0); hi.GetYaxis().SetLabelSize(0.060)

hi.SetLineColor(ROOT.kBlack); hi.SetLineWidth(3); hi.SetLineStyle(1)
hi.SetMarkerColor(ROOT.kBlack); hi.SetMarkerStyle(20); hi.SetMarkerSize(0.8)
hp.SetLineColor(COLOR_ZPLUS); hp.SetLineWidth(3); hp.SetLineStyle(2)
hp.SetMarkerColor(COLOR_ZPLUS); hp.SetMarkerStyle(22); hp.SetMarkerSize(0.9)
hm.SetLineColor(COLOR_ZMINUS); hm.SetLineWidth(3); hm.SetLineStyle(7)
hm.SetMarkerColor(COLOR_ZMINUS); hm.SetMarkerStyle(21); hm.SetMarkerSize(0.8)

hi.Draw("HIST")
hp.Draw("HIST SAME")
hm.Draw("HIST SAME")
hi.Draw("E1 X0 SAME")
hp.Draw("E1 X0 SAME")
hm.Draw("E1 X0 SAME")

# Legend & Info
leg = CMS.cmsLeg(0.65, 0.67, 0.94, 0.86, textSize=0.055)
leg.AddEntry(hi, "inclusive", "lp")
leg.AddEntry(hp, "PPS +z", "lp")
leg.AddEntry(hm, "PPS -z", "lp")
leg.Draw()

ltx = ROOT.TLatex()
ltx.SetNDC()
ltx.SetTextFont(42)
ltx.SetTextSize(0.050)
ltx.DrawLatex(0.18, 0.62, cfg["label"])
if DOUBLE_ARM_POLICY == "exclude":
    ltx.SetTextSize(0.038)
    ltx.DrawLatex(0.18, 0.57, "Double-arm tags excluded")

CMS.CMS_lumi(pad1, iPosX=11)

# Bottom Pad (Ratio)
pad2.cd()
r_p.GetYaxis().SetTitle("tag. / incl. (norm.)")
r_p.GetYaxis().CenterTitle(True)
r_p.GetYaxis().SetRangeUser(ratio_y_min, ratio_y_max)
r_p.GetYaxis().SetNdivisions(ratio_major_intervals, False)
r_p.GetXaxis().SetTitle(cfg['x_title']); r_p.GetXaxis().CenterTitle(True)
r_p.GetXaxis().SetTitleSize(0.13); r_p.GetXaxis().SetLabelSize(0.11); r_p.GetXaxis().SetTitleOffset(1.0)
r_p.GetYaxis().SetTitleSize(0.095); r_p.GetYaxis().SetLabelSize(0.090); r_p.GetYaxis().SetTitleOffset(0.72)

r_p.SetMarkerStyle(22); r_p.SetMarkerColor(COLOR_ZPLUS); r_p.SetLineColor(COLOR_ZPLUS)
r_m.SetMarkerStyle(21); r_m.SetMarkerColor(COLOR_ZMINUS); r_m.SetLineColor(COLOR_ZMINUS)
r_p.Draw("AXIS")

unity = ROOT.TLine(cfg["x_range"][0], 1.0, cfg["x_range"][1], 1.0)
unity.SetLineColor(ROOT.kBlack)
unity.SetLineStyle(2)
unity.Draw("SAME")
r_p.Draw("E1 X0 SAME")
r_m.Draw("E1 X0 SAME")

pad1.RedrawAxis()
pad2.RedrawAxis()

# Save
os.makedirs(OUT_DIR, exist_ok=True)
double_arm_suffix = (
    "" if DOUBLE_ARM_POLICY == "keep" else "_doubleArm-exclude"
)
out_stem = os.path.join(
    OUT_DIR,
    f"correlation_{STREAM.lower()}_{cfg['var_name']}_{run_str}"
    f"{double_arm_suffix}"
)

for extension in ("pdf", "png"):
    c1.SaveAs(f"{out_stem}.{extension}")
print(f"--> Correlation plots saved to {out_stem}.pdf and {out_stem}.png")
