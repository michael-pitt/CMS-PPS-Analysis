import ROOT
import os

# 1. Performance & Batch Setup
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT()

# --- SMART CONFIGURATION ---
MODE = 'JET' # Toggle between 'MUON' and 'JET'

config = {
    'MUON': {
        'label': 'Muon',
        'path_key': 'mu',
        'hlt': "HLT_Mu15",
        'filter': "nano_nLeptons == 1",
        'var_name': "lep_eta",
        'var_expr': "nano_lep_eta[0]",
        'x_range': [-2.5, 2.5],
        'x_title': "Leading Lepton #eta"
    },
    'JET': {
        'label': 'JetMET',
        'path_key': 'mj',
        'hlt': "HLT_PFJet80_L1Jet60",
        'filter': "nano_nJets >= 2",
        'var_name': "y_jj",
        'var_expr': "nano_yJets",
        'x_range': [-4.7, 4.7],
        'x_title': "Dijet Rapidity y_{jj}"
    }
}

cfg = config[MODE]
BASE_DIR = "/eos/cms/store/group/phys_diffraction/CMSLowPU2026/PilotRun/ntuples/"

# 2. Local File Discovery
def get_local_files(cfg):
    input_path = os.path.join(BASE_DIR, cfg['path_key'])
    files = []
    if not os.path.exists(input_path):
        print(f"X Warning: Path {input_path} does not exist.")
        return files
    for root, _, filenames in os.walk(input_path):
        for f in filenames:
            if f.endswith(".root"):
                files.append(os.path.join(root, f))
    return files

all_files = get_local_files(cfg)
if not all_files: exit(f"X ERROR: No files found for {cfg['label']}")

# Extract run range from file paths for the label
run_str = "Run2026B" 

# 3. RDataFrame & Physics Selection
df = ROOT.RDataFrame("Events", all_files)

df_phys = df.Filter(f"{cfg['hlt']} == 1").Filter(cfg['filter']) \
            .Define("plot_var", cfg['var_expr']) \
            .Define("is_Zplus",  "Sum(nano_pps_rpid == 3 || nano_pps_rpid == 23) > 0") \
            .Define("is_Zminus", "Sum(nano_pps_rpid == 103 || nano_pps_rpid == 123) > 0")

# 4. Book Histograms
h_model = ROOT.RDF.TH1DModel("h", f";{cfg['x_title']};Events", 40, *cfg['x_range'])
h_inc = df_phys.Histo1D(h_model, "plot_var")
h_zp  = df_phys.Filter("is_Zplus && !is_Zminus").Histo1D(h_model, "plot_var")
h_zm  = df_phys.Filter("is_Zminus && !is_Zplus").Histo1D(h_model, "plot_var")

# 5. Execution & Styling
hi = h_inc.GetValue(); hp = h_zp.GetValue(); hm = h_zm.GetValue()
for h in [hi, hp, hm]:
    h.SetStats(0)
    h.SetTitle("")

# 6. Plotting Logic
c1 = ROOT.TCanvas("c1", "", 900, 1000)
c1.Divide(1, 2); pad1 = c1.cd(1); pad2 = c1.cd(2)
pad1.SetPad(0, 0.35, 1, 1); pad2.SetPad(0, 0, 1, 0.35)
pad1.SetLogy(1); pad2.SetGridy()
pad1.SetBottomMargin(0.02); pad1.SetTopMargin(0.08); pad1.SetLeftMargin(0.14)
pad2.SetTopMargin(0.02); pad2.SetBottomMargin(0.3); pad2.SetLeftMargin(0.14)

pad1.cd()
hi.SetMaximum(hi.GetMaximum() * 500)
hi.GetXaxis().SetLabelSize(0); hi.GetXaxis().SetTitleSize(0)
hi.GetYaxis().SetTitle("Events")
hi.GetYaxis().SetTitleSize(0.06); hi.GetYaxis().SetTitleOffset(1.1); hi.GetYaxis().SetLabelSize(0.05)

hi.SetLineColor(ROOT.kBlack); hi.SetLineWidth(2); hi.Draw("HIST")
hp.SetLineColor(ROOT.kRed);   hp.SetLineWidth(2); hp.Draw("HIST SAME")
hm.SetLineColor(ROOT.kBlue);  hm.SetLineWidth(2); hm.Draw("HIST SAME")

# Legend & Info
leg = ROOT.TLegend(0.18, 0.65, 0.45, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04)
leg.AddEntry(hi, "Inclusive", "l"); leg.AddEntry(hp, "PPS Tag Z+", "l"); leg.AddEntry(hm, "PPS Tag Z-", "l")
leg.Draw()

ltx = ROOT.TLatex(); ltx.SetNDC(); ltx.SetTextSize(0.035)
ltx.DrawLatex(0.55, 0.86, f"#bf{{Stream:}} {cfg['label']}")
ltx.DrawLatex(0.55, 0.81, f"#bf{{Trigger:}} {cfg['hlt']}")
ltx.SetTextColor(ROOT.kBlack); ltx.DrawLatex(0.55, 0.72, f"Incl: {int(hi.GetEntries()):,}")
ltx.SetTextColor(ROOT.kRed);   ltx.DrawLatex(0.55, 0.67, f"Z+: {int(hp.GetEntries()):,}")
ltx.SetTextColor(ROOT.kBlue);  ltx.DrawLatex(0.55, 0.62, f"Z-: {int(hm.GetEntries()):,}")

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
r_p.GetYaxis().CenterTitle(True); r_p.GetYaxis().SetRangeUser(0.5, 1.5)
r_p.GetXaxis().SetTitle(cfg['x_title']); r_p.GetXaxis().CenterTitle(True)
r_p.GetXaxis().SetTitleSize(0.12); r_p.GetXaxis().SetLabelSize(0.1); r_p.GetXaxis().SetTitleOffset(0.9)
r_p.GetYaxis().SetTitleSize(0.08); r_p.GetYaxis().SetLabelSize(0.08); r_p.GetYaxis().SetTitleOffset(0.8)

r_p.SetMarkerStyle(20); r_p.SetMarkerColor(ROOT.kRed); r_p.SetLineColor(ROOT.kRed); r_p.Draw("P")
r_m.SetMarkerStyle(24); r_m.SetMarkerColor(ROOT.kBlue); r_m.SetLineColor(ROOT.kBlue); r_m.Draw("P SAME")

# Save
if not os.path.exists("plots"): os.makedirs("plots")
out_path = f"plots/correlation_{cfg['label'].lower()}_{cfg['var_name']}.png"
c1.SaveAs(out_path)
print(f"--> Correlation plot saved to {out_path}")