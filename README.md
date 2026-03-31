# CMS_PPS_Analysis

A comprehensive framework for analysis in CMS that use PPS, covering Run 2 and Run 3 data. 

## Analysis Tools
1. **Proton Asymmetry**: Correlation studies between protons and the central system (SD events)
2. **CEP Jets**: Central Exclusive Production of Jets ($M_{jj}$ vs $M_{pp}$ matching).
3. **CEP Muons**: Central Exclusive Production of Muons in High-PU runs (Run2 and Run3 data).

## Repository Structure
The repository is organized by analysis topic. Each directory contains its own NanoAOD modules, configuration data, and plotting scripts.
- `ProtonAsymmetry/`:
	- `python/ProtonAsymModule.py`: NanoAODTools module for SD correlation studies.
	- `scripts/` RDataFrame-based plotting and asymmetry fitting.
- `CEP_jets/`: 
	- `python/ProtonAsymModule.py`: Module for $M_{jj}$ vs $M_{pp}$ matching logic.
- `CEP_muons/`: 
	- `python/JetCEPModule.py`: Logic for exclusive muon production.


## Setup

**Note**: This framework requires specific track-multiplicity variables (`nTrkPV05`, `nTrkPV09`) currently available in customized NanoAOD. See the merge request here: [LINK_TO_CMSSW_PULL_REQUEST_PLACEHOLDER]

## Setup
```bash
cmsrel CMSSW_16_0_4
cd CMSSW_16_0_4/src
cmsenv
git clone https://github.com/michael-pitt/CMS_PPS_Analysis.git
scram b -j
```

## How to Run Production
### Local test

Before submitting to the batch, verify the module logic with a single file:

```bash
nano_postproc.py ./output_dir /eos/cms/store/group/phys_diffraction/CMSLowPU2026/PilotRun/Muon2/Run401866/Muon2_401866_2f4bd762-a9e9-4661-95c7-69d08dec9a17.root \
    -I CMS_PPS_Analysis.ProtonAsymmetry.ProtonAsymModule asymmetry_mu \
    -c "HLT_Mu15 == 1" \
    --bi $CMSSW_BASE/src/CMS_PPS_Analysis/ProtonAsymmetry/data/keep_data_in.txt \
    --bo $CMSSW_BASE/src/CMS_PPS_Analysis/ProtonAsymmetry/data/keep_and_drop_Data_out.txt \
    -s _Asym
```

### Batch submission (HTCondor)

PLACEHOLDER FOR THE INSTRUCTIONS
