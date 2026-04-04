# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: stepMiniAODPPS --python_filename step5_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --filein file:miniAOD.root --fileout file:miniAOD_PPS.root --step NONE --conditions 150X_mcRun3_2026_realistic_v4 --customise SimPPS/Configuration/Utils.setupPPSDirectSimMiniAOD --era Run3_2026 --mc -n -1
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Run3_2026_cff import Run3_2026
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing("analysis")
options.parseArguments()

process = cms.Process("PPSRECO", Run3_2026)

# Basic services
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(8),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary=cms.untracked.bool(False)
)

# MC 2026 GT
process.GlobalTag = GlobalTag(process.GlobalTag, "150X_mcRun3_2026_realistic_v4", "")

process.MINIAODSIMoutput = cms.OutputModule(
    "PoolOutputModule",
    compressionAlgorithm=cms.untracked.string("LZMA"),
    compressionLevel=cms.untracked.int32(4),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string("MINIAODSIM"),
        filterName=cms.untracked.string("")
    ),
    dropMetaData=cms.untracked.string("ALL"),
    eventAutoFlushCompressedSize=cms.untracked.int32(-900),
    fastCloning=cms.untracked.bool(False),
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel=cms.untracked.VPSet(
        cms.untracked.PSet(
            branch=cms.untracked.string("patPackedCandidates_packedPFCandidates__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("recoGenParticles_prunedGenParticles__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("patTriggerObjectStandAlones_slimmedPatTrigger__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("patPackedGenParticles_packedGenParticles__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("patJets_slimmedJets__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("recoVertexs_offlineSlimmedPrimaryVertices__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("recoCaloClusters_reducedEgamma_reducedESClusters_*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("recoGenJets_slimmedGenJets__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("patJets_slimmedJetsPuppi__*"),
            splitLevel=cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch=cms.untracked.string("EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*"),
            splitLevel=cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels=cms.untracked.bool(True),
    splitLevel=cms.untracked.int32(0)
)

process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# start with output-only schedule, same as cmsDriver
process.schedule = cms.Schedule(process.MINIAODSIMoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# apply the PPS simulation
from SimPPS.Configuration.Utils import setupPPSDirectSimMiniAOD
process = setupPPSDirectSimMiniAOD(process)

# add pileup protons
process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
    cms.InputTag("genPUProtons"),
    cms.InputTag("prunedGenParticles")
)

# override PPS geometry
process.load("Geometry.VeryForwardGeometry.geometryRPFromDD_2025_cfi")

# optional: 2026 geometry / alignment / optics overrides can be added below
# once the baseline job is stable

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)