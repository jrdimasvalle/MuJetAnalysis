import FWCore.ParameterSet.Config as cms

process = cms.Process("HLTBending")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')
process.load("MuJetAnalysis.HLTBendingAngle.HLTBendingAngle_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
      'file:/fdata/hepx/store/user/jrdv009/DarkSUSY_mH_125_mGammaD_20000_ctau0_14TeV_madgraph452_bridge224_LHE_pythia8_GEN_SIM_80k_v3/DarkSUSY_mH_125_mGammaD_20000_ctau0_14TeV_madgraph452_bridge224_LHE_pythia8_GEN_SIM_80k_v3/e007db91b5aa2bed7c0fc47ce82a4274/out_sim_23_1_Qio.root',
      'file:/fdata/hepx/store/user/jrdv009/DarkSUSY_mH_125_mGammaD_20000_ctau0_14TeV_madgraph452_bridge224_LHE_pythia8_GEN_SIM_80k_v3/DarkSUSY_mH_125_mGammaD_20000_ctau0_14TeV_madgraph452_bridge224_LHE_pythia8_GEN_SIM_80k_v3/e007db91b5aa2bed7c0fc47ce82a4274/out_sim_44_1_Jlk.root'
   )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_ana.root"),
	closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.HLTBendingAngle)



## messages
print
#print 'Input files:'
print '----------------------------------------'
print process.source.fileNames
print
print 'Output file:'
print '----------------------------------------'
print process.TFileService.fileName
print
