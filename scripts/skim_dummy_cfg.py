import FWCore.ParameterSet.Config as cms
process = cms.Process("HIGHPTSKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
XXXXX
)
)

# =============== Other Statements =====================
## test process.MessageLogger = cms.Service("MessageLogger",
## test                                               destinations   = cms.untracked.vstring('detailedInfo'),
## test                                               categories      = cms.untracked.vstring('eventNumber'),
## test                                               detailedInfo    = cms.untracked.PSet(
## test                                                                         eventNumber = cms.untracked.PSet(
## test                                                                         reportEvery = cms.untracked.int32(1000)
## test                                                                                                 )
## test                                                                                                                                  ),
## test )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#Trigger Selection (done in the loaded filters)

#Further Lepton Selection
process.load('Configuration.Skimming.PbPb_ZMMSkim_cff')
process.load('Configuration.Skimming.PbPb_ZEESkim_cff')
process.load('Configuration.Skimming.PbPb_EMuSkim_cff')

process.eventFilter_dimu = cms.Path(
    process.zMMSkimSequence
    )

process.eventFilter_diele = cms.Path(
    process.zEESkimSequence
    )

process.eventFilter_emu = cms.Path(
    process.emuSkimSequence
    )

process.output_dilep = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('YYYYY'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_dimu','eventFilter_diele','eventFilter_emu')),
    dataset = cms.untracked.PSet(
    dataTier = cms.untracked.string('RECO'),
    filterName = cms.untracked.string('hiHighPt'))
)

process.output_step_dilep = cms.EndPath(process.output_dilep)

process.schedule = cms.Schedule(
    process.eventFilter_dimu,
    process.eventFilter_diele,
    process.eventFilter_emu,
    process.output_step_dilep,
)
