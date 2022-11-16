// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

// added by lupz begin
#ifndef HomogeneousField

#define HomogeneousField

#endif

// include KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

// includes O2
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
// includes O2Physics
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "TableHelper.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include <CCDB/BasicCCDBManager.h>
// added by lupz end

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::track; // added by lupz

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;

constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;

constexpr static int pairTypeEE = VarManager::kJpsiToEE;
constexpr static int pairTypeMuMu = VarManager::kJpsiToMuMu;
constexpr static int pairTypeEMu = VarManager::kElectronMuon;

// added by lupz begin
namespace o2::aod
{
DECLARE_SOA_TABLE(AmbiguousTracksMid, "AOD", "AMBIGUOUSTRACK", //! Table for tracks which are not uniquely associated with a collision
                  o2::soa::Index<>, o2::aod::ambiguous::TrackId, o2::aod::ambiguous::BCIdSlice, o2::soa::Marker<2>);
} // namespace o2::aod
constexpr static uint32_t gkTrackFillMapWithAmbi = VarManager::ObjTypes::Track | VarManager::ObjTypes::AmbiTrack;
// added by lupz end

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a coma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString mixVarsString = fConfigMixingVariables.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
      }
    }
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& event)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<TEventFillMap>(event);
    // TODO: make this condition at compile time
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
    }
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }
      eventSel(1);
    } else {
      eventSel(0);
    }

    if (fMixHandler != nullptr) {
      int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
      hash(hh);
    }
  }

  void processSkimmed(MyEvents::iterator const& event)
  {
    runEventSelection<gkEventFillMap>(event);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", false);
  // TODO: Add process functions subscribing to Framework Collision
};

struct AnalysisTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  // The list of cuts should contain all the track cuts needed later in analysis, including
  //  for candidate electron selection (+ eventual prefilter cuts) and other needs like quarkonium - hadron correlations
  // The user must ensure using them properly in the tasks downstream
  // NOTE: For now, the candidate electron cuts must be provided first, then followed by any other needed selections
  Configurable<string> fConfigCuts{"cfgTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        histDirNames += Form("TrackBarrel_%s;", cut.GetName());
      }

      DefineHistograms(fHistMan, histDirNames.Data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runTrackSelection(TEvent const& event, TTracks const& tracks)
  {
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    // fill event information which might be needed in histograms/cuts that combine track and event properties
    VarManager::FillEvent<TEventFillMap>(event);

    trackSel.reserve(tracks.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& track : tracks) {
      filterMap = 0;
      VarManager::FillTrack<TTrackFillMap>(track);
      if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      iCut = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      }
      trackSel(static_cast<int>(filterMap));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(event, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed tracks", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", false);
};

struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fMuonCuts;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        histDirNames += Form("TrackMuon_%s;", cut.GetName());
      }

      DefineHistograms(fHistMan, histDirNames.Data()); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runMuonSelection(TEvent const& event, TMuons const& muons)
  {
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    VarManager::FillEvent<TEventFillMap>(event);

    muonSel.reserve(muons.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& muon : muons) {
      filterMap = 0;
      VarManager::FillTrack<TMuonFillMap>(muon);
      if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }

      iCut = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut).GetName()), VarManager::fgValues);
          }
        }
      }
      muonSel(static_cast<int>(filterMap));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyMuonTracks const& muons)
  {
    runMuonSelection<gkEventFillMap, gkMuonFillMap>(event, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", false);
};

struct AnalysisEventMixing {
  OutputObj<THashList> fOutputList{"output"};
  // Here one should provide the list of electron and muon candidate cuts in the same order as specified in the above
  // single particle selection tasks to preserve the correspondence between the track cut name and its
  //  bit position in the cuts bitmap
  // TODO: Create a configurable to specify exactly on which of the bits one should run the event mixing
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  Filter filterTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::dqanalysisflags::isMuonSelected > 0;

  HistogramManager* fHistMan;
  // NOTE: The bit mask is required to run pairing just based on the desired electron/muon candidate cuts
  uint32_t fTwoTrackFilterMask = 0;
  uint32_t fTwoMuonFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fTrackMuonHistNames;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  void init(o2::framework::InitContext& context)
  {
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    if (context.mOptions.get<bool>("processBarrelSkimmed") || context.mOptions.get<bool>("processBarrelVnSkimmed")) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsBarrelMEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelMEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelMEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fTrackHistNames.push_back(names);
          fTwoTrackFilterMask |= (uint32_t(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processMuonSkimmed") || context.mOptions.get<bool>("processMuonVnSkimmed")) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsMuonMEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonMEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonMEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fMuonHistNames.push_back(names);
          fTwoMuonFilterMask |= (uint32_t(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processBarrelMuonSkimmed")) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) { // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) {
            std::vector<TString> names = {
              Form("PairsEleMuMEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuMEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuMEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            fTrackMuonHistNames.push_back(names);
            fTwoTrackFilterMask |= (uint32_t(1) << icut);
            fTwoMuonFilterMask |= (uint32_t(1) << icut);
          }
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  template <int TPairType, typename TTracks1, typename TTracks2>
  void runMixedPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;
    if constexpr (TPairType == pairTypeMuMu) {
      ncuts = fMuonHistNames.size();
      histNames = fMuonHistNames;
    }
    if constexpr (TPairType == pairTypeEMu) {
      ncuts = fTrackMuonHistNames.size();
      histNames = fTrackMuonHistNames;
    }

    uint32_t twoTrackFilter = 0;
    for (auto& track1 : tracks1) {
      for (auto& track2 : tracks2) {
        if constexpr (TPairType == VarManager::kJpsiToEE) {
          twoTrackFilter = uint32_t(track1.isBarrelSelected()) & uint32_t(track2.isBarrelSelected()) & fTwoTrackFilterMask;
        }
        if constexpr (TPairType == VarManager::kJpsiToMuMu) {
          twoTrackFilter = uint32_t(track1.isMuonSelected()) & uint32_t(track2.isMuonSelected()) & fTwoMuonFilterMask;
        }
        if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = uint32_t(track1.isBarrelSelected()) & uint32_t(track2.isMuonSelected()) & fTwoTrackFilterMask;
        }

        if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
          continue;
        }
        VarManager::FillPairME<TPairType>(track1, track2);

        constexpr bool eventHasQvector = (VarManager::ObjTypes::ReducedEventQvector > 0);
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<TPairType>(track1, track2);
        }

        for (unsigned int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (uint32_t(1) << icut)) {
            if (track1.sign() * track2.sign() < 0) {
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
            } else {
              if (track1.sign() > 0) {
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
              }
            }
          } // end if (filter bits)
        }   // end for (cuts)
      }     // end for (track2)
    }       // end for (track1)
  }

  // barrel-barrel and muon-muon event mixing
  template <int TPairType, uint32_t TEventFillMap, typename TEvents, typename TTracks>
  void runSameSide(TEvents& events, TTracks const& tracks)
  {
    events.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    GroupSlicer slicerTracks(events, tracksTuple);
    for (auto& [event1, event2] : selfCombinations(hashBin, 100, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);
      auto it1 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event1.index()) {
          it1 = slice;
          break;
        }
      }
      auto tracks1 = std::get<TTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&events);
      auto it2 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event2.index()) {
          it2 = slice;
          break;
        }
      }
      auto tracks2 = std::get<TTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&events);

      runMixedPairing<TPairType>(tracks1, tracks2);
    } // end event loop
  }

  // barrel-muon event mixing
  template <uint32_t TEventFillMap, typename TEvents, typename TTracks, typename TMuons>
  void runBarrelMuon(TEvents& events, TTracks const& tracks, TMuons const& muons)
  {
    events.bindExternalIndices(&muons);
    auto tracksTuple = std::make_tuple(tracks);
    auto muonsTuple = std::make_tuple(muons);
    GroupSlicer slicerTracks(events, tracksTuple);
    GroupSlicer slicerMuons(events, muonsTuple);
    for (auto& [event1, event2] : selfCombinations(hashBin, 100, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);
      auto it1 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event1.index()) {
          it1 = slice;
          break;
        }
      }
      auto tracks1 = std::get<TTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&events);
      auto it2 = slicerMuons.begin();
      for (auto& slice : slicerMuons) {
        if (slice.groupingElement().index() == event2.index()) {
          it2 = slice;
          break;
        }
      }
      auto muons2 = std::get<TMuons>(it2.associatedTables());
      muons2.bindExternalIndices(&events);

      runMixedPairing<pairTypeEMu>(tracks1, muons2);
    } // end event loop
  }

  void processBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMap>(events, tracks);
  }
  void processMuonSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMap>(events, muons);
  }
  void processBarrelMuonSkimmed(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runBarrelMuon<gkEventFillMap>(events, tracks, muons);
  }
  void processBarrelVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    runSameSide<pairTypeEE, gkEventFillMapWithQvector>(events, tracks);
  }
  void processMuonVnSkimmed(soa::Filtered<MyEventsHashSelectedQvector>& events, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    runSameSide<pairTypeMuMu, gkEventFillMapWithQvector>(events, muons);
  }
  // TODO: This is a dummy process function for the case when the user does not want to run any of the process functions (no event mixing)
  //    If there is no process function enabled, the workflow hangs
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventMixing, processBarrelSkimmed, "Run barrel-barrel mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonSkimmed, "Run muon-muon mixing on skimmed muons", false);
  PROCESS_SWITCH(AnalysisEventMixing, processBarrelMuonSkimmed, "Run barrel-muon mixing on skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisEventMixing, processBarrelVnSkimmed, "Run barrel-barrel vn mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processMuonVnSkimmed, "Run muon-muon vn mixing on skimmed tracks", false);
  PROCESS_SWITCH(AnalysisEventMixing, processDummy, "Dummy function", false);
};

struct AnalysisSameEventPairing {

  Produces<aod::Dileptons> dileptonList;
  Produces<aod::DileptonsExtra> dileptonExtraList;
  Produces<aod::DileptonFlow> dileptonFlowList;

  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // added by lupz begin
  Configurable<float> magneticField{"d_bz", -999, "magnetic field"};                                                   // added by lupz
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};                                 // added by lupz
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"}; // added by lupz
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  float d_bz;     // added by lupz
  int mRunNumber; // added by lupz
  o2::base::MatLayerCylSet* lut = nullptr;
  KFParticle KFPV;
  bool flagKF = false;
  // added by lupz begin

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::dqanalysisflags::isMuonSelected > 0;

  HistogramManager* fHistMan;

  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  // TODO: provide as Configurable the list and names of the cuts which should be used in pairing
  uint32_t fTwoTrackFilterMask = 0;
  uint32_t fTwoMuonFilterMask = 0;
  std::vector<std::vector<TString>> fTrackHistNames;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<std::vector<TString>> fTrackMuonHistNames;

  void init(o2::framework::InitContext& context)
  {
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // added by lupz begin
    if (context.mOptions.get<bool>("processJpsiToEESkimmedKFParticle")) // added by lupz
    {
      LOGF(info, "It is running Jpsi->ee with KF Particle"); // added by lupz
      mRunNumber = 0;                                        // added by lupz
      d_bz = 0;                                              // added by lupz
      flagKF = true;
    }
    // added by lupz end

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    if (context.mOptions.get<bool>("processJpsiToEESkimmed") || context.mOptions.get<bool>("processVnJpsiToEESkimmed") || context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processJpsiToEESkimmedKFParticle")) {
      TString cutNames = fConfigTrackCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fTrackHistNames.push_back(names);
          fTwoTrackFilterMask |= (uint32_t(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processJpsiToMuMuSkimmed") || context.mOptions.get<bool>("processJpsiToMuMuVertexingSkimmed") || context.mOptions.get<bool>("processVnJpsiToMuMuSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
          std::vector<TString> names = {
            Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fMuonHistNames.push_back(names);
          fTwoMuonFilterMask |= (uint32_t(1) << icut);
        }
      }
    }
    if (context.mOptions.get<bool>("processElectronMuonSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNamesBarrel = fConfigTrackCuts.value;
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) { // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) {
            std::vector<TString> names = {
              Form("PairsEleMuSEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            fTrackMuonHistNames.push_back(names);
            fTwoTrackFilterMask |= (uint32_t(1) << icut);
            fTwoMuonFilterMask |= (uint32_t(1) << icut);
          }
        }
      }
    }

    // Usage example of ccdb
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // ccdb->setCreatedNotAfter(nolaterthan.value);
    ccdb->setFatalWhenNull(false);                        // added by lupz
    if (!o2::base::GeometryManager::isGeometryLoaded()) { // added by lupz
      ccdb->get<TGeoManager>(geoPath);                    // added by lupz
    }                                                     // added by lupz

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  // added by lupz begin
  float getMagneticField(uint64_t timestamp)
  {
    o2::base::MatLayerCylSet* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    static o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, timestamp);
    if (grpo != nullptr) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    } else {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return 0;
    }
    LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
    float bz = std::lround(5.f * grpo->getL3Current() / 30000.f); // in kG
    LOGF(info, "magnetig field = %f kG", bz);
    return bz;
  }

  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (abs(magneticField) > 99) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = magneticField;
      }
      mRunNumber = lRunNumber;
    }
  }
  // added by lupz end

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2, typename TAmbiTracks>
  void runSameEventPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2, TAmbiTracks const& ambiTracksMid)
  {

    unsigned int ncuts = fTrackHistNames.size();
    std::vector<std::vector<TString>> histNames = fTrackHistNames;
    if constexpr (TPairType == pairTypeMuMu) {
      ncuts = fMuonHistNames.size();
      histNames = fMuonHistNames;
    }
    if constexpr (TPairType == pairTypeEMu) {
      ncuts = fTrackMuonHistNames.size();
      histNames = fTrackMuonHistNames;
    }

    uint32_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    uint32_t dileptonMcDecision = 0; // placeholder, copy of the dqEfficiency.cxx one
    dileptonList.reserve(1);
    dileptonExtraList.reserve(1);

    // added by lupz begin
    // PV information
    float PVNContributors = -999.;
    int PVNDF = -999;
    float PVParameters[8] = {-999.};
    float PVCovariance[36] = {-999.};
    // tracks information
    int trk0IsAmbiguous = 0;
    int trk1IsAmbiguous = 0;
    char trk0Charge;
    char trk1Charge;
    float trk0Parameters[8] = {-999.};
    float trk1Parameters[8] = {-999.};
    // only Geometrical fitting
    float pairMassKFGeo = -999.;        // added by lupz
    float pairMassErrKFGeo = -999.;     // added by lupz
    float pairChi2OverNDFKFGeo = -999.; // added by lupz
    int pairNDFKFGeo = -999;
    float pairDecayLengthKFGeo = -999.;
    float pairDecayLengthOverErrKFGeo = -999.;
    float pairPseudoProperDecayTimeKFGeo = -999.;
    float pairPseudoProperDecayLengthManuallyGeo = -999.;
    float dcaTrk0KFGeo = -999.;
    float dcaTrk1KFGeo = -999.;
    float dcaTrksMaxKFGeo = -999.;
    float dcaBetweenTrksKFGeo = -999.;
    float pairParametersGeo[8] = {-999.};
    float pairCovarianceGeo[36] = {-999.};
    // Geometrical fitting and topological constraint
    float pairMassKFGeoTop = -999.;        // added by lupz
    float pairMassErrKFGeoTop = -999.;     // added by lupz
    float pairChi2OverNDFKFGeoTop = -999.; // added by lupz
    int pairNDFKFGeoTop = -999;
    float pairDecayLengthKFGeoTop = -999.;
    float pairDecayLengthOverErrKFGeoTop = -999.;
    float pairPseudoProperDecayTimeKFGeoTop = -999.;
    float pairPseudoProperDecayLengthManuallyGeoTop = -999.;
    float dcaTrk0KFGeoTop = -999.;
    float dcaTrk1KFGeoTop = -999.;
    float dcaBetweenTrksKFGeoTop = -999.;
    float dcaTrksMaxKFGeoTop = -999.;
    float pairParametersGeoTop[8] = {-999.};
    float pairCovarianceGeoTop[36] = {-999.};
    // Geometrical fitting and mass constraint
    float pairMassKFGeoMass = -999.;        // added by lupz
    float pairMassErrKFGeoMass = -999.;     // added by lupz
    float pairChi2OverNDFKFGeoMass = -999.; // added by lupz
    int pairNDFKFGeoMass = -999;
    float pairDecayLengthKFGeoMass = -999.;
    float pairDecayLengthOverErrKFGeoMass = -999.;
    float pairPseudoProperDecayTimeKFGeoMass = -999.;
    float pairPseudoProperDecayLengthManuallyGeoMass = -999.;
    float dcaTrk0KFGeoMass = -999.;
    float dcaTrk1KFGeoMass = -999.;
    float dcaTrksMaxKFGeoMass = -999.;
    float dcaBetweenTrksKFGeoMass = -999.;
    float pairParametersGeoMass[8] = {-999.};
    float pairCovarianceGeoMass[36] = {-999.};
    // Geometrical fitting + mass and topological constraints
    float pairMassKFGeoMassTop = -999.;        // added by lupz
    float pairMassErrKFGeoMassTop = -999.;     // added by lupz
    float pairChi2OverNDFKFGeoMassTop = -999.; // added by lupz
    int pairNDFKFGeoMassTop = -999;
    float pairDecayLengthKFGeoMassTop = -999.;
    float pairDecayLengthOverErrKFGeoMassTop = -999.;
    float pairPseudoProperDecayTimeKFGeoMassTop = -999.;
    float pairPseudoProperDecayLengthManuallyGeoMassTop = -999.;
    float dcaTrk0KFGeoMassTop = -999.;
    float dcaTrk1KFGeoMassTop = -999.;
    float dcaTrksMaxKFGeoMassTop = -999.;
    float dcaBetweenTrksKFGeoMassTop = -999.;
    float pairParametersGeoMassTop[8] = {-999.};
    float pairCovarianceGeoMassTop[36] = {-99.};

    if (flagKF) {
      CheckAndUpdate(event.runNumber(), event.timestamp());
      LOG(info) << "~~~~~~~~~~~~~~~~Run: " << event.runNumber() << " with magnetic field of " << d_bz << " kZG";
      KFParticle::SetField(d_bz);

      if constexpr ((TPairType == pairTypeEE) && ((TEventFillMap & VarManager::ObjTypes::ReducedEventVtxCov) > 0)) {
        KFPVertex kfpVertex;
        kfpVertex.SetXYZ(event.posX(), event.posY(), event.posZ());
        //kfpVertex.SetCovarianceMatrix(event.covXX(), event.covXY(), event.covYY(), event.covXZ(), event.covYZ(), event.covZZ()); // this is the right one, but the covariance YY and XZ were swaped in run3 data, MC and run2 converted
        kfpVertex.SetCovarianceMatrix(event.covXX(), event.covXY(), event.covXZ(), event.covYY(), event.covYZ(), event.covZZ());
        kfpVertex.SetChi2(event.chi2());
        kfpVertex.SetNDF(2*event.numContrib() - 3); // added on 2022/11/16
        kfpVertex.SetNContributors(event.numContrib());
        PVNContributors = kfpVertex.GetNContributors();
        PVNDF = kfpVertex.GetNDF();
        KFPV = KFParticle(kfpVertex);
      }
    }
    // added by lupz end

    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      // added by lupz begin
      if (flagKF) {
        /*
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::AmbiTrack) > 0) {
          trk0IsAmbiguous = 0;
          trk1IsAmbiguous = 0;
          for (auto& ambiTrackMid : ambiTracksMid) {
            auto ambiTrack = ambiTrackMid.template track_as<MyBarrelTracks>();
            auto ambiTrackWithCov = ambiTrackMid.template track_as<MyBarrelTracksWithCov>();
            if (ambiTrack.globalIndex() == t1.globalIndex() || ambiTrackWithCov.globalIndex() == t1.globalIndex()) {
              trk0IsAmbiguous = 1;
              break;
            }
          }
          for (auto& ambiTrackMid : ambiTracksMid) {
            auto ambiTrack = ambiTrackMid.template track_as<MyBarrelTracks>();
            auto ambiTrackWithCov = ambiTrackMid.template track_as<MyBarrelTracksWithCov>();
            if (ambiTrack.globalIndex() == t2.globalIndex() || ambiTrackWithCov.globalIndex() == t2.globalIndex()) {
              trk1IsAmbiguous = 1;
              break;
            }
          }
        }
        */
       trk0IsAmbiguous = t1.isAmbiguous();
       trk1IsAmbiguous = t2.isAmbiguous();
        if constexpr ((TPairType == pairTypeEE) && (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0) {
          // dauther0;
          std::array<float, 3> trk0ParPos;
          std::array<float, 3> trk0ParMom;
          std::array<float, 21> trk0Cov;

          auto track0ParCov = getTrackParCov(t1);
          track0ParCov.getXYZGlo(trk0ParPos);
          track0ParCov.getPxPyPzGlo(trk0ParMom);
          track0ParCov.getCovXYZPxPyPzGlo(trk0Cov);

          float trk0ParKF[6] = {trk0ParPos[0], trk0ParPos[1], trk0ParPos[2], trk0ParMom[0], trk0ParMom[1], trk0ParMom[2]};
          float trk0CovKF[21];
          for (int i = 0; i < 21; i++) {
            trk0CovKF[i] = trk0Cov[i];
          }

          KFPTrack kfpTrack0;
          kfpTrack0.SetParameters(trk0ParKF);
          kfpTrack0.SetCovarianceMatrix(trk0CovKF);
          kfpTrack0.SetCharge(t1.sign());
          kfpTrack0.SetNDF(1); //added on 2022/11/16
          //kfpTrack0.SetChi2(...); //added on 2022/11/16, do not have this information in AO2D

          int pdgTrack0 = 0;
          if (t1.sign() < 0)
            pdgTrack0 = 11; // e-
          if (t1.sign() > 0)
            pdgTrack0 = -11; // e+
          KFParticle trk0KF(kfpTrack0, pdgTrack0);

          // dauther1;
          std::array<float, 3> trk1ParPos;
          std::array<float, 3> trk1ParMom;
          std::array<float, 21> trk1Cov;

          auto track1ParCov = getTrackParCov(t2);
          track1ParCov.getXYZGlo(trk1ParPos);
          track1ParCov.getPxPyPzGlo(trk1ParMom);
          track1ParCov.getCovXYZPxPyPzGlo(trk1Cov);

          float trk1ParKF[6] = {trk1ParPos[0], trk1ParPos[1], trk1ParPos[2], trk1ParMom[0], trk1ParMom[1], trk1ParMom[2]};
          float trk1CovKF[21];
          for (int i = 0; i < 21; i++) {
            trk1CovKF[i] = trk1Cov[i];
          }

          KFPTrack kfpTrack1;
          kfpTrack1.SetParameters(trk1ParKF);
          kfpTrack1.SetCovarianceMatrix(trk1CovKF);
          kfpTrack1.SetCharge(t2.sign());
          kfpTrack1.SetNDF(1); //added on 2022/11/16
          //kfpTrack1.SetChi2(...); //added on 2022/11/16, do not have this information in AO2D

          int pdgTrack1 = 0;
          if (t2.sign() < 0)
            pdgTrack1 = 11; // e-
          if (t2.sign() > 0)
            pdgTrack1 = -11; // e+
          KFParticle trk1KF(kfpTrack1, pdgTrack1);

          // get PV information
          for (int i = 0; i < 8; i++) {
            PVParameters[i] = KFPV.GetParameter(i);
          }
          for (int i = 0; i < 36; i++) {
            PVCovariance[i] = KFPV.GetCovariance(i);
          }

          // reconstruct Jpsi via KF
          KFParticle JpsiGeo;
          JpsiGeo.SetConstructMethod(2);
          JpsiGeo.AddDaughter(trk0KF);
          JpsiGeo.AddDaughter(trk1KF);

          KFParticle JpsiGeoTop = JpsiGeo;
          // trk0KF.SubtractFromVertex(KFPV);
          // trk1KF.SubtractFromVertex(KFPV);
          // KFPV += JpsiGeoTop;
          JpsiGeoTop.SetProductionVertex(KFPV);
          // trk0KF.SetProductionVertex(JpsiGeoTop);
          // trk1KF.SetProductionVertex(JpsiGeoTop);

          KFParticle JpsiGeoMass = JpsiGeo;
          JpsiGeoMass.SetNonlinearMassConstraint(RecoDecay::getMassPDG(443));

          KFParticle JpsiGeoMassTop = JpsiGeo;
          // trk0KF.SubtractFromVertex(KFPV);
          // trk1KF.SubtractFromVertex(KFPV);
          // KFPV += JpsiGeoMassTop;
          JpsiGeoMassTop.SetProductionVertex(KFPV);
          // trk0KF.SetProductionVertex(JpsiGeoMassTop);
          // trk1KF.SetProductionVertex(JpsiGeoMassTop);
          JpsiGeoMassTop.SetNonlinearMassConstraint(RecoDecay::getMassPDG(443));

          // get daughters information
          trk0Charge = trk0KF.GetQ();
          trk1Charge = trk1KF.GetQ();
          for (int i = 0; i < 8; i++) {
            trk0Parameters[i] = trk0KF.GetParameter(i);
            trk1Parameters[i] = trk1KF.GetParameter(i);
          }

          pairChi2OverNDFKFGeo = JpsiGeo.GetChi2() / JpsiGeo.GetNDF();
          pairNDFKFGeo = JpsiGeo.GetNDF();
          pairMassKFGeo = JpsiGeo.GetMass();
          pairMassErrKFGeo = JpsiGeo.GetErrMass();
          pairDecayLengthKFGeo = JpsiGeo.GetDecayLength();
          pairDecayLengthOverErrKFGeo = JpsiGeo.GetDecayLength() / JpsiGeo.GetErrDecayLength();
          pairPseudoProperDecayTimeKFGeo = JpsiGeo.GetPseudoProperDecayTime(KFPV, pairMassKFGeo);
          pairPseudoProperDecayLengthManuallyGeo = JpsiGeo.GetDecayLengthXY() * (JpsiGeo.GetMass() / JpsiGeo.GetPt());
          // dcaTrk0KFGeo = trk0KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          // dcaTrk1KFGeo = trk1KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          dcaTrk0KFGeo = trk0KF.GetDistanceFromVertex(KFPV);
          dcaTrk1KFGeo = trk1KF.GetDistanceFromVertex(KFPV);
          if (dcaTrk0KFGeo > dcaTrk1KFGeo)
            dcaTrksMaxKFGeo = dcaTrk0KFGeo;
          else
            dcaTrksMaxKFGeo = dcaTrk1KFGeo;
          dcaBetweenTrksKFGeo = trk0KF.GetDistanceFromParticle(trk1KF);
          for (int i = 0; i < 8; i++) {
            pairParametersGeo[i] = JpsiGeo.GetParameter(i);
          }
          for (int i = 0; i < 36; i++) {
            pairCovarianceGeo[i] = JpsiGeo.GetCovariance(i);
          }

          pairChi2OverNDFKFGeoTop = JpsiGeoTop.GetChi2() / JpsiGeoTop.GetNDF();
          pairNDFKFGeoTop = JpsiGeoTop.GetNDF();
          pairMassKFGeoTop = JpsiGeoTop.GetMass();
          pairMassErrKFGeoTop = JpsiGeoTop.GetErrMass();
          pairDecayLengthKFGeoTop = JpsiGeoTop.GetDecayLength();
          pairDecayLengthOverErrKFGeoTop = JpsiGeoTop.GetDecayLength() / JpsiGeoTop.GetErrDecayLength();
          pairPseudoProperDecayTimeKFGeoTop = JpsiGeoTop.GetPseudoProperDecayTime(KFPV, pairMassKFGeoTop);
          pairPseudoProperDecayLengthManuallyGeoTop = JpsiGeoTop.GetDecayLengthXY() * (JpsiGeoTop.GetMass() / JpsiGeoTop.GetPt());
          // dcaTrk0KFGeoTop = trk0KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          // dcaTrk1KFGeoTop = trk1KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          dcaTrk0KFGeoTop = trk0KF.GetDistanceFromVertex(KFPV);
          dcaTrk1KFGeoTop = trk1KF.GetDistanceFromVertex(KFPV);
          if (dcaTrk0KFGeoTop > dcaTrk1KFGeoTop)
            dcaTrksMaxKFGeoTop = dcaTrk0KFGeoTop;
          else
            dcaTrksMaxKFGeoTop = dcaTrk1KFGeoTop;
          dcaBetweenTrksKFGeoTop = trk0KF.GetDistanceFromParticle(trk1KF);
          for (int i = 0; i < 8; i++) {
            pairParametersGeoTop[i] = JpsiGeoTop.GetParameter(i);
          }
          for (int i = 0; i < 36; i++) {
            pairCovarianceGeoTop[i] = JpsiGeoTop.GetCovariance(i);
          }

          pairChi2OverNDFKFGeoMass = JpsiGeoMass.GetChi2() / JpsiGeoMass.GetNDF();
          pairNDFKFGeoMass = JpsiGeoMass.GetNDF();
          pairMassKFGeoMass = JpsiGeoMass.GetMass();
          pairMassErrKFGeoMass = JpsiGeoMass.GetErrMass();
          pairDecayLengthKFGeoMass = JpsiGeoMass.GetDecayLength();
          pairDecayLengthOverErrKFGeoMass = JpsiGeoMass.GetDecayLength() / JpsiGeoMass.GetErrDecayLength();
          pairPseudoProperDecayTimeKFGeoMass = JpsiGeoMass.GetPseudoProperDecayTime(KFPV, pairMassKFGeoMass);
          pairPseudoProperDecayLengthManuallyGeoMass = JpsiGeoMass.GetDecayLengthXY() * (JpsiGeoMass.GetMass() / JpsiGeoMass.GetPt());
          // dcaTrk0KFGeoMass = trk0KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          // dcaTrk1KFGeoMass = trk1KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          dcaTrk0KFGeoMass = trk0KF.GetDistanceFromVertex(KFPV);
          dcaTrk1KFGeoMass = trk1KF.GetDistanceFromVertex(KFPV);
          if (dcaTrk0KFGeoMass > dcaTrk1KFGeoMass)
            dcaTrksMaxKFGeoMass = dcaTrk0KFGeoMass;
          else
            dcaTrksMaxKFGeoMass = dcaTrk1KFGeoMass;
          dcaBetweenTrksKFGeoMass = trk0KF.GetDistanceFromParticle(trk1KF);
          for (int i = 0; i < 8; i++) {
            pairParametersGeoMass[i] = JpsiGeoMass.GetParameter(i);
          }
          for (int i = 0; i < 36; i++) {
            pairCovarianceGeoMass[i] = JpsiGeoMass.GetCovariance(i);
          }

          pairChi2OverNDFKFGeoMassTop = JpsiGeoMassTop.GetChi2() / JpsiGeoMassTop.GetNDF();
          pairNDFKFGeoMassTop = JpsiGeoMassTop.GetNDF();
          pairMassKFGeoMassTop = JpsiGeoMassTop.GetMass();
          pairMassErrKFGeoMassTop = JpsiGeoMassTop.GetErrMass();
          pairDecayLengthKFGeoMassTop = JpsiGeoMassTop.GetDecayLength();
          pairDecayLengthOverErrKFGeoMassTop = JpsiGeoMassTop.GetDecayLength() / JpsiGeoMassTop.GetErrDecayLength();
          pairPseudoProperDecayTimeKFGeoMassTop = JpsiGeoMassTop.GetPseudoProperDecayTime(KFPV, pairMassKFGeoMassTop);
          pairPseudoProperDecayLengthManuallyGeoMassTop = JpsiGeoMassTop.GetDecayLengthXY() * (JpsiGeoMassTop.GetMass() / JpsiGeoMassTop.GetPt());
          // dcaTrk0KFGeoMassTop = trk0KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          // dcaTrk1KFGeoMassTop = trk1KF.GetR(); // distance to the origin of the coordinate system {0,0,0}
          dcaTrk0KFGeoMassTop = trk0KF.GetDistanceFromVertex(KFPV);
          dcaTrk1KFGeoMassTop = trk1KF.GetDistanceFromVertex(KFPV);
          if (dcaTrk0KFGeoMassTop > dcaTrk1KFGeoMassTop)
            dcaTrksMaxKFGeoMassTop = dcaTrk0KFGeoMassTop;
          else
            dcaTrksMaxKFGeoMassTop = dcaTrk1KFGeoMassTop;
          dcaBetweenTrksKFGeoMassTop = trk0KF.GetDistanceFromParticle(trk1KF);
          for (int i = 0; i < 8; i++) {
            pairParametersGeoMassTop[i] = JpsiGeoMassTop.GetParameter(i);
          }
          for (int i = 0; i < 36; i++) {
            pairCovarianceGeoMassTop[i] = JpsiGeoMassTop.GetCovariance(i);
          }
        }
      }
      // added by lupz end
      if constexpr (TPairType == VarManager::kJpsiToEE) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isBarrelSelected()) & fTwoTrackFilterMask;
      }
      if constexpr (TPairType == VarManager::kJpsiToMuMu) {
        twoTrackFilter = uint32_t(t1.isMuonSelected()) & uint32_t(t2.isMuonSelected()) & fTwoMuonFilterMask;
      }
      if constexpr (TPairType == VarManager::kElectronMuon) {
        twoTrackFilter = uint32_t(t1.isBarrelSelected()) & uint32_t(t2.isMuonSelected()) & fTwoTrackFilterMask;
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
      // TODO: FillPair functions need to provide a template argument to discriminate between cases when cov matrix is available or not
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      if constexpr ((TPairType == pairTypeEE) || (TPairType == pairTypeMuMu)) { // call this just for ee or mumu pairs
        VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2);
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<TPairType>(t1, t2);
        }
      }

      // TODO: provide the type of pair to the dilepton table (e.g. ee, mumu, emu...)
      dileptonFilterMap = twoTrackFilter;
      // added by lupz begin
      if constexpr ((TPairType == pairTypeEE) && (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelPID) > 0) {
        dileptonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision, t1.pt(),
                     t1.eta(), t1.phi(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(), t2.pt(), t2.eta(), t2.phi(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                     trk0IsAmbiguous, trk1IsAmbiguous, trk0Parameters,trk1Parameters,// trk0Charge,trk1Charge,
                     pairMassKFGeo, pairChi2OverNDFKFGeo, pairNDFKFGeo, pairDecayLengthKFGeo, pairDecayLengthOverErrKFGeo, pairPseudoProperDecayTimeKFGeo, pairPseudoProperDecayLengthManuallyGeo, dcaTrk0KFGeo, dcaTrk1KFGeo, dcaTrksMaxKFGeo, dcaBetweenTrksKFGeo, pairParametersGeo, pairCovarianceGeo,
                     pairMassKFGeoTop, pairChi2OverNDFKFGeoTop, pairNDFKFGeoTop, pairDecayLengthKFGeoTop, pairDecayLengthOverErrKFGeoTop, pairPseudoProperDecayTimeKFGeoTop, pairPseudoProperDecayLengthManuallyGeoTop, dcaTrk0KFGeoTop, dcaTrk1KFGeoTop, dcaTrksMaxKFGeoTop, dcaBetweenTrksKFGeoTop, pairParametersGeoTop, pairCovarianceGeoTop,
                     pairMassKFGeoMass, pairChi2OverNDFKFGeoMass, pairNDFKFGeoMass, pairDecayLengthKFGeoMass, pairDecayLengthOverErrKFGeoMass, pairPseudoProperDecayTimeKFGeoMass, pairPseudoProperDecayLengthManuallyGeoMass, dcaTrk0KFGeoMass, dcaTrk1KFGeoMass, dcaTrksMaxKFGeoMass, dcaBetweenTrksKFGeoMass, pairParametersGeoMass, pairCovarianceGeoMass,
                     pairMassKFGeoMassTop, pairChi2OverNDFKFGeoMassTop, pairNDFKFGeoMassTop, pairDecayLengthKFGeoMassTop, pairDecayLengthOverErrKFGeoMassTop, pairPseudoProperDecayTimeKFGeoMassTop, pairPseudoProperDecayLengthManuallyGeoMassTop, dcaTrk0KFGeoMassTop, dcaTrk1KFGeoMassTop, dcaTrksMaxKFGeoMassTop, dcaBetweenTrksKFGeoMassTop, pairParametersGeoMassTop, pairCovarianceGeoMassTop,
                     PVParameters, PVCovariance, PVNContributors, PVNDF,
                     VarManager::fgValues[VarManager::kMCVtxX], VarManager::fgValues[VarManager::kMCVtxY], VarManager::fgValues[VarManager::kMCVtxZ],
                     VarManager::fgValues[VarManager::kMCPdgCode], VarManager::fgValues[VarManager::kMCVx], VarManager::fgValues[VarManager::kMCVy], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCPx], VarManager::fgValues[VarManager::kMCPy], VarManager::fgValues[VarManager::kMCPz],
                     event.globalIndex(), t1.reducedeventId(), t2.reducedeventId(),0,0,0); // added by lupz
      }
      // added by lupz end

      // dileptonList(event, VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap, dileptonMcDecision);

      constexpr bool muonHasCov = ((TTrackFillMap & VarManager::ObjTypes::MuonCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedMuonCov) > 0);
      if constexpr ((TPairType == pairTypeMuMu) && muonHasCov) {
        dileptonExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
      }

      if constexpr (eventHasQvector) {
        dileptonFlowList(VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi]);
      }

      for (unsigned int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint32_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
            }
          }
        } // end if (filter bits)
      }   // end for (cuts)
    }     // end loop over pairs
  }

  void processJpsiToEESkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks, nullptr);
  }
  void processJpsiToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, muons, nullptr);
  }
  void processJpsiToMuMuVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelectedWithCov> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, muons, muons, nullptr);
  }
  void processVnJpsiToEESkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToEE, gkEventFillMapWithCovQvector, gkTrackFillMap>(event, tracks, tracks, nullptr);
  }
  void processVnJpsiToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelectedQvector>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCovQvector>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToMuMu, gkEventFillMapWithCovQvector, gkMuonFillMap>(event, muons, muons, nullptr);
  }
  void processElectronMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons, nullptr);
  }
  void processAllSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToEE, gkEventFillMap, gkTrackFillMap>(event, tracks, tracks, nullptr);
    runSameEventPairing<VarManager::kJpsiToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, muons, nullptr);
    runSameEventPairing<VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons, nullptr);
  }
  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }
  // added by lupz begin
  void processJpsiToEESkimmedKFParticle(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelectedWithCov> const& tracks, aod::AmbiguousTracksMid const& ambiTracksMid)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMapWithCov>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kJpsiToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, tracks, ambiTracksMid);
  }
  // added by lupz end

  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToEESkimmed, "Run electron-electron pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToMuMuSkimmed, "Run muon-muon pairing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToMuMuVertexingSkimmed, "Run muon-muon pairing and vertexing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnJpsiToEESkimmed, "Run electron-electron pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processVnJpsiToMuMuSkimmed, "Run muon-muon pairing, with skimmed tracks for vn", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processElectronMuonSkimmed, "Run electron-muon pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processJpsiToEESkimmedKFParticle, "Run Jpsi reconstruction using KFparticle tools, with skimmed tracks", false); // added by lupz
};

struct AnalysisDileptonHadron {
  //
  // This task combines dilepton candidates with a track and could be used for example
  //  in analyses with the dilepton as one of the decay products of a higher mass resonance (e.g. B -> Jpsi + K)
  //    or in dilepton + hadron correlations, etc.
  //
  //  The barrel and muon track filtering tasks can produce multiple parallel decisions, which are used to produce
  //   dileptons which inherit the logical intersection of the track level decisions (see the AnalysisSameEventPairing task).
  //  This can be used also in the dilepton-hadron correlation analysis. However, in this model of the task, we use all the dileptons produced in the
  //    lepton pairing task to combine them with the hadrons selected by the barrel track selection.
  //  To be modified/adapted if new requirements appear

  OutputObj<THashList> fOutputList{"output"};
  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "", "Comma separated list of barrel track cuts"};
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::mass > 2.92f && aod::reducedpair::mass < 3.16f && aod::reducedpair::sign == 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  HistogramManager* fHistMan;

  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expressions
  int fNHadronCutBit;

  void init(o2::framework::InitContext& context)
  {
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // TODO: Create separate histogram directories for each selection used in the creation of the dileptons
    // TODO: Implement possibly multiple selections for the associated track ?
    if (context.mOptions.get<bool>("processSkimmed")) {
      DefineHistograms(fHistMan, "DileptonsSelected;DileptonHadronInvMass;DileptonHadronCorrelation"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString configCutNamesStr = fConfigTrackCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries();
    } else {
      fNHadronCutBit = 0;
    }
  }

  // Template function to run pair - hadron combinations
  template <uint32_t TEventFillMap, typename TEvent, typename TTracks>
  void runDileptonHadron(TEvent const& event, TTracks const& tracks, soa::Filtered<aod::Dileptons> const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    // loop once over dileptons for QA purposes
    for (auto dilepton : dileptons) {
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      fHistMan->FillHistClass("DileptonsSelected", fValuesDilepton);
      // loop over hadrons
      for (auto& hadron : tracks) {
        // TODO: Replace this with a Filter expression
        if (!(uint32_t(hadron.isBarrelSelected()) & (uint32_t(1) << fNHadronCutBit))) {
          continue;
        }
        // TODO: Check whether this hadron is one of the dilepton daughters!
        VarManager::FillDileptonHadron(dilepton, hadron, fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronInvMass", fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronCorrelation", fValuesHadron);
      }
    }
  }

  void processSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelected const& tracks, soa::Filtered<aod::Dileptons> const& dileptons)
  {
    runDileptonHadron<gkEventFillMap>(event, tracks, dileptons);
  }
  // TODO: Add process functions which use cov matrices for secondary vertexing (e.g. B->Jpsi + K)
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonHadron, processSkimmed, "Run dilepton-hadron pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonHadron, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisEventMixing>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonHadron>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,cent,muon");
    }

    if (classStr.Contains("Track")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid");
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
      }
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel", "vertexing-barrel,flow-barrel");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_dimuon", "vertexing-forward,flow-dimuon");
      }
      if (classStr.Contains("EleMu")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_electronmuon");
      }
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "kine");
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }
  } // end loop over histogram classes
}
