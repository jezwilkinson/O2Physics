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

///
/// \file   qaSimpleITSTPCMatch.cxx
/// \author Jeremy Wilkinson jeremy.wilkinson@cern.ch
/// \brief  Task to analyse data for ITS-TPC matching efficiency with simple cluster cuts
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

// ROOT includes
#include "TPDGCode.h"
#include "TEfficiency.h"
#include "THashList.h"

using namespace o2::framework;

struct qaSimpleITSTPCMatch {


  Configurable<int> minITSCl{"minITSCl", 6, "minimum # ITS clusters for matching"};
  Configurable<int> minTPCCl{"minTPCCl", 100, "minimum # TPC clusters for matching"};
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 10.f}, "Pt binning"};
  ConfigurableAxis phiBins{"phiBins", {54, 0.f, 6.284f}, "Phi binning"};
  ConfigurableAxis etaBins{"etaBins", {200, -1.f, 1.f}, "Eta binning"};

  HistogramRegistry histos{"Histos",{},OutputObjHandlingPolicy::AnalysisObject};

  //2D efficiencies
  TEfficiency* effITSTPCMatchingVsPtVsPhi = nullptr;
  TEfficiency* effTPCITSMatchingVsPtVsPhi = nullptr;

  TEfficiency* effITSTPCMatchingVsEtaVsPhi = nullptr;
  TEfficiency* effTPCITSMatchingVsEtaVsPhi = nullptr;

  void init(InitContext&)
  {
    const AxisSpec axisPt{ptBins,"#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{phiBins, "#phi"};
    const AxisSpec axisEta{etaBins, "#eta"};

    histos.add("Data/trackLength", "Track length;Track length (cm)", kTH1F, {{2000, -1000, 1000}});
    
    
    histos.add("Data/pos/etaphi/its_tpc", "ITS-TPC Positive ", kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc", "ITS-TPC Negative ", kTH2D, {axisEta, axisPhi});

    histos.add("Data/pos/ptphi/its_tpc", "ITS-TPC Positive " , kTH2D, {axisPt, axisPhi});
    histos.add("Data/neg/ptphi/its_tpc", "ITS-TPC Negative " , kTH2D, {axisPt, axisPhi});

    histos.add("Data/pos/etaphi/its", "ITS Positive ", kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its", "ITS Negative ", kTH2D, {axisEta, axisPhi});

    histos.add("Data/pos/ptphi/its", "ITS Positive " , kTH2D, {axisPt, axisPhi});
    histos.add("Data/neg/ptphi/its", "ITS Negative " , kTH2D, {axisPt, axisPhi});

    histos.add("Data/pos/etaphi/tpc", "TPC Positive ", kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/tpc", "TPC Negative ", kTH2D, {axisEta, axisPhi});

    histos.add("Data/pos/ptphi/tpc", "TPC Positive " , kTH2D, {axisPt, axisPhi});
    histos.add("Data/neg/ptphi/tpc", "TPC Negative " , kTH2D, {axisPt, axisPhi});

    histos.add("Data/pos/etaphi/all", "ITS-TPC Positive ", kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/all", "ITS-TPC Negative ", kTH2D, {axisEta, axisPhi});

    histos.add("Data/pos/ptphi/all", "ITS-TPC Positive " , kTH2D, {axisPt, axisPhi});
    histos.add("Data/neg/ptphi/all", "ITS-TPC Negative " , kTH2D, {axisPt, axisPhi});



    histos.add("eventSelection", "Event Selection", kTH1D, {{10, 0.5, 10.5, "Selection"}});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel. (no ev. sel)");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

    auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY, TEfficiency*& eff) {
      TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
      TAxis* axisY = histos.get<TH1>(templateHistoY)->GetYaxis();
      if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
        eff = new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray());
      } else {
        eff = new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax());
      }
      listEfficiencyData->Add(eff);
    };

  }

  template <bool doFillHistograms, typename CollType>
  bool isCollisionSelected(const CollType& collision)
  {
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 1);
    }
    if (applyEvSel == 1 && !collision.sel7()) {
      return false;
    } else if (applyEvSel == 2 && !collision.sel8()) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 2);
    }
    if (collision.numContrib() < nMinNumberOfContributors) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 3);
    }
    if ((collision.posZ() < vertexZMin || collision.posZ() > vertexZMax)) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 4);
    }
    return true;
  }

  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
                   TrackCandidates const& tracks)
  {

    if (!isCollisionSelected<true>(collision)) {
      return;
    }
    bool passedITS = false;
    bool passedTPC = false;

    for (const auto& track : tracks) {
      //if (!isTrackSelected<false>(track, HIST("Data/trackSelection"))) {
      //  continue;
      //}
      passedITS = false;
      passedTPC = false;

      if (track.pt() < 0.4) {  // skip low-pT tracks
        continue;
      }
      if (track.itsNCls() >= minITSCl)
        passedITS = true;
      if (track.tpcNClsFound() >= minTPCCl)
        passedTPC = true;
      histos.fill(HIST("Data/trackLength"), track.length());

      if (track.sign() > 0) {
        histos.fill(HIST("Data/pos/ptphi/all"), track.pt(),track.phi());
        histos.fill(HIST("Data/pos/etaphi/all"), track.eta(), track.phi());
      } else {
        histos.fill(HIST("Data/neg/ptphi/all"), track.pt(),track.phi());
        histos.fill(HIST("Data/neg/etaphi/all"), track.eta(), track.phi());
      }

      if (passedITS) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/ptphi/its"), track.pt(),track.phi());
          histos.fill(HIST("Data/pos/etaphi/its"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/ptphi/its"), track.pt(),track.phi());
          histos.fill(HIST("Data/neg/etaphi/its"), track.eta(), track.phi());
        }
      }

      if (passedTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/ptphi/tpc"), track.pt(),track.phi());
          histos.fill(HIST("Data/pos/etaphi/tpc"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/ptphi/tpc"), track.pt(),track.phi());
          histos.fill(HIST("Data/neg/etaphi/tpc"), track.eta(), track.phi());
        }
      }

      if (passedITS && passedTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/ptphi/its_tpc"), track.pt(),track.phi());
          histos.fill(HIST("Data/pos/etaphi/its_tpc"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/ptphi/its_tpc"), track.pt(),track.phi());
          histos.fill(HIST("Data/neg/etaphi/its_tpc"), track.eta(), track.phi());
        }
      }
      
    }
  }
}; // struct qaSimpleITSTPCMatch

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaSimpleITSTPCMatch>(cfgc)};
}
