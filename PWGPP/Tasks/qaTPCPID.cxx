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
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>, GSI Darmstadt

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{         //runtime customisation goes here
  };
  std::swap(workflowOptions, options);
     
}

#include "Framework/runDataProcessing.h"

using namespace o2::dataformats;


struct QaTPCPID {
  //Configurables
  Configurable<float> cutTOF{"cutTOF", 3.f, "TOF nsigma cut for TPC-TOF PID"};
  Configurable<int> pBins{"pBins", 400, "Number of momentum bins"};
  Configurable<float> pMin{"pMin", 0.f, "Lower limit in momentum"};
  Configurable<float> pMax{"pMax", 20.f, "Upper limit in momentum"};
  Configurable<int> nBinsNSigma{"nBinsNSigma",200, "Number of bins for TPC nSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Lower limit for TPC nSigma"};
  Configurable<float> maxNSigma{"maxNSigma",  10.f, "Upper limit for TPC nSigma"};
  
  HistogramRegistry hists{"HistogramsTPCPIDQA"};
  
  static constexpr int Np = 9;
  static constexpr std::string_view hnsigmaTPC[Np] = {"nsigmaTPC/El", "nsigmaTPC/Mu", "nsigmaTPC/Pi",
                                                   "nsigmaTPC/Ka", "nsigmaTPC/Pr", "nsigmaTPC/De",
                                                   "nsigmaTPC/Tr", "nsigmaTPC/He", "nsigmaTPC/Al"};
  static constexpr std::string_view hnsigmaTPCTOF[Np] = {"nsigmaTPCTOF/El", "nsigmaTPCTOF/Mu", "nsigmaTPCTOF/Pi",
                                                   "nsigmaTPCTOF/Ka", "nsigmaTPCTOF/Pr", "nsigmaTPCTOF/De",
                                                   "nsigmaTPCTOF/Tr", "nsigmaTPCTOF/He", "nsigmaTPCTOF/Al"};
  static constexpr const char* partName[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};

  template <uint8_t i>
  void addParticleHistos()
  {
    AxisSpec pAxis{pBins, pMin, pMax, "#it{p} [GeV/#it{c}]"};
    AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, "TPC n_{#sigma}"};
    
    hists.add(hnsigmaTPC[i].data(), Form("TPC signal (%s) without TOF cut",partName[i]), kTH2F, {pAxis, nSigmaAxis}); 
    hists.add(hnsigmaTPCTOF[i].data(), Form("TPC signal (%s) after %.2f#sigma TOF cut",partName[i],cutTOF), kTH2F, {pAxis, nSigmaAxis}); 
    
  }//addParticleHistos
  

  void init(InitContext&)
  {
    for (int i = 0; i < Np; i++) {
      addParticleHistos<i>;
    }//for
  }//init
  
  template <uint8_t i, typename T>
  void fillParticleHistos(const T& t, const float mom, const float tofNSigma, const float tpcNSigma)
  {
    hists.fill(HIST(hnsigmaTPC[i]),mom, tpcNSigma);
    if (tofNSigma < cutTOF) hists.fill(HIST(hnsigmaTPCTOF[i]),mom, tpcNSigma);
    
  }//fillParticleHistos
  
  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra,
                                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                          aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                                          aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                                                          aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                                                          aod::TrackSelection> const& tracks)
  {
    for (auto t : tracks) {
      const float mom = t.tpcInnerParam();
      fillParticleHistos<0>(t, mom, t.tofNSigmaEl(), t.tpcNSigmaEl());
      fillParticleHistos<1>(t, mom, t.tofNSigmaMu(), t.tpcNSigmaMu());
      fillParticleHistos<2>(t, mom, t.tofNSigmaPi(), t.tpcNSigmaPi());
      fillParticleHistos<3>(t, mom, t.tofNSigmaKa(), t.tpcNSigmaKa());
      fillParticleHistos<4>(t, mom, t.tofNSigmaPr(), t.tpcNSigmaPr());
      fillParticleHistos<5>(t, mom, t.tofNSigmaDe(), t.tpcNSigmaDe());
      fillParticleHistos<6>(t, mom, t.tofNSigmaTr(), t.tpcNSigmaTr());
      fillParticleHistos<7>(t, mom, t.tofNSigmaHe(), t.tpcNSigmaHe());
      fillParticleHistos<8>(t, mom, t.tofNSigmaAl(), t.tpcNSigmaAl());
      
    }//for
  }//process
}//struct QaTPCPID


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  Workflowpec w;
  w.push_back(adaptAnalysisTask<QaTPCPID>(cfgc));
  return w;
}