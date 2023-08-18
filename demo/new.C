#pragma once 

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

//#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TVector3.h"

using namespace ana;

//need a variable

const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
	auto const& trk = slc->reco.pfp.at(i).trk;
        if(trk.bestplane == -1) continue;
	////////        if(trk.bestplane == -1 || slc->reco.pfp.at(i).trackScore < 0.5) continue; // rm the track Score because seems like the clearcosmic slices don't go through track/shower fit 

	// First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

	// We require that the distance of the track from the slice is less than
	// 10 cm and that the parent of the track has been marked as the primary.
	const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);

        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

	const bool Contained = ( !isnan(trk.end.x) &&
				 ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
				   ( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
				 !isnan(trk.end.y) &&
				 ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
				 !isnan(trk.end.z) &&
				 ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );
	              
	const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
	const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );

	if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
	  {
	    Longest = trk.len;
	    PTrackInd = i;

	  }
      }
    return PTrackInd;
  });

const Var cosThNuMI( [](const caf::SRSliceProxy* slc) {
    float cosTh(-5.f);

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
	TVector3 muDir( trk.dir.x, trk.dir.y, trk.dir.z );
	muDir = muDir.Unit();
	double thNuMI = muDir.Angle(NuDirection_NuMI);
	cosTh = TMath::Cos(thNuMI);
      }
    return cosTh;
  });

///what cuts to apply
const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
	       ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	     !isnan(slc->vertex.y) &&
	     ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	     !isnan(slc->vertex.z) &&
	     ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 ) );
  });

void new(){
  const std::string dataFile = "";//add file
  
  SpectrumLoader loader(dataFile);

  //double kPOTnumi = 6E20;

  const Binning kTimeBinning = Binning::Simple(100,-10.,10.);
}
