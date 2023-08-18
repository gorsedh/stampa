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
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202208.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202208.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TVector3.h"
//#include ""
#include <string>
#include <fstream>

using namespace ana;

//////////////////////////////////////////////////////
/// Here you can define several variables,         ///
//////////////////////////////////////////////////////

/////++++++++ Slice Variables ++++++++/////

const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
	auto const& trk = slc->reco.pfp.at(i).trk;
	if(trk.bestplane < 0 || trk.bestplane > 3) continue; // used to be == -1
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

const Var kmuon_momentum([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
        if(Contained) p = trk.rangeP.p_muon;
        else p = trk.mcsP.fwdP_muon;
      }
    return p;
  });

const Var kmomentum([](const caf::SRSliceProxy* slc) -> float {
    float p = sqrt(pow(slc->truth.momentum.x,2) + pow(slc->truth.momentum.y,2) + pow(slc->truth.momentum.z,2));

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
        if(Contained) p = trk.rangeP.p_muon;
        else p = trk.mcsP.fwdP_muon;
      }
    return p;
  });


// TVector3 for NuMI beam: Here we'll get the cosine of the angle between the track direction and the NuMI beam direction
// First attempt
TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

// See get_beam_Dir.C
TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

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

const Var kHasTruthMatch(
			 [](const caf::SRSliceProxy *slc) -> double
			 {
			   return ( slc->truth.index != -1);
			 }
			 );

const Var kTruthEnergy(
		       [](const caf::SRSliceProxy *slc) -> double
		       {
			 return ( kHasTruthMatch(slc) ? (float)slc->truth.E : -5.f );
		       });

const Var slcvtxX( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.x;
  });

const Var slcvtxY( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.y;
  });

const Var slcvtxZ( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.z;
  });

const Var TrueslcvtxX( [](const caf::SRSliceProxy* slc) {
    return slc->truth.position.x;
  });

const Var TrueslcvtxY( [](const caf::SRSliceProxy* slc) {
    return slc->truth.position.y;
  });

const Var TrueslcvtxZ( [](const caf::SRSliceProxy* slc) {
    return slc->truth.position.z;
  });

const Var kEndTrkX([](const caf::SRSliceProxy *slc) -> double{
    float endx(-100000.f);
    for(auto& pfp: slc -> reco.pfp) {
      const auto& trk = pfp.trk;
      endx = trk.end.x;
    }
    return endx;
  });

const Var kEndTrkY([](const caf::SRSliceProxy *slc) -> double{
    float endy(-100000.f);
    for(auto& pfp: slc -> reco.pfp) {
      const auto& trk = pfp.trk;
      endy = trk.end.y;
    }
    return endy;
  });

const Var kEndTrkZ([](const caf::SRSliceProxy *slc) -> double{
    float endz(-100000.f);
    for(auto& pfp: slc -> reco.pfp) {
      const auto& trk = pfp.trk;
      endz = trk.end.z;
    }
    return endz;
  });

const Var recophi( [](const caf::SRSliceProxy* slc) {
    return slc->reco.pfp[0].trk.phi;
  });

const Var longestCR( [](const caf::SRSliceProxy* slc) {
    return slc->nuid.crlongtrkdiry;
  });


/////++++++++ Spill Variables ++++++++/////

const SpillMultiVar spillvarTopCRTHitTime([](const caf::SRSpillProxy *sr) ///Here is the CRT time
					  {
					    std::vector<double> rets;
					    for(const auto& hit : sr->crt_hits){
					      if(hit.plane>=30&&hit.plane<=39){
						double this_crttime = -(sr->hdr.ismc ? hit.t0 : hit.t1);
						rets.push_back(this_crttime);
					      }
					    }
					    return rets;
					  });

const SpillMultiVar spillvarOpFlashTime([](const caf::SRSpillProxy *sr)  ///Here is the Optical flash time [for data works weell BUT for MC should substract sr->hdr.triggerinfo.trigger_within_gate]
					{
					  std::vector<double> rets;
					  for(const auto& opflash : sr->opflashes){
					    double this_oft = opflash.firsttime - sr->hdr.triggerinfo.trigger_within_gate;
					    rets.push_back(this_oft);
					  }
					  return rets;
					});

const SpillMultiVar NuMag([](const caf::SRSpillProxy *sr)
			  {
			    std::vector<double> rets;
			    for(const auto& nu : sr->mc.nu){
			      double mag = sqrt(pow(nu.position.x,2) + pow(nu.position.y,2) + pow(nu.position.z,2));
			      rets.push_back(mag);
			    }
			    return rets;
			  });

const SpillMultiVar NuMagMinX([](const caf::SRSpillProxy *sr)
			      {
				std::vector<double> rets;
				double mag_min = DBL_MAX;
				double minX = 0;
				for(const auto& nu : sr->mc.nu){
				  double mag = sqrt(pow(nu.position.x,2) + pow(nu.position.y,2) + pow(nu.position.z,2));
				  if (mag_min > mag) {
				    mag_min = mag;
				    minX = nu.position.x;
				  }
				}
				rets.push_back(minX);
				return rets;
			      });

const SpillMultiVar NuMagMinY([](const caf::SRSpillProxy *sr)
			      {
				std::vector<double> rets;
				double mag_min = DBL_MAX;
				double minY = 0;
				for(const auto& nu : sr->mc.nu){
				  double mag = sqrt(pow(nu.position.x,2) + pow(nu.position.y,2) + pow(nu.position.z,2));
				  if (mag_min > mag) {
				    mag_min = mag;
				    minY = nu.position.y;
				  }

				}
				rets.push_back(minY);
				return rets;
			      });

const SpillMultiVar NuMagMinZ([](const caf::SRSpillProxy *sr)
			      {
				std::vector<double> rets;
				double mag_min = DBL_MAX;
				double minZ = 0;
				for(const auto& nu : sr->mc.nu){
				  double mag = sqrt(pow(nu.position.x,2) + pow(nu.position.y,2) + pow(nu.position.z,2));
				  if (mag_min > mag) {
				    mag_min = mag;
				    minZ = nu.position.z;
				  }

				}
				rets.push_back(minZ);
				return rets;
			      });

const SpillMultiVar NuE([](const caf::SRSpillProxy *sr)
			{
			  std::vector<double> rets;
			  for(const auto& nu : sr->mc.nu){
			    double this_oft = nu.E;
			    rets.push_back(this_oft);
			  }
			  return rets;
			});

const SpillMultiVar NuX([](const caf::SRSpillProxy *sr)
			{
			  std::vector<double> rets;
			  for(const auto& nu : sr->mc.nu){
			    double this_oft = nu.position.x;
			    rets.push_back(this_oft);
			  }
			  return rets;
			});

const SpillMultiVar NuY([](const caf::SRSpillProxy *sr)
			{
			  std::vector<double> rets;
			  for(const auto& nu : sr->mc.nu){
			    double this_oft = nu.position.y;
			    rets.push_back(this_oft);
			  }
			  return rets;
			});

const SpillMultiVar NuZ([](const caf::SRSpillProxy *sr)
			{
			  std::vector<double> rets;
			  for(const auto& nu : sr->mc.nu){
			    double this_oft = nu.position.z;
			    rets.push_back(this_oft);
			  }
			  return rets;
			});



//////////////////////////////////////////////////////
/// here you can define multiple selection cuts,   ///
//////////////////////////////////////////////////////

const Cut kmuon_contained([](const caf::SRSliceProxy* slc) {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
        auto const& trk = slc->reco.pfp.at(i).trk;
        if(trk.bestplane < 0 || trk.bestplane > 3) continue; // used to be == -1

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

        if (Chi2Proton > 60 && Chi2Muon < 30) {
          return true;
        }
      }
    return false;
  });

const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
	       ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	     !isnan(slc->vertex.y) &&
	     ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	     !isnan(slc->vertex.z) &&
	     ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 )     );
  });

const Cut kXFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
	       ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) );
  });

const Cut kYFiducial([](const caf::SRSliceProxy* slc) {
    return (
	    !isnan(slc->vertex.y) &&
	    ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 )  );
  });

const Cut kZFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.z) &&
	     ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 ) );
  });

const Cut kXYFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
	       ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	     !isnan(slc->vertex.y) &&
             ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 )
	     );
    //&& (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    // && kmuon_contained(slc));
  });

const Cut kXZFiducial([](const caf::SRSliceProxy* slc) {
    return( !isnan(slc->vertex.x) &&
	    ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
	      ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	    !isnan(slc->vertex.z) &&
            ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 )
            );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    // && kmuon_contained(slc));
  });

const Cut kYZFiducial([](const caf::SRSliceProxy* slc) {
    return (!isnan(slc->vertex.y) &&
	    ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	    !isnan(slc->vertex.z) &&
            ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 )
            );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kXYNotFiducial([](const caf::SRSliceProxy* slc) {
    return (( (!isnan(slc->vertex.x) &&
	       ( ( slc->vertex.x > -61.94 - 25 && slc->vertex.x < 61.94 + 25 ) ||
		 ( slc->vertex.x < -358.49 + 25) || ( slc->vertex.x > 358.49 - 25 ) )) ||
	      (!isnan(slc->vertex.y) &&
	       ( slc->vertex.y < -181.86 + 25 || slc->vertex.y > 134.96 - 25 )))
            );
    //&& (slc->nuid.crlongtrkdiry > -0.2));
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kXZNotFiducial([](const caf::SRSliceProxy* slc) {
    return(( (!isnan(slc->vertex.x) &&
	      ( ( slc->vertex.x > -61.94 - 25 && slc->vertex.x < 61.94 + 25 ) ||
		( slc->vertex.x < -358.49 + 25) || ( slc->vertex.x > 358.49 - 25 ) )) ||
	     (!isnan(slc->vertex.z) &&
	      ( slc->vertex.z < -894.951 + 30 || slc->vertex.z > 894.951 - 50 )))
           );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kYZNotFiducial([](const caf::SRSliceProxy* slc) {
    return (((!isnan(slc->vertex.y) &&
	      ( slc->vertex.y < -181.86 + 25 || slc->vertex.y > 134.96 - 25 )) ||
	     (!isnan(slc->vertex.z) &&
	      ( slc->vertex.z < -894.951 + 30 || slc->vertex.z > 894.951 - 50 )))
            );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kXYVol([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->truth.position.x) &&
	     ( ( slc->truth.position.x < -71.94 - 5 && slc->truth.position.x > -358.49 + 5 ) ||
	       ( slc->truth.position.x > 71.94 + 5 && slc->truth.position.x < 358.49 - 5 ) ) &&
	     !isnan(slc->truth.position.y) &&
             ( slc->truth.position.y > -181.86 + 5 && slc->truth.position.y < 134.96 - 5 )
	     // );
	     //&& (slc->nuid.crlongtrkdiry > -0.2)
             && kHasNu(slc) );
    // && kmuon_contained(slc));
  });

const Cut kXZVol([](const caf::SRSliceProxy* slc) {
    return( !isnan(slc->truth.position.x) &&
	    ( ( slc->truth.position.x < -71.94 - 5 && slc->truth.position.x > -358.49 + 5 ) ||
	      ( slc->truth.position.x > 71.94 + 5 && slc->truth.position.x < 358.49 - 5 ) ) &&
	    !isnan(slc->truth.position.z) &&
            ( slc->truth.position.z > -894.951 + 5 && slc->truth.position.z < 894.951 - 5 )
            );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    // && kmuon_contained(slc));
  });

const Cut kYZVol([](const caf::SRSliceProxy* slc) {
    return (!isnan(slc->truth.position.y) &&
	    ( slc->truth.position.y > -181.86 + 5 && slc->truth.position.y < 134.96 - 5 ) &&
	    !isnan(slc->truth.position.z) &&
            ( slc->truth.position.z > -894.951 + 3 && slc->truth.position.z < 894.951 - 5 )
	    // );
	    //&&  (slc->nuid.crlongtrkdiry > -0.2)
            && kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kXYNotVol([](const caf::SRSliceProxy* slc) {
    return (( (!isnan(slc->truth.position.x) &&
	       ( ( slc->truth.position.x > -71.94 - 5 && slc->truth.position.x < 71.94 + 5 ) ||
		 ( slc->truth.position.x < -358.49 + 5) || ( slc->truth.position.x > 358.49 - 5 ) )) ||
	      (!isnan(slc->truth.position.y) &&
	       ( slc->truth.position.y < -181.86 + 5 || slc->truth.position.y > 134.96 - 5 )))
            //&& (slc->nuid.crlongtrkdiry > -0.2));
            && kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kXZNotVol([](const caf::SRSliceProxy* slc) {
    return(( (!isnan(slc->truth.position.x) &&
	      ( ( slc->truth.position.x > -71.94 - 5 && slc->truth.position.x < 71.94 + 5 ) ||
		( slc->truth.position.x < -358.49 + 5) || ( slc->truth.position.x > 358.49 - 5 ) )) ||
	     (!isnan(slc->truth.position.z) &&
	      ( slc->truth.position.z < -894.951 + 5 || slc->truth.position.z > 894.951 - 5 )))
	   );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });

const Cut kYZNotVol([](const caf::SRSliceProxy* slc) {
    return (((!isnan(slc->truth.position.y) &&
	      ( slc->truth.position.y < -181.86 + 5 || slc->truth.position.y > 134.96 - 5 )) ||
	     (!isnan(slc->vertex.z) &&
	      ( slc->truth.position.z < -894.951 + 30 || slc->truth.position.z > 894.951 - 50 )))
	    );
    //&&  (slc->nuid.crlongtrkdiry > -0.2)
    //&& kHasNu(slc) );
    //&& kmuon_contained(slc));
  });


const Cut kSlcTrkDiry([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryNu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kHasNu(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryNoNu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && !kHasNu(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryCosth([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && cosThNuMI(slc) > .95) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryCosthMu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && cosThNuMI(slc) > .95 && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryMu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kMunNu([](const caf::SRSliceProxy* slc) -> double {
    if (kHasNu(slc) && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kNotRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x > -61.94 - 25 && slc->vertex.x < -358.49 + 25 ) ||
	       ( slc->vertex.x < 61.94 + 25 && slc->vertex.x > 358.49 - 25 ) ) &&
	     !isnan(slc->vertex.y) &&
	     ( slc->vertex.y < -181.86 + 25 && slc->vertex.y > 134.96 - 25 ) &&
	     !isnan(slc->vertex.z) &&
	     ( slc->vertex.z < -894.951 + 30 && slc->vertex.z > 894.951 - 50 ) );
  });


const SpillCut kcrtHit_cut([](const caf::SRSpillProxy* sr){
    for (auto const& crtHit: sr->crt_hits){
      auto thistime = crtHit.time;
      if (thistime > -0.1 && thistime < 9.7)
	return false;
    }
    return true;
  }
  );

const SpillCut OpFlash_cut([](const caf::SRSpillProxy* sr){
    for (const auto& opflash : sr->opflashes){
      auto thistime =  opflash.firsttime;
      if (thistime > -0.1 && thistime < 9.7) {
	return true;
      } else {
	return false;
      }
    }
    return false;
  }
  );

const SpillCut OpFlashandCRT_cut([](const caf::SRSpillProxy* sr){
    for (const auto& opflash : sr->opflashes){
      auto thistime =  opflash.firsttime;
      if (thistime < -0.1 || thistime > 9.7) {
	return false;
      }
    }
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=30&&hit.plane<=39){
	double this_crttime = -(sr->hdr.ismc ? hit.t0 : hit.t1);
	if(this_crttime > -.01 && this_crttime < 9.7) {
	  return true;
	}
      }
    }
    return false;
  }
  );

const SpillCut NuXYFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return ( !isnan(nu.position.x) &&
	       ( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
		 ( nu.position.x > 61.94 + 25 && nu.position.x < 358.49 - 25 ) ) &&
	       !isnan(nu.position.y) &&
	       ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) );
    }
    return false;
  });

const SpillCut NuXZFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return( !isnan(nu.position.x) &&
	      ( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
		( nu.position.x > 61.94 + 25 && nu.position.x < 358.49 - 25 ) ) &&
              !isnan(nu.position.z) &&
	      ( nu.position.z > -894.951 + 30 && nu.position.z < 894.951 - 50 ) );
    }
    return false;
  });

const SpillCut NuYZFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return (!isnan(nu.position.y) &&
	      ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) &&
	      !isnan(nu.position.z) &&
	      ( nu.position.z > -894.951 + 30 && nu.position.z < 894.951 - 50 ));
    }
    return false;
  });

const SpillCut NuXYNotFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return (( (!isnan(nu.position.x) &&
		 ( ( nu.position.x > -61.94 - 25 && nu.position.x < 61.94 + 25 ) ||
		   ( nu.position.x < -358.49 + 25) || ( nu.position.x > 358.49 - 25 ) )) ||
		(!isnan(nu.position.y) &&
		 ( nu.position.y < -181.86 + 25 || nu.position.y > 134.96 - 25 ))));
    }
    return false;
  });

const SpillCut NuXZNotFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return(( (nu.position.x) &&
	       ( ( nu.position.x > -61.94 - 25 && nu.position.x < 61.94 + 25 ) ||
		 ( nu.position.x < -358.49 + 25) || ( nu.position.x > 358.49 - 25 ) )) ||
	     (!isnan(nu.position.z) &&
	      ( nu.position.z < -894.951 + 30 || nu.position.z > 894.951 - 50 )));
    }
    return false;
  });

const SpillCut NuYZNotFiducial([](const caf::SRSpillProxy* sr) {
    for(const auto& nu : sr->mc.nu){
      return (((!isnan(nu.position.y) &&
		( nu.position.y < -181.86 + 25 || nu.position.y > 134.96 - 25 )) ||
	       (!isnan(nu.position.z) &&
		( nu.position.z < -894.951 + 30 || nu.position.z > 894.951 - 50 ))));
    }
    return false;
  });

void My2DselectRock(){

  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader.t
  const std::string fnewdata_MAJORITY = "/pnfs/icarus/scratch/users/gputnam/DMCP2023G/majority-3t1p/57386892_[3,4,5,7,8,9][0,1,3]*/data*Prescaled*.root";
  const std::string fdata_MAJORITY = "/pnfs/sbn/data/sbn_fd/poms_production/data/Run2/reconstructed/icaruscode_v09_72_00_03/numimajority/flatcaf_unblind/0*/[0,1,2]*/data*.flat.caf*.root"; // Data ICARUS run 2
  const std::string fMC_MAJORITY = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_dirt_plus_cosmics/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/0*/[0,1,2,3,4]*/detsim*.flat.caf*.root";
  const std::string fMC_NOROCK = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_Nu_Phase1_sample/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/0*/[0,1]*/detsim*.flat.caf*.root";
  // Source of events
  const std::string fMC_COSMICS = "/pnfs/icarus/scratch/users/aheggest/crt/mc/NuMI_intime_cosmics/v09_75_01/caf/NuMI_intime_cosmics_246files_merged.root";
  const std::string fdata_OFFBEAM = "/pnfs/icarus/scratch/users/aheggest/crt/data/run_9726/offbeamNuMImajority/v09_75_01/caf/data_offbeamNuMI_crtpmt_Blind_193files_merged.OKTOLOOK.flat.caf.root";
  SpectrumLoader loader(fMC_NOROCK);

  double kPOTnumi = 9.63683E17 * 10E2; // IN data we don't need this. This value is from prescaled data.
  double knorm_norock = .406594;// 2.6759;// .406594; First for rock, second for nonrock

  // Binning
  const Binning kCosThNumiBinning = Binning::Simple(50,-1.,1.);
  const Binning kLengthBinning    = Binning::Simple(40,0.,400);
  const Binning kTimeBinning      = Binning::Simple(60,-25.,25.);
  const Binning kEnergyBinning = Binning::Simple(30,0,12);


  const Binning kPositionX        = Binning::Simple(120,-500.,500.);
  const Binning kPositionY        = Binning::Simple(120,-500.,500.);
  const Binning kPositionZ        = Binning::Simple(120,-1000.,1000.);

  const Binning kMagBinning        = Binning::Simple(100,0,3000);

  // vector<Binning> xy_bins = {kPositionX,kPositionY};
  // vector<Binning> xz_bins = {kPositionX,kPositionZ};
  // vector<Binning> yz_bins = {kPositionY,kPositionZ};

  // Hist Axis
  const HistAxis axtkl   ("trklength",      kLengthBinning,         kLongestTrackLength);
  const HistAxis axCRdir   ("trklength",      kCosThNumiBinning,         longestCR);
  const HistAxis axcos   ("costh",          kCosThNumiBinning,      cosThNuMI);
  const HistAxis axphi   ("costh",          kCosThNumiBinning,      recophi);
  const HistAxis axmumom   ("mom",          kEnergyBinning,      kmuon_momentum);
  const HistAxis axallmom   ("mom",          kEnergyBinning,      kmomentum);
  const HistAxis axenergy   ("energy",          kEnergyBinning,      kTruthEnergy);

  const HistAxis axvtxX  ("vtxx",           kPositionX,             slcvtxX); //kSlcVtxX, kTruthVtxX
  const HistAxis axvtxY  ("vtxy",           kPositionY,             slcvtxY); //kSlcVtxY
  const HistAxis axvtxZ  ("vtxz",           kPositionZ,             slcvtxZ);

  const HistAxis axendX ("endX",            kPositionX,             kEndTrkX);
  const HistAxis axendY ("endY",            kPositionY,             kEndTrkY);
  const HistAxis axendZ ("endZ",            kPositionZ,             kEndTrkZ);

  // const HistAxis axvtxXY  ("vtxxy",           xy_bins,             xy_vars); //Currently unused
  // const HistAxis axvtxXZ  ("vtxxz",           xz_bins,             xz_vars); //Currently unused
  // const HistAxis axvtxYZ  ("vtxyz",           yz_bins,             yz_vars); //Currently unused


  // Let's define the spectrums

  //For Slices
  Spectrum scos       (loader, axcos, kNoCut);
  Spectrum scos_nofid (loader, axcos, !kRFiducial);
  Spectrum scos_Fcuts (loader, axcos,   OpFlash_cut, kNoCut);
  Spectrum scos_FCcuts (loader, axcos,  OpFlash_cut, kNoCut);
  Spectrum scos_Dirycuts (loader, axcos,kNoSpillCut,  kSlcTrkDiry);

  Spectrum sphi       (loader, axphi, kNoCut);

  Spectrum stkl       (loader, axtkl, kNoCut);
  Spectrum stkl_nofid (loader, axtkl, !kRFiducial);
  Spectrum stkl_muon (loader, axtkl, OpFlash_cut, kNoCut);
  Spectrum stkl_neutrino (loader, axtkl, kMunNu);
  Spectrum stkl_Diry (loader, axtkl, kSlcTrkDiry);

  Spectrum sCRdir       (loader, axCRdir, kNoCut);

  Spectrum sMuMom (loader, axmumom, kNoCut);
  Spectrum sNuMom (loader, axmumom, kHasNu);
  Spectrum sAllMom (loader, axallmom, kNoCut);
  Spectrum sDiryMom (loader, axallmom, kSlcTrkDiry);
  Spectrum sEnergy (loader, axenergy, kNoCut);
  Spectrum sNuEnergy (loader, axenergy, kHasNu);
  Spectrum sNuMuEnergy (loader, axenergy, kMunNu);

  Spectrum vtxX       (loader, axvtxX, kNoCut);
  Spectrum vtxX_NFV   (loader, axvtxX, !kXFiducial);
  Spectrum vtxX_FV    (loader, axvtxX, kXFiducial);
  Spectrum vtxY       (loader, axvtxY, kNoCut);
  Spectrum vtxZ       (loader, axvtxZ, kNoCut);

  Spectrum sendX       (loader, axendX, kSlcTrkDiryNu);
  Spectrum sendY       (loader, axendY, kMunNu);
  Spectrum sendZ       (loader, axendZ, kNoCut);


  //spectra for 2D histogram
  Spectrum vtxXY       (loader, axvtxX, axvtxY, kNoSpillCut, kNoCut, kNoShift, kUnweighted);
  Spectrum vtxXZ       (loader, axvtxX, axvtxZ, kNoSpillCut, kNoCut, kNoShift, kUnweighted);
  Spectrum vtxYZ       (loader, axvtxY, axvtxZ, kNoSpillCut, kNoCut, kNoShift, kUnweighted);

  Spectrum vtxXY_NFV       (loader, axvtxX, axvtxY, kNoSpillCut, kXYNotFiducial, kNoShift, kUnweighted);
  Spectrum vtxXZ_NFV       (loader, axvtxX, axvtxZ, kNoSpillCut, kXZNotFiducial, kNoShift, kUnweighted);
  Spectrum vtxYZ_NFV       (loader, axvtxY, axvtxZ, kNoSpillCut, kYZNotFiducial, kNoShift, kUnweighted);

  Spectrum vtxXY_FV       (loader, axvtxX, axvtxY, kNoSpillCut, kXYFiducial, kNoShift, kUnweighted);
  Spectrum vtxXZ_FV       (loader, axvtxX, axvtxZ, kNoSpillCut, kXZFiducial, kNoShift, kUnweighted);
  Spectrum vtxYZ_FV       (loader, axvtxY, axvtxZ, kNoSpillCut, kYZFiducial, kNoShift, kUnweighted);

  Spectrum sendXY       (loader, axendX, axendY, kNoSpillCut, kSlcTrkDiryNu, kNoShift, kUnweighted);
  Spectrum sendXZ       (loader, axendX, axendZ, kNoSpillCut, kSlcTrkDiryNu, kNoShift, kUnweighted);
  Spectrum sendYZ       (loader, axendY, axendZ, kNoSpillCut, kSlcTrkDiryNu, kNoShift, kUnweighted);


  //For Spills
  Spectrum *sTopCRTHitTime = new Spectrum("TopCRTHitTime", kTimeBinning, loader, spillvarTopCRTHitTime, kNoSpillCut);

  Spectrum *sOpFlashTime   = new Spectrum("OpFlashTime", kTimeBinning, loader, spillvarOpFlashTime, kNoSpillCut);

  Spectrum *sNuMag   = new Spectrum("NuMag", kMagBinning, loader, NuMag, kNoSpillCut);

  // Spectrum *sNuEnergy = new Spectrum("NuEnergy",kEnergyBinning, loader, NuE, kNoSpillCut);

  Spectrum *sNuX   = new Spectrum("NuX", kPositionX, loader, NuX, kNoSpillCut);
  Spectrum *sNuY   = new Spectrum("NuY", kPositionY, loader, NuY, kNoSpillCut);
  Spectrum *sNuZ   = new Spectrum("NuZ", kPositionZ, loader, NuZ, kNoSpillCut);


  Spectrum *sNuXY   = new Spectrum("NuXY", loader, kPositionX, NuX, kPositionY, NuY, kNoSpillCut);
  Spectrum *sNuXZ   = new Spectrum("NuXZ", loader, kPositionX, NuX, kPositionZ, NuZ, kNoSpillCut);
  Spectrum *sNuYZ   = new Spectrum("NuYZ", loader, kPositionY, NuY, kPositionZ, NuZ, kNoSpillCut);

  Spectrum *sNuXY_FV   = new Spectrum("NuXY", loader, kPositionX, NuX, kPositionY, NuY, NuXYFiducial);
  Spectrum *sNuXZ_FV   = new Spectrum("NuXZ", loader, kPositionX, NuX, kPositionZ, NuZ, NuXZFiducial);
  Spectrum *sNuYZ_FV   = new Spectrum("NuYZ", loader, kPositionY, NuY, kPositionZ, NuZ, NuYZFiducial);

  Spectrum *sNuXY_NFV   = new Spectrum("NuXY", loader, kPositionX, NuX, kPositionY, NuY, NuXYNotFiducial);
  Spectrum *sNuXZ_NFV   = new Spectrum("NuXZ", loader, kPositionX, NuX, kPositionZ, NuZ, NuXZNotFiducial);
  Spectrum *sNuYZ_NFV   = new Spectrum("NuYZ", loader, kPositionY, NuY, kPositionZ, NuZ, NuYZNotFiducial);


  // actually make the spectra
  loader.Go();

  // For data
  // scos.OverridePOT(1);
  // scos_nofid.OverridePOT(1);
  // scos_FCcuts.OverridePOT(1);
  // scos_Fcuts.OverridePOT(1);
  // scos_Dirycuts.OverridePOT(1);
  // sphi.OverridePOT(1);
  // stkl.OverridePOT(1);
  // stkl_nofid.OverridePOT(1);
  // stkl_muon.OverridePOT(1);
  // stkl_Diry.OverridePOT(1);
  // stkl_neutrino.OverridePOT(1);
  // sCRdir.OverridePOT(1);
  // sAllMom.OverridePOT(1);
  // sDiryMom.OverridePOT(1);
  // sMuMom.OverridePOT(1);
  // sNuMom.OverridePOT(1);
  // sEnergy.OverridePOT(1);
  // sNuEnergy.OverridePOT(1);
  // sNuMuEnergy.OverridePOT(1);

  // vtxX.OverridePOT(1);
  // vtxX_NFV.OverridePOT(1);
  // vtxX_FV.OverridePOT(1);
  // vtxY.OverridePOT(1);
  // vtxZ.OverridePOT(1);
  // sendX.OverridePOT(1);
  // sendY.OverridePOT(1);
  // sendZ.OverridePOT(1);

  // vtxXY.OverridePOT(1);
  // vtxXZ.OverridePOT(1);
  // vtxYZ.OverridePOT(1);
  // vtxXY_NFV.OverridePOT(1);
  // vtxXZ_NFV.OverridePOT(1);
  // vtxYZ_NFV.OverridePOT(1);
  // vtxXY_FV.OverridePOT(1);
  // vtxXZ_FV.OverridePOT(1);
  // vtxYZ_FV.OverridePOT(1);
  // sendXY.OverridePOT(1);
  // sendXZ.OverridePOT(1);
  // sendYZ.OverridePOT(1);


  // sTopCRTHitTime -> OverridePOT(1);
  // sOpFlashTime -> OverridePOT(1);
  // sNuX -> OverridePOT(1);
  // sNuY -> OverridePOT(1);
  // sNuZ -> OverridePOT(1);

  // sNuXY -> OverridePOT(1);
  // sNuXZ -> OverridePOT(1);
  // sNuYZ -> OverridePOT(1);

  // sNuXY_FV -> OverridePOT(1);
  // sNuXZ_FV -> OverridePOT(1);
  // sNuYZ_FV -> OverridePOT(1);

  // sNuXY_NFV -> OverridePOT(1);
  // sNuXZ_NFV -> OverridePOT(1);
  // sNuYZ_NFV -> OverridePOT(1);

  // sNuMag -> OverridePOT(1);




  // Histograms
  TH1* hcos = scos.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hcos_nofid = scos_nofid.ToTH1(kPOTnumi/knorm_norock, kGreen+2);
  TH1* hcos_FCcuts = scos_FCcuts.ToTH1(kPOTnumi/knorm_norock, kRed);
  TH1* hcos_Fcuts = scos_Fcuts.ToTH1(kPOTnumi/knorm_norock, kBlue);
  TH1* hcos_Dirycuts = scos_Dirycuts.ToTH1(kPOTnumi/knorm_norock, kGreen);
  TH1* hphi = sphi.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* htkl = stkl.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* htkl_nofid = stkl_nofid.ToTH1(kPOTnumi/knorm_norock, kGreen+2);
  TH1* htkl_Diry = stkl_Diry.ToTH1(kPOTnumi/knorm_norock, kGreen);
  TH1* htkl_muon = stkl_muon.ToTH1(kPOTnumi/knorm_norock, kBlue);
  TH1* htkl_neutrino = stkl_neutrino.ToTH1(kPOTnumi/knorm_norock, kRed);
  TH1* hCRdir = sCRdir.ToTH1(kPOTnumi/knorm_norock, kBlack);

  TH1* hMuMom = sMuMom.ToTH1(kPOTnumi/knorm_norock, kBlue);
  TH1* hNuMom = sNuMom.ToTH1(kPOTnumi/knorm_norock, kRed);
  TH1* hAllMom = sAllMom.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hDiryMom = sDiryMom.ToTH1(kPOTnumi/knorm_norock, kGreen);
  TH1* hEnergy = sEnergy.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hNuEnergy = sNuEnergy.ToTH1(kPOTnumi/knorm_norock, kRed);
  TH1* hNuMuEnergy = sNuMuEnergy.ToTH1(kPOTnumi/knorm_norock, kBlue);

  TH1* vtxX_hist = vtxX.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* vtxX_NFV_hist = vtxX_NFV.ToTH1(kPOTnumi/knorm_norock, kGreen+2);
  TH1* vtxX_FV_hist = vtxX_FV.ToTH1(kPOTnumi/knorm_norock, kRed);
  TH1* vtxY_hist = vtxY.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* vtxZ_hist = vtxZ.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hendX = sendX.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hendY = sendY.ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1* hendZ = sendZ.ToTH1(kPOTnumi/knorm_norock, kBlack);

  TH2* vtxXY_hist = vtxXY.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxXZ_hist = vtxXZ.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxYZ_hist = vtxYZ.ToTH2(kPOTnumi/knorm_norock);

  TH2* vtxXY_NFV_hist = vtxXY_NFV.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxXZ_NFV_hist = vtxXZ_NFV.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxYZ_NFV_hist = vtxYZ_NFV.ToTH2(kPOTnumi/knorm_norock);

  TH2* vtxXY_FV_hist = vtxXY_FV.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxXZ_FV_hist = vtxXZ_FV.ToTH2(kPOTnumi/knorm_norock);
  TH2* vtxYZ_FV_hist = vtxYZ_FV.ToTH2(kPOTnumi/knorm_norock);

  TH2* hendXY = sendXY.ToTH2(kPOTnumi/knorm_norock);
  TH2* hendXZ = sendXZ.ToTH2(kPOTnumi/knorm_norock);
  TH2* hendYZ = sendYZ.ToTH2(kPOTnumi/knorm_norock);

  TH1D* hTopCRTHitTime = sTopCRTHitTime -> ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1D* hOpFlashTime = sOpFlashTime -> ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1D* hNuMag = sNuMag -> ToTH1(kPOTnumi/knorm_norock, kBlack);


  TH1D* hNuX = sNuX -> ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1D* hNuY = sNuY -> ToTH1(kPOTnumi/knorm_norock, kBlack);
  TH1D* hNuZ = sNuZ -> ToTH1(kPOTnumi/knorm_norock, kBlack);

  TH2* hNuXY = sNuXY -> ToTH2(kPOTnumi/knorm_norock);
  TH2* hNuXZ = sNuXZ -> ToTH2(kPOTnumi/knorm_norock);
  TH2* hNuYZ = sNuYZ -> ToTH2(kPOTnumi/knorm_norock);


  TH2* hNuXY_FV = sNuXY_FV -> ToTH2(kPOTnumi);
  TH2* hNuXZ_FV = sNuXZ_FV -> ToTH2(kPOTnumi);
  TH2* hNuYZ_FV = sNuYZ_FV -> ToTH2(kPOTnumi);

  TH2* hNuXY_NFV = sNuXY_NFV -> ToTH2(kPOTnumi);
  TH2* hNuXZ_NFV = sNuXZ_NFV -> ToTH2(kPOTnumi);
  TH2* hNuYZ_NFV = sNuYZ_NFV -> ToTH2(kPOTnumi);

  //////////// You can save your histograms in a root file and draw them later...
  TFile f("histos_selectRock.root", "RECREATE");

  hcos -> Write("hcos");
  hcos_nofid -> Write("hcos_nofid");
  hcos_FCcuts -> Write("hcos_FCcuts");
  hcos_Fcuts -> Write("hcos_Fcuts");
  hcos_Dirycuts -> Write("hcos_Dirycuts");
  hphi -> Write("hphi");
  htkl -> Write("htkl");
  htkl_nofid -> Write("htkl_nofid");
  htkl_muon -> Write("htkl_muon");
  htkl_Diry -> Write("htkl_Diry");
  htkl_neutrino -> Write("htkl_neutrino");
  hCRdir -> Write("hCRdir");

  vtxX_hist -> Write("VtxX");
  vtxX_NFV_hist -> Write("VtxX_NFV");
  vtxX_FV_hist -> Write("VtxX_FV");
  vtxY_hist -> Write("VtxY");
  vtxZ_hist -> Write("VtxZ");
  hendX -> Write("hendX");
  hendY -> Write("hendY");
  hendZ -> Write("hendZ");

  hTopCRTHitTime -> Write("hTopCRTHitTime");
  hOpFlashTime -> Write("hOpFlashTime");
  hNuMag -> Write("hNuMag");
  hEnergy -> Write("hEnergy");
  hNuEnergy -> Write("hNuEnergy");
  hNuMuEnergy -> Write("hNuMuEnergy");
  hMuMom -> Write("hMuMom");
  hNuMom -> Write("hNuMom");
  hAllMom -> Write("hAllMom");
  hDiryMom -> Write("hDiryMom");

  hNuX -> Write("hNuX");
  hNuY -> Write("hNuY");
  hNuZ -> Write("hNuZ");

  hNuXY -> Write("hNuXY");
  hNuXZ -> Write("hNuXZ");
  hNuYZ -> Write("hNuYZ");

  hNuXY_FV -> Write("hNuXY_FV");
  hNuXZ_FV -> Write("hNuXZ_FV");
  hNuYZ_FV -> Write("hNuYZ_FV");

  hNuXY_NFV -> Write("hNuXY_NFV");
  hNuXZ_NFV -> Write("hNuXZ_NFV");
  hNuYZ_NFV -> Write("hNuYZ_NFV");

  vtxXY_hist -> Write("VtxXY");
  vtxXZ_hist -> Write("VtxXZ");
  vtxYZ_hist -> Write("VtxYZ");

  vtxXY_NFV_hist -> Write("VtxXY_NFV");
  vtxXZ_NFV_hist -> Write("VtxXZ_NFV");
  vtxYZ_NFV_hist -> Write("VtxYZ_NFV");

  vtxXY_FV_hist -> Write("VtxXY_FV");
  vtxXZ_FV_hist -> Write("VtxXZ_FV");
  vtxYZ_FV_hist -> Write("VtxYZ_FV");

  hendXY -> Write("hendXY");
  hendXZ -> Write("hendXZ");
  hendYZ -> Write("hendYZ");

  //////////// ... or draw them now

  hcos -> SetTitle("");
  hcos -> GetYaxis() -> SetTitle("Slices");
  hcos -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcos -> GetYaxis() -> CenterTitle();
  hcos -> GetXaxis() -> CenterTitle();

  hphi -> SetTitle("");
  hphi -> GetYaxis() -> SetTitle("Slices");
  hphi -> GetXaxis() -> SetTitle("Reco cos(#phi)");
  hphi -> GetYaxis() -> CenterTitle();
  hphi -> GetXaxis() -> CenterTitle();

  hCRdir -> SetTitle("");
  hCRdir -> GetYaxis() -> SetTitle("Slices");
  hCRdir -> GetXaxis() -> SetTitle("Reco #theta_{LongestCR}");
  hCRdir -> GetYaxis() -> CenterTitle();
  hCRdir -> GetXaxis() -> CenterTitle();

  hcos_Fcuts -> SetTitle("");
  hcos_Fcuts -> GetYaxis() -> SetTitle("Slices");
  hcos_Fcuts -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcos_Fcuts -> GetYaxis() -> CenterTitle();
  hcos_Fcuts -> GetXaxis() -> CenterTitle();

  htkl -> SetTitle("");
  htkl -> GetYaxis() -> SetTitle("Slices");
  htkl -> GetXaxis() -> SetTitle("Track Length (cm)");
  htkl -> GetYaxis() -> CenterTitle();
  htkl -> GetXaxis() -> CenterTitle();

  hNuMag -> SetTitle("");
  hNuMag -> GetYaxis() -> SetTitle("Neutrinos");
  hNuMag -> GetXaxis() -> SetTitle("Distance from [0,0,0] (cm cubed)");
  hNuMag -> GetYaxis() -> CenterTitle();
  hNuMag -> GetXaxis() -> CenterTitle();

  hNuEnergy -> SetTitle("");
  hNuEnergy -> GetYaxis() -> SetTitle("Slices");
  hNuEnergy -> GetXaxis() -> SetTitle("Energy (GeV)");
  hNuEnergy -> GetYaxis() -> CenterTitle();
  hNuEnergy -> GetXaxis() -> CenterTitle();

  hNuEnergy -> SetTitle("");
  hNuEnergy -> GetYaxis() -> SetTitle("Neutrinos");
  hNuEnergy -> GetXaxis() -> SetTitle("Energy (GeV)");
  hNuEnergy -> GetYaxis() -> CenterTitle();
  hNuEnergy -> GetXaxis() -> CenterTitle();

  hMuMom -> SetTitle("");
  hMuMom -> GetYaxis() -> SetTitle("Muons");
  hMuMom -> GetXaxis() -> SetTitle("Momentum");
  hMuMom -> GetYaxis() -> CenterTitle();
  hMuMom -> GetXaxis() -> CenterTitle();

  hAllMom -> SetTitle("");
  hAllMom -> GetYaxis() -> SetTitle("Slices");
  hAllMom -> GetXaxis() -> SetTitle("Momentum");
  hAllMom -> GetYaxis() -> CenterTitle();
  hAllMom -> GetXaxis() -> CenterTitle();

  vtxX_hist -> SetTitle("");
  vtxX_hist -> GetYaxis() -> SetTitle("Slices");
  vtxX_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxX_hist -> GetYaxis() -> CenterTitle();
  vtxX_hist -> GetXaxis() -> CenterTitle();

  vtxY_hist -> SetTitle("");
  vtxY_hist -> GetYaxis() -> SetTitle("Slices");
  vtxY_hist -> GetXaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxY_hist -> GetYaxis() -> CenterTitle();
  vtxY_hist -> GetXaxis() -> CenterTitle();

  vtxZ_hist -> SetTitle("");
  vtxZ_hist -> GetYaxis() -> SetTitle("Slices");
  vtxZ_hist -> GetXaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxZ_hist -> GetYaxis() -> CenterTitle();
  vtxZ_hist -> GetXaxis() -> CenterTitle();

  hendX -> SetTitle("");
  hendX -> GetYaxis() -> SetTitle("Slices");
  hendX -> GetXaxis() -> SetTitle("X End Location (cm)");
  hendX -> GetYaxis() -> CenterTitle();
  hendX -> GetXaxis() -> CenterTitle();

  hendY -> SetTitle("");
  hendY -> GetYaxis() -> SetTitle("Slices");
  hendY -> GetXaxis() -> SetTitle("Y End Location (cm)");
  hendY -> GetYaxis() -> CenterTitle();
  hendY -> GetXaxis() -> CenterTitle();

  hendZ -> SetTitle("");
  hendZ -> GetYaxis() -> SetTitle("Slices");
  hendZ -> GetXaxis() -> SetTitle("Z End Location (cm)");
  hendZ -> GetYaxis() -> CenterTitle();
  hendZ -> GetXaxis() -> CenterTitle();

  vtxXY_hist -> SetTitle("");
  vtxXY_hist -> GetYaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxXY_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXY_hist -> GetYaxis() -> CenterTitle();
  vtxXY_hist -> GetXaxis() -> CenterTitle();

  vtxXZ_hist -> SetTitle("");
  vtxXZ_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxXZ_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXZ_hist -> GetYaxis() -> CenterTitle();
  vtxXZ_hist -> GetXaxis() -> CenterTitle();

  vtxYZ_hist -> SetTitle("");
  vtxYZ_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxYZ_hist -> GetXaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxYZ_hist -> GetYaxis() -> CenterTitle();
  vtxYZ_hist -> GetXaxis() -> CenterTitle();

  vtxXY_NFV_hist -> SetTitle("");
  vtxXY_NFV_hist -> GetYaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxXY_NFV_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXY_NFV_hist -> GetYaxis() -> CenterTitle();
  vtxXY_NFV_hist -> GetXaxis() -> CenterTitle();

  vtxXZ_NFV_hist -> SetTitle("");
  vtxXZ_NFV_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxXZ_NFV_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXZ_NFV_hist -> GetYaxis() -> CenterTitle();
  vtxXZ_NFV_hist -> GetXaxis() -> CenterTitle();

  vtxYZ_NFV_hist -> SetTitle("");
  vtxYZ_NFV_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxYZ_NFV_hist -> GetXaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxYZ_NFV_hist -> GetYaxis() -> CenterTitle();
  vtxYZ_NFV_hist -> GetXaxis() -> CenterTitle();

  vtxXY_FV_hist -> SetTitle("");
  vtxXY_FV_hist -> GetYaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxXY_FV_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXY_FV_hist -> GetYaxis() -> CenterTitle();
  vtxXY_FV_hist -> GetXaxis() -> CenterTitle();

  vtxXZ_FV_hist -> SetTitle("");
  vtxXZ_FV_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxXZ_FV_hist -> GetXaxis() -> SetTitle("X Vertex Location (cm)");
  vtxXZ_FV_hist -> GetYaxis() -> CenterTitle();
  vtxXZ_FV_hist -> GetXaxis() -> CenterTitle();

  vtxYZ_FV_hist -> SetTitle("");
  vtxYZ_FV_hist -> GetYaxis() -> SetTitle("Z Vertex Location (cm)");
  vtxYZ_FV_hist -> GetXaxis() -> SetTitle("Y Vertex Location (cm)");
  vtxYZ_FV_hist -> GetYaxis() -> CenterTitle();
  vtxYZ_FV_hist -> GetXaxis() -> CenterTitle();

  hendXY -> SetTitle("");
  hendXY -> GetYaxis() -> SetTitle("Y End Location (cm)");
  hendXY -> GetXaxis() -> SetTitle("X End Location (cm)");
  hendXY -> GetYaxis() -> CenterTitle();
  hendXY -> GetXaxis() -> CenterTitle();

  hendXZ -> SetTitle("");
  hendXZ -> GetYaxis() -> SetTitle("Z End Location (cm)");
  hendXZ -> GetXaxis() -> SetTitle("X End Location (cm)");
  hendXZ -> GetYaxis() -> CenterTitle();
  hendXZ -> GetXaxis() -> CenterTitle();

  hendYZ -> SetTitle("");
  hendYZ -> GetYaxis() -> SetTitle("Z End Location (cm)");
  hendYZ -> GetXaxis() -> SetTitle("Y End Location (cm)");
  hendYZ -> GetYaxis() -> CenterTitle();
  hendYZ -> GetXaxis() -> CenterTitle();

  hTopCRTHitTime -> SetTitle("");
  hTopCRTHitTime -> GetYaxis() -> SetTitle("Spills");
  hTopCRTHitTime -> GetXaxis() -> SetTitle("CRT hits Time (#mus)");
  hTopCRTHitTime -> GetYaxis() -> CenterTitle();
  hTopCRTHitTime -> GetXaxis() -> CenterTitle();

  hOpFlashTime -> SetTitle("");
  hOpFlashTime -> GetYaxis() -> SetTitle("Spills");
  hOpFlashTime -> GetXaxis() -> SetTitle("OpFlash Time (#mus)");
  hOpFlashTime -> GetYaxis() -> CenterTitle();
  hOpFlashTime -> GetXaxis() -> CenterTitle();

  hNuX -> SetTitle("");
  hNuX -> GetYaxis() -> SetTitle("Spills");
  hNuX -> GetXaxis() -> SetTitle("Nu position (cm)");
  hNuX -> GetYaxis() -> CenterTitle();
  hNuX -> GetXaxis() -> CenterTitle();

  hNuY -> SetTitle("");
  hNuY -> GetYaxis() -> SetTitle("Spills");
  hNuY -> GetXaxis() -> SetTitle("Nu position (cm)");
  hNuY -> GetYaxis() -> CenterTitle();
  hNuY -> GetXaxis() -> CenterTitle();

  hNuZ -> SetTitle("");
  hNuZ -> GetYaxis() -> SetTitle("Spills");
  hNuZ -> GetXaxis() -> SetTitle("Nu position (cm)");
  hNuZ -> GetYaxis() -> CenterTitle();
  hNuZ -> GetXaxis() -> CenterTitle();

  hNuXY -> SetTitle("");
  hNuXY -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXY -> GetYaxis() -> SetTitle("Y position (cm)");
  hNuXY -> GetXaxis() -> CenterTitle();
  hNuXY -> GetYaxis() -> CenterTitle();

  hNuXZ -> SetTitle("");
  hNuXZ -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXZ -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuXZ -> GetXaxis() -> CenterTitle();
  hNuXZ -> GetYaxis() -> CenterTitle();

  hNuYZ -> SetTitle("");
  hNuYZ -> GetXaxis() -> SetTitle("Y position (cm)");
  hNuYZ -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuYZ -> GetXaxis() -> CenterTitle();
  hNuYZ -> GetYaxis() -> CenterTitle();

  hNuXY_FV -> SetTitle("");
  hNuXY_FV -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXY_FV -> GetYaxis() -> SetTitle("Y position (cm)");
  hNuXY_FV -> GetXaxis() -> CenterTitle();
  hNuXY_FV -> GetYaxis() -> CenterTitle();

  hNuXZ_FV -> SetTitle("");
  hNuXZ_FV -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXZ_FV -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuXZ_FV -> GetXaxis() -> CenterTitle();
  hNuXZ_FV -> GetYaxis() -> CenterTitle();

  hNuYZ_FV -> SetTitle("");
  hNuYZ_FV -> GetXaxis() -> SetTitle("Y position (cm)");
  hNuYZ_FV -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuYZ_FV -> GetXaxis() -> CenterTitle();
  hNuYZ_FV -> GetYaxis() -> CenterTitle();

  hNuXY_NFV -> SetTitle("");
  hNuXY_NFV -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXY_NFV -> GetYaxis() -> SetTitle("Y position (cm)");
  hNuXY_NFV -> GetXaxis() -> CenterTitle();
  hNuXY_NFV -> GetYaxis() -> CenterTitle();

  hNuXZ_NFV -> SetTitle("");
  hNuXZ_NFV -> GetXaxis() -> SetTitle("X position (cm)");
  hNuXZ_NFV -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuXZ_NFV -> GetXaxis() -> CenterTitle();
  hNuXZ_NFV -> GetYaxis() -> CenterTitle();

  hNuYZ_NFV -> SetTitle("");
  hNuYZ_NFV -> GetXaxis() -> SetTitle("Y position (cm)");
  hNuYZ_NFV -> GetYaxis() -> SetTitle("Z position (cm)");
  hNuYZ_NFV -> GetXaxis() -> CenterTitle();
  hNuYZ_NFV -> GetYaxis() -> CenterTitle();

  TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg -> SetFillStyle(0);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(0);
  leg -> AddEntry(hcos,          "w/o any cut",       "l");
  leg -> AddEntry(hcos_nofid,    "Are not in the FV", "l");
  //leg -> AddEntry(hcos_Fcuts,    "Neutrinos", "l");
  leg -> AddEntry(hcos_FCcuts,    "OpFlash cut", "l");
  leg -> AddEntry(hcos_Dirycuts, "Longest CR Diry cut", "l");

  TLegend *leg2 = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg2 -> SetFillStyle(0);
  leg2 -> SetTextSize(0.05);
  leg2 -> SetBorderSize(0);
  leg2 -> AddEntry(hEnergy,            "All",       "l");
  leg2 -> AddEntry(hNuEnergy,          "Neutrino",       "l");
  leg2 -> AddEntry(hNuMuEnergy,    "Neutrino and Muon", "l");

  TLegend *leg3 = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg3 -> SetFillStyle(0);
  leg3 -> SetTextSize(0.05);
  leg3 -> SetBorderSize(0);
  leg3 -> AddEntry(hcos,          "w/o any cut",       "l");
  leg3 -> AddEntry(hcos_nofid,    "Are not in the FV", "l");
  //leg3 -> AddEntry(hcos_Fcuts,    "Slice Contains Neutrino", "l");
  leg3 -> AddEntry(htkl_muon,    "OpFlash cut", "l");
  leg3 -> AddEntry(htkl_Diry, "Longest CR Diry slice cut", "l");

  TLegend *leg4 = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg4 -> SetFillStyle(0);
  leg4 -> SetTextSize(0.05);
  leg4 -> SetBorderSize(0);
  leg4 -> AddEntry(hAllMom,            "All",       "l");
  leg4 -> AddEntry(hMuMom,            "Muon",       "l");
  leg4 -> AddEntry(hNuMom,            "Muon & Neutrino",       "l");
  leg4 -> AddEntry(hDiryMom, "Longest CR Diry slice cut", "l");


  // Draw section

  ofstream outfile;
  outfile.open ("integrals.txt");
  outfile << "POT = " << scos.POT() << endl;
  outfile << "POT_norm =  " << kPOTnumi / scos.POT() << endl;
  outfile << endl;

  TCanvas* ccos = new TCanvas("ccos", "");
  outfile << "costh:" << endl;
  hcos  -> GetYaxis()->SetRangeUser(0.,1.3 * hcos -> GetMaximum());
  hcos  -> Draw("hist");
  outfile << "no cut = " << hcos->Integral() << endl;
  hcos_nofid -> Draw("hist same");
  outfile << "nonfiducial = " << hcos_nofid->Integral() << endl;
  hcos_FCcuts -> Draw("hist same");
  outfile << "OpFlash = " << hcos_FCcuts -> Integral() << endl;
  //hcos_Fcuts -> Draw("hist same");
  //outfile << "Neutrinos  = " << hcos_Fcuts -> Integral() << endl;
  hcos_Dirycuts -> Draw("hist same");
  outfile << "Longest CR Diry cut  = " << hcos_Dirycuts ->Integral() << endl << endl;
  leg -> Draw("same");
  ccos -> Print("cos_numi.png");

  TCanvas* ccos_2 = new TCanvas("ccos_2", "");
  hcos_Dirycuts  -> GetYaxis()->SetRangeUser(0.,1.3 * hcos_Dirycuts -> GetMaximum());
  //hcos_Fcuts  -> Draw("hist");
  //hcos_FCcuts -> Draw("hist");
  hcos_Dirycuts -> Draw("hist");
  leg -> Draw("same");
  ccos_2 -> Print("cos_numi_zoom.png");

  TCanvas* cphi = new TCanvas("cphi", "");
  hphi  -> GetYaxis()->SetRangeUser(0.,1.3 * hphi -> GetMaximum());
  hphi -> Draw("hist");
  cphi -> Print("reco_phi.png");

  TCanvas* cdiry = new TCanvas("cdiry", "");
  hCRdir  -> GetYaxis()->SetRangeUser(0.,1.3 * hCRdir -> GetMaximum());
  hCRdir -> Draw("hist");
  cdiry -> Print("CRdiry.png");

  TCanvas* ctkl = new TCanvas("ctkl", "");
  outfile << "track length" << endl;
  htkl  -> GetYaxis()->SetRangeUser(0.,1.3 * htkl -> GetMaximum());
  htkl  -> Draw("hist");
  outfile << "no cuts = " << htkl->Integral() << endl;
  htkl_nofid -> Draw("hist same");
  outfile << "nonfiducial = " << htkl_nofid->Integral() << endl;
  htkl_muon -> Draw("hist same");
  outfile << "OpFlash cut = " << htkl_muon->Integral() << endl;
  //htkl_neutrino -> Draw("hist same");
  //outfile << "muon & neutrino = " << htkl_neutrino->Integral() << endl;
  htkl_Diry -> Draw("hist same");
  outfile << "LongestCRDiry = " << htkl_Diry->Integral() << endl << endl;
  leg3 -> Draw("same");
  ctkl -> Print("track_length.png");

  TCanvas* ctkl2 = new TCanvas("ctkl", "");
  htkl_muon -> GetYaxis()->SetRangeUser(0.,1.3 * htkl_Diry -> GetMaximum());
  htkl_muon -> Draw("hist");
  //htkl_neutrino -> Draw("hist same");
  htkl_Diry -> Draw("hist same");
  leg3 -> Draw("same");
  ctkl2 -> Print("track_length_zoom.png");

  TCanvas* cvtxX = new TCanvas("cvtxX", "");
  vtxX_hist  -> GetYaxis()->SetRangeUser(0.,1.3 * vtxX_hist -> GetMaximum());
  vtxX_hist  -> Draw("hist");
  vtxX_FV_hist  -> Draw("hist same");
  vtxX_NFV_hist  -> Draw("hist same");
  cvtxX -> Print("combined_VtxX.pdf");

  TCanvas* cvtxY = new TCanvas("cvtxY", "");
  vtxY_hist  -> GetYaxis()->SetRangeUser(0.,1.3 * vtxY_hist -> GetMaximum());
  vtxY_hist  -> Draw("hist");
  cvtxY -> Print("nocut_VtxY.pdf");

  TCanvas* cNuMag = new TCanvas("cNuMag", "");
  hNuMag  -> GetYaxis()->SetRangeUser(0.,1.3 * hNuMag -> GetMaximum());
  hNuMag  -> Draw("hist");
  cNuMag -> Print("hNuMag.png");

  TCanvas* cEnergy = new TCanvas("cNuEnergy", "");
  outfile << "Energy" << endl;
  hEnergy  -> GetYaxis()->SetRangeUser(0.,1.3 * hEnergy -> GetMaximum());
  hEnergy -> Draw("hist");
  outfile << "All = " << hEnergy->Integral() << endl;
  hNuEnergy  -> Draw("hist same");
  outfile << "Neutrino = " << hNuEnergy->Integral() << endl;
  hNuMuEnergy -> Draw("hist same");
  outfile << "Neutrino & Muon = " << hNuMuEnergy->Integral() << endl << endl;
  leg2 -> Draw("same");
  cEnergy -> Print("hEnergy.png");

  TCanvas* cMom = new TCanvas("cMuMom", "");
  outfile << "Momentum" << endl;
  hAllMom  -> GetYaxis()->SetRangeUser(0.,1.3 * hAllMom -> GetMaximum());
  hAllMom -> Draw("hist");
  outfile << "All = " << hAllMom->Integral() << endl;
  hMuMom  -> Draw("hist same");
  outfile << "Muon = " << hMuMom->Integral() << endl;
  hNuMom -> Draw("hist same");
  outfile << "Neutrino = " << hNuMom ->Integral() << endl;
  hDiryMom -> Draw("hist same");
  outfile << "Longest CR Diry cut = " << hDiryMom->Integral() << endl << endl;;
  leg4 -> Draw("same");
  cMom -> Print("hMuMom.png");

  TCanvas* cvtxZ = new TCanvas("cvtxZ", "");
  vtxZ_hist  -> GetYaxis()->SetRangeUser(0.,1.3 * vtxZ_hist -> GetMaximum());
  vtxZ_hist  -> Draw("hist");
  cvtxZ -> Print("nocut_VtxZ.pdf");


  TCanvas* XY = new TCanvas("XY", "");
  vtxXY_hist  -> Draw("colz");
  //vtxXY_hist  -> Scale(1./vtxXY_hist -> Integral());
  XY->Print("XY.png");
  outfile << "XY: " << vtxXY_hist->Integral();
  TCanvas* XZ = new TCanvas("XZ", "");
  vtxXZ_hist  -> Draw("colz");
  //vtxXZ_hist  -> Scale(1./vtxXZ_hist -> Integral());
  XZ->Print("XZ.png");
  outfile << " XZ: " << vtxXZ_hist->Integral();
  TCanvas* YZ = new TCanvas("YZ", "");
  vtxYZ_hist  -> Draw("colz");
  //vtxYZ_hist  -> Scale(1./vtxYZ_hist -> Integral());
  YZ->Print("YZ.png");
  outfile << " YZ: " << vtxYZ_hist->Integral() << endl;
  TCanvas* XY_NFV = new TCanvas("XY_NFV", "");
  vtxXY_NFV_hist  -> Draw("colz");
  //vtxXY_NFV_hist  -> Scale(1./vtxXY_NFV_hist -> Integral());
  XY_NFV->Print("XY_NFV.png");
  outfile << "XY_NFV:" << vtxXY_NFV_hist->Integral();
  TCanvas* XZ_NFV = new TCanvas("XZ_NFV", "");
  vtxXZ_NFV_hist  -> Draw("colz");
  //vtxXZ_NFV_hist  -> Scale(1./vtxXZ_NFV_hist -> Integral());
  XZ_NFV->Print("XZ_NFV.png");
  outfile << " XZ_NFV:" << vtxXZ_NFV_hist->Integral();
  TCanvas* YZ_NFV = new TCanvas("YZ_NFV", "");
  vtxYZ_NFV_hist  -> Draw("colz");
  //vtxYZ_NFV_hist  -> Scale(1./vtxYZ_NFV_hist -> Integral());
  YZ_NFV->Print("YZ_NFV.png");
  outfile << " YZ_NFV:" << vtxYZ_NFV_hist->Integral() << endl;
  TCanvas* XY_FV = new TCanvas("XY_FV", "");
  vtxXY_FV_hist  -> Draw("colz");
  //txXY_FV_hist  -> Scale(1./vtxXY_FV_hist -> Integral());
  XY_FV->Print("XY_FV.png");
  outfile << "XY_FV:" << vtxXY_FV_hist->Integral();
  TCanvas* XZ_FV = new TCanvas("XZ_FV", "");
  vtxXZ_FV_hist  -> Draw("colz");
  //vtxXZ_FV_hist  -> Scale(1./vtxXZ_FV_hist -> Integral());
  XZ_FV->Print("XZ_FV.png");
  outfile << " XZ_FV:" << vtxXZ_FV_hist->Integral();
  TCanvas* YZ_FV = new TCanvas("YZ_FV", "");
  vtxYZ_FV_hist  -> Draw("colz");
  //vtxYZ_FV_hist  -> Scale(1./vtxYZ_FV_hist -> Integral());
  YZ_FV->Print("YZ_FV.png");
  outfile << " YZ_NFV:" << vtxYZ_FV_hist->Integral() << endl << endl;

  TCanvas* cendXY = new TCanvas("cendXY", "");
  hendXY  -> Draw("colz");
  cendXY->Print("endXY.png");

  TCanvas* cendXZ = new TCanvas("cendXZ", "");
  hendXZ  -> Draw("colz");
  cendXZ->Print("endXZ.png");

  TCanvas* cendYZ = new TCanvas("cendYZ", "");
  hendYZ  -> Draw("colz");
  cendYZ->Print("endYZ.png");

  TCanvas* cendX = new TCanvas("cendX", "");
  hendX  -> Draw("hist");
  cendX->Print("endX.png");
  outfile << "CRDiry, OPflash, Neutrino & muon: " << hendX->Integral() << endl;

  TCanvas* cendY = new TCanvas("cendY", "");
  hendY  -> Draw("hist");
  cendY->Print("endY.png");
  outfile << "Neutrino & muon: " << hendY->Integral() << endl;

  TCanvas* cendZ = new TCanvas("cendZ", "");
  hendZ  -> Draw("hist");
  cendZ->Print("endZ.png");
  outfile << "No Cuts: " << hendZ->Integral() << endl << endl;


  TCanvas* ccrtT1 = new TCanvas("ccrtT1", "");
  hTopCRTHitTime -> GetYaxis() -> SetRangeUser(0.,1.3 * hTopCRTHitTime -> GetMaximum());
  hTopCRTHitTime -> Draw("hist");
  leg -> Draw("same");
  ccrtT1 -> Print("crt_hit_time_top.pdf");

  TCanvas* copft = new TCanvas("copft", "");
  hOpFlashTime -> GetYaxis() -> SetRangeUser(0.,1.3 * hOpFlashTime -> GetMaximum());
  hOpFlashTime -> Draw("hist");
  leg -> Draw("same");
  copft -> Print("opFlash_time.pdf");

  TCanvas* cNuX = new TCanvas("cNuX", "");
  hNuX -> GetYaxis() -> SetRangeUser(0.,1.3 * hNuX -> GetMaximum());
  hNuX -> Draw("hist");
  cNuX -> Print("NuX.png");

  TCanvas* cNuY = new TCanvas("cNuY", "");
  hNuY -> GetYaxis() -> SetRangeUser(0.,1.3 * hNuY -> GetMaximum());
  hNuY -> Draw("hist");
  cNuY -> Print("NuY.png");

  TCanvas* cNuZ = new TCanvas("cNuZ", "");
  hNuZ -> GetYaxis() -> SetRangeUser(0.,1.3 * hNuZ -> GetMaximum());
  hNuZ -> Draw("hist");
  cNuZ -> Print("NuZ.png");

  outfile << " " << endl;
  TCanvas* cNuXY = new TCanvas("NuXY", "");
  hNuXY  -> Draw("colz");
  //hNuXY -> Scale(1./hNuXY -> Integral());
  cNuXY->Print("NuXY.png");
  outfile << "NuXY:" << hNuXY->Integral();
  TCanvas* cNuXY_FV = new TCanvas("NuXY_FV", "");
  hNuXY_FV  -> Draw("colz");
  //hNuXY_FV -> Scale(1./hNuXY_FV -> Integral());
  cNuXY_FV->Print("NuXY_FV.png");
  outfile << " NuXY_FV:" << hNuXY_FV->Integral();
  TCanvas* cNuXY_NFV = new TCanvas("NuXY_NFV", "");
  hNuXY_NFV  -> Draw("colz");
  //hNuXY_NFV -> Scale(1./hNuXY_NFV -> Integral());
  cNuXY_NFV->Print("NuXY_NFV.png");
  outfile <<  " NuXY_NFV:" << hNuXY_NFV->Integral() << endl;

  TCanvas* cNuXZ = new TCanvas("NuXZ", "");
  hNuXZ  -> Draw("colz");
  //hNuXZ -> Scale(1./hNuXZ -> Integral());
  cNuXZ->Print("NuXZ.png");
  outfile << "NuXZ:" << hNuXZ->Integral();
  TCanvas* cNuXZ_FV = new TCanvas("NuXZ_FV", "");
  hNuXZ_FV  -> Draw("colz");
  //hNuXZ_FV -> Scale(1./hNuXZ_FV -> Integral());
  cNuXZ_FV->Print("NuXZ_FV.png");
  outfile << " NuXZ_FV:" << hNuXZ_FV->Integral();
  TCanvas* cNuXZ_NFV = new TCanvas("NuXZ_NFV", "");
  hNuXZ_NFV  -> Draw("colz");
  //hNuXZ_NFV -> Scale(1./hNuXZ_NFV -> Integral());
  cNuXZ_NFV->Print("NuXZ_NFV.png");
  outfile << " NuXZ_NFV:" << hNuXZ_NFV->Integral() << endl;

  TCanvas* cNuYZ = new TCanvas("NuYZ", "");
  hNuYZ  -> Draw("colz");
  //hNuYZ -> Scale(1./hNuYZ -> Integral());
  cNuYZ->Print("NuYZ.png");
  outfile << "NuYZ:" << hNuYZ->Integral();
  TCanvas* cNuYZ_FV = new TCanvas("NuYZ_FV", "");
  hNuYZ_FV  -> Draw("colz");
  //hNuYZ_FV -> Scale(1./hNuYZ_FV -> Integral());
  cNuYZ_FV->Print("NuYZ_FV.png");
  outfile << " NuYZ_FV:" << hNuYZ_FV->Integral();
  TCanvas* cNuYZ_NFV = new TCanvas("NuYZ_NFV", "");
  hNuYZ_NFV  -> Draw("colz");
  //hNuYZ_NFV -> Scale(1./hNuYZ_NFV -> Integral());
  cNuYZ_NFV->Print("NuYZ_NFV.png");
  outfile << " NuYZ_NFV:" << hNuYZ_NFV->Integral() << endl;
  outfile.close();
}
