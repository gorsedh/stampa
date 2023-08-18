#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
//#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
//#include "sbnana/SBNAna/Cuts/NueCuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
//#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TVector3.h"

#include <vector>

#include "TStyle.h" ///// GM

using namespace ana;

Color_t color_nue   = kBlue-7;
Color_t color_numu  = kOrange+2;
Color_t color_nc    = kGreen+2;
Color_t color_cos   = kGray+2;
Style_t line_nue    = kSolid;
Style_t line_numu   = kSolid;
Style_t line_nc     = kDashed;
Style_t line_cos    = kSolid;


// Define all the variables and cuts, including binnings

// ---- VARS -----
// A Var returns a number per slice, a.k.a. variables to plot
//const Var kNTrk = SIMPLEVAR(reco.ntrk);
//const Var kNShw = SIMPLEVAR(reco.nshw);

const Var kSlcFlashTime([](const caf::SRSliceProxy *slc) -> double{
    return ((bool)kSlcHasFlash(slc) ? (float)slc->fmatch.time : -10000.f);
    //          return ((bool)kSlcHasFlash(slc) ? (float)slc->fmatch.time_beam : -10000.f);
  });

const Var kSlcTrkDirY([](const caf::SRSliceProxy *slc) -> double
		      {
			return slc->nuid.crlongtrkdiry;
		      });
/*
  const Var kTrkStrY([](const caf::SRSliceProxy *slc) -> double
       {
         return slc -> reco.trk.start.y;
       });

  const Var kTrkEndY([](const caf::SRSliceProxy *slc) -> double
       {
         return slc -> reco.trk.end.y;
       });
*/
/*  const Var kDifTrkY([](const caf::SRSliceProxy *slc) -> double{
            float dify(-5.f);
            for (auto const& trk : slc->reco.trk){
                dify = trk.start.y-trk.end.y;
            }
         return dify;
       });
*/

const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
        auto const& trk = slc->reco.pfp.at(i).trk;
        //if(trk.bestplane == -1 || slc->reco.pfp.at(i).trackScore < 0.5) continue;
        if(trk.bestplane < 0 || trk.bestplane > 3) continue;

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
	//      const bool MaybeMuonContained = ( Contained && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
	  {
	    Longest = trk.len;
	    PTrackInd = i;
	  }
      }
    return PTrackInd;
  });

const Var kVertexX([](const caf::SRSliceProxy* slc){
    float xvertex = slc->vertex.x;
    return xvertex;
  });

const Var kVertexY([](const caf::SRSliceProxy* slc){
    return slc->vertex.y;
  });
//turns out these are already defined in newer CAFAna, so no need to do this here.
/* const Var kHasTruthMatch( */
/*       [](const caf::SRSliceProxy *slc) -> double */
/*       { */
/*         return ( slc->truth.index != -1); */
/*       } */
/*       ); */

/* const Var kTruthEnergy( */
/*        [](const caf::SRSliceProxy *slc) -> double */
/*        { */
/*          return ( slc->truth.index != -1 ? (float)slc->truth.E : -5.f ); */
/*        }); */

// BASICALLY FROM NumuVars but removing the chi2pid stuff
const Var kPTrackInd_ParedDown([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Longest(0);
    bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
        auto const& trk = slc->reco.pfp.at(i).trk;
        // First we calculate the distance of each track to the slice vertex.
        Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
                      pow( slc->vertex.y - trk.start.y, 2.0 ) +
                      pow( slc->vertex.z - trk.start.z, 2.0 ) );
        // We require that the distance of the track from the slice is less than
        // 10 cm and (BH: REMOVED THIS NEXT BIT) that the parent of the track has been marked as the primary.
        AtSlice = ( Atslc < 10.0 );

        Contained = ( !isnan(trk.end.x) &&
		      ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
			( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
                      !isnan(trk.end.y) &&
                      ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
                      !isnan(trk.end.z) &&
                      ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );
        MaybeMuonExiting = ( !Contained && trk.len > 100);
        MaybeMuonContained = ( Contained && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
          {
            Longest = trk.len;
            PTrackInd = i;
          }
      }
    return PTrackInd;
  });
/*
// Momenta:
// Bascially from NumuVars but only the length part of recoP

const Var kRecoP_Length([](const caf::SRSliceProxy* slc) {
    float p(-5.f);
    bool Contained(false);

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackInd(slc));
        p = trk.rangeP.p_muon;
      }
    return p;
  });

const Var kRecoP_MCS([](const caf::SRSliceProxy* slc) {
    float p(-5.f);
    bool Contained(false);

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackInd(slc));
        p = trk.mcsP.fwdP_muon;
      }
    return p;
  });
*/
const Var kRecoMuMom([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
				 ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
				   ( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
				 !isnan(trk.end.y) &&
				 ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
				 !isnan(trk.end.z) &&
				 ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );
        if(Contained) p = trk.rangeP.p_muon; // Muon range momentum for contained tracks
        else p = trk.mcsP.fwdP_muon;         // And MCS for exiting tracks
      }
    return p;
  });
/*
  const Var kTrueMuonP([](const caf::SRSliceProxy* slc) -> float {
      float p(-5.f);

      if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.trk.at(kPTrackInd(slc));
        p = std::hypot(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
      }
      return p;
    });

const Var kPrimTrkCosth([](const caf::SRSliceProxy *slc)
                        {
                          int muIdx = kPrimMuonIdx(slc);
                          if( muIdx < 0 ) return -5.;

                          return (double)slc->reco.trk[muIdx].costh;
                        });
*/
// TVector3 for NuMI beam

// First attempt
TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

// See get_beam_Dir.C
TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

const Var cosThNuMI( [](const caf::SRSliceProxy* slc) {
    float cosTh(-5.f);

    //std::cout << "rFromNuMI: " << rFromNuMI.X() << " " << rFromNuMI.Y() << " " << rFromNuMI.Z() << std::endl;
    //std::cout << "NuDirection_NuMI: " << NuDirection_NuMI.X() << " " << NuDirection_NuMI.Y() << " " << NuDirection_NuMI.Z() << std::endl;

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

const SpillMultiVar spillvarOpFlashTime([](const caf::SRSpillProxy *sr)
					{
					  std::vector<double> rets;
					  for(const auto& opflash : sr->opflashes){
					    double this_oft = opflash.firsttime;
					    rets.push_back(this_oft);
					  }
					  return rets;
					});

  ///////////
  ///////////
  // CUTS////
  //// DIRECTLY FROM Numu Cuts


const Cut kSlcFlashMatchDataCut([](const caf::SRSliceProxy *slc)
                                {
                                  return (kSlcHasFlashMatch(slc) && slc->fmatch.score>0. && slc->fmatch.score<12.);
                                });
/*
const Cut kFlashTimeCut([](const caf::SRSliceProxy *slc)
                        {
                          return (slc->fmatch.time>-0.2 && slc->fmatch.time<9.9);
//                        return (slc->fmatch.time_beam>-0.2 && slc->fmatch.time_beam<9.9);
                        });
*/
const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
             ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
               ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
             !isnan(slc->vertex.y) &&
             ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
             !isnan(slc->vertex.z) &&
             ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 ) );
  });
/*
  const Cut kFContained([](const caf::SRSliceProxy* slc) -> double {
      bool Contained;
    int PTrackCon(-1);
    for (std::size_t i(0); i < slc->reco.trk.size(); ++i)
      {
        auto const& trk = slc->reco.trk.at(i);
        Contained = ( !isnan(trk.end.x) &&
                    ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
                      ( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
                      !isnan(trk.end.y) &&
                      ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
                      !isnan(trk.end.z) &&
                      ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );

        if ( Contained )
          {
            PTrackCon = i;
          }
      }
    return PTrackCon;

  });
*/
const Cut kPTrack([](const caf::SRSliceProxy* slc) {
    return ( kPTrackInd(slc) >= 0 );
  });

const Cut kPTrackContained([](const caf::SRSliceProxy* slc) {
    int Ind = kPTrackInd(slc);
    bool Contained(false);
    if ( Ind >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        Contained = ( !isnan(trk.end.x) &&
		      ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
			( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
                      !isnan(trk.end.y) &&
                      ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
                      !isnan(trk.end.z) &&
                      ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );
      }
    return Contained;
  });

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

const Var dedxVar( [](const caf::SRSliceProxy* slc) {
    float dedxV = -1.f;
    if ( kPTrackInd(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        float biggestdedx = -1.;
        for(Size_t idx = 0; idx < trk.calo[2].points.size(); idx++){
          float tempor = trk.calo[2].points.at(idx).dedx;
          if (tempor>biggestdedx){
            biggestdedx = tempor;
          }
        }
        dedxV = biggestdedx;
      }
    if (dedxV>=0){
      //std::cout << dedxV << "\n";
      return dedxV;
    } else {
             return -1.f;
           }
  });

const Cut kPTrackExiting([](const caf::SRSliceProxy* slc) {
    int Ind = kPTrackInd(slc);
    bool Exiting(false);
    if ( Ind >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
        Exiting = !( !isnan(trk.end.x) &&
		     ( ( trk.end.x < -61.94 - 10. && trk.end.x > -358.49 + 10. ) ||
		       ( trk.end.x > 61.94 + 10. && trk.end.x < 358.49 - 10. ) ) &&
		     !isnan(trk.end.y) &&
		     ( trk.end.y > -181.86 + 10. && trk.end.y < 134.96 - 10. ) &&
		     !isnan(trk.end.z) &&
		     ( trk.end.z > -894.951 + 10. && trk.end.z < 894.951 - 10. ) );
      }
    return Exiting;
  });

/*
// ---- CUTS -----

  const Cut kNueCC        = kIsNue && !kIsNC;
  const Cut kNumuCC       = kIsNumu && !kIsNC;
const Cut kNC           = kIsNC;
const Cut kTotal        = kNoCut;
const Cut kOther        = kIsNue || kIsNC;
const Cut kThisCosmic   = !kHasNu;
const Cut kThisNu       = kHasNu;

const SpillCut kNueCCSpill      = kIsNueSpill && !kIsNCSpill;
const SpillCut kNumuCCSpill     = kIsNumuSpill && !kIsNCSpill;
const SpillCut kNCSpill         = kIsNCSpill;
const SpillCut kTotalSpill      = kNoSpillCut;
const SpillCut kThisCosmicSpill = kIsCosmicSpill;
*/
const Cut kNucrlongtrkdiryCut([](const caf::SRSliceProxy* slc) {
    return slc->nuid.crlongtrkdiry > -0.7;
  });
const Cut kDiryCut([](const caf::SRSliceProxy* slc) {
    return slc->nuid.crlongtrkdiry > -0.2;
  });
/*
  const SpillCut kCRTHitVetoFD_numi([](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          auto thistime = crtHit.time - 1600.; // manually shift to bring beam spill start to zero
          if (thistime > -0.1 && thistime < 9.7 && crtHit.pe > 100)
            return false;
        }
        return true;
      }
  );

  const SpillCut kCRTHitVetoFD_numi_top(
      [](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          auto thistime = crtHit.time - 1600.; // manually shift to bring beam spill start to zero
          if (thistime > -0.1 && thistime < 9.7 && crtHit.position.y < 425 && crtHit.pe > 100)
            return false;
        }
        return true;
      }
  );


  const Cut kCryo0Vol([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
             ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) &&
             !isnan(slc->vertex.y) &&
             ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
             !isnan(slc->vertex.z) &&
             ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 ) );
  });

  const Cut kCryo1Vol([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
             ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) &&
             !isnan(slc->vertex.y) &&
             ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
             !isnan(slc->vertex.z) &&
             ( slc->vertex.z > -894.951 + 30 && slc->vertex.z < 894.951 - 50 ) );
  });
*/


const Cut kClearCosmic([](const caf::SRSliceProxy* slc){
    return (slc -> is_clear_cosmic == 0);
  });
/*
  const Cut kNotClearCosmic   = !kClearCosmic;


  const Cut kcosThnumiCut([](const caf::SRSliceProxy* slc) {
        return cosThNuMI(slc)>0.6f;
        });
*/
// Full Cuts
const Cut kNumuCC_Selection = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack;
//  const Cut kNumuCC_Selection = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack && kFlashTimeCut;
//  const Cut kNumuCC_Selection = kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack;
const Cut kNumuCC_Sel_Cont  = kNumuCC_Selection && kPTrackContained;
const Cut kNumuCC_Sel_Exit  = kNumuCC_Selection && kPTrackExiting;

const Cut kN1FlashScore = kClearCosmic && kRFiducial && kNucrlongtrkdiryCut && kPTrack;
const Cut kN2FlashScore = kRFiducial && kNucrlongtrkdiryCut && kPTrack;

/*
  const Cut kNumuCC_Sel_Cryo0 = kNumuCC_Selection && kCryo0Vol;
  const Cut kNumuCC_Sel_Cryo1 = kNumuCC_Selection && kCryo1Vol;


// N-1 cuts: apply all except one cut
// // AllCuts = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack;

  const Cut kN1ClearCos   = kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack; //&& kFlashTimeCut;
  const Cut kN1Contained  = kClearCosmic && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack; //&& kFlashTimeCut;
//  const Cut kN1Contained  = kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack;
  const Cut kN1FlashScore = kClearCosmic && kRFiducial && kNucrlongtrkdiryCut && kPTrack; // && kFlashTimeCut;
  const Cut kN1FlashTime  = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack;
  const Cut kN1NuTdirY    = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kPTrack; // && kFlashTimeCut;
  const Cut kN1NuScore    = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kPTrack; // && kFlashTimeCut;
//  const Cut kN1TrkLen     = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut && kPTrack; // && kFlashTimeCut;
  const Cut kN1TrkLen     = kClearCosmic && kRFiducial && kSlcFlashMatchDataCut && kNucrlongtrkdiryCut; // && kFlashTimeCut;


    // SpillCuts
  const SpillCut kRFiducialSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//          if (slc.tmatch.index < 0) continue;
          if (kRFiducial(&slc) && kNumuCC(&slc)) {
            std::cout << " -- fmScore " << slc.fmatch.score  << std::endl;
          ++counter;
          }
        }
      return counter > 0;
  });

  const SpillCut kSlcFlashMatchDataSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//          if (slc.tmatch.index < 0) continue;
          if (kSlcFlashMatchDataCut(&slc)) ++counter;
        }
      return counter > 0;
  });

  const SpillCut kFlashTimeSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//          if (slc.tmatch.index < 0) continue;
          if (kFlashTimeCut(&slc)) ++counter;
        }
      return counter > 0;
  });


  const SpillCut kSlcNuScoreSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//          if (slc.tmatch.index < 0) continue;
          if (kNucrlongtrkdiryCut(&slc)) ++counter;
        }
      return counter > 0;
  });
*/
/*  const SpillCut kTrkLenSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//          if (slc.tmatch.index < 0) continue;
if (kTrkLenCut(&slc)) ++counter;
}
return counter > 0;
});
*/
const SpillCut kNumuCC_SelectionSpillCut([](const caf::SRSpillProxy* sr) {
    unsigned int counter(0);
    for (auto const& slc : sr->slc) {
      if (kNumuCC_Selection(&slc)) ++counter;
    }
    return counter > 0;
  });

/*  const SpillCut kNumuCC_Sel_Cryo0_SpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
          if (kNumuCC_Sel_Cryo0(&slc)) ++counter;
        }
      return counter > 0;
  });

  const SpillCut kNumuCC_Sel_Cryo1_SpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
          if (kNumuCC_Sel_Cryo1(&slc)) ++counter;
        }
      return counter > 0;
  });


  const SpillCut kNumuCC_SelectionSpillCut([](const caf::SRSpillProxy* sr) {
      unsigned int counter(0);
        for (auto const& slc : sr->slc) {
//        for (auto const& trk : slc.reco.trk) {
//        for (auto const& chi2 : slc.reco.trk.chi2pid) {
          if (kNumuCC_Selection(&slc))
          {
//            std::cout << "Event " << sr->hdr.run << " " << sr->hdr.evt
//            std::cout << "Event " << sr->hdr.run << " " << sr->hdr.subrun << " " << sr->hdr.evt
//                      << " -- vtx " << slc.vertex.x << ", " << slc.vertex.y << ", " << slc.vertex.z
//                      << " -- nuScore " << slc.nu_score  << " -- trkLength " << trk.len
//                      << " -- trkend " << trk.end.x << ", " << trk.end.y << ", " << trk.end.z
//                      << " -- chi2 muon " << trk.chi2pid.chi2_muon  << " -- chi2 proton " << trk.chi2pid.2.chi2_proton
//                      << " -- fmTime " << slc.fmatch.time << " -- fmScore " << slc.fmatch.score  << std::endl;
//                      << " -- fmScore " << slc.fmatch.score
//                      << std::endl;
            ++counter;
          }
        }//}//}
      return counter > 0;
  });


*/
// Binning

const Binning kEvtBinning       = Binning::Simple(1000000,0.f,1000000.f);
// const Binning kSliceBinning     = Binning::Simple(3,0.f,3.f);
  const Binning kPositionX        = Binning::Simple(120,-600.,600.);
  const Binning kPositionY        = Binning::Simple(120,-600.,600.);
  const Binning kPositionZ        = Binning::Simple(300,-1500.,1500.);
/*//  const Binning kTimeBinning      = Binning::Simple(60,-15.,15.);
  const Binning kTimeBinning      = Binning::Simple(50,-50.,50.);
//  const Binning kTimeBinning      = Binning::Simple(500,-20.,30.);
  const Binning kNucltdiryBinning = Binning::Simple(50,-1.,1.f);
  const Binning kNuScoBinning     = Binning::Simple(50,0.,1.f);*/
const Binning kLengthBinning    = Binning::Simple(40,0.,800);
const Binning kMomBinning    = Binning::Simple(40,0.,8);
const Binning kEnergyBinning    = Binning::Simple(30.,0.,12);
const Binning kCosThNumiBinning  = Binning::Simple(50,-1.,1.);
const Binning kCRTPMTTimeBinning = Binning::Simple(3000, -0.15, 0.15);
const Binning kMatchIDBinning    = Binning::Simple(10, 0, 10);
const Binning dedxBinning       = Binning::Simple(30, 0., 45.);//NEW
//const Binning dedxBinningFid    = Binning::Simple(200, 0., 50.);
//const Binning dedxBinningAllPar = Binning::Simple(200, 0., 50.);
/*  const Binning kdEdxBinning      = Binning::Simple(40,0.f,10.f);
  const Binning kFlashBinning     = Binning::Simple(70,0.f,70.f);
  const Binning kShowereBinning   = Binning::Simple(25,0.f,5.f);
  const Binning kPEBinning        = Binning::Simple(60,0.,600.);
  const Binning kMuPBinning       = Binning::Simple(25,0.,5.);
  const Binning kDifTYBinning     = Binning::Simple(120,-600.,600.);
*/

const Binning kTimeBinning   = Binning::Simple(100,-50.,50.); //attention may duplicate
const Binning kTimeBinningExt   = Binning::Simple(100,-200.,200.); //attention may duplicate

// ---- HistAxis ----
//const HistAxis axslice ("count",          kSliceBinning,          kCounting);
  const HistAxis axvtxX  ("vtxx",           kPositionX,             kSlcVtxX); //kSlcVtxX, kTruthVtxX
  const HistAxis axvtxY  ("vtxy",           kPositionY,             kSlcVtxY); //kSlcVtxY
  const HistAxis axvtxZ  ("vtxz",           kPositionZ,             kSlcVtxZ);
/*  const HistAxis axfms   ("flashScore",     kFlashBinning,          kSlcFlashScore);
  const HistAxis axfmt   ("flashTime",      kTimeBinning,           kSlcFlashTime);
  const HistAxis axtdy   ("cltdiry",        kNucltdiryBinning,      kSlcTrkDirY);
  const HistAxis axnus   ("nuScore",        kNuScoBinning,          kSlcNuScore);*/
const HistAxis axtkl   ("trklength",      kLengthBinning,         kLongestTrackLength);
const HistAxis axcos   ("costh",          kCosThNumiBinning,      cosThNuMI);
const HistAxis axenergy   ("energy",      kEnergyBinning,         kTruthEnergy);
const HistAxis axmom ("momentum", kMomBinning, kRecoMuMom);
const HistAxis axdedx  ("dedx",      dedxBinning,       dedxVar);
//  const HistAxis axdEdx  ("dEdx",           kdEdxBinning,           kRecoShower_BestdEdx);
/*  const HistAxis axtnue  ("truthenergy",    kLowEnergyGeVBinning,   kTruthEnergy);
  const HistAxis axmup   ("RecoMuMom",      kMuPBinning,            kRecoMuMom);
  const HistAxis axmpt   ("TrueMuMom",      kMuPBinning,            kTrueMuonP);
  const HistAxis axpct   ("RecMomcosth",    kCosThNumiBinning,      kPrimTrkCosth);
//  const HistAxis axdty   ("DifseTrkY",      kDifTYBinning,          kDifTrkY);
*/
