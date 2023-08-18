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
#include "TRatioPlot.h"

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

const Var slcvtxX( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.x;
  });

const Var slcvtxY( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.y;
  });

const Var slcvtxZ( [](const caf::SRSliceProxy* slc) {
    return slc->vertex.z;
  });

/*const Var dedxVar( [](const caf::SRSliceProxy* slc) {
    for(Size_t idx = 0; idx < slc->reco.pfp.size(); idx++){
      auto const& trk = slc->reco.pfp.at(idx).trk;//.calo[2].points[0].dedx;
      for(Size_t idx2 = 0; idx <  slc->reco.pfp.at(idx).trk.calo[2].points.size(); idx2++){
      
      }
    }
    //trk.calo[2].points[idxPoint].dedx;
    std::cout << "k\n";
    return slc->vertex.z;
  });//HEREHERE
*/

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
    if (dedxV>0){
      //std::cout << dedxV << "\n";
    return dedxV;
    } else {return -1.f;}
  });


const Var dedxVarAllPart( [](const caf::SRSliceProxy* slc) {
    float dedxV = -1.f;
    int vecSize = slc->reco.pfp.size();
    for ( int i=0; i< vecSize;i++)
      {
        auto const& trk = slc->reco.pfp.at(i).trk;
        float biggestdedx = -1.;
        for(Size_t idx = 0; idx < trk.calo[2].points.size(); idx++){
          auto tempor =trk.calo[2].points.at(idx).dedx;
          if (tempor>biggestdedx){
            biggestdedx = tempor;
          }
        }
        dedxV = biggestdedx;
      }
    if (dedxV>0){
      return dedxV;
    } else {return -1.f;}
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
                                            double this_oft = opflash.firsttime;
                                            rets.push_back(this_oft);
                                          }
                                          return rets;
                                        });



//////////////////////////////////////////////////////
/// here you can define multiple selection cuts,   ///
//////////////////////////////////////////////////////

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



const Cut kNotRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->vertex.x) &&
             ( ( slc->vertex.x > -61.94 - 25 && slc->vertex.x < -358.49 + 25 ) ||
               ( slc->vertex.x < 61.94 + 25 && slc->vertex.x > 358.49 - 25 ) ) &&
             !isnan(slc->vertex.y) &&
             ( slc->vertex.y < -181.86 + 25 && slc->vertex.y > 134.96 - 25 ) &&
             !isnan(slc->vertex.z) &&
             ( slc->vertex.z < -894.951 + 30 && slc->vertex.z > 894.951 - 50 ) );
  });

// const SpillCut kcrtHit_cut([](const caf::SRSpillProxy* sr){
//         for (auto const& crtHit: sr->crt_hits){
//           auto thistime = crtHit.time;
//           if (thistime > -0.1 && thistime < 9.7)
//             return false;
//         }
//         return true;
//       }
//   );

// const SpillCut kOpFlash_cut([](const caf::SRSpillProxy* sr){
//         for (const auto& opflash : sr->opflashes){
//           auto thistime =  opflash.firsttime;
//           if (thistime > -0.1 && thistime < 9.7)
//             return false;
//         }
//         return true;
//       }
//   );

void selectRock2(){

  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader.
  const std::string fdata_MAJORITY = "/pnfs/sbn/data/sbn_fd/poms_production/data/Run2/reconstructed/icaruscode_v09_72_00_03/numimajority/flatcaf_unblind/0*/[0,1,2]*/data*.flat.caf*.root"; // Data ICARUS run 2
  //const std::string fdata_MAJORITY = "/pnfs/sbn/data/sbn_fd/poms_production/data/Run2/reconstructed/icaruscode_v09_72_00_03/numimajority/flatcaf_unblind/0*/0*/data*.flat.caf*.root";
  //const std::string fMC = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_dirt_plus_cosmics/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/0*/0*/detsim*.flat.caf*.root";
  const std::string fMC = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_dirt_plus_cosmics/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/0*/0*/detsim*.flat.caf*.root";



  // Source of events
  SpectrumLoader loader_Data(fdata_MAJORITY);
  SpectrumLoader loader_MC(fMC);

  double kPOTnumi = 6E20; //In Data we don't need this

  // Binning
  /*
  const Binning kCosThNumiBinning = Binning::Simple(50,-1.,1.);
  const Binning kLengthBinning    = Binning::Simple(40,0.,400);
  const Binning kTimeBinning      = Binning::Simple(60,-25.,25.);

  const Binning kPositionX        = Binning::Simple(120,-600.,600.);
  const Binning kPositionY        = Binning::Simple(120,-600.,600.);
  const Binning kPositionZ        = Binning::Simple(300,-1500.,1500.);
  */
  const Binning dedxBinning       = Binning::Simple(200, 0., 50.);//NEW
  const Binning dedxBinningFid    = Binning::Simple(200, 0., 50.);
  const Binning dedxBinningAllPar = Binning::Simple(200, 0., 50.);

  const HistAxis axdedx  ("dedx",      dedxBinning,       dedxVar);
  //const HistAxis axdedxFid  ("dedx",      dedxBinningFid,       dedxVar);
  //const HistAxis axdedxAllPar  ("dedx",      dedxBinningAllPar, dedxVarAllPart);

  Spectrum dedx_all  (loader_Data, axdedx, kNoCut);
  //Spectrum dedx_allFid  (loader_Data, axdedxFid, kRFiducial);
  //Spectrum dedx_allAllPar  (loader_Data, axdedxAllPar, kNoCut);
  
  Spectrum dedx_allMC (loader_MC, axdedx, kNoCut);
  

  // Hist Axis
  /*
  const HistAxis axtkl   ("trklength", kLengthBinning,    kLongestTrackLength);
  const HistAxis axcos   ("costh",     kCosThNumiBinning, cosThNuMI);

  const HistAxis axvtxX  ("vtxx",      kPositionX,        slcvtxX); //kSlcVtxX, kTruthVtxX
  const HistAxis axvtxY  ("vtxy",      kPositionY,        slcvtxY); //kSlcVtxY
  const HistAxis axvtxZ  ("vtxz",      kPositionZ,        slcvtxZ);

  const HistAxis axdedx  ("dedx",      dedxBinning,       dedxVar);
  const HistAxis axdedx  ("dedx",      dedxBinningAllPar, dedxVar);
  */
  // Let's define the spectrums

  //For Slices
  /*
  Spectrum scos       (loader, axcos, kNoCut);
  Spectrum scos_nofid (loader, axcos, !kRFiducial);

  Spectrum stkl       (loader, axtkl, kNoCut);
  Spectrum stkl_nofid (loader, axtkl, !kRFiducial);

  Spectrum vtxX       (loader, axvtxX, kNoCut);
  Spectrum vtxX_NFV   (loader, axvtxX, !kXFiducial);
  Spectrum vtxX_FV    (loader, axvtxX, kXFiducial);
  Spectrum vtxY       (loader, axvtxY, kNoCut);
  Spectrum vtxZ       (loader, axvtxZ, kNoCut);

  Spectrum dedx_all  (loader, axdedx, kNoCut);
  */
  //For Spills
  //Spectrum *sTopCRTHitTime = new Spectrum("TopCRTHitTime", kTimeBinning, loader, spillvarTopCRTHitTime, kNoSpillCut);

  //Spectrum *sOpFlashTime   = new Spectrum("OpFlashTime", kTimeBinning, loader, spillvarOpFlashTime, kNoSpillCut);


  // actually make the spectra
  loader_Data.Go();
  loader_MC.Go();
  
  // For data
  /*
  scos.OverridePOT(1);
  scos_nofid.OverridePOT(1);
  stkl.OverridePOT(1);
  stkl_nofid.OverridePOT(1);
  vtxX.OverridePOT(1);
  vtxX_NFV.OverridePOT(1);
  vtxX_FV.OverridePOT(1);
  vtxY.OverridePOT(1);
  vtxZ.OverridePOT(1);*/
  //dedx_all.OverridePOT(1);
  //dedx_allFid.OverridePOT(1);
  //dedx_allAllPar.OverridePOT(1);

  //dedx_allMC.OverridePOT(1);
  //sTopCRTHitTime -> OverridePOT(1);
  //sOpFlashTime -> OverridePOT(1);

  dedx_all.OverridePOT(1);
  dedx_allMC.OverridePOT(1);
  TH1* hdedx_all = dedx_all.ToTH1(dedx_all.POT(), kBlack);
  TH1* hdedx_allMC = dedx_allMC.ToTH1(dedx_all.POT(), kBlue);

  //TH1* hdedx_allFid = dedx_allFid.ToTH1(dedx_allFid.POT(), kRed);
  
  hdedx_all ->SetTitle("");
  hdedx_all -> GetXaxis() -> SetTitle("dE/dX (MeV/cm)");
  TCanvas* c1 = new TCanvas("c1", "");
  hdedx_all->Draw("histo");
  hdedx_allMC->Draw("hist same");

  //TLegend * leg

  TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9, NULL,"brNDC");
  leg -> SetFillStyle(0);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(0);
  leg -> AddEntry(hdedx_all,  "data","l");
  leg -> AddEntry(hdedx_allMC, "MC", "l");
  leg ->Draw("same");
  //hdedx_allFid->Draw("same");
  c1->Print("highestdedxAllMuons.pdf");

  TCanvas* c2 = new TCanvas("c2", "");
  hdedx_all->Draw("histo");
  hdedx_allMC->Scale(hdedx_all->GetEntries()/hdedx_allMC->GetEntries()); //ATTENTION FROM THIS POINT HISTOMC HAS BEEN SCALED
  hdedx_allMC->Draw("hist same");
  leg ->Draw("same");
  c1->Print("highestdedxAllMuonsRescaled.pdf");

  TCanvas* c3 = new TCanvas("c3", "");
  //TH1F *hpull = new TH1F("hpull","hpull",200,0.,50.); //const Binning dedxBinning       = Binning::Simple(200, 0., 50.);
  auto rp = new TRatioPlot(hdedx_allMC, hdedx_all, "diffsig");
  rp->Draw();
  //hpull ->Add(hdedx_all,hdedx_allMC,1,-1);
  //hpull->Draw("histo");
  c3->Print("residuals.pdf");

  /*
  TH1* hdedx_allAllPart = dedx_allAllPar.ToTH1(dedx_allAllPar.POT(), kBlack);
  hdedx_allAllPart ->SetTitle("");
  hdedx_allAllPart -> GetXaxis() -> SetTitle("dE/dX (MeV/cm)");
  TCanvas* c2 = new TCanvas("c2", "");
  hdedx_allAllPart->Draw("histo");
  c2->Print("highestdedxAllParticles.pdf");
  */
  // Histograms
  /*
  TH1* hcos = scos.ToTH1(scos.POT(), kBlack);
  TH1* hcos_nofid = scos_nofid.ToTH1(scos_nofid.POT(), kGreen+2);
  TH1* htkl = stkl.ToTH1(stkl.POT(), kBlack);
  TH1* htkl_nofid = stkl_nofid.ToTH1(stkl_nofid.POT(), kGreen+2);

  TH1* vtxX_hist = vtxX.ToTH1(vtxX.POT(), kBlack);
  TH1* vtxX_NFV_hist = vtxX_NFV.ToTH1(vtxX_NFV.POT(), kGreen+2);
  TH1* vtxX_FV_hist = vtxX_FV.ToTH1(vtxX_FV.POT(), kRed);
  TH1* vtxY_hist = vtxY.ToTH1(vtxY.POT(), kBlack);
  TH1* vtxZ_hist = vtxZ.ToTH1(vtxZ.POT(), kBlack);

  TH1D* hTopCRTHitTime = sTopCRTHitTime -> ToTH1(sTopCRTHitTime -> POT(), kBlack);
  TH1D* hOpFlashTime = sOpFlashTime -> ToTH1(sOpFlashTime -> POT(), kBlack);

  //////////// You can save your histograms in a root file and draw them later...
  TFile f("histos_selectRock.root", "RECREATE");

  hcos -> Write("hcos");
  hcos_nofid -> Write("hcos_nofid");
  htkl -> Write("htkl");
  htkl_nofid -> Write("htkl_nofid");

  vtxX_hist -> Write("VtxX");
  vtxX_NFV_hist -> Write("VtxX");
  vtxX_FV_hist -> Write("VtxX");
  vtxY_hist -> Write("VtxY");
  vtxZ_hist -> Write("VtxZ");

  hTopCRTHitTime -> Write("hTopCRTHitTime");
  hOpFlashTime -> Write("hOpFlashTime");

  //////////// ... or draw them now

  hcos -> SetTitle("");
  hcos -> GetYaxis() -> SetTitle("Slices");
  hcos -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcos -> GetYaxis() -> CenterTitle();
  hcos -> GetXaxis() -> CenterTitle();

  htkl -> SetTitle("");
  htkl -> GetYaxis() -> SetTitle("Slices");
  htkl -> GetXaxis() -> SetTitle("Track Length (cm)");
  htkl -> GetYaxis() -> CenterTitle();
  htkl -> GetXaxis() -> CenterTitle();

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

  TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg -> SetFillStyle(0);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(0);
  leg -> AddEntry(hcos,          "w/o any cut",       "l");
  leg -> AddEntry(hcos_nofid,    "Are not in the FV", "l");

  // Draw section

  TCanvas* ccos = new TCanvas("ccos", "");
  hcos  -> GetYaxis()->SetRangeUser(0.,1.3 * hcos -> GetMaximum());
  hcos  -> Draw("hist");
  hcos_nofid -> Draw("hist same");
  leg -> Draw("same");
  ccos -> Print("cos_numi.pdf");

  TCanvas* ctkl = new TCanvas("ctkl", "");
  htkl  -> GetYaxis()->SetRangeUser(0.,1.3 * htkl -> GetMaximum());
  htkl  -> Draw("hist");
  htkl_nofid -> Draw("hist same");
  leg -> Draw("same");
  ctkl -> Print("track_length.pdf");

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

  TCanvas* cvtxZ = new TCanvas("cvtxZ", "");
  vtxZ_hist  -> GetYaxis()->SetRangeUser(0.,1.3 * vtxZ_hist -> GetMaximum());
  vtxZ_hist  -> Draw("hist");
  cvtxZ -> Print("nocut_VtxZ.pdf");

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

    */

}
