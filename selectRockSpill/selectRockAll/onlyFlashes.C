#pragma once

#include "TFile.h"
#include "TH1.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
// NuMINumuXSec
//#include "ICARUSCRTPMTMatching.h"
#include "helper_numuCCSelCuts.h"

#include "TCanvas.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "iostream"
#include "fstream"

using namespace ana;

////////////////////////////////////
//            VARS               // see helper file
///////////////////////////////////

////////////////////////////////////
//             Cuts              //
///////////////////////////////////

const SpillCut OpFlash_cut([](const caf::SRSpillProxy* sr){
    for (const auto& opflash : sr->opflashes){
      auto thistime =  opflash.firsttime;
      if (thistime > -0.1 && thistime < 9.7) {
	return true;
      } 
    }
    return false;
  }
  );

const SpillCut OpFlash_MC([](const caf::SRSpillProxy* sr){
    for (const auto& opflash : sr->opflashes){
      auto thistime =  opflash.firsttime -  sr->hdr.triggerinfo.trigger_within_gate;
      if (thistime > -0.1 && thistime < 9.7) {
	return true;
      } 
    }
    return false;
  }
  );

const Cut kSlcTrkDiry([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2) {
      return true;
    } else {
      return false;
    }
  });

const Cut kContained([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kHasMunNu([](const caf::SRSliceProxy* slc) -> double {
    if (kHasNu(slc) && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kSlcTrkDiryTkl([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kLongestTrackLength(slc) > 15) {
      return true;
    } else {
      return false;
    }
  });

const Cut kDiryNuTkl([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kHasNu(slc) && kLongestTrackLength(slc) > 100) {
      return true;
    } else {
      return false;
    }
  });

const Cut kDiryNoNuTkl([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && !kHasNu(slc) && kLongestTrackLength(slc) > 100) {
      return true;
    } else {
      return false;
    }
  });

const Cut kNuContained([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kHasNu(slc) && kLongestTrackLength(slc) > 100 && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kNoNuContained([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && !kHasNu(slc) && kLongestTrackLength(slc) > 100 && kmuon_contained(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kDiryNu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && kHasNu(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kDiryNoNu([](const caf::SRSliceProxy* slc) -> double {
    if (slc->nuid.crlongtrkdiry > -0.2 && !kHasNu(slc)) {
      return true;
    } else {
      return false;
    }
  });

const Cut kTkl([](const caf::SRSliceProxy* slc) -> double {
    if (kLongestTrackLength(slc) > 100) {
      return true;
    } else {
      return false;
    }
  });

////////////////////////////////////
//       MAIN EVENT LOOP          //
///////////////////////////////////

void onlyFlashes(){

  //cpmt.SetGateType(NUMI);

  //  const std::string fnus = "/pnfs/icarus/persistent/users/jskim/data/run_8515/flatcaf/v09_63_00_02/221212_UseOldFlashMatching_FixCAFT1/NUMIMAJORITY/flatcaf_0.root"; //numi Beam ON Majority
  //const std::string fcos = "/pnfs/icarus/persistent/users/jskim/data/run_8515/flatcaf/v09_63_00_02/221212_UseOldFlashMatching_FixCAFT1/OffBeamNUMIMAJORITY/flatcaf_0.root"; //numi Beam OFF Majority
  //const std::string fcosmics_data = "/pnfs/icarus/scratch/users/aheggest/crt/mc/NuMI_intime_cosmics/v09_75_01/caf/NuMI_intime_cosmics_246files_merged.root"; //Cosmics
  //const std::string foffbeam = "/pnfs/icarus/scratch/users/aheggest/crt/data/run_9726/offbeamNuMImajority/v09_75_01/caf/data_offbeamNuMI_crtpmt_Blind_193files_merged.OKTOLOOK.flat.caf.root"; //foffbeam
  const std::string fMC_majority_old = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_Nu_Phase1/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/[0]*/[2,3,4]*/detsim*.flat.caf*.root"; // old MC sample //it was flatcaf/[0]*/[0,1,2,3,4,5,6,7,8,9]*/detsim* //FULL PATH//"/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_Nu_Phase1/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/*/*/detsim*.flat.caf*.root"
  ///pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_Nu_Phase1/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/
  const std::string fdirt_mc_old = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_dirt_plus_cosmics/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/[0,1,2,3]*/[0,1,2,3,4,5,6,7,8,9]*/detsim*.flat.caf*.root"; //[0,1,2,3]*/[0,1,2,3,4,5,6,7,8,9]* //full path // /[0,1,2,3]*/[0,1,2,3,4,5,6,7,8,9]*/
  const std::string fdata_majority_old = "/pnfs/icarus/scratch/users/gputnam/DMCP2023G/majority-3t1p/57386892_[2,3,4,5,6][0,1,2,3]*/data*Prescaled*.root";///57386892_[2,3,4,5][0,1,3]*/ //full path "/pnfs/icarus/scratch/users/gputnam/DMCP2023G/majority-3t1p/57386892_*/data*Prescaled*.root"
  const std::string fdirt_mc = "/pnfs/icarus/scratch/users/aheggest/crt/mc/NuMI_dirt_plus_cosmics/v09_75_01/caf/out/detsim_crtpmt_109files_mergedflat.caf.root"; // dirt MC
  const std::string fdata_majority = "/pnfs/icarus/scratch/users/aheggest/crt/data/run2/NuMImajority/v09_75_01/caf/data_NuMImajority_crtpmt_77files_merged_Blind.OKTOLOOK.flat.caf.root"; //beam on data
  const std::string fMC_majority = "/pnfs/icarus/scratch/users/aheggest/crt/mc/NuMI_nu_plus_cosmics/v09_75_01/caf/out/2619881_*/detsim_2d_icarus_detsim_stage0_stage1*flat.caf.root";
  //const std::string fMC_intimecosmic = "/pnfs/sbn/data/sbn_fd/poms_production/2023A_ICARUS_NuMI_MC_intime_cosmics/pretuned_signal_shape/mc/reconstructed/icaruscode_v09_72_00_03/flatcaf/[0]*/[0]*/detsim*.flat.caf*.root";

  // Source of events
  SpectrumLoader loader_nus(fdirt_mc);
  //SpectrumLoader loader_intime(fMC_intimecosmic);
  SpectrumLoader loader_maj(fMC_majority_old);
  SpectrumLoader loader_cos(fdata_majority_old);


  Spectrum scos_cos (loader_cos,     axcos,   kNoCut);

  Spectrum sxAxis_cos (loader_cos, axvtxX, kNoSpillCut, kSlcTrkDiryTkl); //no spill cut
  Spectrum sxAxisFlashCut_cos (loader_cos, axvtxX, OpFlash_cut, kSlcTrkDiryTkl);
  Spectrum sxAxisFlashMcCut_cos (loader_cos, axvtxX, OpFlash_MC, kSlcTrkDiryTkl);
  Spectrum sxAxis_maj (loader_maj, axvtxX, kNoSpillCut, kSlcTrkDiryTkl); //no spill cut
  Spectrum sxAxisFlashCut_maj (loader_maj, axvtxX, OpFlash_MC, kSlcTrkDiryTkl);
  Spectrum sxAxis_nus (loader_nus, axvtxX, kNoSpillCut, kSlcTrkDiryTkl); //no spill cut
  Spectrum sxAxisFlashCut_nus (loader_nus,axvtxX, OpFlash_MC, kSlcTrkDiryTkl);


  //loader, asse, spillcut, slicecut
  //Spectrum for the vertex variable coi due tagli dopo pescali  
  Spectrum *sOpFlashTimeCos = new Spectrum("OpFlashTime", kTimeBinning, loader_cos, spillvarOpFlashTime, kNoSpillCut);
  Spectrum *sOpFlashTimeMaj = new Spectrum("OpFlashTime", kTimeBinning, loader_maj, spillvarOpFlashTime, kNoSpillCut);
  Spectrum *sOpFlashTimeNus = new Spectrum("OpFlashTime", kTimeBinning, loader_nus, spillvarOpFlashTime, kNoSpillCut);
  Spectrum *sOpFlashTimeCosCUT = new Spectrum("OpFlashTime", kTimeBinning, loader_cos, spillvarOpFlashTime, OpFlash_cut);
  Spectrum *sOpFlashTimeMajCUT = new Spectrum("OpFlashTime", kTimeBinning, loader_maj, spillvarOpFlashTime, OpFlash_cut);
  Spectrum *sOpFlashTimeNusCUT = new Spectrum("OpFlashTime", kTimeBinning, loader_nus, spillvarOpFlashTime, OpFlash_cut);
  
  Spectrum *sOpFlashTimeMajMcCUT = new Spectrum("OpFlashTime", kTimeBinning, loader_maj, spillvarOpFlashTime, OpFlash_MC);
  Spectrum *sOpFlashTimeNusMcCUT = new Spectrum("OpFlashTime", kTimeBinning, loader_nus, spillvarOpFlashTime, OpFlash_MC);
  //maybe OpFlash_MC
  
  //Da fixare //sOpFlashTime[Nus,Maj,Cos]+eventually CUT


  double POT = scos_cos.POT(); //2.7E12;//6.0E20; //NuMI

  loader_nus.Go();
  loader_maj.Go();
  loader_cos.Go();

  //  sOpFlashTimeCos->OverridePOT(1); //tried it, does not work

  ////sOpFlashTime[Nus,Maj,Cos]+eventually CUT
  TH1* hOpFlashTimeNus = sOpFlashTimeNus->ToTH1(scos_cos.POT(),kBlack);
  TH1* hOpFlashTimeNusCUT = sOpFlashTimeNusCUT->ToTH1(scos_cos.POT(),kRed);
  TH1* hOpFlashTimeNusMcCUT = sOpFlashTimeNusMcCUT->ToTH1(scos_cos.POT(),kGreen);

  TH1* hOpFlashTimeCos = sOpFlashTimeCos->ToTH1(scos_cos.POT(),kGreen);
  TH1* hOpFlashTimeCosCUT = sOpFlashTimeCosCUT->ToTH1(scos_cos.POT(),kRed);

  TH1* hOpFlashTimeMaj = sOpFlashTimeMaj->ToTH1(scos_cos.POT(),kBlack);
  TH1* hOpFlashTimeMajCUT = sOpFlashTimeMajCUT->ToTH1(scos_cos.POT(),kRed);
  TH1* hOpFlashTimeMajMcCUT = sOpFlashTimeMajMcCUT->ToTH1(scos_cos.POT(),kGreen);

  TH1* hxAxisFlashCut_cos = sxAxisFlashCut_cos.ToTH1(scos_cos.POT(),kBlue);
  TH1* hxAxisFlashMcCut_cos = sxAxisFlashMcCut_cos.ToTH1(scos_cos.POT(),kGreen);
  TH1* hxAxis_cos = sxAxis_cos.ToTH1(scos_cos.POT(),kRed);
  TH1* hxAxisFlashCut_maj = sxAxisFlashCut_maj.ToTH1(scos_cos.POT(),kBlue);
  TH1* hxAxis_maj = sxAxis_maj.ToTH1(scos_cos.POT(),kRed);
  TH1* hxAxisFlashCut_nus = sxAxisFlashCut_nus.ToTH1(scos_cos.POT(),kBlue);
  TH1* hxAxis_nus = sxAxis_nus.ToTH1(scos_cos.POT(),kRed);

  hxAxis_cos->SetTitle("Data");
  hxAxis_maj->SetTitle("MC");
  hxAxis_nus->SetTitle("Dirt MC");
  //sxAxisFlashCut
  //QUIIIII TCanvas

  /////////////////////////////
  //RIVEDIAMO TIME CONSTRAINT//
  /////////////////////////////  
  TCanvas* cFlashNus = new TCanvas("cFlashNus","");
  hOpFlashTimeNus -> Draw("hist");
  hOpFlashTimeNusCUT -> Draw("hist same");
  hOpFlashTimeNusMcCUT -> Draw("hist same");
  TLegend* legNus = new TLegend(0.75, 0.75, 0.9, 0.9, NULL,"legNus");
  legNus->AddEntry(hOpFlashTimeNus, "No OpFlash cut", "l");
  legNus->AddEntry(hOpFlashTimeNusCUT, "OpFlash cut", "l");
  legNus->Draw("same");
  cFlashNus -> SetTitle("Dirt MC");
  cFlashNus -> Print("flashNus.pdf");

  TCanvas* cFlashMaj = new TCanvas("cFlashMaj","");
  hOpFlashTimeMaj->SetTitle("MC");
  hOpFlashTimeMaj -> Draw("hist");
  hOpFlashTimeMajCUT -> Draw("hist same");
  hOpFlashTimeMajMcCUT -> Draw("hist same");
  TLegend* legMaj = new TLegend(0.6, 0.6, 0.9, 0.9, NULL,"legMaj");
  legMaj -> SetFillStyle(0);
  //legMaj -> SetTextSize(0.2);
  //legMaj -> SetBorderSize(0);
  legMaj->AddEntry(hOpFlashTimeMaj, "No OpFlash cut", "l");
  legMaj->AddEntry(hOpFlashTimeMajCUT, "OpFlash cut", "l");
  //legMaj->Draw();
  //cFlashMaj -> SetTitle("MC");
  cFlashMaj -> Print("flashMaj.pdf");

  TCanvas* cFlashCos = new TCanvas("cFlashCos","");
  hOpFlashTimeCos -> Draw("hist");
  hOpFlashTimeCosCUT -> Draw("hist same");
  TLegend* legCos = new TLegend(0.75, 0.75, 0.9, 0.9, NULL,"legMaj");
  legCos->AddEntry(hOpFlashTimeCos, "No OpFlash cut", "l");
  legCos->AddEntry(hOpFlashTimeCosCUT, "OpFlash cut", "l");
  legCos->Draw("same");
  cFlashCos -> SetTitle("data");
  cFlashCos -> Print("flashCos.pdf");

  //POSIZIONE VERTICE X//
  TLegend* legxAxisFlashNus = new TLegend(0.45, 0.75, 0.7, 0.9, NULL,"legxAxisNus");
  legxAxisFlashNus->AddEntry(hxAxis_nus, "No OpFlash cut", "l");
  legxAxisFlashNus->AddEntry(hxAxisFlashCut_nus, "OpFlash cut", "l");//sxAxisFlashCut_nus
  TCanvas* cxAxisFlashNus = new TCanvas("cxAxisFlashNus","");
  hxAxis_nus-> Draw("hist");
  hxAxisFlashCut_nus-> Draw("hist same");
  legxAxisFlashNus->Draw("same");
  cxAxisFlashNus->Print("xAxisFlashCut_nus.pdf");

  
  TLegend* legxAxisFlashCos = new TLegend(0.45, 0.75, 0.7, 0.9, NULL,"legxAxisCos");
  legxAxisFlashCos->AddEntry(hxAxis_cos, "No OpFlash cut", "l");
  legxAxisFlashCos->AddEntry(hxAxisFlashCut_cos, "OpFlash cut", "l");
  TCanvas* cxAxisFlashCos = new TCanvas("cxAxisFlashCos","");
  hxAxis_cos-> Draw("hist");
  hxAxisFlashCut_cos-> Draw("hist same");
  hxAxisFlashMcCut_cos->Draw("hist same");
  legxAxisFlashCos->Draw("same");
  cxAxisFlashCos->Print("xAxisFlashCut_cos.pdf");

  
  TLegend* legxAxisFlashMaj = new TLegend(0.45, 0.75, 0.7, 0.9, NULL,"legxAxisNus");
  legxAxisFlashMaj->AddEntry(hxAxis_maj, "No OpFlash cut", "l");
  legxAxisFlashMaj->AddEntry(hxAxisFlashCut_maj, "OpFlash cut", "l");
  TCanvas* cxAxisFlashMaj = new TCanvas("cxAxisFlashMaj","");
  hxAxis_maj-> Draw("hist");
  hxAxisFlashCut_maj-> Draw("hist same");
  legxAxisFlashMaj->Draw("same");
  cxAxisFlashMaj->Print("xAxisFlashCut_maj.pdf");


  std::ofstream histooutfile("flashEntriesIndedx.txt");
  if (histooutfile.is_open()) {
    histooutfile << "All entries: " << 0  << " w/spill cut\n";
    histooutfile << "All entries: " << 0 << " no spill cut\n";
    histooutfile << "MC  entries: " << 0 << "\n";
    histooutfile << "MC  entries: " << 0 << " no spill cut\n";
    histooutfile.close();
    std::cout << "Successfully wrote to file." << std::endl;
  }



  ////////////////////////////////////
  //             .PNG              //
  ///////////////////////////////////



  ////////////////////////////////////
  //            INTEGRALS           //
  ///////////////////////////////////



}
