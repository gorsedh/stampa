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

TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

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

const Var lenghtNuMI([](const caf::SRSliceProxy* slc){
  //rec.reco.pfp.trk.len
  float len = -5;
  if ( kPTrackInd(slc)>=0)
    {
      auto const& trk =slc->reco.pfp.at(kPTrackInd(slc)).trk;
      len = trk.len;
    }
  return len;
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

void second(){
  const std::string dataFile = "/pnfs/sbn/data/sbn_fd/poms_production/data/Run2/reconstructed/icaruscode_v09_72_00_03/numimajority/flatcaf_unblind/00/[0,1,2]*/data*.flat.caf*.root"; // Data ICARUS run 2
  SpectrumLoader loader(dataFile);

  //double kPOTnumi = 6E20;

  const Binning kTimeBinning = Binning::Simple(100,-10.,10.);
  const Binning kCosThNumiBinning = Binning::Simple(50,-1.,1.);
  const Binning kLenght = Binning::Simple(100,0,800);

  const HistAxis axcos   ("costh",  kCosThNumiBinning, cosThNuMI  );
  const HistAxis axlen   ("length", kLenght, lenghtNuMI );

  /*Spectrum scos       (loader, axcos, kNoCut);
  Spectrum scos_nofid (loader, axcos, !kRFiducial);
  Spectrum scos_fid   (loader, axcos, kRFiducial);
  */
  Spectrum slen (loader, axlen, kNoCut);
  Spectrum slenFid(loader, axlen, kRFiducial);

  loader.Go();
  slen.OverridePOT(1);
  slenFid.OverridePOT(1);

  TH1* hlen = slen.ToTH1(slen.POT(), kBlack);
  TH1* hlenFid = slenFid.ToTH1(slen.POT(), kGreen);
  TCanvas *c = new TCanvas("c","");
  hlen -> GetYaxis()->SetRangeUser(0.,1.3 * hlen -> GetMaximum());
  hlen->Draw();
  hlenFid->Draw("same");
  c->Print("c.pdf");


  // For data
  /*
  scos.OverridePOT(1);
  scos_nofid.OverridePOT(1);
  scos_fid.OverridePOT(1);

  TH1* hcos = scos.ToTH1(scos.POT(), kBlack);
  TH1* hcosNC = scos_nofid.ToTH1(scos_nofid.POT(), kRed);//no fiducial -> not contained
  TH1* hcosC = scos_fid.ToTH1(scos_fid.POT(), kBlue);

  hcos -> SetTitle("");
  hcos -> GetYaxis() -> SetTitle("Slices");
  hcos -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcos -> GetYaxis() -> CenterTitle();
  hcos -> GetXaxis() -> CenterTitle();

  hcosC -> SetTitle("");
  hcosC -> GetYaxis() -> SetTitle("Slices");
  hcosC -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcosC -> GetYaxis() -> CenterTitle();
  hcosC -> GetXaxis() -> CenterTitle();

  hcosNC -> SetTitle("");
  hcosNC -> GetYaxis() -> SetTitle("Slices");
  hcosNC -> GetXaxis() -> SetTitle("Reco cos(#theta_{NuMI})");
  hcosNC -> GetYaxis() -> CenterTitle();
  hcosNC -> GetXaxis() -> CenterTitle();

  TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.9, NULL,"brNDC");
  leg -> SetFillStyle(0);
  leg -> SetTextSize(0.05);
  leg -> SetBorderSize(0);
  leg -> AddEntry(hcos, "w/o any cut", "l");
  leg -> AddEntry(hcosC,"fiducial", "l");
  leg -> AddEntry(hcosNC,"outside fiducial", "l");
  */
  /* TCanvas* ccos = new TCanvas("ccos", "");
  hcos  -> GetYaxis()->SetRangeUser(0.,1.3 * hcos -> GetMaximum());
  hcos  -> Draw("hist");
  leg -> Draw("same");*/
  /*
  TCanvas* ccosc = new TCanvas("ccosh", "");
  hcosC  -> GetYaxis()->SetRangeUser(0.,1.3 * hcosC -> GetMaximum());
  hcosC  -> Draw("hist");
  leg -> Draw("same");


  TCanvas* ccoscc = new TCanvas("ccoshc", "");
  hcosNC  -> GetYaxis()->SetRangeUser(0.,1.3 * hcosNC -> GetMaximum());
  hcosNC  -> Draw("hist");
  leg -> Draw("same");

  TCanvas* ccos = new TCanvas("ccos", "");
  hcos  -> GetYaxis()->SetRangeUser(0.,1.3 * hcos -> GetMaximum());
  hcos  -> Draw("hist");
  hcosC -> Draw("same");
  hcosNC -> Draw("same");
  leg -> Draw("same");

  ccos -> Print("cos_numi.pdf");
  ccos -> Print("cos_numi.png");

  ccosc -> Print("cos_numiC.pdf");
  ccosc -> Print("cos_numiC.png");

  ccoscc -> Print("cos_numiNC.pdf");
  ccoscc -> Print("cos_numiNC.png");*/
}
