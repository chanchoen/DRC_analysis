#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float truth = std::stof(argv[3]);
  float high = std::stof(argv[4]);
  TString outputname = argv[5];

  gStyle->SetOptFit(1);

  std::pair<double, double> fCalib = std::make_pair(1470./20, 22650./20);  // (copper1.25m) 1929./20., 26660./20.); // <ceren, scint>
  //std::cout<<fCalib.first<<fCalib.second<<std::endl;

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename), false);
  drInterface->GetChain("DRsim");

  TH1F* tEdep = new TH1F("totEdep","Total Energy deposit;MeV;Evt",100,low*1000.,high*1000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tHit_C = new TH1F("Hit_C","# of p.e. of Cerenkov ch.;# of p.e.;Evt",200,0,3000*(truth/20));
  tHit_C->Sumw2(); tHit_C->SetLineColor(kBlue); tHit_C->SetLineWidth(2);
  TH1F* tHit_S = new TH1F("Hit_S","# of p.e. of Scintillation ch.;# of p.e.;Evt",200,0,40000*(truth/20));
  tHit_S->Sumw2(); tHit_S->SetLineColor(kRed); tHit_S->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak","Momentum leak;MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu","Neutrino energy leak;MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);

  TH1F* tE_S = new TH1F("E_S","Scintillation Energy;GeV;Evt",100,low,high);
  tE_S->Sumw2(); tE_S->SetLineColor(kRed); tE_S->SetLineWidth(2);
  TH1F* tE_C = new TH1F("E_C","Cerenkov Energy;GeV;Evt",100,low,high);
  tE_C->Sumw2(); tE_C->SetLineColor(kBlue); tE_C->SetLineWidth(2);
  TH1F* tE_T = new TH1F("E_C+S", "Total Energy;GeV;Evt", 100, low*2, high*2);
  tE_T->Sumw2(); tE_T->SetLineColor(kBlack); tE_T->SetLineWidth(2);

  TH1F* tT_C = new TH1F("time_C","Cerenkov time;ns;p.e.",600,10.,70.);
  tT_C->Sumw2(); tT_C->SetLineColor(kBlue); tT_C->SetLineWidth(2);
  TH1F* tT_S = new TH1F("time_S","Scint time;ns;p.e.",600,10.,70.);
  tT_S->Sumw2(); tT_S->SetLineColor(kRed); tT_S->SetLineWidth(2);
  TH1F* tWav_S = new TH1F("wavlen_S","Scint wavelength;nm;p.e.",120,300.,900.);
  tWav_S->Sumw2(); tWav_S->SetLineColor(kRed); tWav_S->SetLineWidth(2);
  TH1F* tWav_C = new TH1F("wavlen_C","Cerenkov wavelength;nm;p.e.",120,300.,900.);
  tWav_C->Sumw2(); tWav_C->SetLineColor(kBlue); tWav_C->SetLineWidth(2);
  TH1F* tNhit_S = new TH1F("nHits_S","Number of Scint p.e./SiPM;p.e.;n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);
  TH1F* tNhit_C = new TH1F("nHits_C","Number of Cerenkov p.e./SiPM;p.e.;n",50,0.,50.);
  tNhit_C->Sumw2(); tNhit_C->SetLineColor(kBlue); tNhit_C->SetLineWidth(2);

  TH2D* t2DhitC = new TH2D("2D Hit C", "", 420, -0.5, 419.5, 420, -0.5, 419.5); t2DhitC->Sumw2(); t2DhitC->SetStats(0);
  TH2D* t2DhitS = new TH2D("2D Hit S", "", 420, -0.5, 419.5, 420, -0.5, 419.5); t2DhitS->Sumw2(); t2DhitS->SetStats(0);

  std::vector<float> E_Ss, E_Cs;

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float Edep = 0.;
    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      Edep += edep.Edep;
    }
    tEdep->Fill(Edep);

    float Pleak = 0.;
    float Eleak_nu = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
      }
    }
    tP_leak->Fill(Pleak);
    tP_leak_nu->Fill(Eleak_nu);

    int nHitC = 0; int nHitS = 0;
    float sEtmp =0; float cEtmp =0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum;
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        int plateNum = sipm->x; int fiberNum = sipm->y; 
        if ( RecoInterface::IsCerenkov(sipm->x,sipm->y) ) {
          tNhit_C->Fill(sipm->count);
          for (const auto timepair : sipm->timeStruct) {
            tT_C->Fill(timepair.first.first+0.05,timepair.second);
            if (timepair.first.first < 35) {
              nHitC += timepair.second;
	      cEtmp = nHitC /fCalib.first;
              t2DhitC->Fill(60*(moduleNum%7)+fiberNum, 60*(moduleNum/7)+plateNum, timepair.second);
            }
          }
          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_C->Fill(wavpair.first.first,wavpair.second);
          }
        } else {
          tNhit_S->Fill(sipm->count);
          nHitS += sipm->count;
          sEtmp = nHitS /fCalib.second;
          t2DhitS->Fill(60*(moduleNum%7)+fiberNum, 60*(moduleNum/7)+plateNum, sipm->count);
          for (const auto timepair : sipm->timeStruct) {
            tT_S->Fill(timepair.first.first+0.05,timepair.second);
          }
          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_S->Fill(wavpair.first.first,wavpair.second);
          }
        }
      }
    }

    E_Cs.push_back(cEtmp);
    E_Ss.push_back(sEtmp);

    tHit_C->Fill(nHitC);
    tHit_S->Fill(nHitS);
    tE_C->Fill(cEtmp);
    tE_S->Fill(sEtmp);
    tE_T->Fill(cEtmp+sEtmp);

  } // event loop
  drInterface->close();

  TCanvas* c = new TCanvas("c","");

  tEdep->Draw("Hist"); c->SaveAs(outputname+"_Edep.png");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs(outputname+"_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs(outputname+"_Pleak_nu.png");
  c->SetLogy(0);

  tHit_C->Draw("Hist"); c->SaveAs(outputname+"_nHitpEventC.png");
  tHit_S->Draw("Hist"); c->SaveAs(outputname+"_nHitpEventS.png");

  c->cd();
  tE_S->Draw("Hist"); c->Update();
  TPaveStats* statsE_S = (TPaveStats*)c->GetPrimitive("stats");
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetY1NDC(.6); statsE_S->SetY2NDC(.8);

  tE_C->Draw("Hist&sames"); c->Update();
  TPaveStats* statsE_C = (TPaveStats*)c->GetPrimitive("stats");
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetY1NDC(.8); statsE_C->SetY2NDC(1.);
  c->SaveAs(outputname+"_EcsHist.png");

  TF1* grE_C = new TF1("Cfit", "gaus", low, high);
  TF1* grE_S = new TF1("Sfit", "gaus", low, high);
  TF1* grE_DR = new TF1("Totalfit", "gaus", low*2, high*2);
  grE_C->SetLineColor(kBlue);
  grE_S->SetLineColor(kRed);
  grE_DR->SetLineColor(kBlack);
  tE_C->SetOption("p"); tE_C->Fit(grE_C, "R+&same");
  tE_S->SetOption("p"); tE_S->Fit(grE_S, "R+&same");

  c->cd();
  tE_S->SetTitle("");
  tE_S->Draw(""); c->Update();
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetX1NDC(.7);
  statsE_S->SetY1NDC(.4); statsE_S->SetY2NDC(.7);

  tE_C->Draw("sames"); c->Update();
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetX1NDC(.7);
  statsE_C->SetY1NDC(.7); statsE_C->SetY2NDC(1.);
  c->SaveAs(outputname+"_Ecs.png");

  tE_T->SetOption("p"); tE_T->Fit(grE_DR, "R+&same");
  tE_T->SetOption("p"); tE_T->Fit(grE_DR, "R+&same");
  tE_T->Draw(""); c->SaveAs(outputname+"_E_T.png");

  t2DhitS->Draw("COLZ"); c->SaveAs(outputname+"_n2DHitS.png");
  t2DhitC->Draw("COLZ"); c->SaveAs(outputname+"_n2DHitC.png");

  TGraph* grSvsC = new TGraph(entries, &(E_Ss[0]), &(E_Cs[0]));
  grSvsC->SetTitle("SvsC;E_S;E_C");
  grSvsC->SetMarkerSize(0.5); grSvsC->SetMarkerStyle(20);
  grSvsC->GetXaxis()->SetLimits(0.,high);
  grSvsC->GetYaxis()->SetRangeUser(0.,high);
  grSvsC->SetMaximum(high);
  grSvsC->SetMinimum(0.);
  grSvsC->Draw("ap"); c->SaveAs(outputname+"_SvsC.png");

  tT_C->Draw("Hist"); c->SaveAs(outputname+"_tC.png");
  tT_S->Draw("Hist"); c->SaveAs(outputname+"_tS.png");
  tWav_C->Draw("Hist"); c->SaveAs(outputname+"_wavC.png");
  tWav_S->Draw("Hist"); c->SaveAs(outputname+"_wavS.png");
  tNhit_C->Draw("Hist"); c->SaveAs(outputname+"_nhitC.png");
  tNhit_S->Draw("Hist"); c->SaveAs(outputname+"_nhitS.png");
}
