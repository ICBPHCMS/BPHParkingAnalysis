/*
g++ -Wall -o analyzer `root-config --cflags --glibs` -lRooFitCore analyzer.cpp
./analyzer --isEle (0,1) --typeSelection (tightCB) --ntupleList (list.txt) --JOBid (1,2..) --outputFolder ("outfolder") --nMaxEvents (-1, N) --testFile ("path")
*/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <TLorentzVector.h>

using namespace RooFit;

const int kBToKstllMax = 100000;
const int KLepMax = 100000;
const int kProbeTracksMax = 100000;
const float ElectronMass = 0.5109989e-3;
const float MuonMass = 0.10565837;
const float KaonMass = 0.493677;
const float PionMass = 0.139570;

int main(int argc, char **argv){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isEleFinalState = -1;
  std::string typeSelection = "-1";
  std::string ntupleList = "-1";
  std::string JOBid = "-1";
  std::string outputFolder = "-1";
  int nMaxEvents = -1;
  std::string testFile = "-1";
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isEle") {
      if (i + 1 < argc) {
	isEleFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isEle option requires one argument " << std::endl;
	return 1;
      }
    }
  } for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--typeSelection") {
            if (i + 1 < argc) {
                typeSelection = argv[i+1];
                break;
            } else {
            std::cerr << " --typeSelection option requires one argument " << std::endl;
            return 1;
            }
        }
    } 
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--ntupleList") {
      if (i + 1 < argc) {
        ntupleList = argv[i+1];
        break;
      } else {
	std::cerr << " --ntupleList option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--JOBid") {
      if (i + 1 < argc) {
        JOBid = argv[i+1];
        break;
      } else {
	std::cerr << " --JOBid option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--outputFolder") {
      if (i + 1 < argc) {
        outputFolder = argv[i+1];
        break;
      } else {
	std::cerr << " --outputFolder option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--nMaxEvents") {
      if (i + 1 < argc) {
        nMaxEvents = atoi(argv[i+1]);
        break;
      } else {
	std::cerr << " --nMaxEvents option requires one argument " << std::endl;
        return 1;
      }
    }
  } for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--testFile") {
      if (i + 1 < argc) {
        testFile = argv[i+1];
        break;
      } else {
	std::cerr << " --testFile option requires one argument " << std::endl;
        return 1;
      }
    }
  }

  

  if(ntupleList != "-1" && (JOBid == "-1" || outputFolder == "-1")){
    std::cout << " configuration ERROR => splitting file based but missing JOBid and output folder " << std::endl;
    return -1;
  }

  if(ntupleList == "-1" && testFile == "-1"){
    std::cout << " configuration ERROR => need a file list or a test file " << std::endl;
    return -1;
  }


  std::cout << " isEleFinalState = " << isEleFinalState << " typeSelection = " << typeSelection << " ntupleList = " << ntupleList << " JOBid = " << JOBid << " outputFolder = " << outputFolder << " nMaxEvents = " << nMaxEvents << std::endl;


  gROOT->Reset();
  gROOT->Macro("./setStyle.C");
  gSystem->Load("libRooFit") ;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain* t1 = new TChain("Events");
  
  if(ntupleList != "-1"){
    std::string rootFileName;
    std::ifstream inFileLong;
    inFileLong.open(ntupleList.c_str(), std::ios::in);
    while(!inFileLong.eof()){
        if(inFileLong >> rootFileName){
            t1->Add(rootFileName.c_str());
            std::cout << " adding " << rootFileName << std::endl;
        }
    }
  }
  else{
    t1->Add(testFile.c_str());
  }

  
  int nEvts = t1->GetEntries();

  std::cout << " #initial n. events: " << nEvts << std::endl;

  ULong64_t event = 0;
  
  uint nBinTree = 0;

  float BToKstll_pt[kBToKstllMax];
  float BToKstll_fit_mass[kBToKstllMax];
  float BToKstll_mass[kBToKstllMax];
  float BToKstll_cos2D[kBToKstllMax];
  float BToKstll_svprob[kBToKstllMax];
  float BToKstll_l_xy[kBToKstllMax];
  float BToKstll_l_xy_unc[kBToKstllMax];
  float BToKstll_mll[kBToKstllMax];
  
  int lep1_Idx[kBToKstllMax];
  int lep2_Idx[kBToKstllMax];
  int kaon_Idx[kBToKstllMax];  

  float ProbeTracks_pt[kProbeTracksMax];
  float ProbeTracks_eta[kProbeTracksMax];
  float ProbeTracks_phi[kProbeTracksMax];
  float ProbeTracks_DCASig[kProbeTracksMax];
  
  float lep_pt[KLepMax];
  float lep_eta[KLepMax];
  float lep_phi[KLepMax];
  bool Electron_isPFoverlap[KLepMax];
  bool Electron_isLowPt[KLepMax];
  bool Electron_isPF[KLepMax];
  bool Electron_convVeto[KLepMax];
  float Electron_mvaId[KLepMax];
  
  
  t1->SetBranchStatus("*", 0);

  t1->SetBranchStatus("event", 1);                      t1->SetBranchAddress("event", &event);

  t1->SetBranchStatus("ProbeTracks_pt", 1);             t1->SetBranchAddress("ProbeTracks_pt", &ProbeTracks_pt);
  t1->SetBranchStatus("ProbeTracks_eta", 1);            t1->SetBranchAddress("ProbeTracks_eta", &ProbeTracks_eta);
  t1->SetBranchStatus("ProbeTracks_phi", 1);            t1->SetBranchAddress("ProbeTracks_phi", &ProbeTracks_phi);
  t1->SetBranchStatus("ProbeTracks_DCASig", 1);         t1->SetBranchAddress("ProbeTracks_DCASig", &ProbeTracks_DCASig);  
  
  
  if(isEleFinalState){    
      
    t1->SetBranchStatus("nBToKEE", 1);                  t1->SetBranchAddress("nBToKEE", &nBinTree);
    t1->SetBranchStatus("BToKEE_l1Idx", 1);             t1->SetBranchAddress("BToKEE_l1Idx", &lep1_Idx);
    t1->SetBranchStatus("BToKEE_l2Idx", 1);             t1->SetBranchAddress("BToKEE_l2Idx", &lep2_Idx);
    t1->SetBranchStatus("BToKEE_kIdx", 1);              t1->SetBranchAddress("BToKEE_kIdx", &kaon_Idx);
    t1->SetBranchStatus("BToKEE_pt", 1);                t1->SetBranchAddress("BToKEE_pt", &BToKstll_pt);
    t1->SetBranchStatus("BToKEE_mass", 1);              t1->SetBranchAddress("BToKEE_mass", &BToKstll_mass);
    t1->SetBranchStatus("BToKEE_fit_mass", 1);          t1->SetBranchAddress("BToKEE_fit_mass", &BToKstll_fit_mass);
    t1->SetBranchStatus("BToKEE_cos2D", 1);             t1->SetBranchAddress("BToKEE_cos2D", &BToKstll_cos2D);
    t1->SetBranchStatus("BToKEE_svprob", 1);            t1->SetBranchAddress("BToKEE_svprob", &BToKstll_svprob);
    t1->SetBranchStatus("BToKEE_l_xy", 1);              t1->SetBranchAddress("BToKEE_l_xy", &BToKstll_l_xy);
    t1->SetBranchStatus("BToKEE_l_xy_unc", 1);          t1->SetBranchAddress("BToKEE_l_xy_unc", &BToKstll_l_xy_unc);
    //t1->SetBranchStatus("BToKEE_mll_raw", 1);         t1->SetBranchAddress("BToKEE_mll_raw", &BToKstll_mll);
    t1->SetBranchStatus("BToKEE_mll_fullfit", 1);       t1->SetBranchAddress("BToKEE_mll_fullfit", &BToKstll_mll);
    t1->SetBranchStatus("Electron_pt", 1);              t1->SetBranchAddress("Electron_pt", &lep_pt);
    t1->SetBranchStatus("Electron_eta", 1);             t1->SetBranchAddress("Electron_eta", &lep_eta);
    t1->SetBranchStatus("Electron_phi", 1);             t1->SetBranchAddress("Electron_phi", &lep_phi);
    t1->SetBranchStatus("Electron_isPFoverlap", 1);     t1->SetBranchAddress("Electron_isPFoverlap", &Electron_isPFoverlap);
    t1->SetBranchStatus("Electron_isLowPt", 1);         t1->SetBranchAddress("Electron_isLowPt", &Electron_isLowPt);
    t1->SetBranchStatus("Electron_isPF", 1);            t1->SetBranchAddress("Electron_isPF", &Electron_isPF);
    t1->SetBranchStatus("Electron_convVeto", 1);        t1->SetBranchAddress("Electron_convVeto", &Electron_convVeto);
    t1->SetBranchStatus("Electron_mvaId", 1);           t1->SetBranchAddress("Electron_mvaId", &Electron_mvaId);
  }
  else{    
      
    t1->SetBranchStatus("nBToKMuMu", 1);                   t1->SetBranchAddress("nBToKMuMu", &nBinTree);
    t1->SetBranchStatus("BToKMuMu_l1Idx", 1);              t1->SetBranchAddress("BToKMuMu_l1Idx", &lep1_Idx);
    t1->SetBranchStatus("BToKMuMu_l2Idx", 1);              t1->SetBranchAddress("BToKMuMu_l2Idx", &lep2_Idx);
    t1->SetBranchStatus("BToKMuMu_kIdx", 1);               t1->SetBranchAddress("BToKMuMu_kIdx", &kaon_Idx);    
    t1->SetBranchStatus("BToKMuMu_pt", 1);                 t1->SetBranchAddress("BToKMuMu_pt", &BToKstll_pt);
    t1->SetBranchStatus("BToKMuMu_mass", 1);               t1->SetBranchAddress("BToKMuMu_mass", &BToKstll_mass);    
    t1->SetBranchStatus("BToKMuMu_fit_mass", 1);           t1->SetBranchAddress("BToKMuMu_fit_mass", &BToKstll_fit_mass);
    t1->SetBranchStatus("BToKMuMu_cos2D", 1);              t1->SetBranchAddress("BToKMuMu_cos2D", &BToKstll_cos2D);
    t1->SetBranchStatus("BToKMuMu_svprob", 1);             t1->SetBranchAddress("BToKMuMu_svprob", &BToKstll_svprob);
    t1->SetBranchStatus("BToKMuMu_l_xy", 1);               t1->SetBranchAddress("BToKMuMu_l_xy", &BToKstll_l_xy);
    t1->SetBranchStatus("BToKMuMu_l_xy_unc", 1);           t1->SetBranchAddress("BToKMuMu_l_xy_unc", &BToKstll_l_xy_unc);
    //t1->SetBranchStatus("BToKMuMu_mll_raw", 1);          t1->SetBranchAddress("BToKMuMu_mll_raw", &BToKstll_mll);
    t1->SetBranchStatus("BToKMuMu_mll_fullfit", 1);        t1->SetBranchAddress("BToKMuMu_mll_fullfit", &BToKstll_mll);
    t1->SetBranchStatus("Muon_pt", 1);                     t1->SetBranchAddress("Muon_pt", &lep_pt);
    t1->SetBranchStatus("Muon_eta", 1);                    t1->SetBranchAddress("Muon_eta", &lep_eta);
    t1->SetBranchStatus("Muon_phi", 1);                    t1->SetBranchAddress("Muon_phi", &lep_phi);
    
  }

  
  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);
  llMassBoundary.push_back(100.);


  std::string outName = "outMassHistos_Kee_";
  if(!isEleFinalState) outName = "outMassHistos_Kmumu_";
  if(JOBid != "-1") outName = outputFolder; // + "/" +outName + "_JOB_"+JOBid;
  else  outName += ".root";
  TFile outMassHistos(outName.c_str(), "recreate");

  ///histos: 1 per bin plus inclusive
  TH1F* hAlpha[7];
  TH1F* hCLVtx[7];
  TH1F* hLxy[7];
  TH1F* hKaonpt[7];
  TH1F* hBpt[7];
  TH1F* hLep1pt[7];
  TH1F* hLep2pt[7];  
  TH1F* hLep1pt_EB[7];
  TH1F* hLep2pt_EB[7];
  TH1F* hLep1pt_EE[7];
  TH1F* hLep2pt_EE[7];  
  TH2F* pTele2_vs_pTele1[7];
  TH1F* hllRefitMass[7];  
  TH2F* hllRefitMass_vs_Bmass[7];
  TH1F* hBfitMass_all[7];
  TH1F* hBmass_all[7];
  TH1F* hBfitMass_bothPF[7];
  TH1F* hBfitMass_bothLowPt[7];
  TH1F* hBfitMass_bothLowPt_noPFoverlap[7];
  TH1F* hBfitMass_mix[7];
  TH1F* hBfitMass_mix_noPFoverlap[7];

  for(int ij=0; ij<7; ++ij){
    hAlpha[ij] = new TH1F(Form("hAlpha_%d", ij), "", 500, 0, 1.1);
    hAlpha[ij]->Sumw2();
    hAlpha[ij]->SetLineColor(kRed);
    hAlpha[ij]->SetLineWidth(2);

    hCLVtx[ij] = new TH1F(Form("hCLVtx_%d", ij), "", 100, 0., 1.);
    hCLVtx[ij]->Sumw2();
    hCLVtx[ij]->SetLineColor(kRed);
    hCLVtx[ij]->SetLineWidth(2);

    hLxy[ij] = new TH1F(Form("hLxy_%d", ij), "", 100, 0., 100.);
    hLxy[ij]->Sumw2();
    hLxy[ij]->SetLineColor(kRed);
    hLxy[ij]->SetLineWidth(2);  

    hKaonpt[ij] = new TH1F(Form("hKaonpt_%d", ij), "", 100, 0., 10.);
    hKaonpt[ij]->Sumw2();
    hKaonpt[ij]->SetLineColor(kRed);
    hKaonpt[ij]->SetLineWidth(2);   
    
    hBpt[ij] = new TH1F(Form("hBpt_%d", ij), "", 200, 0., 100.);
    hBpt[ij]->Sumw2();
    hBpt[ij]->SetLineColor(kRed);
    hBpt[ij]->SetLineWidth(2);    
    
    hLep1pt[ij] = new TH1F(Form("hLep1pt_%d", ij), "", 100, 0., 10.);
    hLep1pt[ij]->Sumw2();
    hLep1pt[ij]->SetLineColor(kRed);
    hLep1pt[ij]->SetLineWidth(2);   

    hLep2pt[ij] = new TH1F(Form("hLep2pt_%d", ij), "", 100, 0., 10.);
    hLep2pt[ij]->Sumw2();
    hLep2pt[ij]->SetLineColor(kRed);
    hLep2pt[ij]->SetLineWidth(2);
    
    hLep1pt_EB[ij] = new TH1F(Form("hLep1pt_EB_%d", ij), "", 100, 0., 10.);
    hLep1pt_EB[ij]->Sumw2();
    hLep1pt_EB[ij]->SetLineColor(kRed);
    hLep1pt_EB[ij]->SetLineWidth(2);

    hLep2pt_EB[ij] = new TH1F(Form("hLep2pt_EB_%d", ij), "", 100, 0., 10.);
    hLep2pt_EB[ij]->Sumw2();
    hLep2pt_EB[ij]->SetLineColor(kRed);
    hLep2pt_EB[ij]->SetLineWidth(2);

    hLep1pt_EE[ij] = new TH1F(Form("hLep1pt_EE_%d", ij), "", 100, 0., 10.);
    hLep1pt_EE[ij]->Sumw2();
    hLep1pt_EE[ij]->SetLineColor(kRed);
    hLep1pt_EE[ij]->SetLineWidth(2);

    hLep2pt_EE[ij] = new TH1F(Form("hLep2pt_EE_%d", ij), "", 100, 0., 10.);
    hLep2pt_EE[ij]->Sumw2();
    hLep2pt_EE[ij]->SetLineColor(kRed);
    hLep2pt_EE[ij]->SetLineWidth(2);
	  
    pTele2_vs_pTele1[ij] = new TH2F(Form("pTele2_vs_pTele1_%d", ij), "", 300, 0., 30., 300, 0., 30.);
    pTele2_vs_pTele1[ij]->Sumw2();
    pTele2_vs_pTele1[ij]->SetMarkerColor(kRed);
    pTele2_vs_pTele1[ij]->SetMarkerStyle(20);
    
    hllRefitMass[ij] = new TH1F(Form("hllRefitMass_%d", ij), "", 750, 0., 15.);
    hllRefitMass[ij]->Sumw2();
    hllRefitMass[ij]->SetLineColor(kRed);
    hllRefitMass[ij]->SetLineWidth(2);

    hllRefitMass_vs_Bmass[ij] = new TH2F(Form("hllRefitMass_vs_Bmass_%d", ij), "", 500, 0., 15., 500, 0., 15.);
    hllRefitMass_vs_Bmass[ij]->Sumw2();
    hllRefitMass_vs_Bmass[ij]->SetMarkerColor(kRed);
    hllRefitMass_vs_Bmass[ij]->SetMarkerStyle(20);    
    
    hBfitMass_all[ij] = new TH1F(Form("BfitMass_all_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_all[ij]->Sumw2();
    hBfitMass_all[ij]->SetLineColor(kRed);
    hBfitMass_all[ij]->SetLineWidth(2);
    
    hBmass_all[ij] = new TH1F(Form("Bmass_all_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_all[ij]->Sumw2();
    hBmass_all[ij]->SetLineColor(kRed);
    hBmass_all[ij]->SetLineWidth(2);    

    hBfitMass_bothPF[ij] = new TH1F(Form("BfitMass_bothPF_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_bothPF[ij]->Sumw2();
    hBfitMass_bothPF[ij]->SetLineColor(kRed);
    hBfitMass_bothPF[ij]->SetLineWidth(2);

    hBfitMass_bothLowPt[ij] = new TH1F(Form("BfitMass_bothLowPt_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_bothLowPt[ij]->Sumw2();
    hBfitMass_bothLowPt[ij]->SetLineColor(kRed);
    hBfitMass_bothLowPt[ij]->SetLineWidth(2);

    hBfitMass_bothLowPt_noPFoverlap[ij] = new TH1F(Form("BfitMass_bothLowPt_noPFoverlap_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_bothLowPt_noPFoverlap[ij]->Sumw2();
    hBfitMass_bothLowPt_noPFoverlap[ij]->SetLineColor(kRed);
    hBfitMass_bothLowPt_noPFoverlap[ij]->SetLineWidth(2);

    hBfitMass_mix[ij] = new TH1F(Form("BfitMass_mix_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_mix[ij]->Sumw2();
    hBfitMass_mix[ij]->SetLineColor(kRed);
    hBfitMass_mix[ij]->SetLineWidth(2);

    hBfitMass_mix_noPFoverlap[ij] = new TH1F(Form("BfitMass_mix_noPFoverlap_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBfitMass_mix_noPFoverlap[ij]->Sumw2();
    hBfitMass_mix_noPFoverlap[ij]->SetLineColor(kRed);
    hBfitMass_mix_noPFoverlap[ij]->SetLineWidth(2);    
  }


  float nEv_selected[7] = {0.};


  if(nMaxEvents == -1) nMaxEvents = nEvts;
  for(int iEvt = 0; iEvt<nMaxEvents; ++iEvt){
    if(iEvt%500000 == 0) std::cout << " >>> processing event " << iEvt << " " << 1.*iEvt/nEvts*100. << std::endl;

    t1->GetEntry(iEvt);
    
    for(unsigned int i_Btree = 0; i_Btree < nBinTree; ++i_Btree){//i_Btree loop	
    
      int l1_Index = lep1_Idx[i_Btree];
      int l2_Index = lep2_Idx[i_Btree];
      int kaon_Index = kaon_Idx[i_Btree];        
        
      bool all = false;
      bool bothPF = false;
      bool bothLowPt = false;
      bool bothLowPt_noPFoverlap = false;
      bool mix = false;
      bool mix_noPFoverlap = false;
      
      //l1 & l2 selections
	  if( !Electron_convVeto[l1_Index] || Electron_mvaId[l1_Index] < 3.96 || lep_pt[l1_Index] < 1.5) continue; 
	  if( !Electron_convVeto[l2_Index] || Electron_mvaId[l2_Index] < 3.96 || lep_pt[l2_Index] < 0.5) continue;      

      if(typeSelection == "tightCB"){
        if( BToKstll_svprob[i_Btree] < 0.1 ) continue;
        if( ProbeTracks_pt[kaon_Index] < 3. || BToKstll_pt[i_Btree] < 3. ) continue;
        if( BToKstll_cos2D[i_Btree] < 0.999 ) continue;
        if( (BToKstll_l_xy[i_Btree]/sqrt(BToKstll_l_xy_unc[i_Btree])) < 6. ) continue;
        if( ProbeTracks_DCASig[kaon_Index] < 2. ) continue;
      }

      all = bool( !Electron_isPFoverlap[l1_Index] && !Electron_isPFoverlap[l2_Index] );        
      bothPF = bool( Electron_isPF[l1_Index] && Electron_isPF[l2_Index] );      
      bothLowPt = bool( Electron_isLowPt[l1_Index] && Electron_isLowPt[l2_Index] );      
      bothLowPt_noPFoverlap = bool( bothLowPt && all );      
      mix = bool( (Electron_isPF[l1_Index] && Electron_isLowPt[l2_Index]) || (Electron_isLowPt[l1_Index] && Electron_isPF[l2_Index]) );
      mix_noPFoverlap = bool( mix && all ); 
    
      
      float llMass = BToKstll_mll[i_Btree];
      int massBin = -1;
      for(unsigned int kl=0; kl<llMassBoundary.size()-1; ++kl){
        if(llMass >= llMassBoundary[kl] && llMass < llMassBoundary[kl+1]){
            massBin = kl;
            break;
        }
      }

      if(massBin != -1) ++nEv_selected[massBin];
      else ++nEv_selected[6];     
          
    
      
      //histograms for each mass bin
      if(massBin != -1){
          
        if(all){
            hAlpha[massBin]->Fill(BToKstll_cos2D[i_Btree]);
            hCLVtx[massBin]->Fill(BToKstll_svprob[i_Btree]);
            hLxy[massBin]->Fill(BToKstll_l_xy[i_Btree]/sqrt(BToKstll_l_xy_unc[i_Btree]));
            hKaonpt[massBin]->Fill(ProbeTracks_pt[kaon_Index]);
            hBpt[massBin]->Fill(BToKstll_pt[i_Btree]);

            hLep1pt[massBin]->Fill(lep_pt[l1_Index]);
            hLep2pt[massBin]->Fill(lep_pt[l2_Index]);
            if(std::abs(lep_eta[l1_Index]) < 1.47) hLep1pt_EB[massBin]->Fill(lep_pt[l1_Index]);
            else hLep1pt_EE[massBin]->Fill(lep_pt[l1_Index]);
            if(std::abs(lep_eta[l2_Index]) < 1.47) hLep2pt_EB[massBin]->Fill(lep_pt[l2_Index]);
            else hLep2pt_EE[massBin]->Fill(lep_pt[l2_Index]);
            pTele2_vs_pTele1[massBin]->Fill(lep_pt[l1_Index], lep_pt[l2_Index]);

            hllRefitMass[massBin]->Fill(llMass);
            hllRefitMass_vs_Bmass[massBin]->Fill(BToKstll_fit_mass[i_Btree], llMass);

            hBfitMass_all[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
            hBmass_all[massBin]->Fill(BToKstll_mass[i_Btree]);
        }  
        
        if( bothPF ) hBfitMass_bothPF[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
        if( bothLowPt ) hBfitMass_bothLowPt[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
        if( bothLowPt_noPFoverlap ) hBfitMass_bothLowPt_noPFoverlap[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
        if( mix ) hBfitMass_mix[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
        if( mix_noPFoverlap ) hBfitMass_mix_noPFoverlap[massBin]->Fill(BToKstll_fit_mass[i_Btree]);
    
      }
      
      
      //histograms inclusive over all m(ll)
      if( all ){
        hAlpha[6]->Fill(BToKstll_cos2D[i_Btree]);
        hCLVtx[6]->Fill(BToKstll_svprob[i_Btree]);
        hLxy[6]->Fill(BToKstll_l_xy[i_Btree]/sqrt(BToKstll_l_xy_unc[i_Btree]));
        hKaonpt[6]->Fill(ProbeTracks_pt[kaon_Index]);
        hBpt[6]->Fill(BToKstll_pt[i_Btree]);

        hLep1pt[6]->Fill(lep_pt[l1_Index]);
        hLep2pt[6]->Fill(lep_pt[l2_Index]);
        if(std::abs(lep_eta[l1_Index]) < 1.47) hLep1pt_EB[6]->Fill(lep_pt[l1_Index]);
        else hLep1pt_EE[6]->Fill(lep_pt[l1_Index]);
        if(std::abs(lep_eta[l2_Index]) < 1.47) hLep2pt_EB[6]->Fill(lep_pt[l2_Index]);
        else hLep2pt_EE[6]->Fill(lep_pt[l2_Index]);
        pTele2_vs_pTele1[6]->Fill(lep_pt[l1_Index], lep_pt[l2_Index]);
      
        hllRefitMass[6]->Fill(llMass);
        hllRefitMass_vs_Bmass[6]->Fill(BToKstll_fit_mass[i_Btree], llMass);
      
        hBfitMass_all[6]->Fill(BToKstll_fit_mass[i_Btree]);
        hBmass_all[6]->Fill(BToKstll_mass[i_Btree]);
      }
      
      if( bothPF ) hBfitMass_bothPF[6]->Fill(BToKstll_fit_mass[i_Btree]);
      if( bothLowPt ) hBfitMass_bothLowPt[6]->Fill(BToKstll_fit_mass[i_Btree]);
      if( bothLowPt_noPFoverlap ) hBfitMass_bothLowPt_noPFoverlap[6]->Fill(BToKstll_fit_mass[i_Btree]);
      if( mix ) hBfitMass_mix[6]->Fill(BToKstll_fit_mass[i_Btree]);
      if( mix_noPFoverlap ) hBfitMass_mix_noPFoverlap[6]->Fill(BToKstll_fit_mass[i_Btree]);
      
    }//i_Btree loop

  }//loop over events

  

  outMassHistos.cd();
  
  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<7; ++ij){
    hAlpha[ij]->Write(hAlpha[ij]->GetName());
    hCLVtx[ij]->Write(hCLVtx[ij]->GetName());
    hLxy[ij]->Write(hLxy[ij]->GetName());
    hKaonpt[ij]->Write(hKaonpt[ij]->GetName());
    hBpt[ij]->Write(hBpt[ij]->GetName());
    hLep1pt[ij]->Write(hLep1pt[ij]->GetName());
    hLep2pt[ij]->Write(hLep2pt[ij]->GetName());
    hLep1pt_EB[ij]->Write(hLep1pt_EB[ij]->GetName());
    hLep1pt_EE[ij]->Write(hLep1pt_EE[ij]->GetName());
    hLep2pt_EB[ij]->Write(hLep2pt_EB[ij]->GetName());
    hLep2pt_EE[ij]->Write(hLep2pt_EE[ij]->GetName());
    pTele2_vs_pTele1[ij]->Write(pTele2_vs_pTele1[ij]->GetName());
    hllRefitMass[ij]->Write(hllRefitMass[ij]->GetName());
    hllRefitMass_vs_Bmass[ij]->Write(hllRefitMass_vs_Bmass[ij]->GetName());

    hBfitMass_all[ij]->Write(hBfitMass_all[ij]->GetName());
    std::cout << " >>> hBfitMass_all[ij]->GetName() = " << hBfitMass_all[ij]->GetName() << " entries = " << hBfitMass_all[ij]->GetEntries() << std::endl;
    hBmass_all[ij]->Write(hBmass_all[ij]->GetName());
    hBfitMass_bothPF[ij]->Write(hBfitMass_bothPF[ij]->GetName());
    hBfitMass_bothLowPt[ij]->Write(hBfitMass_bothLowPt[ij]->GetName());
    hBfitMass_bothLowPt_noPFoverlap[ij]->Write(hBfitMass_bothLowPt_noPFoverlap[ij]->GetName());
    hBfitMass_mix[ij]->Write(hBfitMass_mix[ij]->GetName());
    hBfitMass_mix_noPFoverlap[ij]->Write(hBfitMass_mix_noPFoverlap[ij]->GetName());
    

    if(ij > 5) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
	      << " \n \t selected Events = " << nEv_selected[ij] << std::endl;

  }
  outMassHistos.Close();
  
}  
