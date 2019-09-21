/*
g++ -Wall -o analyzer_CommonCode `root-config --cflags --glibs` -lRooFitCore analyzer_CommonCode.cpp

./analyzer_CommonCode --isEle (0,1) --dataset (-1, runA, runB, runD, MC) --run (1,2,3,..., -1) --ntupleList (list.txt) --JOBid (1,2..) --outputFolder ("outfolder") --nMaxEvents (-1, N) --saveSelectedNTU (1,0) --outSelectedNTU (path for selected ntuples) --testFile ("path") --nMaxTriplets (1,2,...)
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
const int kMuonMax = 100000;
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
  std::string dataset = "-1";
  std::string BPHRun = "-1";
  std::string ntupleList = "-1";
  std::string JOBid = "-1";
  std::string outputFolder = "-1";
  int nMaxEvents = -1;
  int saveOUTntu = 0;
  std::string outSelectedNTU = "-1";
  std::string testFile = "-1";
  unsigned int nMaxTriplets = -1;
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
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--dataset") {
      if (i + 1 < argc) {
        dataset = argv[i+1];
        break;
      } else {
	std::cerr << " --dataset option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--run") {
      if (i + 1 < argc) {
        BPHRun = argv[i+1];
        break;
      } else {
	std::cerr << " --run option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
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
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveSelectedNTU") {
      if (i + 1 < argc) {
        saveOUTntu = atoi(argv[i+1]);
        break;
      } else {
	std::cerr << " --saveSelectedNTU option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--outSelectedNTU") {
      if (i + 1 < argc) {
        outSelectedNTU = argv[i+1];
        break;
      } else {
	std::cerr << " --outSelectedNTU option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
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
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--nMaxTriplets") {
      if (i + 1 < argc) {
	nMaxTriplets = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --nMaxTriplets option requires one argument " << std::endl;
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

  if(saveOUTntu != 0 && outSelectedNTU == "-1"){
    std::cout << " configuration ERROR => missing output folder for final trees " << std::endl;
    return -1;
  }

  std::cout << " isEleFinalState = " << isEleFinalState << " dataset = " << dataset << " BPHRun = " << BPHRun << " ntupleList = " << ntupleList << " JOBid = " << JOBid << " outputFolder = " << outputFolder
	    << " nMaxEvents = " << nMaxEvents << " saveOUTntu = " << saveOUTntu << " outSelectedNTU = " << outSelectedNTU << " nMaxTriplets " << nMaxTriplets << std::endl;


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

  std::string outNtuName = "selectedEvents_Kee_"+dataset+"_BPHRun"+BPHRun;
  if(!isEleFinalState) outNtuName = "selectedEvents_Kmumu_"+dataset+"_BPHRun"+BPHRun;
  if(JOBid != "-1") outNtuName = outSelectedNTU+"/"+outNtuName+"_JOB_"+JOBid;
  outNtuName += ".root";
  gROOT->cd();
  TFile* newFile;
  TTree* newT;

  Float_t Kstll_pt;
  Float_t Kstll_mass;
  Float_t Kstll_cosAlpha;
  Float_t Kstll_CL_vtx;
  Float_t Kstll_Lxy;

  Float_t Kstll_l1_pt;
  Float_t Kstll_l1_eta;
  Float_t Kstll_l1_phi;
  
  Float_t Kstll_l2_pt;
  Float_t Kstll_l2_eta;
  Float_t Kstll_l2_phi;
  
  Float_t Kstll_llRefitmass;

  Float_t Kstll_k_pt;
  Float_t Kstll_k_eta;
  Float_t Kstll_k_phi;

  if(saveOUTntu){
    newFile = new TFile(outNtuName.c_str(),"recreate");
    newT = new TTree("selectedEvents", "");

    newT->Branch("Kstll_pt", &Kstll_pt, "Kstll_pt/F");         
    newT->Branch("Kstll_mass", &Kstll_mass, "Kstll_mass/F");    
    newT->Branch("Kstll_cosAlpha", &Kstll_cosAlpha, "Kstll_cosAlpha/F");
    newT->Branch("Kstll_CL_vtx", &Kstll_CL_vtx, "Kstll_CL_vtx/F");
    newT->Branch("Kstll_Lxy", &Kstll_Lxy, "Kstll_Lxy/F");

    newT->Branch("Kstll_l1_pt", &Kstll_l1_pt, "Kstll_l1_pt/F");
    newT->Branch("Kstll_l1_eta", &Kstll_l1_eta, "Kstll_l1_eta/F");
    newT->Branch("Kstll_l1_phi", &Kstll_l1_phi, "Kstll_l1_phi/F");
    
    newT->Branch("Kstll_l2_pt", &Kstll_l2_pt, "Kstll_l2_pt/F");
    newT->Branch("Kstll_l2_eta", &Kstll_l2_eta, "Kstll_l2_eta/F");
    newT->Branch("Kstll_l2_phi", &Kstll_l2_phi, "Kstll_l2_phi/F");
    
    newT->Branch("Kstll_llRefitmass", &Kstll_llRefitmass, "Kstll_llRefitmass/F");
    
    newT->Branch("Kstll_k_pt", &Kstll_k_pt, "Kstll_k_pt/F");
    newT->Branch("Kstll_k_eta", &Kstll_k_eta, "Kstll_k_eta/F");
    newT->Branch("Kstll_k_phi", &Kstll_k_phi, "Kstll_k_phi/F");
  }

  std::cout << " >>> so far so good " << std::endl;

  UInt_t run = 0;
  UInt_t lumi = 0;
  ULong64_t event = 0;
  
  //float nnBMX = -1;
    
  std::vector<int>* BToKstll_order_index = 0;
  
  std::vector<int>* Muon_tag_index = 0;
  int Muon_probe_index = -1;

  float BToKstll_pt[kBToKstllMax];
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
  
  float lep_pt[kMuonMax];
  float lep_eta[kMuonMax];
  float lep_phi[kMuonMax];
  
  int   BToKstll_gen_index = -1;  
  
  t1->SetBranchStatus("*", 0);


  t1->SetBranchStatus("run", 1);                        t1->SetBranchAddress("run", &run);
  t1->SetBranchStatus("luminosityBlock", 1);            t1->SetBranchAddress("luminosityBlock", &lumi);
  t1->SetBranchStatus("event", 1);                      t1->SetBranchAddress("event", &event);

  t1->SetBranchStatus("ProbeTracks_pt", 1);             t1->SetBranchAddress("ProbeTracks_pt", &ProbeTracks_pt);
  t1->SetBranchStatus("ProbeTracks_eta", 1);            t1->SetBranchAddress("ProbeTracks_eta", &ProbeTracks_eta);
  t1->SetBranchStatus("ProbeTracks_phi", 1);            t1->SetBranchAddress("ProbeTracks_phi", &ProbeTracks_phi);  
  
  if(dataset == "MC"){ 
    t1->SetBranchStatus("Muon_probe_index", 1);         t1->SetBranchAddress("Muon_probe_index", &Muon_probe_index);
    t1->SetBranchStatus("BToKstll_gen_index", 1);       t1->SetBranchAddress("BToKstll_gen_index", &BToKstll_gen_index);
  }
  
  
  if(isEleFinalState){    
      
    t1->SetBranchStatus("BToKstll_order_index_KEE", 1); t1->SetBranchAddress("BToKstll_order_index_KEE", &BToKstll_order_index);
    t1->SetBranchStatus("Muon_tag_index_KEE", 1);       t1->SetBranchAddress("Muon_tag_index_KEE", &Muon_tag_index);
    t1->SetBranchStatus("BToKEE_l1Idx", 1);             t1->SetBranchAddress("BToKEE_l1Idx", &lep1_Idx);
    t1->SetBranchStatus("BToKEE_l2Idx", 1);             t1->SetBranchAddress("BToKEE_l2Idx", &lep2_Idx);
    t1->SetBranchStatus("BToKEE_kIdx", 1);              t1->SetBranchAddress("BToKEE_kIdx", &kaon_Idx);
    t1->SetBranchStatus("BToKEE_pt", 1);                t1->SetBranchAddress("BToKEE_pt", &BToKstll_pt);
    t1->SetBranchStatus("BToKEE_mass", 1);              t1->SetBranchAddress("BToKEE_mass", &BToKstll_mass);
    t1->SetBranchStatus("BToKEE_cos2D", 1);             t1->SetBranchAddress("BToKEE_cos2D", &BToKstll_cos2D);
    t1->SetBranchStatus("BToKEE_svprob", 1);            t1->SetBranchAddress("BToKEE_svprob", &BToKstll_svprob);
    t1->SetBranchStatus("BToKEE_l_xy", 1);              t1->SetBranchAddress("BToKEE_l_xy", &BToKstll_l_xy);
    t1->SetBranchStatus("BToKEE_l_xy_unc", 1);          t1->SetBranchAddress("BToKEE_l_xy_unc", &BToKstll_l_xy_unc);
    t1->SetBranchStatus("BToKEE_mll_raw", 1);           t1->SetBranchAddress("BToKEE_mll_raw", &BToKstll_mll);    
    t1->SetBranchStatus("Electron_pt", 1);              t1->SetBranchAddress("Electron_pt", &lep_pt);
    t1->SetBranchStatus("Electron_eta", 1);             t1->SetBranchAddress("Electron_eta", &lep_eta);
    t1->SetBranchStatus("Electron_phi", 1);             t1->SetBranchAddress("Electron_phi", &lep_phi);
  }
  else{    
      
    t1->SetBranchStatus("BToKstll_order_index_KMuMu", 1);  t1->SetBranchAddress("BToKstll_order_index_KMuMu", &BToKstll_order_index);
    t1->SetBranchStatus("Muon_tag_index_KMuMu", 1);        t1->SetBranchAddress("Muon_tag_index_KMuMu", &Muon_tag_index);
    t1->SetBranchStatus("BToKMuMu_l1Idx", 1);              t1->SetBranchAddress("BToKMuMu_l1Idx", &lep1_Idx);
    t1->SetBranchStatus("BToKMuMu_l2Idx", 1);              t1->SetBranchAddress("BToKMuMu_l2Idx", &lep2_Idx);
    t1->SetBranchStatus("BToKMuMu_kIdx", 1);               t1->SetBranchAddress("BToKMuMu_kIdx", &kaon_Idx);    
    t1->SetBranchStatus("BToKMuMu_pt", 1);                 t1->SetBranchAddress("BToKMuMu_pt", &BToKstll_pt);
    t1->SetBranchStatus("BToKMuMu_mass", 1);               t1->SetBranchAddress("BToKMuMu_mass", &BToKstll_mass);
    t1->SetBranchStatus("BToKMuMu_cos2D", 1);              t1->SetBranchAddress("BToKMuMu_cos2D", &BToKstll_cos2D);
    t1->SetBranchStatus("BToKMuMu_svprob", 1);             t1->SetBranchAddress("BToKMuMu_svprob", &BToKstll_svprob);
    t1->SetBranchStatus("BToKMuMu_l_xy", 1);               t1->SetBranchAddress("BToKMuMu_l_xy", &BToKstll_l_xy);
    t1->SetBranchStatus("BToKMuMu_l_xy_unc", 1);           t1->SetBranchAddress("BToKMuMu_l_xy_unc", &BToKstll_l_xy_unc);
    t1->SetBranchStatus("BToKMuMu_mll_raw", 1);            t1->SetBranchAddress("BToKMuMu_mll_raw", &BToKstll_mll);
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


  std::string outName = "outMassHistos_Kee_"+dataset+"_BPHRun"+BPHRun;
  if(!isEleFinalState) outName = "outMassHistos_Kmumu_"+dataset+"_BPHRun"+BPHRun;
  if(JOBid != "-1") outName = outputFolder; // + "/" +outName + "_JOB_"+JOBid;
  else  outName += ".root";
  TFile outMassHistos(outName.c_str(), "recreate");

  ///histos: 1 per bin plus inclusive
  TH1F* hAlpha[7];
  TH1F* hCLVtx[7];
  TH1F* hLxy[7];
  TH1F* hKaonpt[7];
  TH1F* hLep1pt[7];
  TH1F* hLep2pt[7];  
  TH1F* hLep1pt_EB[7];
  TH1F* hLep2pt_EB[7];
  TH1F* hLep1pt_EE[7];
  TH1F* hLep2pt_EE[7];
  TH1F* hBpt[7];
  TH1F* hllRefitMass[7];  
  TH2F* hllRefitMass_vs_Bmass[7];
  TH1F* hBmass[7];
  TH2F* pTele2_vs_pTele1[7];

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

    hllRefitMass[ij] = new TH1F(Form("hllRefitMass_%d", ij), "", 750, 0., 15.);
    hllRefitMass[ij]->Sumw2();
    hllRefitMass[ij]->SetLineColor(kRed);
    hllRefitMass[ij]->SetLineWidth(2);

    hllRefitMass_vs_Bmass[ij] = new TH2F(Form("hllRefitMass_vs_Bmass_%d", ij), "", 500, 0., 15., 500, 0., 15.);
    hllRefitMass_vs_Bmass[ij]->Sumw2();
    hllRefitMass_vs_Bmass[ij]->SetMarkerColor(kRed);
    hllRefitMass_vs_Bmass[ij]->SetMarkerStyle(20);

    hBmass[ij] = new TH1F(Form("Bmass_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass[ij]->Sumw2();
    hBmass[ij]->SetLineColor(kRed);
    hBmass[ij]->SetLineWidth(2);  

    hKaonpt[ij] = new TH1F(Form("hKaonpt_%d", ij), "", 100, 0., 10.);
    hKaonpt[ij]->Sumw2();
    hKaonpt[ij]->SetLineColor(kRed);
    hKaonpt[ij]->SetLineWidth(2);    
    
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

    hBpt[ij] = new TH1F(Form("hBpt_%d", ij), "", 200, 0., 100.);
    hBpt[ij]->Sumw2();
    hBpt[ij]->SetLineColor(kRed);
    hBpt[ij]->SetLineWidth(2);
	  
    pTele2_vs_pTele1[ij] = new TH2F(Form("pTele2_vs_pTele1_%d", ij), "", 300, 0., 30., 300, 0., 30.);
    pTele2_vs_pTele1[ij]->Sumw2();
    pTele2_vs_pTele1[ij]->SetMarkerColor(kRed);
    pTele2_vs_pTele1[ij]->SetMarkerStyle(20);
  }


  float nEv_muonTag[1] = {0.};
  float nEv_recoCand[1] = {0.};
  float nEv_chargeSel[1] = {0.};
  float nEv_selected[7] = {0.};


  if(nMaxEvents == -1) nMaxEvents = nEvts;
  for(int iEvt = 0; iEvt<nMaxEvents; ++iEvt){
    if(iEvt%500000 == 0) std::cout << " >>> processing event " << iEvt << " " << 1.*iEvt/nEvts*100. << std::endl;

    t1->GetEntry(iEvt);
    
    int muon_tag_index_event = -1;
    int triplet_sel_index = -1;
    
    
    std::vector<int> triplets_Idx_vec;
    bool atLeastOneTriplet = false;
    

    for(unsigned int iP = 0; iP < BToKstll_order_index->size(); ++iP){//iP loop	
	
      triplet_sel_index = BToKstll_order_index->at(iP);
      muon_tag_index_event = Muon_tag_index->at(triplet_sel_index);
	
      if(muon_tag_index_event == -1 || triplet_sel_index == -1) continue;
      if(dataset == "MC" && Muon_probe_index == -1) continue;
      if(dataset == "MC" && triplet_sel_index != BToKstll_gen_index) continue;
    
      atLeastOneTriplet = true;
    
      triplets_Idx_vec.push_back(triplet_sel_index);
      if( triplets_Idx_vec.size() == nMaxTriplets ){
        break;
      }
    
    }//iP loop

      
    if(!atLeastOneTriplet) continue;

    
    
    for(unsigned int i_tr = 0; i_tr < triplets_Idx_vec.size(); ++i_tr){//goodTripletsFound
        
      int tr_Idx = triplets_Idx_vec.at(i_tr); 
      int l1_Idx = lep1_Idx[tr_Idx];
      int l2_Idx = lep2_Idx[tr_Idx];
      int k_Idx = kaon_Idx[tr_Idx];

      //misleading: not refitted for Kmumu and Kee
      float llInvRefitMass = BToKstll_mll[tr_Idx];
      int massBin = -1;
      for(unsigned int kl=0; kl<llMassBoundary.size()-1; ++kl){
	if(llInvRefitMass >= llMassBoundary[kl] && llInvRefitMass < llMassBoundary[kl+1]){
	  massBin = kl;
	  break;
	}
      }

      if(massBin != -1) ++nEv_selected[massBin];
      else ++nEv_selected[6];

    
      //save output ntuple with selected events fro final plots
      if(saveOUTntu){
        
	Kstll_pt = BToKstll_pt[tr_Idx];     
	Kstll_mass = BToKstll_mass[tr_Idx];
	Kstll_cosAlpha = BToKstll_cos2D[tr_Idx];
	Kstll_CL_vtx = BToKstll_svprob[tr_Idx];
	Kstll_Lxy = BToKstll_l_xy[tr_Idx]/sqrt(BToKstll_l_xy_unc[tr_Idx]);
    
	Kstll_l1_pt = lep_pt[l1_Idx];	
	Kstll_l1_eta = lep_eta[l1_Idx];
	Kstll_l1_phi = lep_phi[l1_Idx];
    
	Kstll_l2_pt = lep_pt[l2_Idx];
	Kstll_l2_eta = lep_eta[l2_Idx];
	Kstll_l2_phi = lep_phi[l2_Idx];
	
	Kstll_llRefitmass = llInvRefitMass;

	Kstll_k_pt = ProbeTracks_pt[k_Idx];
	Kstll_k_eta = ProbeTracks_eta[k_Idx];
	Kstll_k_phi = ProbeTracks_phi[k_Idx];

	newT->Fill();
      }
      //end output Nutuple

    
      //histograms for each mass bin
      if(massBin != -1){
	hAlpha[massBin]->Fill(BToKstll_cos2D[tr_Idx]);
	hCLVtx[massBin]->Fill(BToKstll_svprob[tr_Idx]);
	hLxy[massBin]->Fill(BToKstll_l_xy[tr_Idx]/sqrt(BToKstll_l_xy_unc[tr_Idx]));
	hKaonpt[massBin]->Fill(ProbeTracks_pt[k_Idx]);
	hBpt[massBin]->Fill(BToKstll_pt[tr_Idx]);

	hLep1pt[massBin]->Fill(lep_pt[l1_Idx]);
	hLep2pt[massBin]->Fill(lep_pt[l2_Idx]);
	if(std::abs(lep_eta[l1_Idx]) < 1.47) hLep1pt_EB[massBin]->Fill(lep_pt[l1_Idx]);
	else hLep1pt_EE[massBin]->Fill(lep_pt[l1_Idx]);
	if(std::abs(lep_eta[l2_Idx]) < 1.47) hLep2pt_EB[massBin]->Fill(lep_pt[l2_Idx]);
	else hLep2pt_EE[massBin]->Fill(lep_pt[l2_Idx]);

	hllRefitMass[massBin]->Fill(llInvRefitMass);
	hllRefitMass_vs_Bmass[massBin]->Fill(BToKstll_mass[tr_Idx], llInvRefitMass);
      
	hBmass[massBin]->Fill(BToKstll_mass[tr_Idx]);
      
	pTele2_vs_pTele1[massBin]->Fill(lep_pt[l1_Idx], lep_pt[l2_Idx]);
      }
    
      //histograms inclusive over all m(ll)
      hllRefitMass[6]->Fill(llInvRefitMass);
      hllRefitMass_vs_Bmass[6]->Fill(BToKstll_mass[tr_Idx], llInvRefitMass);
      hBmass[6]->Fill(BToKstll_mass[tr_Idx]);

      hAlpha[6]->Fill(BToKstll_cos2D[tr_Idx]);
      hCLVtx[6]->Fill(BToKstll_svprob[tr_Idx]);
      hLxy[6]->Fill(BToKstll_l_xy[tr_Idx]/sqrt(BToKstll_l_xy_unc[tr_Idx]));
      hKaonpt[6]->Fill(ProbeTracks_pt[k_Idx]);
      hBpt[6]->Fill(BToKstll_pt[tr_Idx]);

      hLep1pt[6]->Fill(lep_pt[l1_Idx]);
      hLep2pt[6]->Fill(lep_pt[l2_Idx]);

      if(std::abs(lep_eta[l1_Idx]) < 1.47) hLep1pt_EB[6]->Fill(lep_pt[l1_Idx]);
      else hLep1pt_EE[6]->Fill(lep_pt[l1_Idx]);
      if(std::abs(lep_eta[l2_Idx]) < 1.47) hLep2pt_EB[6]->Fill(lep_pt[l2_Idx]);
      else hLep2pt_EE[6]->Fill(lep_pt[l2_Idx]);
    
      pTele2_vs_pTele1[6]->Fill(lep_pt[l1_Idx], lep_pt[l2_Idx]);
      
    }//goodTripletsFound
  }//loop over events

  
  if(saveOUTntu){
    newFile->cd();
    newT->Write();
    newFile->Close(); 
  }

  outMassHistos.cd();
  
  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<7; ++ij){
    hAlpha[ij]->Write(hAlpha[ij]->GetName());
    hCLVtx[ij]->Write(hCLVtx[ij]->GetName());
    hLxy[ij]->Write(hLxy[ij]->GetName());
    hKaonpt[ij]->Write(hKaonpt[ij]->GetName());
    hLep1pt[ij]->Write(hLep1pt[ij]->GetName());
    hLep2pt[ij]->Write(hLep2pt[ij]->GetName());
    hLep1pt_EB[ij]->Write(hLep1pt_EB[ij]->GetName());
    hLep1pt_EE[ij]->Write(hLep1pt_EE[ij]->GetName());
    hLep2pt_EB[ij]->Write(hLep2pt_EB[ij]->GetName());
    hLep2pt_EE[ij]->Write(hLep2pt_EE[ij]->GetName());

    hBpt[ij]->Write(hBpt[ij]->GetName());
    std::cout << " >>> hBpt[ij]->GetName() = " << hBpt[ij]->GetName() << " entries = " << hBpt[ij]->GetEntries() << std::endl;

    hllRefitMass[ij]->Write(hllRefitMass[ij]->GetName());
    hllRefitMass_vs_Bmass[ij]->Write(hllRefitMass_vs_Bmass[ij]->GetName());
    hBmass[ij]->Write(hBmass[ij]->GetName());

    pTele2_vs_pTele1[ij]->Write(pTele2_vs_pTele1[ij]->GetName());

    if(ij > 5) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
	      << " \n \t muonTag = " << nEv_muonTag[0] << " \t recoEvts = " << nEv_recoCand[0]
              << " \n \t chargeSel = " << nEv_chargeSel[0] << " \t selected Events = " << nEv_selected[ij] << std::endl;

  }
  outMassHistos.Close();
  
}  
