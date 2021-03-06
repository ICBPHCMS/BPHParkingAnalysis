//g++ -Wall -o analyzeCharged_fastDATA_Kstll `root-config --cflags --glibs` -lRooFitCore analyzeCharged_fastDATA_Kstll.cpp

//./analyzeCharged_fastDATA_Kstll --isEle (0,1) --dataset (-1, runA, runB, runD, MC) --run (1,2,3,..., -1) --typeSelection (tightCB) --ntupleList (list.txt) --JOBid (1,2..) --outputFolder ("outfolder") --nMaxEvents (-1, N) --saveSelectedNTU (1,0) --outSelectedNTU (path for selected ntuples) --testFile ("path")


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

const int kBToKstllMax = 50;
const int kMuonMax = 100;
const int kPFCandMax = 10000;
const int kGenPartMax = 10000;
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
  std::string typeSelection = "-1";
  std::string ntupleList = "-1";
  std::string JOBid = "-1";
  std::string outputFolder = "-1";
  int nMaxEvents = -1;
  int saveOUTntu = 0;
  std::string outSelectedNTU = "-1";
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
    if(std::string(argv[i]) == "--typeSelection") {
      if (i + 1 < argc) {
        typeSelection = argv[i+1];
        break;
      } else {
	std::cerr << " --typeSelection option requires one argument " << std::endl;
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

  std::cout << " isEleFinalState = " << isEleFinalState << " dataset = " << dataset << " BPHRun = " << BPHRun << " typeSelection = " << typeSelection
	    << " ntupleList = " << ntupleList << " JOBid = " << JOBid << " outputFolder = " << outputFolder
	    << " nMaxEvents = " << nMaxEvents << " saveOUTntu = " << saveOUTntu << " outSelectedNTU = " << outSelectedNTU << std::endl;


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
  Float_t Kstll_ctxy;

  Float_t Kstll_l1_pt;
  Float_t Kstll_l1_eta;
  Float_t Kstll_l1_phi;
  Float_t Kstll_l1_preFpt;
  Float_t Kstll_l1_preFeta;
  Float_t Kstll_l1_preFphi;
  Float_t Kstll_l1_pfRelIso03;
  Float_t Kstll_l1_pfRelIso04;
  Float_t Kstll_l1_dxy;
  Float_t Kstll_l1_dz;
  
  Float_t Kstll_l2_pt;
  Float_t Kstll_l2_eta;
  Float_t Kstll_l2_phi;
  Float_t Kstll_l2_preFpt;
  Float_t Kstll_l2_preFeta;
  Float_t Kstll_l2_preFphi;
  Float_t Kstll_l2_pfRelIso03;
  Float_t Kstll_l2_pfRelIso04;
  Float_t Kstll_l2_dxy;
  Float_t Kstll_l2_dz;
  
  Float_t Kstll_llRefitmass;

  Int_t   Kstll_k_charge;
  Float_t Kstll_k_pt;
  Float_t Kstll_k_eta;
  Float_t Kstll_k_phi;
  Float_t Kstll_k_DCASig;
  Float_t Kstll_k_dxy;
  Float_t Kstll_k_dz;

  if(saveOUTntu){
    newFile = new TFile(outNtuName.c_str(),"recreate");
    newT = new TTree("selectedEvents", "");

    newT->Branch("Kstll_pt", &Kstll_pt, "Kstll_pt/F");         
    newT->Branch("Kstll_mass", &Kstll_mass, "Kstll_mass/F");    
    newT->Branch("Kstll_cosAlpha", &Kstll_cosAlpha, "Kstll_cosAlpha/F");
    newT->Branch("Kstll_CL_vtx", &Kstll_CL_vtx, "Kstll_CL_vtx/F");
    newT->Branch("Kstll_Lxy", &Kstll_Lxy, "Kstll_Lxy/F");
    newT->Branch("Kstll_ctxy", &Kstll_ctxy, "Kstll_ctxy/F");

    newT->Branch("Kstll_l1_pt", &Kstll_l1_pt, "Kstll_l1_pt/F");
    newT->Branch("Kstll_l1_eta", &Kstll_l1_eta, "Kstll_l1_eta/F");
    newT->Branch("Kstll_l1_phi", &Kstll_l1_phi, "Kstll_l1_phi/F");
    newT->Branch("Kstll_l1_preFpt", &Kstll_l1_preFpt, "Kstll_l1_preFpt/F");
    newT->Branch("Kstll_l1_preFeta", &Kstll_l1_preFeta, "Kstll_l1_preFeta/F");
    newT->Branch("Kstll_l1_preFphi", &Kstll_l1_preFphi, "Kstll_l1_preFphi/F");
    newT->Branch("Kstll_l1_pfRelIso03", &Kstll_l1_pfRelIso03, "Kstll_l1_pfRelIso03/F");
    newT->Branch("Kstll_l1_pfRelIso04", &Kstll_l1_pfRelIso04, "Kstll_l1_pfRelIso04/F");
    newT->Branch("Kstll_l1_dxy", &Kstll_l1_dxy, "Kstll_l1_dxy/F");
    newT->Branch("Kstll_l1_dz", &Kstll_l1_dz, "Kstll_l1_dz/F");
    
    newT->Branch("Kstll_l2_pt", &Kstll_l2_pt, "Kstll_l2_pt/F");
    newT->Branch("Kstll_l2_eta", &Kstll_l2_eta, "Kstll_l2_eta/F");
    newT->Branch("Kstll_l2_phi", &Kstll_l2_phi, "Kstll_l2_phi/F");
    newT->Branch("Kstll_l2_preFpt", &Kstll_l2_preFpt, "Kstll_l2_preFpt/F");
    newT->Branch("Kstll_l2_preFeta", &Kstll_l2_preFeta, "Kstll_l2_preFeta/F");
    newT->Branch("Kstll_l2_preFphi", &Kstll_l2_preFphi, "Kstll_l2_preFphi/F");
    newT->Branch("Kstll_l2_pfRelIso03", &Kstll_l2_pfRelIso03, "Kstll_l2_pfRelIso03/F");
    newT->Branch("Kstll_l2_pfRelIso04", &Kstll_l2_pfRelIso04, "Kstll_l2_pfRelIso04/F");
    newT->Branch("Kstll_l2_dxy", &Kstll_l2_dxy, "Kstll_l2_dxy/F");
    newT->Branch("Kstll_l2_dz", &Kstll_l2_dz, "Kstll_l2_dz/F");
    
    newT->Branch("Kstll_llRefitmass", &Kstll_llRefitmass, "Kstll_llRefitmass/F");

    newT->Branch("Kstll_k_charge", &Kstll_k_charge, "Kstll_k_charge/I");    
    newT->Branch("Kstll_k_pt", &Kstll_k_pt, "Kstll_k_pt/F");
    newT->Branch("Kstll_k_eta", &Kstll_k_eta, "Kstll_k_eta/F");
    newT->Branch("Kstll_k_phi", &Kstll_k_phi, "Kstll_k_phi/F");
    newT->Branch("Kstll_k_DCASig", &Kstll_k_DCASig, "Kstll_k_DCASig/F");    
    newT->Branch("Kstll_k_dxy", &Kstll_k_dxy, "Kstll_k_dxy/F");
    newT->Branch("Kstll_k_dz", &Kstll_k_dz, "Kstll_k_dz/F");
  }

  std::cout << " >>> so far so good " << std::endl;

  UInt_t run = 0;
  UInt_t lumi = 0;
  ULong64_t event = 0;
  
  //float nnBMX = -1;
    
  int BToKstll_sel_index = -1;
  int BToKstll_llsel_index = -1;
  std::vector<int>* BToKstll_order_index = 0;
  
  std::vector<int>* Muon_tag_index = 0;
  int Muon_sel_index = -1;
  int Muon_probe_index = -1;

  int BToKstll_lep1_isLowPt[kBToKstllMax];
  int BToKstll_lep2_isLowPt[kBToKstllMax];
  int BToKstll_lep1_isPFLep[kBToKstllMax];
  int BToKstll_lep2_isPFLep[kBToKstllMax];
  int BToKstll_isLowPtEle[kBToKstllMax];
  int BToKstll_lep1_isPFCand[kBToKstllMax];
  int BToKstll_lep2_isPFCand[kBToKstllMax];
  int BToKstll_lep1_isLT[kBToKstllMax];
  int BToKstll_lep2_isLT[kBToKstllMax];
  
  float BToKstll_lep1_seedBDT_unbiased[kBToKstllMax];
  float BToKstll_lep1_seedBDT_ptbiased[kBToKstllMax];
  float BToKstll_lep2_seedBDT_unbiased[kBToKstllMax];
  float BToKstll_lep2_seedBDT_ptbiased[kBToKstllMax];

  float BToKstll_B_pt[kBToKstllMax];
  float BToKstll_B_mass[kBToKstllMax];
  float BToKstll_B_cosAlpha[kBToKstllMax];
  float BToKstll_B_CL_vtx[kBToKstllMax];
  float BToKstll_B_Lxy[kBToKstllMax];
  float BToKstll_B_ctxy[kBToKstllMax];
  
  int   BToKstll_lep1_charge[kBToKstllMax];
  int   BToKstll_lep1_index[kBToKstllMax];
  float BToKstll_lep1_pt[kBToKstllMax];
  float BToKstll_lep1_eta[kBToKstllMax];
  float BToKstll_lep1_phi[kBToKstllMax];
  
  int   BToKstll_lep2_charge[kBToKstllMax];
  int   BToKstll_lep2_index[kBToKstllMax];  
  float BToKstll_lep2_pt[kBToKstllMax];
  float BToKstll_lep2_eta[kBToKstllMax];
  float BToKstll_lep2_phi[kBToKstllMax];
  
  float BToKstll_ll_mass[kBToKstllMax];
  
  int   BToKstll_kaon_charge[kBToKstllMax];
  int   BToKstll_kaon_index[kBToKstllMax];
  float BToKstll_kaon_pt[kBToKstllMax];
  float BToKstll_kaon_eta[kBToKstllMax];
  float BToKstll_kaon_phi[kBToKstllMax];
  float BToKstll_kaon_DCASig[kBToKstllMax];

  float PFCand_dxy[kPFCandMax];
  float PFCand_dz[kPFCandMax];
  
  int   BToKstll_gen_index = -1;
  float BToKstll_gendR_lep1FromB;
  float BToKstll_gendR_lep2FromB;
  float BToKstll_gendR_KFromB;
  
  float Lepton_pfRelIso03[kMuonMax];
  float Lepton_pfRelIso04[kMuonMax];
  float Lepton_pt[kMuonMax];
  float Lepton_eta[kMuonMax];
  float Lepton_phi[kMuonMax];
  float Lepton_dxy[kMuonMax];
  float Lepton_dz[kMuonMax];
  
  
  t1->SetBranchStatus("*", 0);


  t1->SetBranchStatus("run", 1);                        t1->SetBranchAddress("run", &run);
  t1->SetBranchStatus("luminosityBlock", 1);            t1->SetBranchAddress("luminosityBlock", &lumi);
  t1->SetBranchStatus("event", 1);                      t1->SetBranchAddress("event", &event);
  
  t1->SetBranchStatus("BToKstll_sel_index", 1);         t1->SetBranchAddress("BToKstll_sel_index", &BToKstll_sel_index);
  t1->SetBranchStatus("BToKstll_llsel_index", 1);       t1->SetBranchAddress("BToKstll_llsel_index", &BToKstll_llsel_index);
  t1->SetBranchStatus("BToKstll_order_index", 1);       t1->SetBranchAddress("BToKstll_order_index", &BToKstll_order_index);

  t1->SetBranchStatus("Muon_tag_index", 1);             t1->SetBranchAddress("Muon_tag_index", &Muon_tag_index);
  t1->SetBranchStatus("Muon_sel_index", 1);             t1->SetBranchAddress("Muon_sel_index", &Muon_sel_index);
  
  t1->SetBranchStatus("BToKstll_lep1_isLowPt", 1);      t1->SetBranchAddress("BToKstll_lep1_isLowPt", &BToKstll_lep1_isLowPt);
  t1->SetBranchStatus("BToKstll_lep2_isLowPt", 1);      t1->SetBranchAddress("BToKstll_lep2_isLowPt", &BToKstll_lep2_isLowPt);
  t1->SetBranchStatus("BToKstll_lep1_isPFLep", 1);      t1->SetBranchAddress("BToKstll_lep1_isPFLep", &BToKstll_lep1_isPFLep);
  t1->SetBranchStatus("BToKstll_lep2_isPFLep", 1);      t1->SetBranchAddress("BToKstll_lep2_isPFLep", &BToKstll_lep2_isPFLep);
  t1->SetBranchStatus("BToKstll_isLowPtEle", 1);        t1->SetBranchAddress("BToKstll_isLowPtEle", &BToKstll_isLowPtEle);
  
  t1->SetBranchStatus("BToKstll_lep1_isPFCand", 1);      t1->SetBranchAddress("BToKstll_lep1_isPFCand", &BToKstll_lep1_isPFCand);
  t1->SetBranchStatus("BToKstll_lep2_isPFCand", 1);      t1->SetBranchAddress("BToKstll_lep2_isPFCand", &BToKstll_lep2_isPFCand);

  t1->SetBranchStatus("BToKstll_lep1_isLT", 1);      t1->SetBranchAddress("BToKstll_lep1_isLT", &BToKstll_lep1_isLT);
  t1->SetBranchStatus("BToKstll_lep2_isLT", 1);      t1->SetBranchAddress("BToKstll_lep2_isLT", &BToKstll_lep2_isLT);
  
  t1->SetBranchStatus("BToKstll_lep1_seedBDT_unbiased", 1);   t1->SetBranchAddress("BToKstll_lep1_seedBDT_unbiased", &BToKstll_lep1_seedBDT_unbiased);
  t1->SetBranchStatus("BToKstll_lep1_seedBDT_ptbiased", 1);   t1->SetBranchAddress("BToKstll_lep1_seedBDT_ptbiased", &BToKstll_lep1_seedBDT_ptbiased);
  t1->SetBranchStatus("BToKstll_lep2_seedBDT_unbiased", 1);   t1->SetBranchAddress("BToKstll_lep2_seedBDT_unbiased", &BToKstll_lep2_seedBDT_unbiased);
  t1->SetBranchStatus("BToKstll_lep2_seedBDT_ptbiased", 1);   t1->SetBranchAddress("BToKstll_lep2_seedBDT_ptbiased", &BToKstll_lep2_seedBDT_ptbiased);

  t1->SetBranchStatus("BToKstll_B_pt", 1);              t1->SetBranchAddress("BToKstll_B_pt", &BToKstll_B_pt);
  t1->SetBranchStatus("BToKstll_B_mass", 1);            t1->SetBranchAddress("BToKstll_B_mass", &BToKstll_B_mass);
  t1->SetBranchStatus("BToKstll_B_cosAlpha", 1);        t1->SetBranchAddress("BToKstll_B_cosAlpha", &BToKstll_B_cosAlpha);
  t1->SetBranchStatus("BToKstll_B_CL_vtx", 1);          t1->SetBranchAddress("BToKstll_B_CL_vtx", &BToKstll_B_CL_vtx);
  t1->SetBranchStatus("BToKstll_B_Lxy", 1);             t1->SetBranchAddress("BToKstll_B_Lxy", &BToKstll_B_Lxy);
  t1->SetBranchStatus("BToKstll_B_ctxy", 1);            t1->SetBranchAddress("BToKstll_B_ctxy", &BToKstll_B_ctxy);
  
  t1->SetBranchStatus("BToKstll_lep1_charge", 1);       t1->SetBranchAddress("BToKstll_lep1_charge", &BToKstll_lep1_charge);
  t1->SetBranchStatus("BToKstll_lep1_index", 1);        t1->SetBranchAddress("BToKstll_lep1_index", &BToKstll_lep1_index);
  t1->SetBranchStatus("BToKstll_lep1_pt", 1);           t1->SetBranchAddress("BToKstll_lep1_pt", &BToKstll_lep1_pt);
  t1->SetBranchStatus("BToKstll_lep1_eta", 1);          t1->SetBranchAddress("BToKstll_lep1_eta", &BToKstll_lep1_eta);
  t1->SetBranchStatus("BToKstll_lep1_phi", 1);          t1->SetBranchAddress("BToKstll_lep1_phi", &BToKstll_lep1_phi);
  
  t1->SetBranchStatus("BToKstll_lep2_charge", 1);       t1->SetBranchAddress("BToKstll_lep2_charge", &BToKstll_lep2_charge);
  t1->SetBranchStatus("BToKstll_lep2_index", 1);        t1->SetBranchAddress("BToKstll_lep2_index", &BToKstll_lep2_index);
  t1->SetBranchStatus("BToKstll_lep2_pt", 1);           t1->SetBranchAddress("BToKstll_lep2_pt", &BToKstll_lep2_pt);
  t1->SetBranchStatus("BToKstll_lep2_eta", 1);          t1->SetBranchAddress("BToKstll_lep2_eta", &BToKstll_lep2_eta);
  t1->SetBranchStatus("BToKstll_lep2_phi", 1);          t1->SetBranchAddress("BToKstll_lep2_phi", &BToKstll_lep2_phi);
  
  t1->SetBranchStatus("BToKstll_ll_mass", 1);           t1->SetBranchAddress("BToKstll_ll_mass", &BToKstll_ll_mass);          

  t1->SetBranchStatus("BToKstll_kaon_charge", 1);       t1->SetBranchAddress("BToKstll_kaon_charge", &BToKstll_kaon_charge);  
  t1->SetBranchStatus("BToKstll_kaon_index", 1);        t1->SetBranchAddress("BToKstll_kaon_index", &BToKstll_kaon_index);  
  t1->SetBranchStatus("BToKstll_kaon_pt", 1);           t1->SetBranchAddress("BToKstll_kaon_pt", &BToKstll_kaon_pt);
  t1->SetBranchStatus("BToKstll_kaon_eta", 1);          t1->SetBranchAddress("BToKstll_kaon_eta", &BToKstll_kaon_eta);
  t1->SetBranchStatus("BToKstll_kaon_phi", 1);          t1->SetBranchAddress("BToKstll_kaon_phi", &BToKstll_kaon_phi);
  t1->SetBranchStatus("BToKstll_kaon_DCASig", 1);       t1->SetBranchAddress("BToKstll_kaon_DCASig", &BToKstll_kaon_DCASig);

  t1->SetBranchStatus("PFCand_dxy", 1);                 t1->SetBranchAddress("PFCand_dxy", &PFCand_dxy);
  t1->SetBranchStatus("PFCand_dz", 1);                  t1->SetBranchAddress("PFCand_dz", &PFCand_dz);

  
  if(dataset == "MC"){ 
    t1->SetBranchStatus("Muon_probe_index", 1);       t1->SetBranchAddress("Muon_probe_index", &Muon_probe_index);
    t1->SetBranchStatus("BToKstll_gen_index", 1);       t1->SetBranchAddress("BToKstll_gen_index", &BToKstll_gen_index);
    t1->SetBranchStatus("BToKstll_gendR_lep1FromB", 1); t1->SetBranchAddress("BToKstll_gendR_lep1FromB", &BToKstll_gendR_lep1FromB);
    t1->SetBranchStatus("BToKstll_gendR_lep2FromB", 1); t1->SetBranchAddress("BToKstll_gendR_lep2FromB", &BToKstll_gendR_lep2FromB);
    t1->SetBranchStatus("BToKstll_gendR_KFromB", 1);    t1->SetBranchAddress("BToKstll_gendR_KFromB", &BToKstll_gendR_KFromB);
  }
  
  
  if(isEleFinalState){    
      
    t1->SetBranchStatus("Electron_pfRelIso03_all", 1);  t1->SetBranchAddress("Electron_pfRelIso03_all", &Lepton_pfRelIso03);
    t1->SetBranchStatus("Electron_pt", 1);              t1->SetBranchAddress("Electron_pt", &Lepton_pt);
    t1->SetBranchStatus("Electron_eta", 1);             t1->SetBranchAddress("Electron_eta", &Lepton_eta);
    t1->SetBranchStatus("Electron_phi", 1);             t1->SetBranchAddress("Electron_phi", &Lepton_phi);
    t1->SetBranchStatus("Electron_dxy", 1);             t1->SetBranchAddress("Electron_dxy", &Lepton_dxy);
    t1->SetBranchStatus("Electron_dz", 1);              t1->SetBranchAddress("Electron_dz", &Lepton_dz);
  }
  else{    
      
    t1->SetBranchStatus("Muon_pfRelIso03_all", 1);      t1->SetBranchAddress("Muon_pfRelIso03_all", &Lepton_pfRelIso03);
    t1->SetBranchStatus("Muon_pfRelIso04_all", 1);      t1->SetBranchAddress("Muon_pfRelIso04_all", &Lepton_pfRelIso04);
    t1->SetBranchStatus("Muon_pt", 1);                  t1->SetBranchAddress("Muon_pt", &Lepton_pt);
    t1->SetBranchStatus("Muon_eta", 1);                 t1->SetBranchAddress("Muon_eta", &Lepton_eta);
    t1->SetBranchStatus("Muon_phi", 1);                 t1->SetBranchAddress("Muon_phi", &Lepton_phi);
    t1->SetBranchStatus("Muon_dxy", 1);                 t1->SetBranchAddress("Muon_dxy", &Lepton_dxy);
    t1->SetBranchStatus("Muon_dz", 1);                  t1->SetBranchAddress("Muon_dz", &Lepton_dz);        
    
    //if(typeSelection == "NN_BkgR" || typeSelection == "NN_SigEff"){
    //    t1->SetBranchStatus("nnBMX", 1);                t1->SetBranchAddress("nnBMX", &nnBMX);
    //}
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
  TH1F* hDCASig[7];
  TH1F* hLxy[7];
  TH1F* hctxy[7];
  TH1F* hKaonpt[7];
  TH1F* hLep1pt[7];
  TH1F* hLep2pt[7];
  TH1F* hLep1pt_PFLep[7];
  TH1F* hLep2pt_PFLep[7];
  TH1F* hLep1pt_l1notPFLep[7];
  TH1F* hLep2pt_l2notPFLep[7];
  TH1F* hLep1pt_l1PFLep_l2notPFLep[7];
  TH1F* hLep2pt_l1PFLep_l2notPFLep[7];
  TH1F* hLep1pt_l1notPFLep_l2notPFLep[7];
  TH1F* hLep2pt_l1notPFLep_l2notPFLep[7];  
  TH1F* hLep1pt_PFCand[7];
  TH1F* hLep2pt_PFCand[7];
  TH1F* hLep1pt_LT[7];
  TH1F* hLep2pt_LT[7];  
  TH1F* hLep1pt_EB[7];
  TH1F* hLep2pt_EB[7];
  TH1F* hLep1pt_EE[7];
  TH1F* hLep2pt_EE[7];
  TH1F* hBpt[7];
  TH1F* hllRefitMass[7];  
  TH2F* hllRefitMass_vs_Bmass[7];
  TH1F* hBmass[7];
  TH1F* hBmass_llt[7];
  TH1F* hBmass_l1PFLep_l2notPFLep[7];
  TH1F* hBmass_l1notPFLep_l2PFLep[7];
  TH1F* hBmass_l1notPFLep_l2notPFLep[7];
  TH1F* hBmass_l1PFCand_l2PFCand[7];
  TH1F* hBmass_l1PFCand_l2notPFCand[7];
  TH1F* hBmass_l1notPFCand_l2PFCand[7];
  TH1F* hBmass_l1LT_l2LT[7];
  TH1F* hBmass_l1LT_l2notLT[7];
  TH1F* hBmass_l1notLT_l2LT[7];
  TH2F* BDTele1_vs_pTele1[7];
  TH2F* BDTele2_vs_pTele2[7];
  TH2F* BDTele2_vs_BDTele1[7];
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
 
    hDCASig[ij] = new TH1F(Form("hDCASig_%d", ij), "", 100, -50., 50.);
    hDCASig[ij]->Sumw2();
    hDCASig[ij]->SetLineColor(kRed);
    hDCASig[ij]->SetLineWidth(2);

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

    hBmass_llt[ij] = new TH1F(Form("Bmass_llt_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_llt[ij]->Sumw2();
    hBmass_llt[ij]->SetLineColor(kRed);
    hBmass_llt[ij]->SetLineWidth(2);
    
    hBmass_l1PFLep_l2notPFLep[ij] = new TH1F(Form("Bmass_l1PFLep_l2notPFLep_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1PFLep_l2notPFLep[ij]->Sumw2();
    hBmass_l1PFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hBmass_l1PFLep_l2notPFLep[ij]->SetLineWidth(2);
    
    hBmass_l1notPFLep_l2PFLep[ij] = new TH1F(Form("Bmass_l1notPFLep_l2PFLep_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1notPFLep_l2PFLep[ij]->Sumw2();
    hBmass_l1notPFLep_l2PFLep[ij]->SetLineColor(kRed);
    hBmass_l1notPFLep_l2PFLep[ij]->SetLineWidth(2);
    
    hBmass_l1notPFLep_l2notPFLep[ij] = new TH1F(Form("Bmass_l1notPFLep_l2notPFLep_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1notPFLep_l2notPFLep[ij]->Sumw2();
    hBmass_l1notPFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hBmass_l1notPFLep_l2notPFLep[ij]->SetLineWidth(2);   
    
    hBmass_l1PFCand_l2PFCand[ij] = new TH1F(Form("Bmass_l1PFCand_l2PFCand_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1PFCand_l2PFCand[ij]->Sumw2();
    hBmass_l1PFCand_l2PFCand[ij]->SetLineColor(kRed);
    hBmass_l1PFCand_l2PFCand[ij]->SetLineWidth(2);
    
    hBmass_l1PFCand_l2notPFCand[ij] = new TH1F(Form("Bmass_l1PFCand_l2notPFCand_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1PFCand_l2notPFCand[ij]->Sumw2();
    hBmass_l1PFCand_l2notPFCand[ij]->SetLineColor(kRed);
    hBmass_l1PFCand_l2notPFCand[ij]->SetLineWidth(2);
    
    hBmass_l1notPFCand_l2PFCand[ij] = new TH1F(Form("Bmass_l1notPFCand_l2PFCand_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1notPFCand_l2PFCand[ij]->Sumw2();
    hBmass_l1notPFCand_l2PFCand[ij]->SetLineColor(kRed);
    hBmass_l1notPFCand_l2PFCand[ij]->SetLineWidth(2);
    
    hBmass_l1LT_l2LT[ij] = new TH1F(Form("Bmass_l1LT_l2LT_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1LT_l2LT[ij]->Sumw2();
    hBmass_l1LT_l2LT[ij]->SetLineColor(kRed);
    hBmass_l1LT_l2LT[ij]->SetLineWidth(2);

    hBmass_l1LT_l2notLT[ij] = new TH1F(Form("Bmass_l1LT_l2notLT_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1LT_l2notLT[ij]->Sumw2();
    hBmass_l1LT_l2notLT[ij]->SetLineColor(kRed);
    hBmass_l1LT_l2notLT[ij]->SetLineWidth(2);
    
    hBmass_l1notLT_l2LT[ij] = new TH1F(Form("Bmass_l1notLT_l2LT_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass_l1notLT_l2LT[ij]->Sumw2();
    hBmass_l1notLT_l2LT[ij]->SetLineColor(kRed);
    hBmass_l1notLT_l2LT[ij]->SetLineWidth(2);   

    hctxy[ij] = new TH1F(Form("hctxy_%d", ij), "", 1000, 0., 10.);
    hctxy[ij]->Sumw2();
    hctxy[ij]->SetLineColor(kRed);
    hctxy[ij]->SetLineWidth(2);

    hKaonpt[ij] = new TH1F(Form("hKaonpt_%d", ij), "", 100, 0., 10.);
    hKaonpt[ij]->Sumw2();
    hKaonpt[ij]->SetLineColor(kRed);
    hKaonpt[ij]->SetLineWidth(2);    
    
    hLep1pt[ij] = new TH1F(Form("hLep1pt_%d", ij), "", 100, 0., 10.);
    hLep1pt[ij]->Sumw2();
    hLep1pt[ij]->SetLineColor(kRed);
    hLep1pt[ij]->SetLineWidth(2);
    
    hLep1pt_PFLep[ij] = new TH1F(Form("hLep1pt_PFLep_%d", ij), "", 100, 0., 10.);
    hLep1pt_PFLep[ij]->Sumw2();
    hLep1pt_PFLep[ij]->SetLineColor(kRed);
    hLep1pt_PFLep[ij]->SetLineWidth(2);

    hLep1pt_l1notPFLep[ij] = new TH1F(Form("hLep1pt_l1notPFLep_%d", ij), "", 100, 0., 10.);
    hLep1pt_l1notPFLep[ij]->Sumw2();
    hLep1pt_l1notPFLep[ij]->SetLineColor(kRed);
    hLep1pt_l1notPFLep[ij]->SetLineWidth(2);
    
    hLep1pt_l1PFLep_l2notPFLep[ij] = new TH1F(Form("hLep1pt_l1PFLep_l2notPFLep_%d", ij), "", 100, 0., 10.);
    hLep1pt_l1PFLep_l2notPFLep[ij]->Sumw2();
    hLep1pt_l1PFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hLep1pt_l1PFLep_l2notPFLep[ij]->SetLineWidth(2);
    
    hLep1pt_l1notPFLep_l2notPFLep[ij] = new TH1F(Form("hLep1pt_l1notPFLep_l2notPFLep_%d", ij), "", 100, 0., 10.);
    hLep1pt_l1notPFLep_l2notPFLep[ij]->Sumw2();
    hLep1pt_l1notPFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hLep1pt_l1notPFLep_l2notPFLep[ij]->SetLineWidth(2);    
    
    hLep1pt_PFCand[ij] = new TH1F(Form("hLep1pt_PFCand_%d", ij), "", 100, 0., 10.);
    hLep1pt_PFCand[ij]->Sumw2();
    hLep1pt_PFCand[ij]->SetLineColor(kRed);
    hLep1pt_PFCand[ij]->SetLineWidth(2);

    hLep1pt_LT[ij] = new TH1F(Form("hLep1pt_LT_%d", ij), "", 100, 0., 10.);
    hLep1pt_LT[ij]->Sumw2();
    hLep1pt_LT[ij]->SetLineColor(kRed);
    hLep1pt_LT[ij]->SetLineWidth(2);    

    hLep2pt[ij] = new TH1F(Form("hLep2pt_%d", ij), "", 100, 0., 10.);
    hLep2pt[ij]->Sumw2();
    hLep2pt[ij]->SetLineColor(kRed);
    hLep2pt[ij]->SetLineWidth(2);

    hLep2pt_PFLep[ij] = new TH1F(Form("hLep2pt_PFLep_%d", ij), "", 100, 0., 10.);
    hLep2pt_PFLep[ij]->Sumw2();
    hLep2pt_PFLep[ij]->SetLineColor(kRed);
    hLep2pt_PFLep[ij]->SetLineWidth(2);
    
    hLep2pt_l2notPFLep[ij] = new TH1F(Form("hLep2pt_l2notPFLep_%d", ij), "", 100, 0., 10.);
    hLep2pt_l2notPFLep[ij]->Sumw2();
    hLep2pt_l2notPFLep[ij]->SetLineColor(kRed);
    hLep2pt_l2notPFLep[ij]->SetLineWidth(2);    

    hLep2pt_l1PFLep_l2notPFLep[ij] = new TH1F(Form("hLep2pt_l1PFLep_l2notPFLep_%d", ij), "", 100, 0., 10.);
    hLep2pt_l1PFLep_l2notPFLep[ij]->Sumw2();
    hLep2pt_l1PFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hLep2pt_l1PFLep_l2notPFLep[ij]->SetLineWidth(2);
    
    hLep2pt_l1notPFLep_l2notPFLep[ij] = new TH1F(Form("hLep2pt_l1notPFLep_l2notPFLep_%d", ij), "", 100, 0., 10.);
    hLep2pt_l1notPFLep_l2notPFLep[ij]->Sumw2();
    hLep2pt_l1notPFLep_l2notPFLep[ij]->SetLineColor(kRed);
    hLep2pt_l1notPFLep_l2notPFLep[ij]->SetLineWidth(2);
    
    hLep2pt_PFCand[ij] = new TH1F(Form("hLep2pt_PFCand_%d", ij), "", 100, 0., 10.);
    hLep2pt_PFCand[ij]->Sumw2();
    hLep2pt_PFCand[ij]->SetLineColor(kRed);
    hLep2pt_PFCand[ij]->SetLineWidth(2);

    hLep2pt_LT[ij] = new TH1F(Form("hLep2pt_LT_%d", ij), "", 100, 0., 10.);
    hLep2pt_LT[ij]->Sumw2();
    hLep2pt_LT[ij]->SetLineColor(kRed);
    hLep2pt_LT[ij]->SetLineWidth(2);

    //
    hLep1pt_EB[ij] = new TH1F(Form("hLep1pt_EB_%d", ij), "", 100, 0., 10.);
    hLep1pt_EB[ij]->Sumw2();
    hLep1pt_EB[ij]->SetLineColor(kRed);
    hLep1pt_EB[ij]->SetLineWidth(2);

    hLep2pt_EB[ij] = new TH1F(Form("hLep2pt_EB_%d", ij), "", 100, 0., 10.);
    hLep2pt_EB[ij]->Sumw2();
    hLep2pt_EB[ij]->SetLineColor(kRed);
    hLep2pt_EB[ij]->SetLineWidth(2);

    //
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
    
    BDTele1_vs_pTele1[ij] = new TH2F(Form("BDTele1_vs_pTele1_%d", ij), "", 300, 0., 30., 60, 0., 15.);
    BDTele1_vs_pTele1[ij]->Sumw2();
    BDTele1_vs_pTele1[ij]->SetMarkerColor(kRed);
    BDTele1_vs_pTele1[ij]->SetMarkerStyle(20);    
    
    BDTele2_vs_pTele2[ij] = new TH2F(Form("BDTele2_vs_pTele2_%d", ij), "", 300, 0., 30., 60, 0., 15.);
    BDTele2_vs_pTele2[ij]->Sumw2();
    BDTele2_vs_pTele2[ij]->SetMarkerColor(kRed);
    BDTele2_vs_pTele2[ij]->SetMarkerStyle(20);
    
    BDTele2_vs_BDTele1[ij] = new TH2F(Form("BDTele2_vs_BDTele1_%d", ij), "", 60, 0., 15., 60, 0., 15.);
    BDTele2_vs_BDTele1[ij]->Sumw2();
    BDTele2_vs_BDTele1[ij]->SetMarkerColor(kRed);
    BDTele2_vs_BDTele1[ij]->SetMarkerStyle(20);
	  
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
    bool isllt = false;         
    bool l1PFCand_l2PFCand = false;
    bool l1PFCand_l2notPFCand = false;
    bool l1notPFCand_l2PFCand = false;
    bool l1LT_l2LT = false;
    bool l1LT_l2notLT = false;
    bool l1notLT_l2LT = false;
    bool l1PFLep_l2notPFLep = false;
    bool l1notPFLep_l2PFLep = false;
    bool l1notPFLep_l2notPFLep = false;
    bool goodTripletFound = false;

    //choose the best vtxCL candidate provided it survives the selections
    //priority is given to the llt triplets in case of LTT
    for(int iL=0; iL<3; ++iL){
      if(goodTripletFound) break;

      for(unsigned int iP = 0; iP < BToKstll_order_index->size(); ++iP){	
	
	triplet_sel_index = BToKstll_order_index->at(iP);
	muon_tag_index_event = Muon_tag_index->at(triplet_sel_index);
	
	if( iL == 0 && (BToKstll_lep1_isPFLep[triplet_sel_index] != 1 || BToKstll_lep2_isPFLep[triplet_sel_index] != 1) ) continue;
	else if( iL == 1 && ( (BToKstll_lep1_isPFLep[triplet_sel_index] == 1 && BToKstll_lep2_isPFLep[triplet_sel_index] == 1) 
                         || (BToKstll_lep1_isPFLep[triplet_sel_index] != 1 && BToKstll_lep2_isPFLep[triplet_sel_index] != 1) ) ) continue;
    else if( iL == 2 && (BToKstll_lep1_isPFLep[triplet_sel_index] == 1 || BToKstll_lep2_isPFLep[triplet_sel_index] == 1) ) continue;
	
	if(muon_tag_index_event == -1 || triplet_sel_index == -1) continue;
	if(dataset == "MC" && Muon_probe_index == -1) continue;
	if(dataset == "MC" && triplet_sel_index != BToKstll_gen_index) continue;
	
	//analysis
	if(BToKstll_lep1_charge[triplet_sel_index]*BToKstll_lep2_charge[triplet_sel_index] > 0.) continue;
	
	//to synch. with Riccardo ~ tight selection                                                           
	if(typeSelection == "tightCB"){
	  if((BToKstll_kaon_pt[triplet_sel_index] < 1.5 || BToKstll_B_pt[triplet_sel_index] < 10.)) continue;
	  if(BToKstll_B_cosAlpha[triplet_sel_index] < 0.999) continue;
	  if(BToKstll_B_CL_vtx[triplet_sel_index] < 0.1) continue;
	  if(BToKstll_B_Lxy[triplet_sel_index] < 6.) continue;
	}

	isllt = (iL == 0) ? true : false;

    l1PFCand_l2PFCand = bool(BToKstll_lep1_isPFCand[triplet_sel_index] == 1 && BToKstll_lep2_isPFCand[triplet_sel_index]== 1);
    l1PFCand_l2notPFCand = bool(BToKstll_lep1_isPFCand[triplet_sel_index]== 1 && BToKstll_lep2_isPFCand[triplet_sel_index]!= 1);
    l1notPFCand_l2PFCand = bool(BToKstll_lep1_isPFCand[triplet_sel_index]!= 1 && BToKstll_lep2_isPFCand[triplet_sel_index]== 1);
    
    l1LT_l2LT = bool(BToKstll_lep1_isLT[triplet_sel_index] == 1 && BToKstll_lep2_isLT[triplet_sel_index]== 1);
    l1LT_l2notLT = bool(BToKstll_lep1_isLT[triplet_sel_index]== 1 && BToKstll_lep2_isLT[triplet_sel_index]!= 1);
    l1notLT_l2LT = bool(BToKstll_lep1_isLT[triplet_sel_index]!= 1 && BToKstll_lep2_isLT[triplet_sel_index]== 1);
    
    l1PFLep_l2notPFLep = bool(BToKstll_lep1_isPFLep[triplet_sel_index] == 1 && BToKstll_lep2_isPFLep[triplet_sel_index] != 1);
    l1notPFLep_l2PFLep = bool(BToKstll_lep1_isPFLep[triplet_sel_index] != 1 && BToKstll_lep2_isPFLep[triplet_sel_index] == 1);
    l1notPFLep_l2notPFLep = bool(BToKstll_lep1_isPFLep[triplet_sel_index] != 1 && BToKstll_lep2_isPFLep[triplet_sel_index] != 1);    
    
	goodTripletFound = true;
	break;
      }
    }

    if(!goodTripletFound) continue;
  
    /*         
    if(muon_tag_index_event == -1) continue;
    if(dataset == "MC" && Muon_probe_index == -1) continue;
    ++nEv_muonTag[0];
        
    if(triplet_sel_index == -1)continue;
    if(dataset == "MC" && triplet_sel_index != BToKstll_gen_index) continue;
    ++nEv_recoCand[0];
    
    //opposite sign leptons
    //expect eff 1 because opposite charge is already required at nanoAOD level
    if(BToKstll_lep1_charge[triplet_sel_index]*BToKstll_lep2_charge[triplet_sel_index] > 0.) continue;
    ++nEv_chargeSel[0];
    //to synch. with Riccardo ~ tight selection
    if(typeSelection == "tightCB"){
      if((BToKstll_kaon_pt[triplet_sel_index] < 1.5 || BToKstll_B_pt[triplet_sel_index] < 10.)) continue;
      if(BToKstll_B_cosAlpha[triplet_sel_index] < 0.999) continue;
      if(BToKstll_B_CL_vtx[triplet_sel_index] < 0.1) continue;
      if(BToKstll_B_Lxy[triplet_sel_index] < 6.) continue;
    }
    */


    //MVA selection to get same background rejection as with the cut base
    //if(typeSelection == "NN_BkgR" && nnBMX < 0.993) continue;

    //MVA selection to get same signal efficiency as with the cut base
    //if(typeSelection == "NN_SigEff" && nnBMX < 0.998) continue;



    //misleading: not refitted for Kmumu and Kee
    float llInvRefitMass = BToKstll_ll_mass[triplet_sel_index];
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
        
	Kstll_pt = BToKstll_B_pt[triplet_sel_index];     
	Kstll_mass = BToKstll_B_mass[triplet_sel_index];
	Kstll_cosAlpha = BToKstll_B_cosAlpha[triplet_sel_index];
	Kstll_CL_vtx = BToKstll_B_CL_vtx[triplet_sel_index];
	Kstll_Lxy = BToKstll_B_Lxy[triplet_sel_index];
	Kstll_ctxy = BToKstll_B_ctxy[triplet_sel_index];
    
	Kstll_l1_pt = BToKstll_lep1_pt[triplet_sel_index];	
	Kstll_l1_eta = BToKstll_lep1_eta[triplet_sel_index];
	Kstll_l1_phi = BToKstll_lep1_phi[triplet_sel_index];
	Kstll_l1_preFpt = Lepton_pt[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l1_preFeta = Lepton_eta[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l1_preFphi = Lepton_phi[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l1_pfRelIso03 = Lepton_pfRelIso03[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l1_dxy = Lepton_dxy[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l1_dz = Lepton_dz[BToKstll_lep1_index[triplet_sel_index]];
    
	Kstll_l2_pt = BToKstll_lep2_pt[triplet_sel_index];
	Kstll_l2_eta = BToKstll_lep2_eta[triplet_sel_index];
	Kstll_l2_phi = BToKstll_lep2_phi[triplet_sel_index];
	Kstll_l2_preFpt = Lepton_pt[BToKstll_lep2_index[triplet_sel_index]];
	Kstll_l2_preFeta = Lepton_eta[BToKstll_lep2_index[triplet_sel_index]];
	Kstll_l2_preFphi = Lepton_phi[BToKstll_lep2_index[triplet_sel_index]];
	Kstll_l2_pfRelIso03 = Lepton_pfRelIso03[BToKstll_lep2_index[triplet_sel_index]];
	Kstll_l2_dxy = Lepton_dxy[BToKstll_lep2_index[triplet_sel_index]];
	Kstll_l2_dz = Lepton_dz[BToKstll_lep2_index[triplet_sel_index]];

	if(isEleFinalState){
	  Kstll_l1_pfRelIso04 = -1.;
	  Kstll_l2_pfRelIso04 = -1.;
	}
	else{
	Kstll_l1_pfRelIso04 = Lepton_pfRelIso04[BToKstll_lep1_index[triplet_sel_index]];
	Kstll_l2_pfRelIso04 = Lepton_pfRelIso04[BToKstll_lep2_index[triplet_sel_index]];
	}
	
	Kstll_llRefitmass = llInvRefitMass;

	Kstll_k_charge = BToKstll_kaon_charge[triplet_sel_index];
	Kstll_k_pt = BToKstll_kaon_pt[triplet_sel_index];
	Kstll_k_eta = BToKstll_kaon_eta[triplet_sel_index];
	Kstll_k_phi = BToKstll_kaon_phi[triplet_sel_index];
	Kstll_k_DCASig = BToKstll_kaon_DCASig[triplet_sel_index];
	Kstll_k_dxy = PFCand_dxy[BToKstll_kaon_index[triplet_sel_index]];
	Kstll_k_dz = PFCand_dz[BToKstll_kaon_index[triplet_sel_index]];

	newT->Fill();
    }
    //end output Nutuple

    
    //histograms for each mass bin
    if(massBin != -1){
      hAlpha[massBin]->Fill(BToKstll_B_cosAlpha[triplet_sel_index]);
      hCLVtx[massBin]->Fill(BToKstll_B_CL_vtx[triplet_sel_index]);
      hDCASig[massBin]->Fill(BToKstll_kaon_DCASig[triplet_sel_index]);
      hLxy[massBin]->Fill(BToKstll_B_Lxy[triplet_sel_index]);
      hctxy[massBin]->Fill(BToKstll_B_ctxy[triplet_sel_index]);
      hKaonpt[massBin]->Fill(BToKstll_kaon_pt[triplet_sel_index]);
      hBpt[massBin]->Fill(BToKstll_B_pt[triplet_sel_index]);

      hLep1pt[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      hLep2pt[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      if(l1notPFLep_l2PFLep || l1notPFLep_l2notPFLep)hLep1pt_l1notPFLep[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      if(l1PFLep_l2notPFLep || l1notPFLep_l2notPFLep)hLep2pt_l2notPFLep[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      if(l1PFCand_l2notPFCand || l1PFCand_l2PFCand) hLep1pt_PFCand[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      if(l1notPFCand_l2PFCand || l1PFCand_l2PFCand) hLep2pt_PFCand[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      if(l1LT_l2notLT || l1LT_l2LT) hLep1pt_LT[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      if(l1notLT_l2LT || l1LT_l2LT) hLep2pt_LT[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      if(std::abs(BToKstll_lep1_eta[triplet_sel_index]) < 1.47) hLep1pt_EB[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      else hLep1pt_EE[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
      if(std::abs(BToKstll_lep2_eta[triplet_sel_index]) < 1.47) hLep2pt_EB[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      else hLep2pt_EE[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);

      hllRefitMass[massBin]->Fill(llInvRefitMass);
      hllRefitMass_vs_Bmass[massBin]->Fill(BToKstll_B_mass[triplet_sel_index], llInvRefitMass);
      
      hBmass[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(isllt){ 
          hBmass_llt[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
          hLep1pt_PFLep[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
          hLep2pt_PFLep[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      }
      if(l1PFLep_l2notPFLep){
          hBmass_l1PFLep_l2notPFLep[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
          hLep1pt_l1PFLep_l2notPFLep[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
          hLep2pt_l1PFLep_l2notPFLep[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      }
      if(l1notPFLep_l2PFLep) hBmass_l1notPFLep_l2PFLep[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(l1notPFLep_l2notPFLep){
          hBmass_l1notPFLep_l2notPFLep[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
          hLep1pt_l1notPFLep_l2notPFLep[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
          hLep2pt_l1notPFLep_l2notPFLep[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
      }
      if(l1PFCand_l2PFCand) hBmass_l1PFCand_l2PFCand[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(l1PFCand_l2notPFCand) hBmass_l1PFCand_l2notPFCand[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(l1notPFCand_l2PFCand) hBmass_l1notPFCand_l2PFCand[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(l1LT_l2LT) hBmass_l1LT_l2LT[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]); 
      if(l1LT_l2notLT) hBmass_l1LT_l2notLT[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      if(l1notLT_l2LT) hBmass_l1notLT_l2LT[massBin]->Fill(BToKstll_B_mass[triplet_sel_index]);
      

      BDTele1_vs_pTele1[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index], BToKstll_lep1_seedBDT_unbiased[triplet_sel_index]);
      BDTele2_vs_pTele2[massBin]->Fill(BToKstll_lep2_pt[triplet_sel_index], BToKstll_lep2_seedBDT_unbiased[triplet_sel_index]);
      BDTele2_vs_BDTele1[massBin]->Fill(BToKstll_lep1_seedBDT_unbiased[triplet_sel_index], BToKstll_lep2_seedBDT_unbiased[triplet_sel_index]);
      pTele2_vs_pTele1[massBin]->Fill(BToKstll_lep1_pt[triplet_sel_index], BToKstll_lep2_pt[triplet_sel_index]);
    }
    
    //histograms inclusive over all m(ll)
    hllRefitMass[6]->Fill(llInvRefitMass);
    hllRefitMass_vs_Bmass[6]->Fill(BToKstll_B_mass[triplet_sel_index], llInvRefitMass);
    hBmass[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(isllt){
        hBmass_llt[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
        hLep1pt_PFLep[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
        hLep2pt_PFLep[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    }
    if(l1PFLep_l2notPFLep){
        hBmass_l1PFLep_l2notPFLep[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
        hLep1pt_l1PFLep_l2notPFLep[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
        hLep2pt_l1PFLep_l2notPFLep[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);        
    }
    if(l1notPFLep_l2PFLep) hBmass_l1notPFLep_l2PFLep[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1notPFLep_l2notPFLep){
        hBmass_l1notPFLep_l2notPFLep[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
        hLep1pt_l1notPFLep_l2notPFLep[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
        hLep2pt_l1notPFLep_l2notPFLep[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);         
    }
    if(l1PFCand_l2PFCand) hBmass_l1PFCand_l2PFCand[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1PFCand_l2notPFCand) hBmass_l1PFCand_l2notPFCand[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1notPFCand_l2PFCand) hBmass_l1notPFCand_l2PFCand[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1LT_l2LT) hBmass_l1LT_l2LT[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1LT_l2notLT) hBmass_l1LT_l2notLT[6]->Fill(BToKstll_B_mass[triplet_sel_index]);
    if(l1notLT_l2LT) hBmass_l1notLT_l2LT[6]->Fill(BToKstll_B_mass[triplet_sel_index]);

    hAlpha[6]->Fill(BToKstll_B_cosAlpha[triplet_sel_index]);
    hCLVtx[6]->Fill(BToKstll_B_CL_vtx[triplet_sel_index]);
    hDCASig[6]->Fill(BToKstll_kaon_DCASig[triplet_sel_index]);
    hLxy[6]->Fill(BToKstll_B_Lxy[triplet_sel_index]);
    hctxy[6]->Fill(BToKstll_B_ctxy[triplet_sel_index]);
    hKaonpt[6]->Fill(BToKstll_kaon_pt[triplet_sel_index]);
    hBpt[6]->Fill(BToKstll_B_pt[triplet_sel_index]);

    hLep1pt[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    hLep2pt[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    if(l1notPFLep_l2PFLep || l1notPFLep_l2notPFLep)hLep1pt_l1notPFLep[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    if(l1PFLep_l2notPFLep || l1notPFLep_l2notPFLep)hLep2pt_l2notPFLep[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    if(l1PFCand_l2notPFCand || l1PFCand_l2PFCand) hLep1pt_PFCand[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    if(l1notPFCand_l2PFCand || l1PFCand_l2PFCand) hLep2pt_PFCand[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    if(l1LT_l2notLT || l1LT_l2LT) hLep1pt_LT[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    if(l1notLT_l2LT || l1LT_l2LT) hLep2pt_LT[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    if(std::abs(BToKstll_lep1_eta[triplet_sel_index]) < 1.47) hLep1pt_EB[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    else hLep1pt_EE[6]->Fill(BToKstll_lep1_pt[triplet_sel_index]);
    if(std::abs(BToKstll_lep2_eta[triplet_sel_index]) < 1.47) hLep2pt_EB[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    else hLep2pt_EE[6]->Fill(BToKstll_lep2_pt[triplet_sel_index]);
    
    BDTele1_vs_pTele1[6]->Fill(BToKstll_lep1_pt[triplet_sel_index], BToKstll_lep1_seedBDT_unbiased[triplet_sel_index]);
    BDTele2_vs_pTele2[6]->Fill(BToKstll_lep2_pt[triplet_sel_index], BToKstll_lep2_seedBDT_unbiased[triplet_sel_index]);
    BDTele2_vs_BDTele1[6]->Fill(BToKstll_lep1_seedBDT_unbiased[triplet_sel_index], BToKstll_lep2_seedBDT_unbiased[triplet_sel_index]);
    pTele2_vs_pTele1[6]->Fill(BToKstll_lep1_pt[triplet_sel_index], BToKstll_lep2_pt[triplet_sel_index]);
  
    /*
    if(massBin == 3 && BToKstll_B_mass[triplet_sel_index] < 6.)
      std::cout << run << " " << lumi << " " << event << std::endl;
    */
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
    hDCASig[ij]->Write(hDCASig[ij]->GetName());
    hLxy[ij]->Write(hLxy[ij]->GetName());
    hctxy[ij]->Write(hctxy[ij]->GetName());
    hKaonpt[ij]->Write(hKaonpt[ij]->GetName());
    hLep1pt[ij]->Write(hLep1pt[ij]->GetName());
    hLep1pt_PFLep[ij]->Write(hLep1pt_PFLep[ij]->GetName());
    hLep1pt_l1notPFLep[ij]->Write(hLep1pt_l1notPFLep[ij]->GetName());
    hLep1pt_l1PFLep_l2notPFLep[ij]->Write(hLep1pt_l1PFLep_l2notPFLep[ij]->GetName());
    hLep1pt_l1notPFLep_l2notPFLep[ij]->Write(hLep1pt_l1notPFLep_l2notPFLep[ij]->GetName());
    hLep1pt_PFCand[ij]->Write(hLep1pt_PFCand[ij]->GetName());
    hLep1pt_LT[ij]->Write(hLep1pt_LT[ij]->GetName());
    hLep2pt[ij]->Write(hLep2pt[ij]->GetName());
    hLep2pt_PFLep[ij]->Write(hLep2pt_PFLep[ij]->GetName());
    hLep2pt_l2notPFLep[ij]->Write(hLep2pt_l2notPFLep[ij]->GetName());
    hLep2pt_l1PFLep_l2notPFLep[ij]->Write(hLep2pt_l1PFLep_l2notPFLep[ij]->GetName());
    hLep2pt_l1notPFLep_l2notPFLep[ij]->Write(hLep2pt_l1notPFLep_l2notPFLep[ij]->GetName());
    hLep2pt_PFCand[ij]->Write(hLep2pt_PFCand[ij]->GetName());
    hLep2pt_LT[ij]->Write(hLep2pt_LT[ij]->GetName());
    hLep1pt_EB[ij]->Write(hLep1pt_EB[ij]->GetName());
    hLep1pt_EE[ij]->Write(hLep1pt_EE[ij]->GetName());
    hLep2pt_EB[ij]->Write(hLep2pt_EB[ij]->GetName());
    hLep2pt_EE[ij]->Write(hLep2pt_EE[ij]->GetName());

    hBpt[ij]->Write(hBpt[ij]->GetName());
    std::cout << " >>> hBpt[ij]->GetName() = " << hBpt[ij]->GetName() << " entries = " << hBpt[ij]->GetEntries() << std::endl;

    hllRefitMass[ij]->Write(hllRefitMass[ij]->GetName());
    hllRefitMass_vs_Bmass[ij]->Write(hllRefitMass_vs_Bmass[ij]->GetName());
    hBmass[ij]->Write(hBmass[ij]->GetName());
    hBmass_llt[ij]->Write(hBmass_llt[ij]->GetName());
    hBmass_l1notPFLep_l2PFLep[ij]->Write(hBmass_l1notPFLep_l2PFLep[ij]->GetName());
    hBmass_l1PFLep_l2notPFLep[ij]->Write(hBmass_l1PFLep_l2notPFLep[ij]->GetName());
    hBmass_l1notPFLep_l2notPFLep[ij]->Write(hBmass_l1notPFLep_l2notPFLep[ij]->GetName());
    hBmass_l1PFCand_l2PFCand[ij]->Write(hBmass_l1PFCand_l2PFCand[ij]->GetName());
    hBmass_l1PFCand_l2notPFCand[ij]->Write(hBmass_l1PFCand_l2notPFCand[ij]->GetName());
    hBmass_l1notPFCand_l2PFCand[ij]->Write(hBmass_l1notPFCand_l2PFCand[ij]->GetName());
    hBmass_l1LT_l2LT[ij]->Write(hBmass_l1LT_l2LT[ij]->GetName());
    hBmass_l1LT_l2notLT[ij]->Write(hBmass_l1LT_l2notLT[ij]->GetName());
    hBmass_l1notLT_l2LT[ij]->Write(hBmass_l1notLT_l2LT[ij]->GetName());

    BDTele1_vs_pTele1[ij]->Write(BDTele1_vs_pTele1[ij]->GetName());
    BDTele2_vs_pTele2[ij]->Write(BDTele2_vs_pTele2[ij]->GetName());
    BDTele2_vs_BDTele1[ij]->Write(BDTele2_vs_BDTele1[ij]->GetName());
    pTele2_vs_pTele1[ij]->Write(pTele2_vs_pTele1[ij]->GetName());

    if(ij > 5) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
	      << " \n \t muonTag = " << nEv_muonTag[0] << " \t recoEvts = " << nEv_recoCand[0]
              << " \n \t chargeSel = " << nEv_chargeSel[0] << " \t selected Events = " << nEv_selected[ij] << std::endl;

  }
  outMassHistos.Close();

}  
