// Author: T. Strebler (IC)
// Date:   31 May 2018
//
// Wrapper for NanoAOD tree
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//
// Create the nanoAODTree object from the pointer to the tree, then access the stored objects from it
// Common TTree functions GetEntry (entry), GetEntries() are implemented



#ifndef NANOAODTREE_H
#define NANOAODTREE_H



#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>

const int kLeptonMax = 100; // set to 2 times nTriplets
const int kBToKstllMax = 100;

const int kMuonMax = 100;
const int kElectronMax = 100;
const int kBToKpipiMax = 1000000;
const int kBToKmumuMax = 50000;
const int kBToKeeMax = 50000;
const int kGenPartMax = 10000;
const int kTrigObjMax = 1000;
const int kLowPtGsfTrackMax = 100;
const int kPFCandMax = 10000;
const int kLostTrackMax = 100;

using namespace std;

class NanoAODTree {
public :
   TChain          *_tree; 

   //Old NanoAOD branches used
   int run;
   int luminosityBlock;
   long event;

   uint nMuon;
   int Muon_charge[kMuonMax];
   float Muon_pt[kMuonMax];
   float Muon_eta[kMuonMax];
   float Muon_phi[kMuonMax];
   float Muon_mass[kMuonMax];
   float Muon_dxy[kMuonMax];
   float Muon_dz[kMuonMax];
   float Muon_pfRelIso04_all[kMuonMax];
   bool Muon_softId[kMuonMax];
   bool Muon_mediumId[kMuonMax];

   int BToKstll_configuration[kBToKstllMax];

   uint nBToKstll;
   float BToKstll_B_CL_vtx[kBToKstllMax];
   
   int BToKstll_lep1_charge[kBToKstllMax];
   float BToKstll_lep1_pt[kBToKstllMax];
   float BToKstll_lep1_eta[kBToKstllMax];
   float BToKstll_lep1_phi[kBToKstllMax];
   int BToKstll_lep1_index[kBToKstllMax];
   int BToKstll_lep1_lowPt_index[kBToKstllMax];
   int BToKstll_lep1_pfCand_index[kBToKstllMax];
   int BToKstll_lep1_lostTrack_index[kBToKstllMax];
   int BToKstll_lep1_isPFLep[kBToKstllMax];
   int BToKstll_lep1_isLowPt[kBToKstllMax];
   int BToKstll_lep1_isPFCand[kBToKstllMax];
   
   int BToKstll_lep2_charge[kBToKstllMax];
   float BToKstll_lep2_pt[kBToKstllMax];
   float BToKstll_lep2_eta[kBToKstllMax];
   float BToKstll_lep2_phi[kBToKstllMax];
   int BToKstll_lep2_index[kBToKstllMax];
   int BToKstll_lep2_lowPt_index[kBToKstllMax];
   int BToKstll_lep2_pfCand_index[kBToKstllMax];
   int BToKstll_lep2_lostTrack_index[kBToKstllMax];   
   int BToKstll_lep2_isPFLep[kBToKstllMax];
   int BToKstll_lep2_isLowPt[kBToKstllMax];
   int BToKstll_lep2_isPFCand[kBToKstllMax];   
   
   //int BToKstll_kaon_charge[kBToKstllMax];
   float BToKstll_kaon_pt[kBToKstllMax];
   float BToKstll_kaon_eta[kBToKstllMax];
   float BToKstll_kaon_phi[kBToKstllMax];
   int BToKstll_muTrg_index[kBToKstllMax];


   uint nElectron;
   int Electron_charge[kElectronMax];
   float Electron_pt[kElectronMax];
   float Electron_eta[kElectronMax];
   float Electron_phi[kElectronMax];
   float Electron_mass[kElectronMax];
   float Electron_dxy[kElectronMax];
   float Electron_dz[kElectronMax];

   uint nBToKpipi;
   float BToKpipi_CL_vtx[kBToKpipiMax];
   float BToKpipi_cosAlpha[kBToKpipiMax];
   float BToKpipi_Lxy[kBToKpipiMax];
   float BToKpipi_Kpi_CL_vtx[kBToKpipiMax];
   float BToKpipi_Kpi_mass[kBToKpipiMax];
   float BToKpipi_Kpi_pt[kBToKpipiMax];
   float BToKpipi_mass[kBToKpipiMax];
   int BToKpipi_piBu_charge[kBToKpipiMax];
   float BToKpipi_piBu_pt[kBToKpipiMax];
   float BToKpipi_piBu_eta[kBToKpipiMax];
   float BToKpipi_piBu_phi[kBToKpipiMax];
   float BToKpipi_kaon_pt[kBToKpipiMax];
   float BToKpipi_kaon_eta[kBToKpipiMax];
   float BToKpipi_kaon_phi[kBToKpipiMax];
   int BToKpipi_kaon_charge[kBToKpipiMax];
   float BToKpipi_kaon_dz[kBToKpipiMax];
   float BToKpipi_piD0_pt[kBToKpipiMax];
   float BToKpipi_piD0_eta[kBToKpipiMax];
   float BToKpipi_piD0_phi[kBToKpipiMax];
   float BToKpipi_piD0_dz[kBToKpipiMax];
 
   uint nBToKmumu;
   float BToKmumu_CL_vtx[kBToKmumuMax];
   float BToKmumu_mumu_CL_vtx[kBToKmumuMax];
   float BToKmumu_mumu_mass[kBToKmumuMax];
   float BToKmumu_mass[kBToKmumuMax];
   int BToKmumu_kaon_charge[kBToKmumuMax];
   float BToKmumu_kaon_pt[kBToKmumuMax];
   float BToKmumu_kaon_eta[kBToKmumuMax];
   float BToKmumu_kaon_phi[kBToKmumuMax];
   float BToKmumu_mu1_pt[kBToKmumuMax];
   int BToKmumu_mu1_charge[kBToKmumuMax];
   float BToKmumu_mu1_eta[kBToKmumuMax];
   float BToKmumu_mu1_phi[kBToKmumuMax];
   int BToKmumu_mu1_index[kBToKmumuMax];
   int BToKmumu_mu2_charge[kBToKmumuMax];
   float BToKmumu_mu2_pt[kBToKmumuMax];
   float BToKmumu_mu2_eta[kBToKmumuMax];
   float BToKmumu_mu2_phi[kBToKmumuMax];
   int BToKmumu_mu2_index[kBToKmumuMax];


   uint nBToKee;
   float BToKee_CL_vtx[kBToKeeMax];
   float BToKee_ee_CL_vtx[kBToKeeMax];
   float BToKee_ee_mass[kBToKeeMax];
   float BToKee_mass[kBToKeeMax];
   int BToKee_kaon_charge[kBToKeeMax];
   float BToKee_kaon_pt[kBToKeeMax];
   float BToKee_kaon_eta[kBToKeeMax];
   float BToKee_kaon_phi[kBToKeeMax];
   int BToKee_ele1_index[kBToKeeMax];
   int BToKee_ele1_charge[kBToKeeMax];
   float BToKee_ele1_pt[kBToKeeMax];
   float BToKee_ele1_eta[kBToKeeMax];
   float BToKee_ele1_phi[kBToKeeMax];
   int BToKee_ele2_index[kBToKeeMax];
   int BToKee_ele2_charge[kBToKeeMax];
   float BToKee_ele2_pt[kBToKeeMax];
   float BToKee_ele2_eta[kBToKeeMax];
   float BToKee_ele2_phi[kBToKeeMax];

   uint nGenPart;
   int GenPart_pdgId[kGenPartMax];
   int GenPart_genPartIdxMother[kGenPartMax];
   float GenPart_pt[kGenPartMax];
   float GenPart_eta[kGenPartMax];
   float GenPart_phi[kGenPartMax];

   bool HLT_Mu8p5_IP3p5_part0;
   bool HLT_Mu8p5_IP3p5_part1;
   bool HLT_Mu8p5_IP3p5_part2;
   bool HLT_Mu8p5_IP3p5_part3;
   bool HLT_Mu8p5_IP3p5_part4;
   bool HLT_Mu8p5_IP3p5_part5;
   bool HLT_Mu10p5_IP3p5_part0;
   bool HLT_Mu10p5_IP3p5_part1;
   bool HLT_Mu10p5_IP3p5_part2;
   bool HLT_Mu10p5_IP3p5_part3;
   bool HLT_Mu10p5_IP3p5_part4;
   bool HLT_Mu10p5_IP3p5_part5;
   bool HLT_Mu9_IP6_part0;
   bool HLT_Mu9_IP6_part1;
   bool HLT_Mu9_IP6_part2;
   bool HLT_Mu9_IP6_part3;
   bool HLT_Mu9_IP6_part4;
   bool HLT_Mu9_IP6_part5;
   bool HLT_Mu8_IP3_part0;
   bool HLT_Mu8_IP3_part1;
   bool HLT_Mu8_IP3_part2;
   bool HLT_Mu8_IP3_part3;
   bool HLT_Mu8_IP3_part4;
   bool HLT_Mu8_IP3_part5;

   //for MC
   bool HLT_Mu8_IP6_part0;
   bool HLT_Mu8_IP6_part1;
   bool HLT_Mu8_IP6_part2;
   bool HLT_Mu8_IP6_part3;
   bool HLT_Mu8_IP6_part4;
   bool HLT_Mu8_IP6_part5;

   bool HLT_Mu8_IP5_part0;
   bool HLT_Mu8_IP5_part1;
   bool HLT_Mu8_IP5_part2;
   bool HLT_Mu8_IP5_part3;
   bool HLT_Mu8_IP5_part4;
   bool HLT_Mu8_IP5_part5;

   bool HLT_Mu9_IP4_part0;
   bool HLT_Mu9_IP4_part1;
   bool HLT_Mu9_IP4_part2;
   bool HLT_Mu9_IP4_part3;
   bool HLT_Mu9_IP4_part4;
   bool HLT_Mu9_IP4_part5;

   bool HLT_Mu7_IP4_part0;
   bool HLT_Mu7_IP4_part1;
   bool HLT_Mu7_IP4_part2;
   bool HLT_Mu7_IP4_part3;
   bool HLT_Mu7_IP4_part4;
   bool HLT_Mu7_IP4_part5;

   bool HLT_Mu9_IP5_part0;
   bool HLT_Mu9_IP5_part1;
   bool HLT_Mu9_IP5_part2;
   bool HLT_Mu9_IP5_part3;
   bool HLT_Mu9_IP5_part4;
   bool HLT_Mu9_IP5_part5;

   bool HLT_Mu12_IP6_part0;
   bool HLT_Mu12_IP6_part1;
   bool HLT_Mu12_IP6_part2;
   bool HLT_Mu12_IP6_part3;
   bool HLT_Mu12_IP6_part4;
   bool HLT_Mu12_IP6_part5;

   uint nTrigObj;
   int TrigObj_id[kTrigObjMax];
   float TrigObj_pt[kTrigObjMax];
   float TrigObj_eta[kTrigObjMax];
   float TrigObj_phi[kTrigObjMax];
   int TrigObj_filterBits[kTrigObjMax];
   
   uint nLowPtGsfTrack;
   float LowPtGsfTrack_pt[kLowPtGsfTrackMax];
   float LowPtGsfTrack_eta[kLowPtGsfTrackMax];
   float LowPtGsfTrack_phi[kLowPtGsfTrackMax];
   float LowPtGsfTrack_charge[kLowPtGsfTrackMax];
   float LowPtGsfTrack_seedBDT_unbiased[kLowPtGsfTrackMax];

   uint nPFCand;
   float PFCand_pt[kPFCandMax];
   float PFCand_eta[kPFCandMax];
   float PFCand_phi[kPFCandMax];
   float PFCand_charge[kPFCandMax];
   float PFCand_mass[kPFCandMax];
   int PFCand_pdgId[kPFCandMax];
   float PFCand_DCASig[kPFCandMax];
   float PFCand_dz[kPFCandMax];

   uint nLostTrack;
   float LostTrack_pt[kLostTrackMax];
   float LostTrack_eta[kLostTrackMax];
   float LostTrack_phi[kLostTrackMax];
   float LostTrack_charge[kLostTrackMax];
   
   // methods
   NanoAODTree (TChain* tree);
   ~NanoAODTree();
   void Init(TChain* tree);
   Int_t GetEntry(int entry);
   Long64_t GetEntries();
   TChain* GetTree();

};



NanoAODTree::NanoAODTree (TChain* tree)
{
  Init(tree);
}


NanoAODTree::~NanoAODTree() {}



void NanoAODTree::Init(TChain* tree)
{

  // Set branch addresses and branch pointers
  if (!tree) return;
  _tree = tree;  
  _tree->SetMakeClass(1); // needed especially when compiling
  _tree->GetEntries(); //ROOT issue with TChain which only works after calling GetEntries... (sigh)
  
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("luminosityBlock",&luminosityBlock);
  _tree->SetBranchAddress("event",&event);

  _tree->SetBranchAddress("nMuon",&nMuon);  
  _tree->SetBranchAddress("Muon_charge",&Muon_charge);  
  _tree->SetBranchAddress("Muon_pt",&Muon_pt);  
  _tree->SetBranchAddress("Muon_eta",&Muon_eta);
  _tree->SetBranchAddress("Muon_phi",&Muon_phi);
  _tree->SetBranchAddress("Muon_mass",&Muon_mass);
  _tree->SetBranchAddress("Muon_dxy",&Muon_dxy);
  _tree->SetBranchAddress("Muon_dz",&Muon_dz);
  _tree->SetBranchAddress("Muon_pfRelIso04_all",&Muon_pfRelIso04_all);
  _tree->SetBranchAddress("Muon_softId",&Muon_softId);
  _tree->SetBranchAddress("Muon_mediumId",&Muon_mediumId);

  _tree->SetBranchAddress("nElectron",&nElectron);  
  _tree->SetBranchAddress("Electron_charge",&Electron_charge);  
  _tree->SetBranchAddress("Electron_pt",&Electron_pt);  
  _tree->SetBranchAddress("Electron_eta",&Electron_eta);
  _tree->SetBranchAddress("Electron_phi",&Electron_phi);
  _tree->SetBranchAddress("Electron_mass",&Electron_mass);
  _tree->SetBranchAddress("Electron_dxy",&Electron_dxy);
  _tree->SetBranchAddress("Electron_dz",&Electron_dz);


  int BToKstll_info = _tree->SetBranchAddress("nBToKstll",&nBToKstll);
  if(BToKstll_info >= 0){
    _tree->SetBranchAddress("BToKstll_configuration",&BToKstll_configuration);
    _tree->SetBranchAddress("BToKstll_B_CL_vtx",&BToKstll_B_CL_vtx);
    _tree->SetBranchAddress("BToKstll_lep1_charge",&BToKstll_lep1_charge);
    _tree->SetBranchAddress("BToKstll_lep1_pt",&BToKstll_lep1_pt);
    _tree->SetBranchAddress("BToKstll_lep1_eta",&BToKstll_lep1_eta);
    _tree->SetBranchAddress("BToKstll_lep1_phi",&BToKstll_lep1_phi);
    _tree->SetBranchAddress("BToKstll_lep1_index",&BToKstll_lep1_index);
    _tree->SetBranchAddress("BToKstll_lep1_lowPt_index",&BToKstll_lep1_lowPt_index);
    _tree->SetBranchAddress("BToKstll_lep1_pfCand_index",&BToKstll_lep1_pfCand_index);
    _tree->SetBranchAddress("BToKstll_lep1_lostTrack_index",&BToKstll_lep1_lostTrack_index);
    _tree->SetBranchAddress("BToKstll_lep1_isPFLep",&BToKstll_lep1_isPFLep);
    _tree->SetBranchAddress("BToKstll_lep1_isLowPt",&BToKstll_lep1_isLowPt);
    _tree->SetBranchAddress("BToKstll_lep1_isPFCand",&BToKstll_lep1_isPFCand);
    _tree->SetBranchAddress("BToKstll_lep2_charge",&BToKstll_lep2_charge);
    _tree->SetBranchAddress("BToKstll_lep2_pt",&BToKstll_lep2_pt);
    _tree->SetBranchAddress("BToKstll_lep2_eta",&BToKstll_lep2_eta);
    _tree->SetBranchAddress("BToKstll_lep2_phi",&BToKstll_lep2_phi);
    _tree->SetBranchAddress("BToKstll_lep2_index",&BToKstll_lep2_index);
    _tree->SetBranchAddress("BToKstll_lep2_lowPt_index",&BToKstll_lep2_lowPt_index);
    _tree->SetBranchAddress("BToKstll_lep2_pfCand_index",&BToKstll_lep2_pfCand_index);
    _tree->SetBranchAddress("BToKstll_lep2_lostTrack_index",&BToKstll_lep2_lostTrack_index);    
    _tree->SetBranchAddress("BToKstll_lep2_isPFLep",&BToKstll_lep2_isPFLep);
    _tree->SetBranchAddress("BToKstll_lep2_isLowPt",&BToKstll_lep2_isLowPt);
    _tree->SetBranchAddress("BToKstll_lep2_isPFCand",&BToKstll_lep2_isPFCand);
    //_tree->SetBranchAddress("BToKstll_kaon_charge",&BToKstll_kaon_charge);
    _tree->SetBranchAddress("BToKstll_kaon_pt",&BToKstll_kaon_pt);
    _tree->SetBranchAddress("BToKstll_kaon_eta",&BToKstll_kaon_eta);
    _tree->SetBranchAddress("BToKstll_kaon_phi",&BToKstll_kaon_phi);
    _tree->SetBranchAddress("BToKstll_muTrg_index",&BToKstll_muTrg_index);
  }


  int BToKpipi_info = _tree->SetBranchAddress("nBToKpipi",&nBToKpipi);
  if(BToKpipi_info>=0){
    _tree->SetBranchAddress("BToKpipi_CL_vtx",&BToKpipi_CL_vtx);
    _tree->SetBranchAddress("BToKpipi_Lxy",&BToKpipi_Lxy);
    _tree->SetBranchAddress("BToKpipi_cosAlpha",&BToKpipi_cosAlpha);
    _tree->SetBranchAddress("BToKpipi_Kpi_CL_vtx",&BToKpipi_Kpi_CL_vtx);
    _tree->SetBranchAddress("BToKpipi_Kpi_mass",&BToKpipi_Kpi_mass);
    _tree->SetBranchAddress("BToKpipi_Kpi_pt",&BToKpipi_Kpi_pt);
    _tree->SetBranchAddress("BToKpipi_mass",&BToKpipi_mass);
    _tree->SetBranchAddress("BToKpipi_piBu_charge",&BToKpipi_piBu_charge);
    _tree->SetBranchAddress("BToKpipi_piBu_pt",&BToKpipi_piBu_pt);
    _tree->SetBranchAddress("BToKpipi_piBu_eta",&BToKpipi_piBu_eta);
    _tree->SetBranchAddress("BToKpipi_piBu_phi",&BToKpipi_piBu_phi);
    _tree->SetBranchAddress("BToKpipi_kaon_charge",&BToKpipi_kaon_charge);
    _tree->SetBranchAddress("BToKpipi_kaon_pt",&BToKpipi_kaon_pt);
    _tree->SetBranchAddress("BToKpipi_kaon_eta",&BToKpipi_kaon_eta);
    _tree->SetBranchAddress("BToKpipi_kaon_phi",&BToKpipi_kaon_phi);
    _tree->SetBranchAddress("BToKpipi_kaon_dz",&BToKpipi_kaon_dz);
    _tree->SetBranchAddress("BToKpipi_piD0_pt",&BToKpipi_piD0_pt);
    _tree->SetBranchAddress("BToKpipi_piD0_eta",&BToKpipi_piD0_eta);
    _tree->SetBranchAddress("BToKpipi_piD0_phi",&BToKpipi_piD0_phi);
    _tree->SetBranchAddress("BToKpipi_piD0_dz",&BToKpipi_piD0_dz);
  }

  int BToKmumu_info = _tree->SetBranchAddress("nBToKmumu",&nBToKmumu);
  if(BToKmumu_info>=0){
    _tree->SetBranchAddress("BToKmumu_CL_vtx",&BToKmumu_CL_vtx);
    _tree->SetBranchAddress("BToKmumu_mumu_CL_vtx",&BToKmumu_mumu_CL_vtx);
    _tree->SetBranchAddress("BToKmumu_mumu_mass",&BToKmumu_mumu_mass);
    _tree->SetBranchAddress("BToKmumu_mass",&BToKmumu_mass);
    _tree->SetBranchAddress("BToKmumu_kaon_charge",&BToKmumu_kaon_charge);
    _tree->SetBranchAddress("BToKmumu_kaon_pt",&BToKmumu_kaon_pt);
    _tree->SetBranchAddress("BToKmumu_kaon_eta",&BToKmumu_kaon_eta);
    _tree->SetBranchAddress("BToKmumu_kaon_phi",&BToKmumu_kaon_phi);
    _tree->SetBranchAddress("BToKmumu_mu1_charge",&BToKmumu_mu1_charge);
    _tree->SetBranchAddress("BToKmumu_mu1_pt",&BToKmumu_mu1_pt);
    _tree->SetBranchAddress("BToKmumu_mu1_eta",&BToKmumu_mu1_eta);
    _tree->SetBranchAddress("BToKmumu_mu1_phi",&BToKmumu_mu1_phi);
    _tree->SetBranchAddress("BToKmumu_mu1_index",&BToKmumu_mu1_index);
    _tree->SetBranchAddress("BToKmumu_mu2_charge",&BToKmumu_mu2_charge);
    _tree->SetBranchAddress("BToKmumu_mu2_pt",&BToKmumu_mu2_pt);
    _tree->SetBranchAddress("BToKmumu_mu2_eta",&BToKmumu_mu2_eta);
    _tree->SetBranchAddress("BToKmumu_mu2_phi",&BToKmumu_mu2_phi);
    _tree->SetBranchAddress("BToKmumu_mu2_index",&BToKmumu_mu2_index);
  }

  int BToKee_info = _tree->SetBranchAddress("nBToKee",&nBToKee);
  if(BToKee_info>=0){
    _tree->SetBranchAddress("BToKee_CL_vtx",&BToKee_CL_vtx);
    _tree->SetBranchAddress("BToKee_ee_CL_vtx",&BToKee_ee_CL_vtx);
    _tree->SetBranchAddress("BToKee_ee_mass",&BToKee_ee_mass);
    _tree->SetBranchAddress("BToKee_mass",&BToKee_mass);
    _tree->SetBranchAddress("BToKee_kaon_charge",&BToKee_kaon_charge);
    _tree->SetBranchAddress("BToKee_kaon_pt",&BToKee_kaon_pt);
    _tree->SetBranchAddress("BToKee_kaon_eta",&BToKee_kaon_eta);
    _tree->SetBranchAddress("BToKee_kaon_phi",&BToKee_kaon_phi);
    _tree->SetBranchAddress("BToKee_ele1_pt",&BToKee_ele1_pt);
    _tree->SetBranchAddress("BToKee_ele1_eta",&BToKee_ele1_eta);
    _tree->SetBranchAddress("BToKee_ele1_phi",&BToKee_ele1_phi);
    _tree->SetBranchAddress("BToKee_ele1_charge",&BToKee_ele1_charge);
    _tree->SetBranchAddress("BToKee_ele1_index",&BToKee_ele1_index);
    _tree->SetBranchAddress("BToKee_ele2_pt",&BToKee_ele2_pt);
    _tree->SetBranchAddress("BToKee_ele2_eta",&BToKee_ele2_eta);
    _tree->SetBranchAddress("BToKee_ele2_phi",&BToKee_ele2_phi);
    _tree->SetBranchAddress("BToKee_ele2_charge",&BToKee_ele2_charge);
    _tree->SetBranchAddress("BToKee_ele2_index",&BToKee_ele2_index);
  }

  int isMC = _tree->SetBranchAddress("nGenPart",&nGenPart);
  if(isMC>=0){
    _tree->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
    _tree->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);
    _tree->SetBranchAddress("GenPart_pt",&GenPart_pt);
    _tree->SetBranchAddress("GenPart_eta",&GenPart_eta);
    _tree->SetBranchAddress("GenPart_phi",&GenPart_phi);
  }

  
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part0",&HLT_Mu8p5_IP3p5_part0); 
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part1",&HLT_Mu8p5_IP3p5_part1);
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part2",&HLT_Mu8p5_IP3p5_part2);
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part3",&HLT_Mu8p5_IP3p5_part3);
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part4",&HLT_Mu8p5_IP3p5_part4);
  _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part5",&HLT_Mu8p5_IP3p5_part5);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part0",&HLT_Mu10p5_IP3p5_part0);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part1",&HLT_Mu10p5_IP3p5_part1);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part2",&HLT_Mu10p5_IP3p5_part2);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part3",&HLT_Mu10p5_IP3p5_part3);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part4",&HLT_Mu10p5_IP3p5_part4);
  _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part5",&HLT_Mu10p5_IP3p5_part5);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part0",&HLT_Mu9_IP6_part0);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part1",&HLT_Mu9_IP6_part1);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part2",&HLT_Mu9_IP6_part2);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part3",&HLT_Mu9_IP6_part3);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part4",&HLT_Mu9_IP6_part4);
  _tree->SetBranchAddress("HLT_Mu9_IP6_part5",&HLT_Mu9_IP6_part5);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part0",&HLT_Mu8_IP3_part0);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part1",&HLT_Mu8_IP3_part1);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part2",&HLT_Mu8_IP3_part2);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part3",&HLT_Mu8_IP3_part3);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part4",&HLT_Mu8_IP3_part4);
  _tree->SetBranchAddress("HLT_Mu8_IP3_part5",&HLT_Mu8_IP3_part5);
  
  _tree->SetBranchAddress("HLT_Mu8_IP6_part0",&HLT_Mu8_IP6_part0);
  _tree->SetBranchAddress("HLT_Mu8_IP6_part1",&HLT_Mu8_IP6_part1);
  _tree->SetBranchAddress("HLT_Mu8_IP6_part2",&HLT_Mu8_IP6_part2);
  _tree->SetBranchAddress("HLT_Mu8_IP6_part3",&HLT_Mu8_IP6_part3);
  _tree->SetBranchAddress("HLT_Mu8_IP6_part4",&HLT_Mu8_IP6_part4);
  _tree->SetBranchAddress("HLT_Mu8_IP6_part5",&HLT_Mu8_IP6_part5);
  
  _tree->SetBranchAddress("HLT_Mu8_IP5_part0",&HLT_Mu8_IP5_part0);
  _tree->SetBranchAddress("HLT_Mu8_IP5_part1",&HLT_Mu8_IP5_part1);
  _tree->SetBranchAddress("HLT_Mu8_IP5_part2",&HLT_Mu8_IP5_part2);
  _tree->SetBranchAddress("HLT_Mu8_IP5_part3",&HLT_Mu8_IP5_part3);
  _tree->SetBranchAddress("HLT_Mu8_IP5_part4",&HLT_Mu8_IP5_part4);
  _tree->SetBranchAddress("HLT_Mu8_IP5_part5",&HLT_Mu8_IP5_part5);

  _tree->SetBranchAddress("HLT_Mu9_IP4_part0",&HLT_Mu9_IP4_part0);
  _tree->SetBranchAddress("HLT_Mu9_IP4_part1",&HLT_Mu9_IP4_part1);
  _tree->SetBranchAddress("HLT_Mu9_IP4_part2",&HLT_Mu9_IP4_part2);
  _tree->SetBranchAddress("HLT_Mu9_IP4_part3",&HLT_Mu9_IP4_part3);
  _tree->SetBranchAddress("HLT_Mu9_IP4_part4",&HLT_Mu9_IP4_part4);
  _tree->SetBranchAddress("HLT_Mu9_IP4_part5",&HLT_Mu9_IP4_part5);
  
  _tree->SetBranchAddress("HLT_Mu7_IP4_part0",&HLT_Mu7_IP4_part0);
  _tree->SetBranchAddress("HLT_Mu7_IP4_part1",&HLT_Mu7_IP4_part1);
  _tree->SetBranchAddress("HLT_Mu7_IP4_part2",&HLT_Mu7_IP4_part2);
  _tree->SetBranchAddress("HLT_Mu7_IP4_part3",&HLT_Mu7_IP4_part3);
  _tree->SetBranchAddress("HLT_Mu7_IP4_part4",&HLT_Mu7_IP4_part4);
  _tree->SetBranchAddress("HLT_Mu7_IP4_part5",&HLT_Mu7_IP4_part5);
  
  _tree->SetBranchAddress("HLT_Mu9_IP5_part0",&HLT_Mu9_IP5_part0);
  _tree->SetBranchAddress("HLT_Mu9_IP5_part1",&HLT_Mu9_IP5_part1);
  _tree->SetBranchAddress("HLT_Mu9_IP5_part2",&HLT_Mu9_IP5_part2);
  _tree->SetBranchAddress("HLT_Mu9_IP5_part3",&HLT_Mu9_IP5_part3);
  _tree->SetBranchAddress("HLT_Mu9_IP5_part4",&HLT_Mu9_IP5_part4);
  _tree->SetBranchAddress("HLT_Mu9_IP5_part5",&HLT_Mu9_IP5_part5);
  
  _tree->SetBranchAddress("HLT_Mu12_IP6_part0",&HLT_Mu12_IP6_part0);
  _tree->SetBranchAddress("HLT_Mu12_IP6_part1",&HLT_Mu12_IP6_part1);
  _tree->SetBranchAddress("HLT_Mu12_IP6_part2",&HLT_Mu12_IP6_part2);
  _tree->SetBranchAddress("HLT_Mu12_IP6_part3",&HLT_Mu12_IP6_part3);
  _tree->SetBranchAddress("HLT_Mu12_IP6_part4",&HLT_Mu12_IP6_part4);
  _tree->SetBranchAddress("HLT_Mu12_IP6_part5",&HLT_Mu12_IP6_part5);
  
  _tree->SetBranchAddress("nTrigObj",&nTrigObj);
  _tree->SetBranchAddress("TrigObj_id",&TrigObj_id);
  _tree->SetBranchAddress("TrigObj_pt",&TrigObj_pt);
  _tree->SetBranchAddress("TrigObj_eta",&TrigObj_eta);
  _tree->SetBranchAddress("TrigObj_phi",&TrigObj_phi);
  _tree->SetBranchAddress("TrigObj_filterBits",&TrigObj_filterBits);

  int LowPtGsfTrack_info = _tree->SetBranchAddress("nLowPtGsfTrack",&nLowPtGsfTrack);
  if(LowPtGsfTrack_info){
    _tree->SetBranchAddress("LowPtGsfTrack_pt",&LowPtGsfTrack_pt);
    _tree->SetBranchAddress("LowPtGsfTrack_eta",&LowPtGsfTrack_eta);
    _tree->SetBranchAddress("LowPtGsfTrack_phi",&LowPtGsfTrack_phi);
    _tree->SetBranchAddress("LowPtGsfTrack_charge",&LowPtGsfTrack_charge);
    _tree->SetBranchAddress("LowPtGsfTrack_seedBDT_unbiased",&LowPtGsfTrack_seedBDT_unbiased);
  }

  int PFCand_info = _tree->SetBranchAddress("nPFCand",&nPFCand);
  if(PFCand_info){
    _tree->SetBranchAddress("PFCand_pt",&PFCand_pt);
    _tree->SetBranchAddress("PFCand_eta",&PFCand_eta);
    _tree->SetBranchAddress("PFCand_phi",&PFCand_phi);
    _tree->SetBranchAddress("PFCand_charge",&PFCand_charge);
    _tree->SetBranchAddress("PFCand_mass",&PFCand_mass);
    _tree->SetBranchAddress("PFCand_pdgId",&PFCand_pdgId);
    _tree->SetBranchAddress("PFCand_DCASig",&PFCand_DCASig);
    _tree->SetBranchAddress("PFCand_dz",&PFCand_dz);
  }
  
  int LostTrack_info = _tree->SetBranchAddress("nLostTrack",&nLostTrack);
  if(LostTrack_info){
    _tree->SetBranchAddress("LostTrack_pt",&LostTrack_pt);
    _tree->SetBranchAddress("LostTrack_eta",&LostTrack_eta);
    _tree->SetBranchAddress("LostTrack_phi",&LostTrack_phi);
    _tree->SetBranchAddress("LostTrack_charge",&LostTrack_charge);
  }
    
}


Int_t NanoAODTree::GetEntry(int entry)
{

  int out = _tree->GetEntry(entry);

  if(nMuon>kMuonMax) return -1;
  if(nElectron>kElectronMax) return -1;
  if(nBToKpipi>kBToKpipiMax) return -1;
  if(nBToKmumu>kBToKmumuMax) return -1;
  if(nBToKee>kBToKeeMax) return -1;
  if(nGenPart>kGenPartMax) return -1;
  if(nTrigObj>kTrigObjMax)  return -1;
  if(nLowPtGsfTrack>kLowPtGsfTrackMax) return -1;
  if(nPFCand>kPFCandMax) return -1;
  if(nLostTrack>kLostTrackMax) return -1;

  return out;

} 

Long64_t NanoAODTree::GetEntries()
{
    return _tree->GetEntries();
}

TChain* NanoAODTree::GetTree()
{
    return _tree;
}

#endif 
