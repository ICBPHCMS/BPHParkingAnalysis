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


const int kMuonMax = 1000000;
const int kElectronMax = 1000000;
const int kBToKEEMax = 1000000;
const int kBToKMuMuMax = 1000000;
const int kGenPartMax = 10000;
const int kTrigObjMax = 1000;
const int kProbeTracksMax = 1000000;
const int kTriggerMuonMax = 1000000;


using namespace std;

class NanoAODTree {
public :
   TChain          *_tree; 

   //Old NanoAOD branches used
   int run;
   int luminosityBlock;
   long event;
   
   uint nBToKEE;
   int BToKEE_l1Idx[kBToKEEMax];
   int BToKEE_l2Idx[kBToKEEMax];
   int BToKEE_kIdx[kBToKEEMax];
   float BToKEE_svprob[kBToKEEMax];
   
   uint nBToKMuMu;
   int BToKMuMu_l1Idx[kBToKMuMuMax];
   int BToKMuMu_l2Idx[kBToKMuMuMax];
   int BToKMuMu_kIdx[kBToKMuMuMax];
   float BToKMuMu_svprob[kBToKMuMuMax];   
   
   uint nProbeTracks;
   float ProbeTracks_pt[kProbeTracksMax];
   float ProbeTracks_eta[kProbeTracksMax];
   float ProbeTracks_phi[kProbeTracksMax];
   float ProbeTracks_mass[kProbeTracksMax];
   int ProbeTracks_pdgId[kProbeTracksMax];
   float ProbeTracks_DCASig[kProbeTracksMax];
   float ProbeTracks_vz[kProbeTracksMax];
   int ProbeTracks_isMatchedToEle[kProbeTracksMax];
   int ProbeTracks_isMatchedToMuon[kProbeTracksMax];
   
   uint nTriggerMuon;   
   float TriggerMuon_pt[kTriggerMuonMax];
   float TriggerMuon_eta[kTriggerMuonMax];
   float TriggerMuon_phi[kTriggerMuonMax];
   float TriggerMuon_mass[kTriggerMuonMax];
   float TriggerMuon_vz[kTriggerMuonMax]; 
      
   uint nElectron;
   int Electron_charge[kElectronMax];
   float Electron_pt[kElectronMax];
   float Electron_eta[kElectronMax];
   float Electron_phi[kElectronMax];
   float Electron_mass[kElectronMax];
   float Electron_dxy[kElectronMax];
   float Electron_vz[kElectronMax];
   int Electron_isPFoverlap[kElectronMax];
   int Electron_isPF[kElectronMax];
   int Electron_isLowPt[kElectronMax];
   
   uint nMuon;
   int Muon_charge[kMuonMax];
   float Muon_pt[kMuonMax];
   float Muon_eta[kMuonMax];
   float Muon_phi[kMuonMax];
   float Muon_mass[kMuonMax];
   float Muon_dxy[kMuonMax];
   float Muon_vz[kMuonMax];
   int Muon_isTriggering[kMuonMax];
   int Muon_isPFcand[kMuonMax];
   float Muon_pfRelIso04_all[kMuonMax];
   bool Muon_softId[kMuonMax];
   bool Muon_mediumId[kMuonMax];
   
   uint nGenPart;
   int GenPart_pdgId[kGenPartMax];
   int GenPart_genPartIdxMother[kGenPartMax];
   float GenPart_pt[kGenPartMax];
   float GenPart_eta[kGenPartMax];
   float GenPart_phi[kGenPartMax];

   uint nTrigObj;
   int TrigObj_id[kTrigObjMax];
   float TrigObj_pt[kTrigObjMax];
   float TrigObj_eta[kTrigObjMax];
   float TrigObj_phi[kTrigObjMax];
   int TrigObj_filterBits[kTrigObjMax];

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
  _tree->SetBranchAddress("Muon_vz",&Muon_vz);
  _tree->SetBranchAddress("Muon_isTriggering",&Muon_isTriggering);
  _tree->SetBranchAddress("Muon_isPFcand",&Muon_isPFcand);
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
  _tree->SetBranchAddress("Electron_vz",&Electron_vz);
  _tree->SetBranchAddress("Electron_isPFoverlap",&Electron_isPFoverlap);
  _tree->SetBranchAddress("Electron_isPF",&Electron_isPF);
  _tree->SetBranchAddress("Electron_isLowPt",&Electron_isLowPt);

  
  int BToKEE_info = _tree->SetBranchAddress("nBToKEE",&nBToKEE);
  if(BToKEE_info>=0){
    _tree->SetBranchAddress("BToKEE_l1Idx",&BToKEE_l1Idx);
    _tree->SetBranchAddress("BToKEE_l2Idx",&BToKEE_l2Idx);
    _tree->SetBranchAddress("BToKEE_kIdx",&BToKEE_kIdx);
    _tree->SetBranchAddress("BToKEE_svprob",&BToKEE_svprob);
  }

  
  int BToKMuMu_info = _tree->SetBranchAddress("nBToKMuMu",&nBToKMuMu);
  if(BToKMuMu_info>=0){
    _tree->SetBranchAddress("BToKMuMu_l1Idx",&BToKMuMu_l1Idx);
    _tree->SetBranchAddress("BToKMuMu_l2Idx",&BToKMuMu_l2Idx);
    _tree->SetBranchAddress("BToKMuMu_kIdx",&BToKMuMu_kIdx);
    _tree->SetBranchAddress("BToKMuMu_svprob",&BToKMuMu_svprob);
  }  

  
  int Tracks_info = _tree->SetBranchAddress("nProbeTracks",&nProbeTracks);
  if(Tracks_info){
    _tree->SetBranchAddress("ProbeTracks_pt",&ProbeTracks_pt);
    _tree->SetBranchAddress("ProbeTracks_eta",&ProbeTracks_eta);
    _tree->SetBranchAddress("ProbeTracks_phi",&ProbeTracks_phi);
    _tree->SetBranchAddress("ProbeTracks_mass",&ProbeTracks_mass);
    _tree->SetBranchAddress("ProbeTracks_pdgId",&ProbeTracks_pdgId);
    _tree->SetBranchAddress("ProbeTracks_DCASig",&ProbeTracks_DCASig);
    _tree->SetBranchAddress("ProbeTracks_vz",&ProbeTracks_vz);
    _tree->SetBranchAddress("ProbeTracks_isMatchedToEle",&ProbeTracks_isMatchedToEle);
    _tree->SetBranchAddress("ProbeTracks_isMatchedToMuon",&ProbeTracks_isMatchedToMuon);
  }  

  
  int TriggerMuon_info = _tree->SetBranchAddress("nTriggerMuon",&nTriggerMuon);
  if(TriggerMuon_info){
    _tree->SetBranchAddress("TriggerMuon_pt",&TriggerMuon_pt);
    _tree->SetBranchAddress("TriggerMuon_eta",&TriggerMuon_eta);
    _tree->SetBranchAddress("TriggerMuon_phi",&TriggerMuon_phi);
    _tree->SetBranchAddress("TriggerMuon_mass",&TriggerMuon_mass);
    _tree->SetBranchAddress("TriggerMuon_vz",&TriggerMuon_vz);
  }  


  int isMC = _tree->SetBranchAddress("nGenPart",&nGenPart);
  if(isMC>=0){
    _tree->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
    _tree->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);
    _tree->SetBranchAddress("GenPart_pt",&GenPart_pt);
    _tree->SetBranchAddress("GenPart_eta",&GenPart_eta);
    _tree->SetBranchAddress("GenPart_phi",&GenPart_phi);
  }

   
  _tree->SetBranchAddress("nTrigObj",&nTrigObj);
  _tree->SetBranchAddress("TrigObj_id",&TrigObj_id);
  _tree->SetBranchAddress("TrigObj_pt",&TrigObj_pt);
  _tree->SetBranchAddress("TrigObj_eta",&TrigObj_eta);
  _tree->SetBranchAddress("TrigObj_phi",&TrigObj_phi);
  _tree->SetBranchAddress("TrigObj_filterBits",&TrigObj_filterBits);
  
}


Int_t NanoAODTree::GetEntry(int entry)
{

  int out = _tree->GetEntry(entry);

  if(nMuon>kMuonMax) return -1;
  if(nElectron>kElectronMax) return -1;
  if(nBToKEE>kBToKEEMax) return -1;
  if(nBToKMuMu>kBToKMuMuMax) return -1;
  if(nGenPart>kGenPartMax) return -1;
  if(nTrigObj>kTrigObjMax)  return -1;
  if(nProbeTracks>kProbeTracksMax) return -1;
  if(nTriggerMuon>kTriggerMuonMax) return -1;

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
