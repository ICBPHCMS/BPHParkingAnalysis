// BToKstllNtupleProducer_BParkingNANO --isMC (0,1,2) --isResonant (0, 1, -1) --isEleFS (0, 1) --isKstFS (0, 1) --output ("outfile") --input ("inputFile")

// new features => gen_index refers to the closest in dR - provided gen B and decays are within acceptance -, 
//                 no minimum dR required
//                 saved branch with dR of gen-reco lep1, lep2, kaon
//              => for gen_tag_muon do not require minimum pT
//                 saved branch with highest pT of extra muon (tag candidate)
//              => matching for reco tag muon is dR < 0.01
//              => isMC == 2 is for old - Thomas' - MC; isMC == 1 is new; isMC == 0 is DATA 
// 

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "NanoAODTree.h"


float JPsiMass_ = 3.0969;
float BuMass_ = 5.279;
float KaonMass_ = 0.493677;
float MuonMass_ = 0.10565837;
float ElectronMass_ = 0.5109989e-3;


float maxEtacceptance_ = 2.6; //2.4
float minPtacceptance_ = 0.2;  //1.

bool comparePairs(const std::pair<int, float>& i, const std::pair<int, float>& j){
  return i.second > j.second;
}


int main(int argc, char **argv){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isMC = -1;
  int isResonant = -1;
  int isEleFinalState = -1;
  int isKstFinalState = -1;
  string output = "";
  string input = "";
  string listFilesTXT = "";
  bool inputTXT = false;
  float genMuPtCut = 5.;
  bool overwrite = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isMC") {
      if (i + 1 < argc) {
	isMC = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isMC option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isResonant") {
      if (i + 1 < argc) {
	isResonant = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isResonant option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isEleFS") {
      if (i + 1 < argc) {
	isEleFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isElFinalState option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isKstFS") {
      if (i + 1 < argc) {
	isKstFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isKstFinalState option requires one argument " << std::endl;
	return 1;
      }
    }
  }
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--output") {
      if (i + 1 < argc) {
	output = argv[i+1];
	break;
      } else {
	std::cerr << "--output option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--input") {
      if (i + 1 < argc) {
        input = argv[i+1];
        break;
      } else {
	std::cerr << "--intput option requires one argument." << std::endl;
        return 1;
      }
    }
    else if(std::string(argv[i]) == "--inputTXT") {
      if (i + 1 < argc) {
        inputTXT = true;
        listFilesTXT = argv[i+1];
        break;
      } else {
	std::cerr << "--intputTXT option requires one argument." << std::endl;
        return 1;
      }
    }
  }

  if(output == ""){
    std::cerr << "--output argument required" << std::endl;
    return 1;
  }
  if(listFilesTXT == "" && input == ""){
    std::cerr << "--input argument required" << std::endl;
    return 1;
  }

  /*
  bool overwrite = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }
  */



  //Always saveFullNanoAOD info because we are skimming
  bool saveFullNanoAOD = true;
  /*for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveFullNanoAOD") {
      saveFullNanoAOD = true;
      break;
    }
    }*/
 


  if(isKstFinalState){
    std::cout << " implementation missing " << std::endl;
    return -10;
  }

  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TChain* oldtree = new TChain("Events");
  if(inputTXT){
    string reader;
    std::ifstream inFileLong;
    inFileLong.open(listFilesTXT.c_str(), std::ios::in);

    while(!inFileLong.eof()){
      inFileLong >> reader;
      if( inFileLong.eof() ) break;
      std::cout << " Adding " << reader << std::endl;
      oldtree->Add(reader.c_str());
    }
  }
  else   oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BToKstllTree","BToKstllTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches
  int _BToKstll_sel_index_KEE = -1;
  int _BToKstll_sel_index_KMuMu = -1;
  std::vector<int> _BToKstll_order_index_KEE;
  std::vector<int> _BToKstll_order_index_KMuMu;
  int _Muon_sel_index_KEE = -1; //Probe muon with selection algo.
  int _Muon_sel_index_KMuMu = -1;
  std::vector<int> _Muon_tag_index_KEE; //Probe muon with selection algo.
  std::vector<int> _Muon_tag_index_KMuMu;
  int _Muon_probe_index = -1; //Probe muon for Acc.xEff. = _Muon_sel_index in data

  tree_new->Branch("BToKstll_sel_index_KEE",&_BToKstll_sel_index_KEE,"BToKstll_sel_index_KEE/I");
  tree_new->Branch("BToKstll_sel_index_KMuMu",&_BToKstll_sel_index_KMuMu,"BToKstll_sel_index_KMuMu/I");
  tree_new->Branch("BToKstll_order_index_KEE", &_BToKstll_order_index_KEE);
  tree_new->Branch("BToKstll_order_index_KMuMu", &_BToKstll_order_index_KMuMu);
  tree_new->Branch("Muon_sel_index_KEE",&_Muon_sel_index_KEE,"Muon_sel_index_KEE/I");
  tree_new->Branch("Muon_sel_index_KMuMu",&_Muon_sel_index_KMuMu,"Muon_sel_index_KMuMu/I");
  tree_new->Branch("Muon_tag_index_KEE",&_Muon_tag_index_KEE);
  tree_new->Branch("Muon_tag_index_KMuMu",&_Muon_tag_index_KMuMu);
  tree_new->Branch("Muon_probe_index",&_Muon_probe_index,"Muon_probe_index/I");


  int _GenPart_BToKstll_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_lep1FromB_index = -1;
  int _GenPart_lep2FromB_index = -1;
  int _BToKstll_gen_index = -1;  //index of reco closest to gen (each daughter within 0.1)
  float _BToKstll_gendR_lep1FromB;
  float _BToKstll_gendR_lep2FromB;
  float _BToKstll_gendR_KFromB;
  float _BToKstll_gen_llMass = -1;
  float _BToKstll_gen_mass = -1;
  float _BToKstll_gen_muonTag_hpT = -1;

  //the 3 following are the index of the reco closest to the gen (dR < 0.1)
  //only for no LT
  int _Lep1_gen_index = -1;
  int _Lep2_gen_index = -1;
  int _KPFCand_gen_index = -1;  

  if(isMC != 0){
    tree_new->Branch("GenPart_BToKstll_index",&_GenPart_BToKstll_index,"GenPart_BToKstll_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_lep1FromB_index",&_GenPart_lep1FromB_index,"GenPart_lep1FromB_index/I");
    tree_new->Branch("GenPart_lep2FromB_index",&_GenPart_lep2FromB_index,"GenPart_lep2FromB_index/I");
    tree_new->Branch("BToKstll_gen_index",&_BToKstll_gen_index,"BToKstll_gen_index/I");
    tree_new->Branch("BToKstll_gendR_lep1FromB",&_BToKstll_gendR_lep1FromB,"BToKstll_gendR_lep1FromB/F");
    tree_new->Branch("BToKstll_gendR_lep2FromB",&_BToKstll_gendR_lep2FromB,"BToKstll_gendR_lep2FromB/F");
    tree_new->Branch("BToKstll_gendR_KFromB",&_BToKstll_gendR_KFromB,"BToKstll_gendR_KFromB/F");
    tree_new->Branch("BToKstll_gen_llMass",&_BToKstll_gen_llMass,"BToKstll_gen_llMass/F");
    tree_new->Branch("BToKstll_gen_mass",&_BToKstll_gen_mass,"BToKstll_gen_mass/F");
    tree_new->Branch("BToKstll_gen_muonTag_hpT",&_BToKstll_gen_muonTag_hpT,"BToKstll_gen_muonTag_hpT/F");
    tree_new->Branch("Lep1_gen_index",&_Lep1_gen_index,"Lep1_gen_index/I");
    tree_new->Branch("Lep2_gen_index",&_Lep2_gen_index,"Lep2_gen_index/I");
    tree_new->Branch("KPFCand_gen_index",&_KPFCand_gen_index,"KPFCand_gen_index/I");
  }



  bool _Muon_isHLT_BPHParking[kMuonMax];
  tree_new->Branch("Muon_isHLT_BPHParking",_Muon_isHLT_BPHParking,"Muon_isHLT_BPHParking[nMuon]/O");
  

  int nentries = tree->GetEntries();
  std::cout << " Nentries = " << nentries << std::endl;
  std::cout << " isMC = " << isMC << std::endl;
  
  for (int iEntry = 0; iEntry < nentries; ++iEntry){
    
    bool debug = false;
    //if(iEntry == 419) debug = true;

    int out = tree->GetEntry(iEntry);
    if(out<0){
      std::cout << " Error retrievieng entry #" << iEntry << std::endl;
      return -1;
    }

    if(iEntry%10000==0) std::cout << " Entry #" << iEntry << " " << int(100*float(iEntry)/nentries) << "%" << std::endl;

    if(debug)    std::cout << " iEntry = " << iEntry << std::endl;

    _BToKstll_sel_index_KEE = -1;
    _BToKstll_sel_index_KMuMu = -1;
    _BToKstll_order_index_KEE.clear();
    _BToKstll_order_index_KMuMu.clear();
    _Muon_sel_index_KEE = -1; //tag muon
    _Muon_sel_index_KMuMu = -1;
    _Muon_tag_index_KEE.clear();
    _Muon_tag_index_KMuMu.clear();
    _Muon_probe_index = -1;

    _GenPart_BToKstll_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_lep1FromB_index = -1;
    _GenPart_lep2FromB_index = -1;
    _BToKstll_gen_index = -1; 
    _BToKstll_gendR_lep1FromB = -1; 
    _BToKstll_gendR_lep2FromB = -1; 
    _BToKstll_gendR_KFromB = -1; 
    _BToKstll_gen_llMass = -1;
    _BToKstll_gen_mass = -1;
    _BToKstll_gen_muonTag_hpT = -1;
    _Lep1_gen_index = -1;
    _Lep2_gen_index = -1;
    _KPFCand_gen_index = -1;

    for( int isEleCh=0; isEleCh<2; isEleCh++ ){//loop over KEE (1) and KMuMu (0)

      //Select the BToKll candidate with reco criteria
      int nBinTree = (isEleCh == 1) ? tree->nBToKEE : tree->nBToKMuMu;
      float best_B_CL_vtx = -1.;;
      std::vector<std::pair<int, float>> B_vtxCL_idx_val;
      _Muon_tag_index_KEE.resize(nBinTree);
      _Muon_tag_index_KMuMu.resize(nBinTree);


      if(debug) std::cout << " >>> nBinTree = " << nBinTree << std::endl;

    

      for(int i_Btree=0; i_Btree<nBinTree; ++i_Btree){//loop over triplets            
      
	int l1_Index = (isEleCh == 1) ? tree->BToKEE_l1Idx[i_Btree] : tree->BToKMuMu_l1Idx[i_Btree];
	int l2_Index = (isEleCh == 1) ? tree->BToKEE_l2Idx[i_Btree] : tree->BToKMuMu_l2Idx[i_Btree];
	int kaon_Index = (isEleCh == 1) ? tree->BToKEE_kIdx[i_Btree] : tree->BToKMuMu_kIdx[i_Btree];
      
	//KEE cleaning
	if( isEleCh == 1 ){
	  if( tree->Electron_isPFoverlap[l1_Index] == 1 || tree->Electron_isPFoverlap[l2_Index] == 1 ) continue;
          if( tree->ProbeTracks_isMatchedToEle[kaon_Index] == 1 ) continue;
	}
      
	if( isEleCh == 0 ){
	  //probe muon different from trigger muon
          if( tree->Muon_isTriggering[l1_Index] == 1 || tree->Muon_isTriggering[l2_Index] == 1 ) continue;
          //KMuMu cleaning
          if( tree->ProbeTracks_isMatchedToMuon[kaon_Index] == 1 ) continue;
	}      
      
	//require lepton-1 charge * lepton-2 charge < 0
	int l1_charge = (isEleCh == 1) ? tree->Electron_charge[l1_Index] : tree->Muon_charge[l1_Index];
	int l2_charge = (isEleCh == 1) ? tree->Electron_charge[l2_Index] : tree->Muon_charge[l2_Index];
	if(l1_charge * l2_charge > 0.) continue;
      
	//consider only triplets with trigger muon
      
	if(isEleCh == 1) _Muon_tag_index_KEE[i_Btree] = -1;
	else _Muon_tag_index_KMuMu[i_Btree] = -1;      
      
	int nTriggerMuon = tree->nTriggerMuon;

      
	for(int i_mu=0; i_mu < nTriggerMuon; i_mu++){//loop over trigger muons
          
	  TLorentzVector lep1_tlv;
	  TLorentzVector lep2_tlv;
	  TLorentzVector kaon_tlv;
	  TLorentzVector TriggerMuon_tlv;
	
	  lep1_tlv.SetPtEtaPhiM( (isEleCh == 1) ? tree->Electron_pt[l1_Index] : tree->Muon_pt[l1_Index],
				 (isEleCh == 1) ? tree->Electron_eta[l1_Index] : tree->Muon_eta[l1_Index],
				 (isEleCh == 1) ? tree->Electron_phi[l1_Index] : tree->Muon_phi[l1_Index],
				 (isEleCh == 1) ? ElectronMass_ : MuonMass_);
	  lep2_tlv.SetPtEtaPhiM( (isEleCh == 1) ? tree->Electron_pt[l2_Index] : tree->Muon_pt[l2_Index],
				 (isEleCh == 1) ? tree->Electron_eta[l2_Index] : tree->Muon_eta[l2_Index],
				 (isEleCh == 1) ? tree->Electron_phi[l2_Index] : tree->Muon_phi[l2_Index],
				 (isEleCh == 1) ? ElectronMass_ : MuonMass_);
	  kaon_tlv.SetPtEtaPhiM(tree->ProbeTracks_pt[kaon_Index],
				tree->ProbeTracks_eta[kaon_Index],
				tree->ProbeTracks_phi[kaon_Index],
				tree->ProbeTracks_mass[kaon_Index]);        
	  TriggerMuon_tlv.SetPtEtaPhiM(tree->TriggerMuon_pt[i_mu],
				       tree->TriggerMuon_eta[i_mu],
				       tree->TriggerMuon_phi[i_mu],
				       MuonMass_);
	
	  float dR_lep1FromTrMu = lep1_tlv.DeltaR(TriggerMuon_tlv);
	  float dR_lep2FromTrMu = lep2_tlv.DeltaR(TriggerMuon_tlv);   
	  float dR_kaonFromTrMu = kaon_tlv.DeltaR(TriggerMuon_tlv);
        
	  float l1_vz = (isEleCh == 1) ? tree->Electron_vz[l1_Index] : tree->Muon_vz[l1_Index];
	  float l2_vz = (isEleCh == 1) ? tree->Electron_vz[l2_Index] : tree->Muon_vz[l2_Index];

	  //On each triplet object require dR>0.4 and vz<1
	  if( dR_lep1FromTrMu < 0.4 || dR_lep2FromTrMu < 0.4 || dR_kaonFromTrMu < 0.4 || 
	      std::abs(tree->TriggerMuon_vz[i_mu] - l1_vz) > 1. || 
	      std::abs(tree->TriggerMuon_vz[i_mu] - l2_vz) > 1. || 
	      std::abs(tree->TriggerMuon_vz[i_mu] - tree->ProbeTracks_vz[kaon_Index]) > 1. ) continue;        

	  if(isEleCh == 1) _Muon_tag_index_KEE[i_Btree] = i_mu;
	  else _Muon_tag_index_KMuMu[i_Btree] = i_mu;
	  
	  break;        
	}//loop over trigger muons
    
    
	//force rank only among triplets with extra trigger muon
	float B_CL_vtx = ( (isEleCh == 1) ? tree->BToKEE_svprob[i_Btree] : tree->BToKMuMu_svprob[i_Btree] ) +
	  ( (isEleCh == 1) ? ((_Muon_tag_index_KEE[i_Btree] == -1) ? -1 : 0) 
	    : ((_Muon_tag_index_KMuMu[i_Btree] == -1) ? -1 : 0) );
      
	B_vtxCL_idx_val.push_back(std::pair<int, float>(i_Btree, B_CL_vtx));
      
	if( best_B_CL_vtx < 0. || B_CL_vtx>best_B_CL_vtx ){
	  best_B_CL_vtx = B_CL_vtx;
	  if(isEleCh == 1) _BToKstll_sel_index_KEE = i_Btree;
	  else _BToKstll_sel_index_KMuMu = i_Btree;
	}
      
      }//loop over triplets



      std::sort(B_vtxCL_idx_val.begin(), B_vtxCL_idx_val.end(), comparePairs);
      int ipCount = 0;
      if(isEleCh == 1) _BToKstll_order_index_KEE.resize(B_vtxCL_idx_val.size());
      else _BToKstll_order_index_KMuMu.resize(B_vtxCL_idx_val.size());
      for(auto ip : B_vtxCL_idx_val){
	if(isEleCh == 1) _BToKstll_order_index_KEE[ipCount] = ip.first;
	else _BToKstll_order_index_KMuMu[ipCount] = ip.first;
	++ipCount;
      }
    
      B_vtxCL_idx_val.clear();
    
      int BToKstll_sel_index = (isEleCh == 1) ? _BToKstll_sel_index_KEE : _BToKstll_sel_index_KMuMu;

      if(debug) std::cout << " isEleCh " << isEleCh << " BToKstll_sel_index = " << BToKstll_sel_index << std::endl;

      if(isMC == 0 && BToKstll_sel_index<0) continue;

      //assign Muon_sel_index for triplet with the best CL(B-vtx)
      if(BToKstll_sel_index>=0){
	if(isEleCh == 1) _Muon_sel_index_KEE = _Muon_tag_index_KEE[BToKstll_sel_index];
	else _Muon_sel_index_KMuMu = _Muon_tag_index_KMuMu[BToKstll_sel_index];
      }
      if(debug) std::cout << " isEleCh " << isEleCh << " iEntry = " << iEntry << " _Muon_sel_index = " << 
		  ((isEleCh == 1) ? _Muon_sel_index_KEE : _Muon_sel_index_KMuMu) << std::endl;

      //!!! Can have selected B->Kstll even without additional tag muons
      //Needs to ask both BToKstll_sel_index>=0 && Muon_sel_index>=0 for analysis on probe side
    
    }//loop over KEE (1) and KMuMu (0)
    
    


    //Select the BToKstll candidate based on gen matching

    if(isMC != 0){
      int nGenPart = tree->nGenPart;

      int leptonID = (isEleFinalState == 1) ? 11 : 13;
      
      if(isResonant){
	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_lep1FromB_index = -1;
	  _GenPart_lep2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){
	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];

	      int lep1_index = -1;
	      int lep2_index = -1;

	      if(abs(pdgId)==443 && mother_index == i_Bu){

		for(int j_gen=0; j_gen<nGenPart; j_gen++){
		  int pdgId = tree->GenPart_pdgId[j_gen];
		  float partPt = tree->GenPart_pt[i_gen];
		  float partEta = tree->GenPart_eta[i_gen];
		  if(partPt < minPtacceptance_) continue;
		  if(std::abs(partEta) > maxEtacceptance_) continue;
		  int mother_index = tree->GenPart_genPartIdxMother[j_gen];
		  if(abs(pdgId)==leptonID && mother_index == i_gen && lep1_index < 0)
		    lep1_index = j_gen;
		  else if(abs(pdgId)==leptonID && mother_index == i_gen)
		    lep2_index = j_gen;
		  if(lep1_index >= 0 && lep2_index >= 0) break;
		}

		if(lep1_index >= 0 && lep2_index >= 0){
		  _GenPart_JPsiFromB_index = i_gen;
		  _GenPart_lep1FromB_index = lep1_index;
		  _GenPart_lep2FromB_index = lep2_index;
		}
		else break;

	      }

	      else if(abs(pdgId)==321 && mother_index == i_Bu) _GenPart_KFromB_index = i_gen;

	      else if(mother_index == i_Bu) break; //Additional B decay products
	    }
	  }

	  if(_GenPart_JPsiFromB_index >= 0 && _GenPart_KFromB_index >=0){
	    _GenPart_BToKstll_index = i_Bu;
	    break;
	  }
	}	
      }//resonant
      else{ 

	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){
	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_lep1FromB_index = -1;
	  _GenPart_lep2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu]) == 521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      float partPt = tree->GenPart_pt[i_gen];
	      float partEta = tree->GenPart_eta[i_gen];
	      if(partPt < minPtacceptance_) continue;
	      if(std::abs(partEta) > maxEtacceptance_) continue;
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	      if(abs(pdgId) == leptonID && mother_index == i_Bu && _GenPart_lep1FromB_index < 0)
		_GenPart_lep1FromB_index = i_gen;
	      else if(abs(pdgId)==leptonID && mother_index == i_Bu)
		_GenPart_lep2FromB_index = i_gen;
	      else if(abs(pdgId)==321 && mother_index == i_Bu)
		_GenPart_KFromB_index = i_gen;
	      else if(mother_index == i_Bu) break;
	    }
	  }//if B

	  if(_GenPart_lep1FromB_index >= 0 && _GenPart_lep2FromB_index >= 0 && _GenPart_KFromB_index >= 0){
	    _GenPart_BToKstll_index = i_Bu;
	    break;
	  }
	} //loop over B
      }// non resonant
    
      if(_GenPart_BToKstll_index >= 0){

	//lep1FromB stored a leading daughter
	if(tree->GenPart_pt[_GenPart_lep2FromB_index] > tree->GenPart_pt[_GenPart_lep1FromB_index]){
	  int i_temp = _GenPart_lep1FromB_index;
	  _GenPart_lep1FromB_index = _GenPart_lep2FromB_index;
	  _GenPart_lep2FromB_index = i_temp;
	}

	TLorentzVector gen_KFromB_tlv;
	TLorentzVector gen_lep1FromB_tlv;
	TLorentzVector gen_lep2FromB_tlv;
	
	gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				    tree->GenPart_eta[_GenPart_KFromB_index],
				    tree->GenPart_phi[_GenPart_KFromB_index],
				    KaonMass_);
	gen_lep1FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_lep1FromB_index],
				      tree->GenPart_eta[_GenPart_lep1FromB_index],
				      tree->GenPart_phi[_GenPart_lep1FromB_index],
				       (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	gen_lep2FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_lep2FromB_index],
				       tree->GenPart_eta[_GenPart_lep2FromB_index],
				       tree->GenPart_phi[_GenPart_lep2FromB_index],
				       (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);

	_BToKstll_gen_llMass = (gen_lep1FromB_tlv+gen_lep2FromB_tlv).Mag();
	_BToKstll_gen_mass = (gen_lep1FromB_tlv+gen_lep2FromB_tlv+gen_KFromB_tlv).Mag();
	
	float best_dR = -1.;
    
	int nBinTree = isEleFinalState ? tree->nBToKEE : tree->nBToKMuMu;
	
	for(int i_Btree=0; i_Btree<nBinTree; ++i_Btree){

	  TLorentzVector kaon_tlv;
	  TLorentzVector lep1_tlv;
	  TLorentzVector lep2_tlv;

	  kaon_tlv.SetPtEtaPhiM(tree->BToKstll_kaon_pt[i_Btree],
				tree->BToKstll_kaon_eta[i_Btree],
				tree->BToKstll_kaon_phi[i_Btree],
				KaonMass_);
	  lep1_tlv.SetPtEtaPhiM(tree->BToKstll_lep1_pt[i_Btree],
				tree->BToKstll_lep1_eta[i_Btree],
				tree->BToKstll_lep1_phi[i_Btree],
				(isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	  lep2_tlv.SetPtEtaPhiM(tree->BToKstll_lep2_pt[i_Btree],
				tree->BToKstll_lep2_eta[i_Btree],
				tree->BToKstll_lep2_phi[i_Btree],
				(isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	  
	  float dR_KFromB = kaon_tlv.DeltaR(gen_KFromB_tlv);
	  float dR_lep1FromB = min(lep1_tlv.DeltaR(gen_lep1FromB_tlv), lep2_tlv.DeltaR(gen_lep1FromB_tlv));
	  float dR_lep2FromB = min(lep1_tlv.DeltaR(gen_lep2FromB_tlv), lep2_tlv.DeltaR(gen_lep2FromB_tlv));
	  //Should check that same objects not selected twice

	  float dR_tot = dR_KFromB + dR_lep1FromB + dR_lep2FromB; //In case several BToKmumu matches, take the closest one in dR_tot

	  //if( dR_lep1FromB <0.1 && dR_lep2FromB <0.1
	  //  && 
	  if((best_dR <0. || dR_tot < best_dR) && dR_tot < 0.1){
	    best_dR = dR_tot;
	    _BToKstll_gen_index = i_Btree;
	    _BToKstll_gendR_lep1FromB = dR_lep1FromB;
	    _BToKstll_gendR_lep2FromB = dR_lep2FromB;
	    _BToKstll_gendR_KFromB = dR_KFromB;
	  }
	}

	/*
	float best_dR_lep1FromB = -1.;
	float best_dR_lep2FromB = -1.;
	float best_dR_KFromB = -1.;
	*/

	//for the moment just check best reco closest to gen
	//for lepton-lepton combination only
	/*
	if(!isLeptonTrack){  
	for(int i_mu=0; i_mu<nMuon; i_mu++){
	  TLorentzVector mu_tlv;
	  mu_tlv.SetPtEtaPhiM(tree->Muon_pt[i_mu],
			      tree->Muon_eta[i_mu],
			      tree->Muon_phi[i_mu],
			      MuonMass_);
	  float dR_mu1FromB = mu_tlv.DeltaR(gen_mu1FromB_tlv);
	  float dR_mu2FromB = mu_tlv.DeltaR(gen_mu2FromB_tlv);
	  if(dR_mu1FromB<0.1 && (best_dR_mu1FromB<0. || dR_mu1FromB<best_dR_mu1FromB)){
	    _Muon_mu1FromB_index = i_mu;
	    best_dR_mu1FromB = dR_mu1FromB;
	  }
	  if(dR_mu2FromB<0.1 && (best_dR_mu2FromB<0. || dR_mu2FromB<best_dR_mu2FromB)){
	    _Muon_mu2FromB_index = i_mu;
	    best_dR_mu2FromB = dR_mu2FromB;
	  }
	}
	}//end of lepton-lepton comb
	int nPFCand = tree->nPFCand;
	for(int i_pf=0;i_pf<nPFCand;i_pf++){
	  if(abs(tree->PFCand_pdgId[i_pf])!=211) continue;
	  TLorentzVector PFCand_tlv;
	  PFCand_tlv.SetPtEtaPhiM(tree->PFCand_pt[i_pf],tree->PFCand_eta[i_pf],tree->PFCand_phi[i_pf],tree->PFCand_mass[i_pf]);
	  float dR_KFromB = PFCand_tlv.DeltaR(gen_KFromB_tlv);
	  if(dR_KFromB<0.1 && (best_dR_KFromB<0. || dR_KFromB<best_dR_KFromB)){
	    _PFCand_genKFromB_index = i_pf;
	    best_dR_KFromB = dR_KFromB;
	  }
	}
	*/

	//Gen muon pt filter on tag
	bool isTagMuonHighPt = false;
	for(int i_gen=0; i_gen<nGenPart; i_gen++){

	  int pdgId = tree->GenPart_pdgId[i_gen];
	  if(abs(pdgId)!=13) continue;

	  //exclude the gen muons that are from the chosen B
	  if(!isEleFinalState && 
	     (i_gen == _GenPart_lep1FromB_index || i_gen == _GenPart_lep2FromB_index)) continue;
	  	  
	  TLorentzVector gen_tagMu_tlv;
	  gen_tagMu_tlv.SetPtEtaPhiM(tree->GenPart_pt[i_gen],
				     tree->GenPart_eta[i_gen],
				     tree->GenPart_phi[i_gen],
				     MuonMass_);

	  //In case there are several copies of same muon (with FSR for instance)
	  if(gen_tagMu_tlv.DeltaR(gen_lep1FromB_tlv) > 0.01 && gen_tagMu_tlv.DeltaR(gen_lep2FromB_tlv) > 0.01 &&
	     gen_tagMu_tlv.Pt() > genMuPtCut && gen_tagMu_tlv.Pt() > _BToKstll_gen_muonTag_hpT){
	    _BToKstll_gen_muonTag_hpT = gen_tagMu_tlv.Pt();
	    _Muon_probe_index = i_gen;
	    isTagMuonHighPt = true;
	    //break;
	  }
	}
	if(!isTagMuonHighPt) continue; //Skip events where the gen filter in MC has been applied to the probe muon
      }


      
      //Require tag muon non matched to gen muons that are from B
      /*    
      if(!isEleFianlState && !isLeptonTrack){
	bool isTagMuonSoftID = false;
	for(int i_mu=0; i_mu<nMuon; ++i_mu){
	  if(i_mu == _Lep_mu1FromB_index || i_mu==_Muon_mu2FromB_index) continue;
	  if(tree->Muon_softId[i_mu] && tree->Muon_pt[i_mu] > 8.){
	    isProbeMuonSoftID = true;
	    _Muon_selgen_index = i_mu;
	    break;
	  }
	  
	}
	if(!isProbeMuonSoftID) continue; //Skip events where there is no probe muon passing the soft ID
      }
      */
    }


    /*
    if(isMC == 0){
      _Muon_probe_index = _Muon_sel_index;
    }
    */
    
    tree_new->Fill();
  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}


 
