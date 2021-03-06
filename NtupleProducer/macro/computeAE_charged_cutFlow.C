//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//electron final state
//root -l computeAE_charged_cutFlow.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKee.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsiee.root", 1)' 
//muon final state
//root -l computeAE_charged_cutFlow.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKmumu.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsimumu.root", 0)' 


//root -l computeAE_charged_cutFlow.C'("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKee_18_08_14.root" , "/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKJPsiee_18_08_14.root", 1)'
//root -l computeAE_charged_cutFlow.C'("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKmumu_18_09_3.root" , "/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKJPsimumu_18_09_3.root", 0)'

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

void computeAE_charged_cutFlow(std::string nonResonantFile, std::string ResonantFile, int isEleFinalState){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add(nonResonantFile.c_str());
  t2->Add(ResonantFile.c_str());


  //muon tag with soft ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  int nMuonTag_nonReso = t1->GetEntries();
  int nMuonTag_Reso = t2->GetEntries();

  std::cout << " #muon tag events: nnreso = " << nMuonTag_nonReso << " resonant = " << nMuonTag_Reso << std::endl;

  //if false => move to reco ee invariant mass bins

  //selections
  std::string cut_muonTag = "Muon_sel_index != -1";
  // should require that muonTag does not match with gen probes?
  std::string cut_genAcc = "GenPart_e2FromB_index != -1 && GenPart_pt[GenPart_e2FromB_index] > 1. && abs(GenPart_eta[GenPart_e2FromB_index]) < 2.4 && GenPart_e1FromB_index != -1 && GenPart_pt[GenPart_e1FromB_index] > 1. && abs(GenPart_eta[GenPart_e1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
  std::string cut_genEff = "BToKee_gen_index != -1 && BToKee_gen_index == BToKee_sel_index";

  std::string cut_chargeEff = "BToKee_sel_index != -1 && BToKee_ele1_charge[BToKee_sel_index]*BToKee_ele2_charge[BToKee_sel_index] < 0. && BToKee_gen_index != -1 && BToKee_gen_index == BToKee_sel_index";
  std::string cut_pTEff = "BToKee_pt[BToKee_sel_index] >= 10. && BToKee_kaon_pt[BToKee_sel_index] >= 1.5";
  std::string cut_alphaEff = "BToKee_cosAlpha[BToKee_sel_index] >= 0.999";
  std::string cut_vtxCLEff = "BToKee_CL_vtx[BToKee_sel_index] >= 0.1";  
  std::string cut_DCAEff  = "BToKee_Lxy[BToKee_sel_index] >= 6";

  if(!isEleFinalState){
    cut_muonTag = "Muon_probe_index != -1";
    cut_genAcc = "GenPart_mu2FromB_index != -1 && GenPart_pt[GenPart_mu2FromB_index] > 1. && abs(GenPart_eta[GenPart_mu2FromB_index]) < 2.4 && GenPart_mu1FromB_index != -1 && GenPart_pt[GenPart_mu1FromB_index] > 1. && abs(GenPart_eta[GenPart_mu1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
    cut_genEff = "BToKmumu_gen_index != -1 && BToKmumu_gen_index == BToKmumu_sel_index";

    cut_chargeEff = "BToKmumu_sel_index != -1 && BToKmumu_mu1_charge[BToKmumu_sel_index]*BToKmumu_mu2_charge[BToKmumu_sel_index] < 0. && BToKmumu_gen_index != -1 && BToKmumu_gen_index == BToKmumu_sel_index";
    cut_pTEff = "BToKmumu_pt[BToKmumu_sel_index] >= 10. && BToKmumu_kaon_pt[BToKmumu_sel_index] >= 1.5";
    cut_alphaEff = "BToKmumu_cosAlpha[BToKmumu_sel_index] >= 0.999";
    cut_vtxCLEff = "BToKmumu_CL_vtx[BToKmumu_sel_index] >= 0.1";
    cut_DCAEff  = "BToKmumu_Lxy[BToKmumu_sel_index] >= 6";
  }

  std::vector<std::string> llMassCut;
  std::vector<std::string> llGenMassCut;

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);
  llMassBoundary.push_back(100.);

  for(int ij=0; ij<6; ++ij){
    std::string cut = Form("BToKee_sel_index != -1 && BToKee_eeKFit_ee_mass[BToKee_sel_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_sel_index] < %.2f",
                           llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    std::string gencut = Form("BToKee_gen_eeMass > %.2f && BToKee_gen_eeMass < %.2f",
                              llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    if(!isEleFinalState){
      cut = Form("BToKmumu_sel_index != -1 && BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] > %.2f && BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] < %.2f",
		 llMassBoundary.at(ij), llMassBoundary.at(ij+1));
      gencut = Form("BToKmumu_gen_mumuMass > %.2f && BToKmumu_gen_mumuMass < %.2f",
		    llMassBoundary.at(ij), llMassBoundary.at(ij+1));

    }

    llMassCut.push_back(cut);
    llGenMassCut.push_back(gencut);
  }

  //non resonant first - JPsi last
  float nEv_muonTag[7] = {0.};
  float nEv_genAcc[7] = {0.};
  float nEv_genEff[7] = {0.};
  float nEv_recoEff[7] = {0.};  /* same as genEff, but computed in reco ee invarinat mass bins*/
  float nEv_chargeEff[7] = {0.};
  float nEv_pTEff[7] = {0.};
  float nEv_alphaEff[7] = {0.};
  float nEv_vtxCLEff[7] = {0.};
  float nEv_DCAEff[7] = {0.};

  for(int ij=0; ij<6; ++ij){
    nEv_muonTag[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_genAcc[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_genEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_recoEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");

    nEv_chargeEff[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_pTEff[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_alphaEff[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_vtxCLEff[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    nEv_DCAEff[ij] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str(), "goff");
  }
  nEv_muonTag[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_genAcc[6] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_genEff[6] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_recoEff[6] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_chargeEff[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_pTEff[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_alphaEff[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_vtxCLEff[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  nEv_DCAEff[6] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str(), "goff");
  


  //upgrate to 2D plot for future or something better than printout
  for(int ij=0; ij<6; ++ij){
    // std::cout << " \n mass bin \t nonResonant = " << llGenMassCut.at(ij) << "\t        \t Resonant (J/Psi) \t   \t \t" << std::endl;
    std::cout << " \n   \t nonResonant ll mass-bin [" << llMassBoundary.at(ij) << "-"<< llMassBoundary.at(ij+1) << "] \t       \t Resonant (J/Psi) \t   \t \t" << std::endl;
    std::cout << " muonTag = \t " << nEv_muonTag[ij] << " \t              \t \t        \t " << nEv_muonTag[6]  << " \t " << std::endl;

    /*
    std::cout << " genAcc = \t " << nEv_genAcc[ij] << " \t /muonTag = " << nEv_genAcc[ij]/nEv_muonTag[ij]
	      << "               \t " << nEv_genAcc[6] << " \t  /muonTag = " << nEv_genAcc[6]/nEv_muonTag[6] 
	      << " \t double ratio = \t" << (nEv_genAcc[ij]/nEv_muonTag[ij])/(nEv_genAcc[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x genEff = \t " << nEv_genEff[ij] << " \t  /muonTag = " << nEv_genEff[ij]/nEv_muonTag[ij]
	      << "                 \t " << nEv_genEff[6] << " \t  /muonTag = " << nEv_genEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_genEff[ij]/nEv_muonTag[ij])/(nEv_genEff[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x recoEff = \t " << nEv_recoEff[ij] << " \t  /muonTag = " << nEv_recoEff[ij]/nEv_muonTag[ij]
	      << "                 \t " << nEv_recoEff[6] << " \t  /muonTag = " << nEv_recoEff[6]/nEv_muonTag[6]
	      << " \t  double ratio = \t" << (nEv_recoEff[ij]/nEv_muonTag[ij])/(nEv_recoEff[6]/nEv_muonTag[6]) << "\n";
    */
    std::cout << " x recoCand = \t " << nEv_chargeEff[ij] << " \t  /muonTag = " << nEv_chargeEff[ij]/nEv_muonTag[ij]
	      << "                 \t " << nEv_chargeEff[6] << " \t  /muonTag = " << nEv_chargeEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_chargeEff[ij]/nEv_muonTag[ij])/(nEv_chargeEff[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x pT = \t " << nEv_pTEff[ij] << " \t  /muonTag = " << nEv_pTEff[ij]/nEv_muonTag[ij]
	      << "               \t " << nEv_pTEff[6] << " \t  /muonTag = " << nEv_pTEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_pTEff[ij]/nEv_muonTag[ij])/(nEv_pTEff[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x alpha = \t " << nEv_alphaEff[ij] << " \t  /muonTag = " << nEv_alphaEff[ij]/nEv_muonTag[ij]
	      << "               \t " << nEv_alphaEff[6] << " \t  /muonTag = " << nEv_alphaEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_alphaEff[ij]/nEv_muonTag[ij])/(nEv_alphaEff[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x vtxCL = \t " << nEv_vtxCLEff[ij] << " \t  /muonTag = " << nEv_vtxCLEff[ij]/nEv_muonTag[ij]
	      << "               \t " << nEv_vtxCLEff[6] << " \t  /muonTag = " << nEv_vtxCLEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_vtxCLEff[ij]/nEv_muonTag[ij])/(nEv_vtxCLEff[6]/nEv_muonTag[6]) << "\n";

    std::cout << " x IPsig = \t " << nEv_DCAEff[ij] << " \t  /muonTag = " << nEv_DCAEff[ij]/nEv_muonTag[ij]
	      << "                 \t " << nEv_DCAEff[6]  << " \t  /muonTag = " << nEv_DCAEff[6]/nEv_muonTag[6]
	      << " \t double ratio = \t" << (nEv_DCAEff[ij]/nEv_muonTag[ij])/(nEv_DCAEff[6]/nEv_muonTag[6]) << std::endl;
  }


}
