#include <TF1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TMath.h>
#include <fstream>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <THStack.h>
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Point3D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > LV;

struct LepInfo{
        float pt;
        float eta;
        float phi;
};

struct FatJetInfo{
	float sdmass;
        float pt;
        float eta;
        float phi;
        float tau21;
        float tau2;
        float tau1;
        float deepMDW;
        float deepW;
	float deepMDZ;
        float deepZ;
        float deepT;
        float deepMDbb;
};

struct JetInfo{
        float pt;
        float eta;
        float phi;
        bool passbloose;
        bool passbmedium;
        bool overlapwithfatjet;
};

struct Event{
        vector<LepInfo> lep;
        vector<FatJetInfo> fatjet;
        vector<JetInfo> jet;
        float mll;
        float ht;
        float st;
        float deltaphill;
	float njet_overlapremoved;
	float nbloose;
	float nbmedium;
	float nbloose_overlapremoved;
	float nbmedium_overlapremoved;
	float mtmax;
	float mtmin;
	float ptll;
	int flavor;
	float met;
};

bool addtion_cut(const Event& evt){
	return evt.lep.size()==2;
}

bool make_BDT(const char* input,const char* output){
	Event evt;
	//input
	TChain *chain=new TChain("t");
        vector<LV> *lepp4=0,*jetp4=0,*fatjetp4=0;
	LV *met=0;
	vector<int> *pdgid=0;
	vector<bool> *jetbloose=0,*jetbmedium=0;
	vector<float> *fatjetsdm=0,*tau21=0,*tau2=0,*tau1=0,*deepMDW=0,*deepW=0,*deepMDZ=0,*deepZ=0,*deepT=0,*deepMDbb=0;
        chain->SetBranchAddress("SS2jet_lep_p4",                &lepp4);
	chain->SetBranchAddress("SS2jet_lep_pdgid",		&pdgid);
        chain->SetBranchAddress("SS2jet_jet_p4",                &jetp4);
        chain->SetBranchAddress("SS2jet_jet_passBloose",        &jetbloose);
        chain->SetBranchAddress("SS2jet_jet_passBmedium",       &jetbmedium);
        chain->SetBranchAddress("SS2jet_fatjet_p4",             &fatjetp4);
        chain->SetBranchAddress("SS2jet_fatjet_msoftdrop",      &fatjetsdm);
        chain->SetBranchAddress("SS2jet_fatjet_tau21",          &tau21);
        chain->SetBranchAddress("SS2jet_fatjet_tau2",           &tau2);
        chain->SetBranchAddress("SS2jet_fatjet_tau1",           &tau1);
        chain->SetBranchAddress("SS2jet_fatjet_deepMD_W",       &deepMDW);
        chain->SetBranchAddress("SS2jet_fatjet_deep_W",         &deepW);
        chain->SetBranchAddress("SS2jet_fatjet_deepMD_Z",       &deepMDZ);
        chain->SetBranchAddress("SS2jet_fatjet_deep_Z",         &deepZ);
        chain->SetBranchAddress("SS2jet_fatjet_deep_T",         &deepT);
        chain->SetBranchAddress("SS2jet_fatjet_deepMD_bb",      &deepMDbb);
        chain->SetBranchAddress("Common_met_p4",		&met);
	chain->Add(input);
	//input
	//output
	TFile *f=TFile::Open(output,"recreate");
	TTree *tr=new TTree("t","t");
	int flavor,njet;
	float lep1pt,lep2pt,mtmax,mtmin,ptll,mll,fatjetpt,fatjetsdmass,fatjettau21,ptmiss,ht,st;
	tr->Branch("lep1pt",&lep1pt);
	tr->Branch("lep2pt",&lep2pt);
	tr->Branch("flavor",&flavor);
        tr->Branch("mtmax",&mtmax);
	tr->Branch("mtmin",&mtmin);
	tr->Branch("ptll",&ptll);
	tr->Branch("mll",&mll);
	tr->Branch("fatjetpt",&fatjetpt);
	tr->Branch("fatjetsdmass",&fatjetsdmass);
	tr->Branch("fatjettau21",&fatjettau21);
	tr->Branch("njet",&njet);
	tr->Branch("met",&ptmiss);
	tr->Branch("ht",&ht);
	tr->Branch("st",&st);


	//output
        for(unsigned int i=0;i<chain->GetEntries();i++){
                chain->GetEntry(i);
		//build the event
		evt.lep.clear();
		evt.fatjet.clear();
		evt.jet.clear();
		evt.ht=0;
		evt.st=0;
		evt.nbloose=0;
		evt.nbmedium=0;
		evt.nbloose_overlapremoved=0;
		evt.nbmedium_overlapremoved=0;
		evt.njet_overlapremoved=0;
		evt.met=met->pt();
		evt.mtmax=lepp4->at(0).mt()>lepp4->at(1).mt()?lepp4->at(0).mt():lepp4->at(1).mt();
                evt.mtmin=lepp4->at(0).mt()<lepp4->at(1).mt()?lepp4->at(0).mt():lepp4->at(1).mt();
		evt.ptll=(lepp4->at(0)+lepp4->at(1)).pt();
		switch(pdgid->at(0)*pdgid->at(1)){
			case 121:evt.flavor=1;break;
			case 169:evt.flavor=2;break;
			case 143:evt.flavor=pdgid->at(0)==11?3:4;break;
			default:evt.flavor=0;break;
		}
		
		for(unsigned int ilep=0;ilep<lepp4->size();ilep++){
			LepInfo temp;
			temp.pt=lepp4->at(ilep).pt();
			temp.eta=lepp4->at(ilep).eta();
			temp.phi=lepp4->at(ilep).phi();

			evt.lep.push_back(temp);

			evt.st+=lepp4->at(ilep).pt();
		}

		for(unsigned int ifatjet=0;ifatjet<fatjetp4->size();ifatjet++){
			FatJetInfo temp;
			temp.sdmass=fatjetsdm->at(ifatjet);
			temp.pt=fatjetp4->at(ifatjet).pt();
			temp.eta=fatjetp4->at(ifatjet).eta();
			temp.phi=fatjetp4->at(ifatjet).phi();
			temp.tau21=tau21->at(ifatjet);
			temp.tau2=tau2->at(ifatjet);
			temp.tau1=tau1->at(ifatjet);
			temp.deepMDW=deepMDW->at(ifatjet);
			temp.deepW=deepW->at(ifatjet);
			temp.deepMDZ=deepZ->at(ifatjet);
			temp.deepZ=deepZ->at(ifatjet);
			temp.deepT=deepT->at(ifatjet);
			temp.deepMDbb=deepMDbb->at(ifatjet);

			evt.fatjet.push_back(temp);

			evt.st+=fatjetp4->at(ifatjet).pt();
			evt.ht+=fatjetp4->at(ifatjet).pt();
		}

		for(unsigned int ijet=0;ijet<jetp4->size();ijet++){
			JetInfo temp;
			temp.pt=jetp4->at(ijet).pt();
			temp.eta=jetp4->at(ijet).eta();
			temp.phi=jetp4->at(ijet).phi();
			temp.passbloose=jetbloose->at(ijet);
			temp.passbmedium=jetbmedium->at(ijet);
			bool overlap=false;
			for(unsigned int ifatjet=0;ifatjet<fatjetp4->size();ifatjet++){
				if(ROOT::Math::VectorUtil::DeltaR(jetp4->at(ijet),fatjetp4->at(ifatjet))<0.8){
					overlap=true;
					break;
				}				
			}
			temp.overlapwithfatjet=overlap;
			evt.jet.push_back(temp);
			if(!overlap){
				evt.st+=jetp4->at(ijet).pt();
				evt.ht+=jetp4->at(ijet).pt();
				evt.njet_overlapremoved++;
			}
			if(jetbloose->at(ijet)){
				evt.nbloose=0;
				if(!overlap) evt.nbloose_overlapremoved++;
			}
			if(jetbmedium->at(ijet)){
				evt.nbmedium++;
				if(!overlap) evt.nbmedium_overlapremoved++;
			}
		}
		evt.st+=met->pt();
		evt.mll=(lepp4->at(0)+lepp4->at(1)).M();
		evt.deltaphill=fabs(ROOT::Math::VectorUtil::DeltaPhi(lepp4->at(0),lepp4->at(1)));
		//build the event
		//fill in output treee
		if(addtion_cut(evt)){//addition cut
			flavor=evt.flavor;
			lep1pt=evt.lep.at(0).pt;
			lep2pt=evt.lep.at(1).pt;
			mtmax=evt.mtmax;
			mtmin=evt.mtmin;
			ptll=evt.ptll;
			mll=evt.mll;
			fatjetpt=evt.fatjet.at(0).pt;
			fatjetsdmass=evt.fatjet.at(0).pt;
			fatjettau21=evt.fatjet.at(0).tau21;
			njet=evt.jet.size();
			ptmiss=evt.met;
			ht=evt.ht;
			st=evt.st;

			tr->Fill();
		}
		//fill in output tree
		
        }
	tr->Write();
	f->Close();
        return true;
}


