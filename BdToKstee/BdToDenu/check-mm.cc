#include <iostream>
#include <assert.h>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGaxis.h"

void calc_angles(TLorentzVector p_K,TLorentzVector p_Pi,TLorentzVector p_Muplus,TLorentzVector p_Muminus,double& thetal, double& thetak, double& phi){
	TLorentzVector p_MuMu,p_KPi,p_B;
    TVector3 B_frame,KPi_frame,MuMu_frame;
    double cosphi,sinphi;

    p_KPi = p_K+p_Pi;
    p_MuMu = p_Muplus+p_Muminus;
    p_B = p_KPi+p_MuMu;

    B_frame = -p_B.BoostVector();

    p_K.Boost(B_frame);
    p_Pi.Boost(B_frame);
    p_Muplus.Boost(B_frame);
    p_Muminus.Boost(B_frame);
    p_KPi.Boost(B_frame);
    p_MuMu.Boost(B_frame);

    cosphi = ((p_Muplus.Vect().Unit()).Cross((p_Muminus.Vect().Unit()))).Dot((p_K.Vect().Unit()).Cross((p_Pi.Vect().Unit())));
    sinphi = ((p_Muplus.Vect().Unit()).Cross((p_Muminus.Vect().Unit()))).Cross((p_K.Vect().Unit()).Cross((p_Pi.Vect().Unit()))).Dot(p_KPi.Vect().Unit());
    phi = TMath::ATan2(sinphi,cosphi);

    MuMu_frame = -p_MuMu.BoostVector();
    KPi_frame = -p_KPi.BoostVector();

    p_Muplus.Boost(MuMu_frame);
    p_K.Boost(KPi_frame);

    thetal = acos((p_Muplus.Vect().Unit()).Dot(p_MuMu.Vect().Unit()));
    thetak = acos((p_K.Vect().Unit()).Dot(p_KPi.Vect().Unit()));

}

int main(int argc, char* argv[]){

	gROOT->ProcessLine(".x lhcbstyle.C");
	gStyle->SetPalette(kBird);
    gStyle->SetPadTopMargin(0.06);
    TGaxis::SetMaxDigits(4); 

    TFile* f = new TFile(argv[1]);
	TTree* t = (TTree*)f->Get("MCDecayTreeTuple/MCDecayTree");

	double K_TRUEP_X,K_TRUEP_Y,K_TRUEP_Z,K_TRUEP_E;
    double Pi_TRUEP_X,Pi_TRUEP_Y,Pi_TRUEP_Z,Pi_TRUEP_E;
    double muplus_TRUEP_X,muplus_TRUEP_Y,muplus_TRUEP_Z,muplus_TRUEP_E;
    double muminus_TRUEP_X,muminus_TRUEP_Y,muminus_TRUEP_Z,muminus_TRUEP_E;

    TLorentzVector true_p_K,true_p_Pi,true_p_muplus,true_p_muminus;
    double true_thetal,true_thetak,true_phi;

	t->SetBranchAddress("Kplus_TRUEP_X",&K_TRUEP_X);
    t->SetBranchAddress("Kplus_TRUEP_Y",&K_TRUEP_Y);
    t->SetBranchAddress("Kplus_TRUEP_Z",&K_TRUEP_Z);
    t->SetBranchAddress("Kplus_TRUEP_E",&K_TRUEP_E);

    t->SetBranchAddress("piminus_TRUEP_X",&Pi_TRUEP_X);
    t->SetBranchAddress("piminus_TRUEP_Y",&Pi_TRUEP_Y);
    t->SetBranchAddress("piminus_TRUEP_Z",&Pi_TRUEP_Z);
    t->SetBranchAddress("piminus_TRUEP_E",&Pi_TRUEP_E);

    t->SetBranchAddress("muplus_TRUEP_X",&muplus_TRUEP_X);
    t->SetBranchAddress("muplus_TRUEP_Y",&muplus_TRUEP_Y);  
    t->SetBranchAddress("muplus_TRUEP_Z",&muplus_TRUEP_Z);
    t->SetBranchAddress("muplus_TRUEP_E",&muplus_TRUEP_E);

    t->SetBranchAddress("muminus_TRUEP_X",&muminus_TRUEP_X);
    t->SetBranchAddress("muminus_TRUEP_Y",&muminus_TRUEP_Y);
    t->SetBranchAddress("muminus_TRUEP_Z",&muminus_TRUEP_Z);
    t->SetBranchAddress("muminus_TRUEP_E",&muminus_TRUEP_E);

    int n_bins = 40;

    TH1F* h_thetal = new TH1F("h_thetal","h_thetal",n_bins,-1,1);
    TH1F* h_thetak = new TH1F("h_thetak","h_thetak",n_bins,-1,1);
    TH1F* h_phi = new TH1F("h_phi","h_phi",n_bins,-TMath::Pi(),TMath::Pi());
    TH1F* h_Kst_M = new TH1F("h_Kst_M","h_Kst_M",n_bins,600,1500);
    TH1F* h_B0_M = new TH1F("h_B0_M","h_B0_M",n_bins,0,6000);

    h_thetal->SetXTitle("cos#it{#theta}_{#it{l}}");
    h_thetak->SetXTitle("cos#it{#theta}_{#it{K}}");
    h_phi->SetXTitle("#it{#phi}");
    h_B0_M->SetXTitle("#it{m}(#it{K^{+}#pi^{#font[122]{-}}e^{+}e^{#font[122]{-}}})");
    h_Kst_M->SetXTitle("#it{m}(#it{K^{+}#pi^{#font[122]{-}}})");

    for(int ev=0;ev<t->GetEntries();ev++){
    	t->GetEntry(ev);

    	true_p_K.SetPxPyPzE(K_TRUEP_X,K_TRUEP_Y,K_TRUEP_Z,K_TRUEP_E);
        true_p_Pi.SetPxPyPzE(Pi_TRUEP_X,Pi_TRUEP_Y,Pi_TRUEP_Z,Pi_TRUEP_E);
        true_p_muplus.SetPxPyPzE(muplus_TRUEP_X,muplus_TRUEP_Y,muplus_TRUEP_Z,muplus_TRUEP_E);
        true_p_muminus.SetPxPyPzE(muminus_TRUEP_X,muminus_TRUEP_Y,muminus_TRUEP_Z,muminus_TRUEP_E);

        calc_angles(true_p_K,true_p_Pi,true_p_muplus,true_p_muminus,true_thetal,true_thetak,true_phi);
        h_thetal->Fill(cos(true_thetal));
        h_thetak->Fill(cos(true_thetak));
		h_phi->Fill(true_phi); 

        h_B0_M->Fill((true_p_K+true_p_Pi+true_p_muplus+true_p_muminus).M());
        h_Kst_M->Fill((true_p_K+true_p_Pi).M());
    }

    TApplication *theApp = new TApplication("app", &argc, argv);
    TQObject::Connect("TCanvas", "Closed()", "TApplication", gApplication, "Terminate()");

    TCanvas* c_thetal = new TCanvas("c_thetal","c_thetal",600,600);
    TCanvas* c_thetak = new TCanvas("c_thetak","c_thetak",600,600);
    TCanvas* c_phi = new TCanvas("c_phi","c_phi",600,600);
    TCanvas* c_B0_M = new TCanvas("c_B0_M","c_B0_M",600,600);
    TCanvas* c_Kst_M = new TCanvas("c_Kst_M","c_Kst_M",600,600);

    c_thetal->cd();
    h_thetal->SetMinimum(0);
    h_thetal->Draw();
    c_thetal->Update();

    c_thetak->cd();
    h_thetak->SetMinimum(0);
    h_thetak->Draw();
    c_thetak->Update();

    c_phi->cd();
    h_phi->SetMinimum(0);
    h_phi->Draw();
    c_phi->Update();

    c_B0_M->cd();
    h_B0_M->SetMinimum(0);
    h_B0_M->Draw();
    c_B0_M->Update();

    c_Kst_M->cd();
    h_Kst_M->SetMinimum(0);
    h_Kst_M->Draw();
    c_Kst_M->Update();

    theApp->Run();
}