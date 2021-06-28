#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>


void Mismatch_pid(){
    TFile f("../build/output.root", "READ");
    TFile *f1 = new TFile("../build/mismatch.root", "RECREATE");
    //TCanvas *c1 = new TCanvas("c1", "Canvas", 1200, 1200);
    TH2F *MISMATCH;
    TH2F *PID_RECO;
    f.GetObject("Mismatch",MISMATCH);
    f.GetObject("PID_reco",PID_RECO);

    TH2F *MM_PID = new TH2F("mismatch_pid","Mismatch;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.);
    TH2F *COR_MM_PID = new TH2F("cor_mismatch_pid","Mismatch w/out big errors;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.4);
    PID_RECO->Sumw2();
    MISMATCH->Sumw2();
    
    MM_PID->Divide(MISMATCH, PID_RECO);
    //COR_MM_PID->Clone("Mismatch");
    Double_t entries_zero = 0.;
    Double_t entries_all = MM_PID->GetEntries();
    std::cout << "\n" << MM_PID->GetEntries() << "\n";
    for(int i = 0; i < MM_PID->GetEntries(); i++)
    {
        Double_t cont_bin = MM_PID->GetBinContent(i);
        if(MM_PID->GetBinError(i) > 0.3*cont_bin)
        {
            COR_MM_PID->SetBinContent(i,0);
            entries_zero = entries_zero + 1.0;
        }
        else {
            COR_MM_PID->SetBinContent(i, cont_bin);
            COR_MM_PID->SetBinError(i, MM_PID->GetBinError(i));
        }
    }   
    COR_MM_PID->SetEntries(entries_all - entries_zero);
    std::cout << MM_PID->GetEntries() << "\n";
    std::cout << COR_MM_PID->GetEntries() << "\n";

    MM_PID->Write();
    COR_MM_PID->Write();
}
