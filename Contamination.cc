#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>


void Contamination(){
    TFile mism("../build/mismatch.root", "READ");
    TFile sec_contam("../build/sec_contam.root", "READ");
    TFile *f1 = new TFile("../build/contam.root", "RECREATE");
    //TCanvas *c1 = new TCanvas("c1", "Canvas", 1200, 1200);
    TH2F *MM_PID;
    TH2F *SEC_CONT;

    mism.GetObject("mismatch_pid",MM_PID);
    sec_contam.GetObject("Sec_Cont",SEC_CONT);

    TH2F *CONT = new TH2F("CONT","Contamination;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.);
    TH2F *COR_CONT = new TH2F("COR__CONT","Contamination w/out big errors;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.4);
    
    CONT->Add(MM_PID, SEC_CONT);
    COR_CONT->Clone("Sec_Cont");
    std::cout << "\n" << CONT->GetEntries() << "\n";
    for(int i = 0; i < CONT->GetEntries(); i++)
    {
        float_t cont_bin = 100*CONT->GetBinContent(i);
        if(CONT->GetBinError(i) > 1000)
        {
            COR_CONT->SetBinContent(i,0);
        }
        else {
            COR_CONT->SetBinContent(i, cont_bin);
            COR_CONT->SetBinError(i, SEC_CONT->GetBinError(i));
        }
    }   
    CONT->Write();
    COR_CONT->Write();
}

