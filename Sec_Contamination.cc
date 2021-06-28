#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>


void Sec_Contamination(){
    TFile f("../build/output.root", "READ");
    TFile *f1 = new TFile("../build/sec_contam.root", "RECREATE");
    
    TH2F *PID_SEC;
    TH2F *MISMATCH;
    TH2F *PID_RECO;
    f.GetObject("PID_sec",PID_SEC);
    f.GetObject("PID_reco",PID_RECO);

    TH2F *SEC_CONT = new TH2F("Sec_Cont","Secondary Contamination;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.);
    TH2F *COR_SEC_CONT = new TH2F("COR_SEC_CONT","Secondary contamination, errors < 0.3;Y;p_{T}", 100, 0., 2.5, 100, 0., 2.4);
    PID_RECO->Sumw2();
    PID_SEC->Sumw2();
    
    SEC_CONT->Divide(PID_SEC, PID_RECO);
    //COR_SEC_CONT->Clone("Sec_Cont");
    //std::cout << "\n" << SEC_CONT->GetEntries() << "\n";
    Double_t entries_zero = 0.;
    Double_t entries_nonzero = SEC_CONT->GetEntries();
    for(int i = 0; i < SEC_CONT->GetEntries(); i++)
    {
        Double_t cont_bin = SEC_CONT->GetBinContent(i);
        if((SEC_CONT->GetBinError(i) > 0.5*cont_bin) && cont_bin != 0.0)
        {
            COR_SEC_CONT->SetBinContent(i,0);
            entries_zero = entries_zero + 1.0;
        }
        else {
            COR_SEC_CONT->SetBinContent(i, cont_bin);
            COR_SEC_CONT->SetBinError(i, SEC_CONT->GetBinError(i));
        }
    }   
    COR_SEC_CONT->SetEntries(entries_nonzero - entries_zero);
    std::cout << SEC_CONT->Integral() << "\n";
    PID_SEC->Write();
    PID_RECO->Write();
    SEC_CONT->Write();
    COR_SEC_CONT->Write();
}


