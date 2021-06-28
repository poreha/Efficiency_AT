#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>


void Efficiency(){
    TFile f("../build/output.root", "READ");
    TFile *f1 = new TFile("../build/EffAcc.root", "RECREATE");
    
    TH2F *PID_PRIM;
    TH2F *PDG_PRIM;

    f.GetObject("PID_prim",PID_PRIM);
    f.GetObject("PDG_prim",PDG_PRIM);

    TH2F *EffAcc = new TH2F("EffAcc","Acceptance(x)Efficiency (pid+reco);Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    TH2F *Cor_EffAcc = new TH2F("Cor_EffAcc","Acceptance(x)Efficiency (pid+reco), errors < 0.3;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    
    PID_PRIM->Sumw2();
    PDG_PRIM->Sumw2();
//  main algo
    EffAcc->Divide(PID_PRIM,PDG_PRIM);
    
    Double_t entries_zero = 0.;
//  corrected  
    for(int i = 0; i < EffAcc->GetEntries(); i++){
        Double_t EffAcc_bin = EffAcc->GetBinContent(i);
        if((EffAcc->GetBinError(i) > 0.3*EffAcc_bin) || (EffAcc_bin + EffAcc->GetBinError(i)) > 1)
        {
            Cor_EffAcc->SetBinContent(i,0.0);
            entries_zero = entries_zero + 1.0;
        }
        else {
            Cor_EffAcc->SetBinContent(i, EffAcc_bin);
            Cor_EffAcc->SetBinError(i, EffAcc->GetBinError(i));
        }
    }
    Double_t entries_nonzero = Cor_EffAcc->GetEntries();
    std::cout << EffAcc->GetEntries() << "\n";
    Cor_EffAcc->SetEntries(entries_nonzero - entries_zero);
    std::cout << Cor_EffAcc->GetEntries() << "\n";
    Cor_EffAcc->Sumw2();
    Cor_EffAcc->SetMinimum(0.2);

    PID_PRIM->Write();
    PDG_PRIM->Write();
    Cor_EffAcc->Write();
    EffAcc->Write();
}


