#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

//make it like func(filename, hist name(s))
void AcceptanceEfficiency(){
    //add asserts
    TFile f("../build/output.root", "READ"); // make a pointer
    TFile *f1 = new TFile("../build/AccEff.root", "RECREATE");
    
    TH2F *PDG_PRIM_RECO;
    TH2F *PDG_PRIM;

    f.GetObject("PDG_prim_RECO",PDG_PRIM_RECO);
    f.GetObject("PDG_prim",PDG_PRIM);

    TH2F *AccEff = new TH2F("AccEffEff","Acceptance(x)Efficiency;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    TH2F *Cor_AccEff = new TH2F("Cor_AccEff","Acceptance(x)Efficiency, errors < 0.3;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    
    PDG_PRIM_RECO->Sumw2();
    PDG_PRIM->Sumw2();
//  main algo
    AccEff->Divide(PDG_PRIM_RECO,PDG_PRIM);
    
    Double_t entries_zero = 0.;
//  corrected  
    for(int i = 0; i < AccEff->GetEntries(); i++){
        Double_t AccEff_bin = AccEff->GetBinContent(i);
        if((AccEff->GetBinError(i) > 0.3*AccEff_bin) || (AccEff_bin + AccEff->GetBinError(i)) > 1)
        {
            Cor_AccEff->SetBinContent(i,0.0);
            entries_zero = entries_zero + 1.0;
        }
        else {
            Cor_AccEff->SetBinContent(i, AccEff_bin);
            Cor_AccEff->SetBinError(i, AccEff->GetBinError(i));
        }
    }
    Double_t entries_nonzero = Cor_AccEff->GetEntries();
    std::cout << AccEff->GetEntries() << "\n";
    Cor_AccEff->SetEntries(entries_nonzero - entries_zero);
    std::cout << Cor_AccEff->GetEntries() << "\n";
    Cor_AccEff->Sumw2();

    PDG_PRIM_RECO->Write();
    PDG_PRIM->Write();
    Cor_AccEff->Write();
    AccEff->Write();
}


