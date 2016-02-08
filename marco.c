{
    TFile *f = new TFile("./samples/2012_combined.root", "UPDATE");
    f->cd();

    TList *list2 = (TList*) f->GetListOfKeys();

    for(int i = 0; i < 134  ; i++) {
        TH1F *h = (TH1F*) list2->At(i);
        TString name = h->GetName();
        name.ReplaceAll(".5_", ".50_");
        name.ReplaceAll(".0_", ".00_");
        TH1F* newH = (TH1F*) h->Clone(name);

        newH->Write();
    }

    f->Close();
}
