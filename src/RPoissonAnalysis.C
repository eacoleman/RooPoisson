/*******************************************************************************
* RPoissonAnalysis.C
* 
*
*******************************************************************************/

// RooPoisson Inclusions
#include "include/RPoissonAnalysis.h"



RPoissonAnalysis::RPoissonAnalysis() {

}

RPoissonAnalysis::~RPoissonAnalysis() {

}


TH1F* RPoissonAnalysis::getTemplHisto(string process, int iProp) {
    //
    char tag[50];
    sprintf(tag,"%s%.2f", process.data(), mcSigTemplVal.at(i));

    //
    cout << "Template mass "<< tag << endl;
    return mcSigTemplHistosScaled[tag];
}


TH1F* RPoissonAnalysis::getBkgHisto(string process) {
    return mcTotalBkgHistosScaled[proess.data()];
}


void RPoissonAnalysis::typeTag(char* nameToTag) {
    char tName[200];

    cout << " - Name to modify is: " << nameToTag << endl;

    for (int j = 0; j < processes.size(); ++j) {
        sprintf(tName,"%s%s", nameToTag, processes.at(j));

        if (j<(processes.size()-1)) sprintf(nameToTag,"%s_", tName);
        else sprintf(nameToTag,"%s", tName);
    }

    cout << " - New name is: " nameToTag << endl;
}


void RPoissonAnalysis::doToys(int nExp, int iTemplate) {
    
    // alert the user what is happening
    if (!systematics || !systematicsPDF) 
         cout << " - template " << sProp << " value "
              << mcSignalTemplMass[iTemplate] << endl;
    else cout << " - pdf " << iTemplate << endl;

    // cleanup from any previous runs of doToys
    if (toy_error != 0) {
        delete toy_mean;
        delete toy_error;
        delete toy_pull;
        delete toy_bias;
    }

    if (toy_LL != 0) delete toy_LL;

    //
    massPoint = mcSignalTemplMass[iTemplate];

    // initialize new toy histograms 
    toyMean   = new TH1F("mean"  ,sProp,100, minPropVal, maxPropVal);
    toyBias   = new TH1F("bias"  ,(sProp+string("bias")),100, -3.5, 3.5);
    toyPull   = new TH1F("pull"  ,"pull",200, -10, 10);
    toyError  = new TH1F("error" ,(sProp+string("uncertainty")),500, 0, 0.4);
    toyLL     = new TH2F("LL"  ,"LL residuals",9, -0.5, 8.5,200,-100,100);

    // histogram styling
    toyMean->GetXaxis()->SetNdivisions(50205);
    toyBias->GetXaxis()->SetNdivisions(50205);

    toyMean->SetFillColor(44);
    toyBias->SetFillColor(44);
    toyPull->SetFillColor(44);
    toyError->SetFillColor(44);

    // begin performing experiments
    int j, stat;
    for (int i=0; i < nExp; i++) {
        int countFailures = 0;
        do {
            // get the result of the experiment
            do { j = generate_toy(iTemplate);} while (j==0);

            stat = fitAll();
            if (stat !=0) {
                ++countFailures;
                cout << " - RooFit failure " << countFailures << endl;
            }

        } while ((stat !=0) && countFailures<10);

        if (stat !=0) {
            cout << "ERROR: TOO MANY CONSECUTIVE FAILURES. \n";
            exit(1);
        }

        nFitFailed += countFailures;

        pair<double,double> result = fitHandler->minFake();
        if (result.second >= 0.) {
            cout << " - fake fit result: " << result.first 
                 << " +/- " << result.second << endl;

            // fill the toy histograms with the result of our pseudoexperiment
            toyMean->Fill(result.first);
            toyBias->Fill(result.first-massPoint);
            toyPull->Fill((massPoint-result.first)/result.second);
            toyError->Fill(result.second);
        } else {
            ++nFitFailed;
        }

        ++nFitTried;
    }

    cout << " - failed fits: " << nFitFailed << " / " << nFitTried << endl;
}


void RPoissonAnalysis::setOverflowBins(TH1F* histo) {
    //
    Int_t _last_bin = histo->GetNbinsX();
    Double_t _overflow = histo->GetBinContent(_last_bin+1);
    Int_t n_entries = (Int_t) histo->GetEntries();

    //
    h->AddBinContent(1,histo->GetBinContent(0));
    h->AddBinContent(_last_bin,_overflow);

    //
    h->SetBinContent(0,0);
    h->SetBinContent(_last_bin+1,0);
    h->SetEntries(n_entries);
}


void RPoissonAnalysis::assembleDatasets() {
    char tName[150];

    // create a dictionary of TH1s to operate on
    cout << " - combining "<< processes.size() <<" categories\n";
    map<string,TH1*> mapToImport;

    // iterate through all the types
    for (int itype = 0; itype < processes.size(); ++itype) {
        sprintf(tName, "%s", processes.at(itype));

        cout << " - process type " << tName
             << " has mean value " << datasets[tName]->GetMean()) << endl;

        mapToImport[tName] = datasets[tName];
    }

    // create dataset histogram and print it out
    data = new RooDataHist("data", "combined data",
                           *propVal, Index(*sample), 
                           Import(mapToImport));
    data->Print();
}


void RPoissonAnalysis::run() {

}

void RPoissonAnalysis::save() {

}