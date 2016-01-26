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
    cout << "Template " << sProp << ": " << tag << endl;
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


double RPoissonFitHandler::fitPoint(int index) {
    cout << "Fit with template " << sProp 
         << ": " << mcSigTemplVal[index] << endl;

    fittedTemplate = index;

    // some cleanup
    if (fitResults)  {
        delete fitResults;

        if (bkgsyst) {
            delete pdffit;
            delete extBackgroundPdf;
        }
    }

    // name of the model we're fitting
    char tName[50];
    sprintf(tName,"model%.2f", mcSigTemplVal[index]);

    //
    pdffit = w->pdf(tName) ;
    pdffit->Print();

    //
    if (useRatio && fixBckg == 3) {
        char tag[50];

        for (int itype = 0; itype < processes.size(); ++itype) {
            //
            sprintf(tag,"%s%.2f",processes.at(itype),mcSigTemplVal[index]);
            sprintf(tName,"ratio%s%.2f", processes.at(itype), mcSigTemplVal[index]);

            // S/(S+B) significance metric
            w->var(tName)->setVal(mcSigTemplHistosScaled[tag]->Integral()/
                                  (mcSigTemplHistosScaled[tag]->Integral()
                                   +mcTotalBkgHistosScaled[processes.at(itype)]->Integral()));
            w->var(tName)->setConstant(1);
        }
    } else if (useRatio && fixBkg == 1) {

        for (int itype = 0; itype < processes.size(); ++itype) {
            sprintf(tName,"ratio%s%.2f", processes.at(itype), mcSigTemplVal[index]);
            w->var(tName)->setVal(1.0);
            w->var(tName)->setConstant(1);
        }
    }

    // fit the pdf to the analysis data
    fitResults = pdffit->fitTo(*data, Save(), PrintLevel(1)) ;
    
    // get the chi-squared likelihood from our fit results
    double fitChi2 = fitResults->minNll();
    cout << " - LL is: " << fitChi2 << endl;
    return fitChi2;
}


int RPoissonFitHandler::fitAll() {
    // cleanup and prepare for the next round of fits
    chiSquared.clear();
    chiSquared.resize(mcSigTemplVal.size()+1);

    // fit a model for each propVal
    for (int i = 0; i <= mcSigTemplVal.size(); ++i) {
        chiSquared.at(i) = fitPoint(i);
        if (std::isinf(chiSquared.at(i)) || std::isnan(chiSquared.at(i))) return 1;
    }

    return 0;
}


pair<double,double> RPoissonFitHandler::minimize(bool fake,
                                                 bool all = true, 
                                                 int  pointsToUse = 2) {
    int    pts   = 0;
    double minLL = std::numeric_limits<double>::infinity();

    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared.at(i) >= 0) {
            ++pts;
            chiSquared.assign(i, min(chiSquared.at(i), minLL));
        }
    }

    TVectorD x(pts), ex(pts),
             y(pts), ey(pts);
    pts = 0;

    int minPt    = 0;
    int failures = 0;

    for (int i = 0; i <= mcSigTemplVal.size(); ++i) {
        if (chiSquared.at(i) >= 0){
            x[pts]  = mcSigTemplVal.at(i);
            y[pts]  = chiSquared.at(i) - minLL;
            ey[pts] = 9.4;
            
            //cout << x[pts] << " " << y[pts] << " " << y[minPt] << endl;
            
            if (y[pts] <= y[minPt]) minPt = pts;
            ++pts;
        } else ++failures;
    }

    int points, min, max;

    if (all) { 
        points = pts; 
        min    = 0;
        max    = pts - 1;
    } else {
        if (pointsToUse < 2)   pointsToUse = 2;
        if (pointsToUse > pts) pointsToUse = 4;
        points = 2 * pointsToUse + 1;
        
        if (minPt < pointsToUse) minPt = pointsToUse;
        if (minPt > maxTemplate - pointsToUse - 1) minPt = maxTemplate
                                                            - pointsToUse - 1;
        min = minPt - pointsToUse;
        max = minPt + pointsToUse;
    }

    // some cleanup
    if(gr != 0) delete gr;

    // prepare a graph for the LL plot
    gr = new TGraphErrors(x, y, ex, ey);
    gr->SetName("gr");
    gr->Fit("pol2", "Q", "", mcSigTemplVal.at(min), mcSigTemplVal.at(max));
    gr->GetXaxis()->SetTitle("mlbwa  m(lb) [GeV]");
    gr->GetYaxis()->SetTitle("-log (L/L_{max})");
    gr->GetYaxis()->SetTitleOffset(1.25);

    // get useful params from fits
    Double_t a = gr->GetFunction("pol2")->GetParameter(2);
    if (!fake)
        fitVal = gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal);

    // announce the minimum-likelihood propVal
    cout << " - minimum " << sProp << ": " 
         << gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal) 
         << " +/- " << 1/sqrt(2*a) << endl;

    if (failures) cout << " - fits failed: " << failures << endl;

    for (int i = 0; i <= mcSigTemplVal.size(); ++i) {
        if (chiSquared.at(i) >= 0){
            //
            toyLL->Fill(i,gr->GetFunction("pol2")->Eval(mcSigTemplVal.at(i))
                                                         -(chiSquared.at(i) - minLL));

            //
            cout << i << " " << chiSquared.at(i) << " " 
                 << (chiSquared.at(i) - minLL) << " " 
                 << gr->GetFunction("pol2")->Eval(mcSigTemplVal.at(i))
                                                  -(chiSquared.at(i) - minLL)
                 << endl;
        }
    }

    if (a<0.) return pair<double,double>(0.,0.);

    return pair<double,double>(gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal),
                               1/sqrt(2*a));
}

void RPoissonFitHandler::runCalibration(int numberOfExps = 1000) {

}

void RPoissonAnalysis::run() {

}

void RPoissonAnalysis::save() {

}