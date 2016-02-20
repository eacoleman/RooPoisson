/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RPoissonAnalysis.C                                                          *
 * Author: Evan Coleman, 2016                                                  *
 *                                                                             *
 * Purpose: Core class for Poisson likelihood analyses in ROOT.                *
 * Usage  : RPoissonAnalysis *r = new RPoissonAnalysis();                      *
 *          r-><Setting to Change> = <Setting>;                                *
 *          r->run();                                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// RooPoisson Inclusions
#include "include/RPoissonAnalysis.h"
#include "TStyleHandler.C"

using namespace RooFit;


RPoissonAnalysis::RPoissonAnalysis() {

}


void RPoissonAnalysis::setup() {
    //
    sort(mcSigTemplVal.begin(), mcSigTemplVal.end());

    //
    TStyleHandler::setTDRStyle();
    c_min = new TCanvas("c_min","", canvWid, canvHei);
    TStyleHandler::initStyle(c_min);

    // initialize some helper variables
    nFitTried      = 0;
    nFitFailed     = 0;
    fittedTempl    = -1;

    nTotSample.resize(processes.size());

    //
    char evHistoName[15]   = "h_nEvents";

    //
    char propValPref[15] = "Reconstructed ";
    propVal = new RooRealVar(sProp.c_str(), strcat(propValPref, sProp.c_str()),
                             minPropVal, maxPropVal);
    sample  = new RooCategory("sample", "sample") ;

    //set types of sample
    for (int itype = 0; itype < processes.size(); ++itype) 
        sample->defineType(processes.at(itype).c_str()); 

    //
    char  name[300], 
         sname[150],
         hname[150], 
            tag[50];

    //
    TFile* theFile;
    TH1::AddDirectory(kFALSE);

    // open data file
    theFile = new TFile (dataFileLoc.c_str());
    for (int ihisto = 0; ihisto < processes.size(); ++ihisto) {
        // get the data histograms
        sprintf(sname, "%s"        , processes.at(ihisto).c_str());
        sprintf(hname, "%s__Data_%s", sAcro.c_str(), processes.at(ihisto).c_str());
        cout << hname << endl;

        datasets[sname] = (TH1F*) gDirectory->Get(hname);

        if (datasets[sname]==0) assert(false);

        cout << " - got dataset " << processes.at(ihisto) << " " 
             << datasets[sname] << endl;
    }
    
    // cleanup
    theFile->Close();
    delete theFile;

    //
    assembleDatasets();

    //
    for (int i = 0; i < processes.size(); ++i) {
        nTotSample.at(i) = datasets[processes.at(i)]->GetEntries();
        cout << " - in "<< processes.at(i) << " " << nTotSample.at(i) << endl;
    }

    //
    toyDataHisto = new TH1F("toyDataHisto", "toyDataHisto", 
                            floor((maxPropVal - minPropVal)/toyDataHistoDiv),
                            minPropVal, maxPropVal);
    toyDataHisto->Reset();

    //
    data->fillHistogram(toyDataHisto, *propVal);
    dataHisto = (TH1F*) toyDataHisto->Clone("dataHisto");

    // open data file
    theFile = new TFile (dataFileLoc.c_str());
    theFile->cd();

    // loop through the interaction types we want to check
    for (int itype = 0; itype < processes.size(); ++itype) {

        // loop through the propVals we want to check
        for (unsigned int iVal = 0; iVal <  mcSigTemplVal.size(); ++iVal) {

            // get the histogram for the interaction and propVal we want
            sprintf(hname, "%s__%s_%.1f_%s", sAcro.c_str(), sigProcess.c_str(), 
                                             mcSigTemplVal.at(iVal), 
                                             processes.at(itype).c_str());
            sprintf(tag  , "%s%.1f"        , processes.at(itype).c_str(), 
                                             mcSigTemplVal.at(iVal));
            cout << " - signal template: " << hname << " " << tag << endl;

            printf("itype %d iVal %d \n", itype, iVal);

            TH1F* histo = (TH1F*) gDirectory->Get(hname);

            if (histo==0) { 
                cout << "Histo does not exist\n";
                histo = new TH1F(hname,hname,100,0,200); 
            }

            histo->Rebin(5);
            histo->SetLineColor(4);

            mcSigTemplHistosScaled[tag] = histo;
        }

    }
    // cleanup
    theFile->Close();
    delete theFile;

    // plot the systematics GEN templates if they are available
    cout << " - retrieved systematics signal GEN template from " 
         << dataFileLoc << endl;
    theFile = new TFile (dataFileLoc.c_str()) ;

    unsigned int maxTemplates;
    maxTemplates = mcSigTemplVal.size();

    // loop through the propVals
    for (unsigned int i = 0; i < maxTemplates; ++i){

        float tVal = mcSigTemplVal.at(i);

        //loop through the interaction types
        for (int itype = 0; itype < processes.size(); ++itype) {
            // get the histogram name
            sprintf(tag,"%s%.1f",processes.at(itype).c_str(),tVal);
            sprintf(hname, "%s__%s_%.1f_%s", sAcro.c_str(), sigProcess.c_str(), 
                                             tVal, processes.at(itype).c_str());
            
            // get the histogram and store in mcSigHistos_gen
            cout << " - signal GEN template "  << processes.at(itype) 
                 << " " <<i<<" "<< tag <<" "<< hname <<endl;

            //
            TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
            histo->Rebin(5);
            histo->SetLineColor(4);
            mcSigTemplHistosScaled_gen[tag] = histo;
        }
    }
    // cleanup
    theFile->Close();
    delete theFile;
    

    // Get the background templates
    cout << " - retrieving background templates\n";

    // open the histogram file
    theFile = new TFile (dataFileLoc.c_str());
    theFile->cd();

    // iterate through types
    for (int itype = 0; itype < processes.size(); ++itype) {

        //iterate through background MC types
        for (unsigned int bkgType = 0; bkgType < mcBkgLabels.size(); ++bkgType) {
            // format the name of the histo we want to get
            sprintf(hname, "%s__%s_%s", sAcro.c_str(), mcBkgLabels.at(bkgType).Data(), 
                    processes.at(itype).c_str());
            
            // 
            TH1F* histo  = (TH1F*) gDirectory->Get(hname);

            //
            if (histo == 0) { 
                cout << " - WARNING: histo " << hname << " does not exist! \n";
                histo = new TH1F(hname, hname, 100, 0, 200); 
            }

            // rebin the histogram and add to the total backgrounds histogram
            histo->Rebin(5);
            histo->SetLineColor(2);
            mcBkgHistosScaled[tag] = histo;
            if (mcTotalBkgHistosScaled.find(processes.at(itype)) == mcTotalBkgHistosScaled.end()) {
                mcTotalBkgHistosScaled[processes.at(itype)] = (TH1F*) histo->Clone("mcTotalBkgHistosScaled");
            } else {
                mcTotalBkgHistosScaled[processes.at(itype)]->Add(histo);
            }
        }

        cout << " _ total mean: " << mcTotalBkgHistosScaled[processes.at(itype)]->GetMean() << "\n"
             << "   entries: "    << mcTotalBkgHistosScaled[processes.at(itype)]->Integral()
             << endl;
    }

    // clean up, report to console
    theFile->Close();
    delete theFile;

    cout << " - found "<< mcBkgHistosScaled.size() << "background histos\n";

    //
    cout << " - retrieving systematics background GEN template from "
         << dataFileLoc << endl;
    theFile = new TFile (dataFileLoc.c_str());

    // iterate through the background MC types
    for (unsigned int i = 0; i < mcBkgLabels.size(); ++i) {

        // iterate through types
        for (int itype = 0; itype < processes.size(); ++itype) {
            //format the name of the histo we want to get
            sprintf(hname, "%s__%s_%s", sAcro.c_str(), mcBkgLabels.at(i).Data(), 
                    processes.at(itype).c_str());
            cout << hname << endl;

            TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
            
            if (histo == 0) { 
                cout << "Histo does not exist\n";
                histo = new TH1F(hname, hname, 100, 0, 200); 
            }

            //rebin the histogram and add to the total backgrounds_GEN histogram
            histo->Rebin(5);

            //
            if (mcTotalBkgHistosScaled_gen.find(processes.at(itype)) 
                  == mcTotalBkgHistosScaled_gen.end())               {
                mcTotalBkgHistosScaled_gen[processes.at(itype)] 
                   = (TH1F*) histo->Clone("mcTotalBkgHistosScaled_gen");
            } else {
                mcTotalBkgHistosScaled_gen[processes.at(itype)]->Add(histo);
            }
        }
    }

    // clean up
    theFile->Close();
    delete theFile;

    // re-initialize some variables, establish the workspace
    gr         = 0; 
    grc        = 0;
    toyMean    = 0; 
    fitResults = 0;

    // establish the RooWorkspace
    workspace = new RooWorkspace("w", "workspace") ;
    workspace->import(*propVal);
    workspace->import(*sample);

    // define the datasets initially
    // loop through the interaction types
    for (int itype = 0; itype < processes.size(); ++itype) {
        sprintf(hname, "histo_bck%s", processes.at(itype).c_str());
        cout << hname << endl;
        
        workspace->import(*(new RooDataHist(hname, hname, *propVal, 
                                            mcTotalBkgHistosScaled[processes.at(itype)])));
        
        sprintf(hname,"HistPdf::background%s(%s,histo_bck%s)", 
                processes.at(itype).c_str(), sProp.c_str(), processes.at(itype).c_str());
        cout << hname << endl;
        
        workspace->factory(hname);

        sprintf(hname,"histo_bck_gen%s",processes.at(itype).c_str());
        workspace->import(*(new RooDataHist(hname, hname, *propVal, 
                                            mcTotalBkgHistosScaled_gen[processes.at(itype)])));

        sprintf(hname,"HistPdf::background_gen%s(%s,histo_bck%s)", 
                processes.at(itype).c_str(), sProp.c_str(),  processes.at(itype).c_str());
        workspace->factory(hname);


        //
        float totalBkg = (float) mcTotalBkgHistosScaled[processes.at(itype)]->Integral();

        //
        sprintf(hname,"nbrBkgEvts%s", processes.at(itype).c_str());
        cout << hname << endl;

        //
        nbrBkgEvts[processes.at(itype)] 
            = new RooRealVar(hname, hname, totalBkg, minNumBkg, maxNumBkg);

        //
        workspace->import(*nbrBkgEvts[processes.at(itype)]);

        // loop through the propVals
        for (unsigned int i = 0; i < mcSigTemplVal.size(); ++i) {
            float tVal = mcSigTemplVal.at(i);
            cout << " - building "  << processes.at(itype) << " pdf for "
                 << sProp << " " << tVal << endl;

            sprintf(tag, "%s%.1f", processes.at(itype).c_str(), tVal);

            //
            sprintf(tag  , "%s%.1f"         , processes.at(itype).c_str(), tVal);
            sprintf(hname, "histo_sgn%s%.1f", processes.at(itype).c_str(), tVal);
            cout << hname << " " << tag << endl;
            
            //
            workspace->import( *(new RooDataHist(hname, hname, *propVal, mcSigTemplHistosScaled[tag]) ));
            sprintf(hname,"HistPdf::signal%s%.1f(%s,histo_sgn%s%.1f)",
                    processes.at(itype).c_str(),tVal,sProp.c_str(),processes.at(itype).c_str(),tVal);
            cout << hname << endl;
            
            //
            workspace->factory(hname);

            //
            sprintf(hname, "histo_sgn_gen%s%.1f", processes.at(itype).c_str(), tVal);
            workspace->import( *(new RooDataHist(hname, hname, *propVal, 
                                                 mcSigTemplHistosScaled_gen[tag])));
            //
            sprintf(hname,"HistPdf::signal_gen%s%.1f(%s,histo_sgn_gen%s%.1f)",
                    processes.at(itype).c_str(), tVal, sProp.c_str(), processes.at(itype).c_str(), tVal);
            workspace->factory(hname);

            //
            sprintf(hname,"SUM::model%s%.1f( ratio%s%.1f[0,1]*signal%s%.1f, background%s )",
                    processes.at(itype).c_str(), tVal, processes.at(itype).c_str(), tVal, 
                    processes.at(itype).c_str(), tVal, processes.at(itype).c_str());

            cout << hname << endl;

            workspace->factory(hname);
        }
    }

    // this is only for the pdf systematics templates:
    char name2[300];
    for (unsigned int i = 0; i < mcSigTemplVal.size(); ++i) {
        float tVal = mcSigTemplVal.at(i);
        cout << " - building simultaneous pdf for " << sProp << " " 
             << tVal << endl;

        sprintf(name2, "SIMUL::model%.1f(sample", tVal);
        for (int j = 0; j < processes.size(); ++j) {
            sprintf(name , "%s,%s=model%s%.1f", name2, processes.at(j).c_str(),
                    processes.at(j).c_str(), tVal);
            sprintf(name2, "%s", name);
        }

        sprintf(name,"%s)",name2);

        workspace->factory(name);
    }

    toyError = 0;
    toyLL    = new TH2F("LL", "LL residuals", toyLLResXBin, toyLLResXMin, toyLLResXMax,                     
                                              toyLLResYBin, toyLLResYMin, toyLLResYMax);
}


TH1F* RPoissonAnalysis::getTemplHisto(string process, int iProp) {
    //
    char tag[50];
    sprintf(tag,"%s%.1f", process.data(), mcSigTemplVal.at(iProp));

    //
    cout << "Template " << sProp << ": " << tag << endl;
    return mcSigTemplHistosScaled[tag];
}


TH1F* RPoissonAnalysis::getBkgHisto(string process) {
    return mcTotalBkgHistosScaled[process.data()];
}


int RPoissonAnalysis::generateToy(int templToUse) {

    // some cleanup
    if (toyDataHisto != 0) {
        delete toyDataHisto;
    }

    toyDataHisto = (TH1F*) dataHisto->Clone("toyDataHisto");
    toyDataHisto->Reset();

    char hname[50], tag[50];
    int n;
    totGenSig = totGenBkg = 0;

    for (int itype = 0; itype < processes.size(); ++itype) {

        delete datasets[processes.at(itype)];
        
        sprintf(tag, "%s%.1f", processes.at(itype).c_str(), mcSigTemplVal[nomTemplIndex]);

        int binLow  = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(minPropVal);
        int binHigh = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(maxPropVal);

        cout << tag << endl;
        double sigProb;

        sigProb = mcSigTemplHistosScaled_gen[tag]->Integral(binLow, binHigh) /
                    (mcSigTemplHistosScaled_gen[tag]->Integral(binLow, binHigh)
                      + mcTotalBkgHistosScaled_gen[processes.at(itype)]->Integral(binLow, binHigh));

        n = _random.Binomial(nTotSample[itype], sigProb);
        cout << "Generate " << n << " " << processes.at(itype) << " signal events of ";

        cout << sProp  << " " << mcSigTemplVal[templToUse];
        
        cout << ", with binomial prob " << sigProb << endl;
        genSig[processes.at(itype)] = n; 
        totGenSig += n;
        genBkg[processes.at(itype)] = nTotSample[itype] - n;

        if (genBkg[processes.at(itype)] < 0.){
            cout << "ERROR: number of background events < 0.\n";
            assert(false);
        }

        totGenBkg += genBkg[processes.at(itype)];
        cout << " - generated " << genBkg[processes.at(itype)] << " " 
             << " background events\n";
    

        sprintf(hname,"signal_gen%s%.1f", processes.at(itype).c_str(),
                                          mcSigTemplVal.at(templToUse));
        
        cout << hname << endl; 
        datasets[processes.at(itype)] 
            = workspace->pdf(hname)->generateBinned(*propVal, 
                                                    genSig[processes.at(itype)])->createHistogram(hname,*propVal);
        cout << hname << endl;

        sprintf(hname,"background_gen%s", processes.at(itype).c_str());
        
        datasets[processes.at(itype)]->Add(workspace->pdf(hname)->generateBinned(*propVal, 
                                                                                 genBkg[processes.at(itype)])->createHistogram(hname,*propVal));
    }

    if (data!=0) {
        delete data;
    }

    assembleDatasets();

    data->fillHistogram(toyDataHisto, *propVal);

    return data->numEntries();
}


void RPoissonAnalysis::doToys(int nExp, int iTemplate) {
    
    // alert the user what is happening
    cout << " - template " << sProp << " value "
         << mcSigTemplVal.at(iTemplate) << endl;

    // cleanup from any previous runs of doToys
    if (toyError != 0) {
        delete toyMean;
        delete toyError;
        delete toyPull;
        delete toyBias;
    }

    if (toyLL != 0) delete toyLL;

    //
    float propPoint = mcSigTemplVal.at(iTemplate);

    // initialize new toy histograms
    char * sPropPref = new char[sProp.length()];
    strcpy(sPropPref, sProp.c_str());
    
    toyMean   = new TH1F("mean"  , sPropPref, toyMeanBins, 
                                   propPoint - toyMeanDelt, propPoint + toyMeanDelt);
    toyBias   = new TH1F("bias"  , strcat(sPropPref, " bias"), toyBiasBins, 
                                   -abs(toyBiasDelt), abs(toyBiasDelt));
    toyPull   = new TH1F("pull"  , "pull", toyPullBins, -abs(toyPullDelt), 
                                   abs(toyPullDelt));
    toyError  = new TH1F("error" , strcat(sPropPref, " uncertainty"), 
                                   toyErrorBins, toyErrorMin, toyErrorMax);
    toyLL     = new TH2F("LL"    , "LL residuals", toyLLResXBin, toyLLResXMin, 
                                   toyLLResXMax, toyLLResYBin, toyLLResYMin, 
                                   toyLLResYMax);


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
            do { j = generateToy(iTemplate); } while ( j == 0 );

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
        
        pair<double,double> result = minimize(true);
        if (result.second >= 0.) {
            cout << " - fake fit result: " << result.first 
                 << " +/- " << result.second << endl;
            
            // fill the toy histograms with the result of our pseudoexperiment
            toyMean->Fill(result.first);
            toyBias->Fill(result.first-propPoint);
            toyPull->Fill((propPoint-result.first)/result.second);
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
    histo->AddBinContent(1,histo->GetBinContent(0));
    histo->AddBinContent(_last_bin,_overflow);

    //
    histo->SetBinContent(0,0);
    histo->SetBinContent(_last_bin+1,0);
    histo->SetEntries(n_entries);
}


void RPoissonAnalysis::assembleDatasets() {
    char tName[150];

    // create a dictionary of TH1s to operate on
    cout << " - combining "<< processes.size() <<" categories\n";
    map<string,TH1*> mapToImport;

    // iterate through all the types
    for (int itype = 0; itype < processes.size(); ++itype) {
        sprintf(tName, "%s", processes.at(itype).c_str());

        cout << " - process type " << tName
             << " has mean value " << datasets[tName]->GetMean() << endl;

        mapToImport[tName] = datasets[tName];
    }

    // create dataset histogram and print it out
    data = new RooDataHist("data", "combined data",
                           *propVal, Index(*sample), 
                           Import(mapToImport));
    data->Print();
}


double RPoissonAnalysis::fitPoint(int index) {
    cout << "Fit with template " << sProp 
         << ": " << mcSigTemplVal.at(index) << endl;

    fittedTempl = index;

    // some cleanup
    if (fitResults) delete fitResults;

    // name of the model we're fitting
    char tName[50];
    sprintf(tName,"model%.1f", mcSigTemplVal.at(index));

    //
    pdffit = workspace->pdf(tName) ;
    pdffit->Print();

    //
    char tag[50];

    //
    for (int itype = 0; itype < processes.size(); ++itype) {
        //
        sprintf(tag,"%s%.1f",processes.at(itype).c_str(),mcSigTemplVal.at(index));
        sprintf(tName,"ratio%s%.1f", processes.at(itype).c_str(), mcSigTemplVal.at(index));

        // S/(S+B) significance metric
        workspace->var(tName)->setVal(mcSigTemplHistosScaled[tag]->Integral()/
                              (mcSigTemplHistosScaled[tag]->Integral()
                               +mcTotalBkgHistosScaled[processes.at(itype)]->Integral()));
        workspace->var(tName)->setConstant(1);
    }

    // fit the pdf to the analysis data
    fitResults = pdffit->fitTo(*data, Save(), PrintLevel(1)) ;
    
    // get the chi-squared likelihood from our fit results
    double fitChi2 = fitResults->minNll();
    cout << " - LL is: " << fitChi2 << endl;
    return fitChi2;
}


int RPoissonAnalysis::fitAll() {
    // cleanup and prepare for the next round of fits
    chiSquared.clear();
    chiSquared.resize(mcSigTemplVal.size()+1);

    // fit a model for each propVal
    for (int i = 0; i < mcSigTemplVal.size(); ++i) {
        chiSquared.at(i) = fitPoint(i);
        if (std::isinf(chiSquared.at(i)) || std::isnan(chiSquared.at(i))) return 1;
    }

    return 0;
}


pair<double,double> RPoissonAnalysis::minimize(bool fake,
                                               bool all = true, 
                                               int  pointsToUse = 2) {
    //
    int    pts   = 0;
    double minLL = std::numeric_limits<double>::infinity();   
    
    //
    for (int i = 0; i < mcSigTemplVal.size(); ++i) {
        if (chiSquared.at(i) >= 0) {
            ++pts;
            minLL = min(chiSquared[i], minLL);
        }
    }
    
    //
    TVectorD x(pts), ex(pts),
             y(pts), ey(pts);
    pts = 0;

    int minPt    = 0;
    int failures = 0;

    //
    for (int i = 0; i < mcSigTemplVal.size(); ++i) {
        if (chiSquared.at(i) >= 0){
            x[pts]  = mcSigTemplVal.at(i);
            y[pts]  = chiSquared.at(i) - minLL;
            ey[pts] = 9.4;                                                                          //HARDCODED REDALERT
            
            //cout << x[pts] << " " << y[pts] << " " << y[minPt] << endl;
            
            if (y[pts] <= y[minPt]) minPt = pts;
            ++pts;
        } else ++failures;
    }

    //
    int points, min, max;

    //
    if (all) { 
        points = pts; 
        min    = 0;
        max    = pts - 1;
    } else {
        if (pointsToUse < 2)   pointsToUse = 2;
        if (pointsToUse > pts) pointsToUse = 4;
        points = 2 * pointsToUse + 1;
        
        if (minPt < pointsToUse) minPt = pointsToUse;
        if (minPt > mcSigTemplVal.size() - pointsToUse - 1) 
            minPt = mcSigTemplVal.size() - pointsToUse - 1;
        min = minPt - pointsToUse;
        max = minPt + pointsToUse;
    }

    // some cleanup
    if(gr != 0) delete gr;

    // prepare a graph for the LL plot
    gr = new TGraphErrors(x, y, ex, ey);
    gr->SetName("gr");
    gr->Fit("pol2", "Q", "", mcSigTemplVal.at(min), mcSigTemplVal.at(max));
    gr->GetXaxis()->SetTitle(TString(sAcro.data()) + TString(sProp.data()) + 
                             TString("[") + TString(sUnit.data()) + TString("]"));
    gr->GetYaxis()->SetTitle("-log (L/L_{max})");
    gr->GetYaxis()->SetTitleOffset(1.25);

    // get useful params from fits
    Double_t a = gr->GetFunction("pol2")->GetParameter(2);
    if (!fake) {
      fitVal = gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal);
      fitErr = 1/sqrt(2*a);
    }

    // announce the minimum-likelihood propVal
    cout << " - minimum " << sProp << ": " 
         << gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal) 
         << " +/- " << 1/sqrt(2*a) << endl;

    // announce any failures
    if (failures) cout << " - fits failed: " << failures << endl;

    // 
    for (int i = 0; i < mcSigTemplVal.size(); ++i) {
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

    //
    if (a < 0.) return pair<double,double>(0., 0.);

    //
    return pair<double,double>(gr->GetFunction("pol2")->GetMinimumX(minPropVal,maxPropVal),
                               1/sqrt(2*a));
}

void RPoissonAnalysis::getCalibration(int numberOfExps = 1000) {
    // capitalize the first letter of sProp and store it in another string
    string sPropFirstCap    = sProp.substr(0);
           sPropFirstCap[0] = toupper(sPropFirstCap[0]);

    // keep track of timing
    time_t start, endLoop, end;
    time   (&start);

    // keep track of our fit success/failure rates
    nFitFailed = 0;
    nFitTried  = 0;

    // initialize arrays for storing fit likelihoods
    int min = 0; 
    int max = mcSigTemplVal.size();
    cout << " min = " << min << " and max = " << max << endl;

    TVectorD x(max-min+1), ex(max-min+1),
             y(max-min+1), ey(max-min+1);

    // instantiate helper variables for loops
    int pts = 0;
    int i;
    char hname[50];

    // open the outfile
    TFile * out = new TFile(calibFileLoc.c_str(), "RECREATE");

    // loop over templates
    for (i = min; i < max; ++i) {
        // get the toy statistics and fit them with a gaussian
        doToys(numberOfExps, i);
        toyMean->Fit("gaus");

        // store the points for bias plot
        x[pts]  = mcSigTemplVal.at(i); 
        y[pts]  = toyMean->GetFunction("gaus")->GetParameter(1);

        // store the point errors for bias plot
        ex[pts] = 0.;
        ey[pts] = toyMean->GetFunction("gaus")->GetParameter(2)/sqrt(numberOfExps);

        cout << " - y pts is " << y[pts] << endl;

        ++pts;

        // write the toy histograms to our output file
        sprintf(hname,"mean%s_%.1f", sPropFirstCap.data(), mcSigTemplVal.at(i));
        toyMean->Clone(hname)->Write();

        sprintf(hname,"bias%s_%.1f", sPropFirstCap.data(), mcSigTemplVal.at(i));
        toyBias->Clone(hname)->Write();

        sprintf(hname,"err%s_%.1f" , sPropFirstCap.data(), mcSigTemplVal.at(i));
        toyError->Clone(hname)->Write();

        sprintf(hname,"pull%s_%.1f", sPropFirstCap.data(), mcSigTemplVal.at(i));
        toyPull->Clone(hname)->Write();

        sprintf(hname,"LL_%.1f"    , mcSigTemplVal.at(i));
        toyLL->Clone(hname)->Write();
    }

    // loop again over templates
    for (i = min, pts = 0; i < max; ++i, ++pts) {
        cout << " _ template " << sProp << ": " << mcSigTemplVal.at(i) << endl;
        cout << "   fit: " << y[pts] << " / " << ey[pts] << endl;
    }

    // some cleanup
    if (grc != 0) {
        delete grc;
    }

    // initialize the bias plot
    grc = new TGraphErrors(x, y, ex, ey);
    grc->SetName("grc");

    // establish a fit for the bias
    TF1* f1 = new TF1("f1", "pol1", minPropVal, maxPropVal);
    grc->Fit("f1");

    // create a line with slope m=0 and intercept 1 for the bias plot
    TF1* f2 = new TF1("f2", "pol1", minPropVal, maxPropVal);
    f2->SetParameter(0,0.);
    f2->SetParameter(1,1.);
    f2->SetLineColor(4);

    // draw our bias plot with fits
    grc->Draw("a*");
    f2->Draw("same");

    // write everything to our outfile
    grc->Write();
    f1->Write();
    f2->Write();
    out->Close();

    // announce failure rate
    cout << " - failed fits: " << nFitFailed << " / " << nFitTried << endl;

    // announce timing result
    time(&end);
    double dif = difftime(end, start);
    printf(" - it took  %.2lf seconds to do the whole thing, %.2lf per loop.\n", 
            dif , dif/numberOfExps);
}


void RPoissonAnalysis::calibrate() {
    char  tag[20] = "";
    float fitUnc  = fitErr;
    float nomPVal = mcSigTemplVal.at(nomTemplIndex);
    string sPropFirstCap    = sProp.substr(0);
           sPropFirstCap[0] = toupper(sPropFirstCap[0]);

    // reset and stylize the canvas
    c_min = new TCanvas("c_min","", canvWid, canvHei);
    TStyleHandler::initStyle(c_min);

    // initialize helper variables
    TFile* theFile = new TFile(calibFileLoc.c_str());

    char  name[200], 
          hname[50];
    int   minP   = (calibLastPts ? 0 : 1),
          maxP   = mcSigTemplVal.size() - (calibLastPts ? 1 : 2);
    int   points = maxP - minP + 1;

    cout << points << " " << minP << " " << maxP << endl;
    float tMinPVal = mcSigTemplVal.at(minP), 
          tMaxPVal = mcSigTemplVal.at(maxP);

    // initialize arrays for statistics
    TVectorD propV(points),  propErrV(points),
             meanV(points),  meanErrV(points),
             biasV(points),  biasErrV(points),
             pullV(points),  pullErrV(points),
             pullWV(points), pullWErrV(points);

    // loop through our min and max points
    for (int i = 0; i <= maxP - minP; ++i) {
        float propVal = mcSigTemplVal.at(i+minP);

        // collect the calibration histograms for this template
        sprintf(hname,"mean%s_%.1f", sPropFirstCap.data(), propVal);
        TH1F* toyMean  = (TH1F*) gDirectory->Get(hname);

        sprintf(hname,"bias%s_%.1f", sPropFirstCap.data(), propVal);
        TH1F* toyBias  = (TH1F*) gDirectory->Get(hname) ;

        sprintf(hname,"err%s_%.1f" , sPropFirstCap.data(), propVal);
        TH1F* toyError = (TH1F*) gDirectory->Get(hname) ;

        sprintf(hname,"pull%s_%.1f", sPropFirstCap.data(), propVal);
        TH1F* toyPull  = (TH1F*) gDirectory->Get(hname) ;

        if(propVal == 0 || toyPull == 0 || toyMean == 0 || toyBias == 0 || toyError == 0) assert(false);

        // fit out plots with proper gaussians
        toyMean->Fit("gaus","Q");
        toyBias->Fit("gaus","Q");
        toyPull->Rebin(5);
        toyPull->Fit("gaus","Q");

        // save statistics in arrays
        propV(i)     = propVal;
        meanV(i)     = toyMean->GetFunction("gaus")->GetParameter(1);
        biasV(i)     = toyBias->GetFunction("gaus")->GetParameter(1);
        pullV(i)     = toyPull->GetFunction("gaus")->GetParameter(1);
        pullWV(i)    = toyPull->GetFunction("gaus")->GetParameter(2);

        propErrV(i)  = 0.;
        meanErrV(i)  = toyBias->GetFunction("gaus")->GetParameter(2)/sqrt(nPseudoexperiments);
        biasErrV(i)  = toyBias->GetFunction("gaus")->GetParameter(2)/sqrt(nPseudoexperiments);
        pullErrV(i)  = toyPull->GetFunction("gaus")->GetParError(1)/sqrt(nPseudoexperiments);
        pullWErrV(i) = toyPull->GetFunction("gaus")->GetParError(2)/sqrt(nPseudoexperiments);

        // report this template information
        cout << " _ Template " << sProp << ": " << propVal << endl;
        cout << "     " << sProp << ": " << meanV[i] << " +/- " << meanErrV[i] << endl;
        cout << "\t  -  unc: "   << toyError->GetMean();
        cout << "\t  - pull: "   << pullV[i] << " / " << pullWV[i] << endl;
    }

    // propVal plot & fit
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TGraphErrors *propGraph = new TGraphErrors(propV, meanV, propErrV, meanErrV);
    propGraph->SetName("propGraph");
    TF1* meanFit = new TF1("meanFit","pol1", tMinPVal-10, tMaxPVal+10);
    meanFit->SetLineStyle(1);
    propGraph->Fit("meanFit");
    TF1* f2 = new TF1("f2","pol1", tMinPVal-10, tMaxPVal+10);
    f2->SetParameter(0,0.);
    f2->SetParameter(1,1.);
    f2->SetLineColor(4);
    propGraph->GetXaxis()->SetNdivisions(50205);
    propGraph->GetXaxis()->SetTitle(TString("Generated ") + TString(sProp.data()) + 
                                    TString(" [") + TString(sUnit.data()) + TString("]"));
    propGraph->GetYaxis()->SetTitleOffset(1.22);
    propGraph->GetYaxis()->SetTitle(TString("Estimated ") + TString(sProp.data()) + 
                                    TString(" [") + TString(sUnit.data()) + TString("]"));
    propGraph->Draw("apz");
    f2->Draw("same");
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_%s_%s.pdf", sProp.c_str(), tag);
    c_min->Print(hname);
    sprintf(hname,"cal_%s_%s.C"  , sProp.c_str(), tag);
    c_min->Print(hname);
    sprintf(hname,"cal_%s_%s.png", sProp.c_str(), tag);
    c_min->Print(hname);

    // Pull plot & fit
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TGraphErrors *pullGraph = new TGraphErrors(propV, pullV, propErrV, pullErrV);
    pullGraph->SetName("pullGraph");
    TF1* pullFit = new TF1("pullFit","pol0", tMinPVal-10, tMaxPVal+10);
    pullGraph->Fit("pullFit");
    TF1* f3 = new TF1("f3","pol1", tMinPVal-10, tMaxPVal+10);
    f3->SetParameter(0,0.);
    f3->SetParameter(1,0.);
    f3->SetLineColor(4);
    pullGraph->SetMinimum(-1);
    pullGraph->SetMaximum(+1);
    pullGraph->GetXaxis()->SetNdivisions(50205);
    pullGraph->GetXaxis()->SetTitle(TString("Generated ") + TString(sProp.data()) + 
                                    TString(" [") + TString(sUnit.data()) + TString("]"));
    pullGraph->GetYaxis()->SetTitleOffset(1.22);
    pullGraph->GetYaxis()->SetTitle(TString("Pull mean") + 
                                    TString(" [") + TString(sUnit.data()) + TString("]"));
    pullGraph->Draw("apz");
    f3->Draw("same");
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_pull_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pull_%s.C"  , tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pull_%s.png", tag);
    c_min->Print(hname);

    TGraphErrors *pullWGraph = new TGraphErrors(propV, pullWV, propErrV, pullWErrV);
    pullWGraph->SetName("pullWGraph");
    TF1* f4 = new TF1("f3","pol1", tMinPVal-10, tMaxPVal+10);
    f4->SetParameter(0,1.);
    f4->SetParameter(1,0.);
    f4->SetLineColor(4);
    f4->SetLineStyle(2);
    TF1* pullWFit = new TF1("pullWFit","pol0", tMinPVal-10, tMaxPVal+10);
    pullWGraph->GetXaxis()->SetNdivisions(50205);
    pullWGraph->GetXaxis()->SetTitle(TString("Generated ") + TString(sProp.data()) + 
                                     TString(" [") + TString(sUnit.data()) + TString("]"));
    pullWGraph->GetYaxis()->SetTitleOffset(1.22);
    pullWGraph->GetYaxis()->SetTitle("Pull width");
    pullWGraph->SetMinimum(0.75);
    pullWGraph->SetMaximum(1.25);
    pullWGraph->Draw("apz");
    f4->Draw("same");
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_pullW_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pullW_%s.C"  , tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pullW_%s.png", tag);
    c_min->Print(hname);

    // bias plot & fit
    TGraphErrors *biasGraph = new TGraphErrors(propV, biasV, propErrV, biasErrV);
    biasGraph->SetName("biasGraph");
    TF1* biasFit = new TF1("biasFit","pol1",tMinPVal-10, tMaxPVal+10);
    biasGraph->Fit("biasFit");
    f3->SetLineStyle(2);
    biasGraph->GetXaxis()->SetNdivisions(50205);
    biasGraph->GetXaxis()->SetTitle(TString("Generated ") + TString(sProp.data()) + 
                                    TString(" [") + TString(sUnit.data()) + TString("]"));
    biasGraph->GetYaxis()->SetTitleOffset(1.22);
    biasGraph->GetYaxis()->SetTitle(TString("Bias") + TString(" [") + 
                                    TString(sUnit.data()) + TString("]"));
    biasGraph->SetMinimum(-0.5);
    biasGraph->SetMaximum(+0.5);
    biasGraph->Draw("apz");
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_bias_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_bias_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_bias_%s.png", tag);
    c_min->Print(hname);

    // TODO: keep this? Make optional?
    // begin printing nominal information
    cout << "Properties at " << nomPVal << ": " << endl;
    
    // save propVal error at nominal
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(kTRUE);
    sprintf(hname,"err%s_%.1f", sPropFirstCap.data(), nomPVal);
    TH1F* toyErr = (TH1F*) gDirectory->Get(hname);
    cout << "Toy err is "<<toyErr<<endl;
    toyErr->GetXaxis()->SetTitle(TString("Uncertainty") + TString(" [") + 
                                 TString(sUnit.data()) + TString("]"));
    toyErr->GetXaxis()->SetNdivisions(505);
    toyErr->GetYaxis()->SetTitleOffset(1.4);
    toyErr->Draw();
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"err%s_%.1f_%s.pdf", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"err%s_%.1f_%s.C"  , sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"err%s_%.1f_%s.png", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    cout << "Mean uncertainty "<< toyErr->GetMean() << endl;

    // save mean propVal at nominal
    gStyle->SetOptFit(1111);
    sprintf(hname,"mean%s_%.1f", sPropFirstCap.data(), nomPVal);
    TH1F* toyMean = (TH1F*) gDirectory->Get(hname) ;
    toyMean->Fit("gaus","Q");
    toyMean->GetXaxis()->SetTitle(TString(sPropFirstCap.data()) + TString(" [") + 
                                  TString(sUnit.data()) + TString("]"));
    toyMean->GetYaxis()->SetTitleOffset(1.2);
    toyMean->GetYaxis()->SetTitle("Events/bin");
    toyMean->Draw();
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"mean%s_%.1f_%s.pdf", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"mean%s_%.1f_%s.C"  , sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"mean%s_%.1f_%s.png", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);

    // save propVal pull at nominal
    gStyle->SetOptFit(1111);
    sprintf(hname,"pull%s_%.1f", sPropFirstCap.data(), nomPVal);
    TH1F* toyPull  = (TH1F*) gDirectory->Get(hname) ;
    toyPull->Fit("gaus","Q");
    toyPull->GetXaxis()->SetTitle("Pull");
    toyPull->GetYaxis()->SetTitle("Events/bin");
    toyPull->GetYaxis()->SetTitleOffset(1.4);
    toyPull->Draw();
    TStyleHandler::CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"pull%s_%.1f_%s.pdf", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"pull%s_%.1f_%s.C"  , sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);
    sprintf(hname,"pull%s_%.1f_%s.png", sPropFirstCap.data(), nomPVal, tag);
    c_min->Print(hname);

    // calculate final results from fit parameters
    cout << "\n\n - performing inversion...";
    float a = meanFit->GetParameter(1);
    float b = meanFit->GetParameter(0);
    float ae = meanFit->GetParError(1);
    float be = meanFit->GetParError(0);
    float err = sqrt(ae*ae*((fitVal-b)/(a*a))*((fitVal-b)/(a*a)) + be*be/(a*a));
    calVal = (fitVal-b)/a;

    // report results 
    cout << a << " " << b << endl;
    cout << " - template "   << sProp << ": " << fitVal << " - Fit: " << a 
                                              << " / " << b << endl;
    cout << " - calibrated " << sProp << ": " << (fitVal-b)/a << endl;
    cout << " - calibrated error: "    << err << endl;
    cout << " - pull at "              << nomPVal << ": " << pullWFit->Eval(nomPVal) << endl;
    cout << " - calibrated stat unc: " << pullWFit->Eval(nomPVal)*fitUnc << endl;

    theFile->Close();
}


void RPoissonAnalysis::run() {
    // TODO: Check everything is in order before running
    // BIG TODO: Format output!
    
    // some sorting and cleanup before we begin
    sort(mcSigTemplVal.begin(), mcSigTemplVal.begin() + mcSigTemplVal.size());
    transform(sProp.begin(), sProp.end(), sProp.begin(), ::tolower);
    transform(sAcro.begin(), sAcro.end(), sAcro.begin(), ::tolower);

    setup();

    fitAll();
    minimize(false);

    getCalibration(nPseudoexperiments);
    calibrate();
}
