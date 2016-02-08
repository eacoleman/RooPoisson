/*******************************************************************************
* RPoissonAnalysis.C
* 
*
*******************************************************************************/

// RooPoisson Inclusions
#include "include/RPoissonAnalysis.h"
#include "TStyleHandler.C"

using namespace RooFit;

RPoissonAnalysis::RPoissonAnalysis() {

}

void RPoissonAnalysis::setup() {

    //
    TStyleHandler::setTDRStyle();

    int W = 800;                                                                                    //HARDCODED
    int H = 600;                                                                                    //HARDCODED

    // Simple example of macro: plot with CMS name and lumi text
    //  (this script does not pretend to work in all configurations)
    // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
    // For instance: 
    //               iPeriod = 3 means: 7 TeV + 8 TeV
    //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
    // Initiated by: Gautier Hamel de Monchenault (Saclay)
    int H_ref = 600;                                                                                //HARDCODED
    int W_ref = 800;                                                                                //HARDCODED

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref; 
    float L = 0.16*W_ref;
    float R = 0.04*W_ref;

    // setup a new canvas
    c_min = new TCanvas("c_min","", 600, 600);                                                      //HARDCODED
    c_min->SetFillColor(0);
    c_min->SetBorderMode(0);
    c_min->SetFrameFillStyle(0);
    c_min->SetFrameBorderMode(0);
    c_min->SetLeftMargin( L/W );
    c_min->SetRightMargin( R/W );
    c_min->SetTopMargin( T/H );
    c_min->SetBottomMargin( B/H );
    c_min->SetTickx(0);
    c_min->SetTicky(0);

    // initialize some helper variables
    nFitTried      = 0;
    nFitFailed     = 0;
    hFitFailed     = new TH1F("hFitFailed","hFitFailed", 11, minPropVal, maxPropVal);               //HARDCODED
    hFitFailedDiff = new TH1F("hFitFailedDiff","hFitFailedDiff", 60, -30., 30.);                    //HARDCODED
    fittedTempl    = -1;

    nTotSample.resize(processes.size());

    //
    char dataHistoName[20] = "PeakMassTree_";
    char evHistoName[15]   = "h_nEvents";

    //
    char propValPref[15] = "Reconstructed ";
    propVal = new RooRealVar(sProp.c_str(), strcat(propValPref, sProp.c_str()),
                             lowerCut, upperCut);
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
        if (fixedSample) nTotSample.at(i) = datasets[processes.at(i)]->GetEntries();
        cout << " - in "<< processes.at(i) << " " << nTotSample.at(i) << endl;
    }

    //
    toyDataHisto = new TH1F("toyDataHisto", "toyDataHisto", 40, 0., 400);                           //HARDCODED
    toyDataHisto->Reset();

    //
    data->fillHistogram(toyDataHisto, *propVal);
    dataHisto = (TH1F*) toyDataHisto->Clone("dataHisto");

    // Get signal propVal templates
    templMean   = new TH1F("tmean", "Mean of the templates", 50, minPropVal, maxPropVal);           //HARDCODED
    templRMS    = new TH1F("trms" , "RMS of the templates" , 50, minPropVal, maxPropVal);           //HARDCODED

    // open data file
    theFile = new TFile (dataFileLoc.c_str());
    theFile->cd();

    // loop through the interaction types we want to check
    for (int itype = 0; itype < processes.size(); ++itype) {

        // loop through the propVals we want to check
        for (unsigned int iVal = 0; iVal <  mcSigTemplVal.size(); ++iVal) {

            // get the histogram for the interaction and propVal we want
            sprintf(hname, "%s__TTbar_%.2f_%s", sAcro.c_str(), mcSigTemplVal.at(iVal), processes.at(itype).c_str());
            sprintf(tag  , "%s%.2f"              , processes.at(itype).c_str()   , mcSigTemplVal.at(iVal));
            cout << " - signal template: " << hname << " " << tag << endl;

            printf("itype %d iVal %d \n", itype, iVal);

            TH1F* histo = (TH1F*) gDirectory->Get(hname);

            if (histo==0) { 
                cout << "Histo does not exist\n";
                histo = new TH1F(hname,hname,100,0,200);                                            //HARDCODED
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
    if (systematics) {
        cout << " - retrieved systematics signal GEN template from " 
             << systFileLoc << endl;
        theFile = new TFile (systFileLoc.c_str()) ;

        unsigned int maxTemplates;
        if (systematicsPDF) maxTemplates = 41;                                                      //HARDCODED
        else maxTemplates = mcSigTemplVal.size();

        // loop through the propVals
        for (unsigned int i = 0; i < maxTemplates; ++i){

            float tVal = 0.;
            if (!systematicsPDF) tVal = mcSigTemplVal.at(i);

            //loop through the interaction types
            for (int itype = 0; itype < processes.size(); ++itype) {
                // get the histogram name
                if (systematicsPDF) {
                    sprintf(tag,"%s_pdf%i",processes.at(itype).c_str(),i);
                    sprintf(hname,"%s__%s_pdf%i", sAcro.c_str(), processes.at(itype).c_str(),i);
                } else {
                    sprintf(tag,"%s%.2f",processes.at(itype).c_str(),tVal);
                    sprintf(hname, "%s__TTbar_%.2f_%s", sAcro.c_str(), tVal, processes.at(itype).c_str());
                }

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
    }

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
    if (systematics) {
        //
        cout << " - retrieving systematics background GEN template from "
             << systFileLoc << endl;
        theFile = new TFile (systFileLoc.c_str());

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
                    histo = new TH1F(hname, hname, 100, 0, 200);                                    //HARDCODED
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
    }

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

        if (systematics) {
            sprintf(hname,"histo_bck_gen%s",processes.at(itype).c_str());
            workspace->import(*(new RooDataHist(hname, hname, *propVal, 
                                                mcTotalBkgHistosScaled_gen[processes.at(itype)])));

            sprintf(hname,"HistPdf::background_gen%s(%s,histo_bck%s)", 
                    processes.at(itype).c_str(), sProp.c_str(),  processes.at(itype).c_str());
            workspace->factory(hname);
        }

        //
        float totalBkg = (float) mcTotalBkgHistosScaled[processes.at(itype)]->Integral();

        //
        sprintf(hname,"nbrBkgEvts%s", processes.at(itype).c_str());
        cout << hname << endl;

        //
        nbrBkgEvts[processes.at(itype)] 
            = new RooRealVar(hname, hname, totalBkg, 0, 10000);                                     //HARDCODED

        //
        workspace->import(*nbrBkgEvts[processes.at(itype)]);

        // loop through the propVals
        for (unsigned int i = 0; i < mcSigTemplVal.size(); ++i) {
            float tVal = mcSigTemplVal.at(i);
            cout << " - building "  << processes.at(itype) << " pdf for "
                 << sProp << " " << tVal << endl;

            sprintf(tag, "%s%.2f", processes.at(itype).c_str(), tVal);

            //
            sprintf(tag  , "%s%.2f"         , processes.at(itype).c_str(), tVal);
            sprintf(hname, "histo_sgn%s%.2f", processes.at(itype).c_str(), tVal);
            cout << hname << " " << tag << endl;
            
            //
            workspace->import( *(new RooDataHist(hname, hname, *propVal, mcSigTemplHistosScaled[tag]) ));
            sprintf(hname,"HistPdf::signal%s%.2f(%s,histo_sgn%s%.2f)",
                    processes.at(itype).c_str(),tVal,sProp.c_str(),processes.at(itype).c_str(),tVal);
            cout << hname << endl;
            
            //
            workspace->factory(hname);

            if (systematics && !systematicsPDF) {
                //
                sprintf(hname, "histo_sgn_gen%s%.2f", processes.at(itype).c_str(), tVal);
                workspace->import( *(new RooDataHist(hname, hname, *propVal, 
                                             mcSigTemplHistosScaled_gen[tag])));

                //
                sprintf(hname,"HistPdf::signal_gen%s%.2f(propVal,histo_sgn_gen%s%.2f)",
                        processes.at(itype).c_str(), tVal, processes.at(itype).c_str(), tVal);
                workspace->factory(hname);
            }

            //
            if (useRatio) {
                sprintf(hname,"SUM::model%s%.2f( ratio%s%.2f[0,1]*signal%s%.2f, background%s )",
                        processes.at(itype).c_str(), tVal, processes.at(itype).c_str(), tVal, 
                        processes.at(itype).c_str(), tVal, processes.at(itype).c_str());
            } else {
                sprintf(hname,"SUM::model%s%.2f( Nsig%s%.2f[0,1000]*signal%s%.2f, nbrBkgEvts%s*background%s )",
                        processes.at(itype).c_str(), tVal, processes.at(itype).c_str(), tVal, 
                        processes.at(itype).c_str(), tVal, processes.at(itype).c_str(), 
                        processes.at(itype).c_str());
            }

            cout << hname << endl;

            workspace->factory(hname);
        }

        if (systematics && systematicsPDF) {
            for (int i = 0; i < 41; ++i){                                                           //HARDCODED
                sprintf(hname,"histo_sgn_gen%s%i",processes.at(itype).c_str(),i);
                sprintf(tag,"%s_pdf%i",processes.at(itype).c_str(),i);
                workspace->import( *(new RooDataHist(hname, hname, *propVal, 
                                             mcSigTemplHistosScaled_gen[tag])));

                sprintf(hname, "HistPdf::signal_gen%s%i(%s,histo_sgn_gen%s%i)",
                        processes.at(itype).c_str(), i, sProp.c_str(), processes.at(itype).c_str(), i);
                cout << hname << endl;

                workspace->factory(hname);
            }
        }

    }

    // this is only for the pdf systematics templates:
    char name2[300];
    for (unsigned int i = 0; i < mcSigTemplVal.size(); ++i) {
        float tVal = mcSigTemplVal.at(i);
        cout << " - building simultaneous pdf for " << sProp << " " 
             << tVal << endl;

        sprintf(name2, "SIMUL::model%.2f(sample", tVal);
        for (int j = 0; j < processes.size(); ++j) {
            sprintf(name , "%s,%s=model%s%.2f", name2, processes.at(j).c_str(),
                    processes.at(j).c_str(), tVal);
            sprintf(name2, "%s", name);
        }

        sprintf(name,"%s)",name2);

        workspace->factory(name);
    }

    if (bkgsyst) {
        char bkgMeanPref[5] = "bkg ";
        bkgMean     = new RooRealVar("bkgmean" , strcat(bkgMeanPref, sProp.c_str()), 180.);
        bkgWidth    = new RooRealVar("bkgwidth", "bkg fit width", 20.);
        bkgHistoPDF = new RooGaussian("background", "background PDF",
                                      *propVal,*bkgMean,*bkgWidth);
    }

    toyError = 0;
    toyLL    = new TH2F("LL", "LL residuals", 9, minPropVal, maxPropVal,                            //HARDCODED
                        200, -100, 100);                                                            //HARDCODED
    
    if (!useRatio) {
        for (int itype = 0 ; itype < processes.size(); ++itype) {
            sprintf(hname,"nbrBkgEvts%s", processes.at(itype).c_str());
            workspace->var(hname)->setVal(mcTotalBkgHistosScaled[processes.at(itype)]->Integral());
            workspace->var(hname)->setConstant(1);
        }
    }
}


TH1F* RPoissonAnalysis::getTemplHisto(string process, int iProp) {
    //
    char tag[50];
    sprintf(tag,"%s%.2f", process.data(), mcSigTemplVal.at(iProp));

    //
    cout << "Template " << sProp << ": " << tag << endl;
    return mcSigTemplHistosScaled[tag];
}


TH1F* RPoissonAnalysis::getBkgHisto(string process) {
    return mcTotalBkgHistosScaled[process.data()];
}


void RPoissonAnalysis::typeTag(char* nameToTag) {
    char tName[200];

    cout << " - Name to modify is: " << nameToTag << endl;

    for (int j = 0; j < processes.size(); ++j) {
        sprintf(tName,"%s%s", nameToTag, processes.at(j).c_str());

        if (j<(processes.size()-1)) sprintf(nameToTag,"%s_", tName);
        else sprintf(nameToTag,"%s", tName);
    }

    cout << " - New name is: " << nameToTag << endl;
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

        if (!fixedSample) {

            if (!systematics || !systematicsPDF) 
                 sprintf(tag, "%s%.2f"  , processes.at(itype).c_str(), mcSigTemplVal[templToUse]);
            else sprintf(tag, "%s_pdf%i", processes.at(itype).c_str(), templToUse);

            int binLow  = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(lowerCut);
            int binHigh = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(upperCut);

            if (systematics) {
                n = _random.Poisson((float) mcSigTemplHistosScaled_gen[tag]->Integral(binLow, binHigh));
            } else {
                n = _random.Poisson((float) mcSigTemplHistosScaled[tag]->Integral(binLow, binHigh));
            }

            cout << " - generating "<< n << " " << processes.at(itype) << " signal events of ";

            if (!systematics || !systematicsPDF) cout << sProp  << " " << mcSigTemplVal[templToUse];
            else                                 cout << "pdf " << templToUse;
            
            cout << ", with poisson mean "
                 << (float) mcSigTemplHistosScaled[tag]->Integral(binLow, binHigh)
                 << endl;

            genSig[processes.at(itype)] = n; 
            totGenSig += n;

            if (systematics) {
                n = _random.Poisson((float) mcTotalBkgHistosScaled_gen[processes.at(itype)]->Integral(binLow, binHigh));
            } else {
                n = _random.Poisson((float) mcTotalBkgHistosScaled[processes.at(itype)]->Integral(binLow, binHigh));
            }

            genBkg[processes.at(itype)] = n; 
            totGenBkg += n;

            cout << " - generating "<< n << " " << " background events"
                 << ", with Poisson mean "
                 << (float) mcTotalBkgHistosScaled[processes.at(itype)]->Integral(binLow, binHigh)
                 << endl;

        } else {
            sprintf(tag, "%s%.2f", processes.at(itype).c_str(), mcSigTemplVal[4]);                  //HARDCODED REDALERT

            int binLow  = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(lowerCut);
            int binHigh = mcSigTemplHistosScaled[tag]->GetXaxis()->FindBin(upperCut);

            cout << tag << endl;
            double sigProb;

            if (systematics) {
                if (systematicsPDF) sprintf(tag, "%s_pdf%i", processes.at(itype).c_str(), templToUse);
                sigProb = mcSigTemplHistosScaled_gen[tag]->Integral(binLow, binHigh) /
                            (mcSigTemplHistosScaled_gen[tag]->Integral(binLow, binHigh)
                              + mcTotalBkgHistosScaled_gen[processes.at(itype)]->Integral(binLow, binHigh));
            } else {
                sigProb = mcSigTemplHistosScaled[tag]->Integral(binLow, binHigh) /
                            (mcSigTemplHistosScaled[tag]->Integral(binLow, binHigh)
                              + mcTotalBkgHistosScaled[processes.at(itype)]->Integral(binLow, binHigh));
            }

            n = _random.Binomial(nTotSample[itype], sigProb);
            cout << "Generate " << n << " " << processes.at(itype) << " signal events of ";

            if (!systematics || !systematicsPDF) cout << sProp  << " " << mcSigTemplVal[templToUse];
            else                                 cout << "pdf " << templToUse;
            
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
        }

        if (systematics) {

            if (!systematicsPDF) sprintf(hname,"signal_gen%s%.2f", 
                                         processes.at(itype).c_str(),
                                         mcSigTemplVal[templToUse]);
            else sprintf(hname,"signal_gen%s%i",
                         processes.at(itype).c_str(), templToUse);
        } else {
            sprintf(hname,"signal%s%.2f", processes.at(itype).c_str(),
                    mcSigTemplVal[templToUse]);
        }

        cout << hname << endl;
        datasets[processes.at(itype)] 
            = workspace->pdf(hname)->generateBinned(*propVal, 
                                                    genSig[processes.at(itype)])->createHistogram(hname,*propVal);
        cout << hname << endl;

        if (systematics) {
            sprintf(hname,"background_gen%s", processes.at(itype).c_str());
        } else {
            sprintf(hname,"background%s", processes.at(itype).c_str());
        }

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
    if (!systematics || !systematicsPDF) 
         cout << " - template " << sProp << " value "
              << mcSigTemplVal.at(iTemplate) << endl;
    else cout << " - pdf " << iTemplate << endl;

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
    
    toyMean   = new TH1F("mean"  , sPropPref, 100, minPropVal, maxPropVal);                         //HARDCODED
    toyBias   = new TH1F("bias"  , strcat(sPropPref, " bias"), 100, -3.5, 3.5);                     //HARDCODED
    toyPull   = new TH1F("pull"  , "pull", 200, -10, 10);                                           //HARDCODED
    toyError  = new TH1F("error" , strcat(sPropPref, " uncertainty"), 500, 0, 0.4);                 //HARDCODED
    toyLL     = new TH2F("LL"    , "LL residuals", 9, -0.5, 8.5, 200, -100, 100);                   //HARDCODED

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
        cout << "EAC683" << endl;
        pair<double,double> result = minimize(true);cout << "EAC684" << endl;
        if (result.second >= 0.) {cout << "EAC685" << endl;
            cout << " - fake fit result: " << result.first 
                 << " +/- " << result.second << endl;
            cout << "EAC688" << endl;
            // fill the toy histograms with the result of our pseudoexperiment
            toyMean->Fill(result.first);cout << "EAC690" << endl;
            toyBias->Fill(result.first-propPoint);cout << "EAC691" << endl;
            toyPull->Fill((propPoint-result.first)/result.second);cout << "EAC692" << endl;
            toyError->Fill(result.second); cout << "EAC693" << endl;
        } else {
            ++nFitFailed;
        }
        cout << "EAC697" << endl;
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
    if (fitResults)  {
        delete fitResults;

        if (bkgsyst) {
            delete pdffit;
        }
    }

    // name of the model we're fitting
    char tName[50];
    sprintf(tName,"model%.2f", mcSigTemplVal.at(index));

    //
    pdffit = workspace->pdf(tName) ;
    pdffit->Print();

    //
    if (useRatio && fixBkg == 3) {
        char tag[50];

        //
        for (int itype = 0; itype < processes.size(); ++itype) {
            //
            sprintf(tag,"%s%.2f",processes.at(itype).c_str(),mcSigTemplVal.at(index));
            sprintf(tName,"ratio%s%.2f", processes.at(itype).c_str(), mcSigTemplVal.at(index));

            // S/(S+B) significance metric
            workspace->var(tName)->setVal(mcSigTemplHistosScaled[tag]->Integral()/
                                  (mcSigTemplHistosScaled[tag]->Integral()
                                   +mcTotalBkgHistosScaled[processes.at(itype)]->Integral()));
            workspace->var(tName)->setConstant(1);
        }
    } else if (useRatio && fixBkg == 1) {
        //
        for (int itype = 0; itype < processes.size(); ++itype) {
            sprintf(tName,"ratio%s%.2f", processes.at(itype).c_str(), mcSigTemplVal.at(index));
            workspace->var(tName)->setVal(1.0);
            workspace->var(tName)->setConstant(1);
        }
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
        chiSquared.at(i) = fitPoint(i); cout << "EAC812" << endl;
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
    cout << "EAC826 " << chiSquared.size() << endl;
    //
    for (int i = 0; i < mcSigTemplVal.size(); ++i) { cout << " EAC" << i << endl;
        if (chiSquared.at(i) >= 0) {
            ++pts;
            minLL = min(chiSquared[i], minLL);
        }
    }
    cout << "834" << endl;
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
            ey[pts] = 9.4;
            
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
    gr->GetXaxis()->SetTitle("mlbwa  m(lb) [GeV]"); //HARDCODED
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

void RPoissonAnalysis::runCalibration(int numberOfExps = 1000) {
    // keep track of timing
    time_t start, end;
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

    // set the name of the outfile (TODO: abstract)
    TString name = TString("calibration_");
    name += lLumi;
    name.ReplaceAll ( " " , "" );
    cout << name << endl;

    // open the outfile
    TFile * out = new TFile(name + ".root", "RECREATE");
    cout << " - opened output file with name " << name << endl;

    // loop over templates
    for (i = min; i < max; ++i) {
        cout << "doing toy " << i << endl;
        // get the toy statistics and fit them with a gaussian
        doToys(numberOfExps, i); cout << "EAC961 " << toyMean->GetEntries() << endl;
        toyMean->Fit("gaus"); cout << "EAC962" << endl;

        // store the points for bias plot
        x[pts]  = mcSigTemplVal.at(i); 
        y[pts]  = toyMean->GetFunction("gaus")->GetParameter(1);cout << "EAC966" << endl;

        // store the point errors for bias plot
        ex[pts] = 0.;
        ey[pts] = toyMean->GetFunction("gaus")->GetParameter(2)/sqrt(numberOfExps);

        cout << " - y pts is " << y[pts] << endl;

        ++pts;

        // write the toy histograms to our output file
        sprintf(hname,"meanMass_%.2f", mcSigTemplVal.at(i));
        toyMean->Clone(hname)->Write();

        sprintf(hname,"biasMass_%.2f", mcSigTemplVal.at(i));
        toyBias->Clone(hname)->Write();

        sprintf(hname,"errMass_%.2f",  mcSigTemplVal.at(i));
        toyError->Clone(hname)->Write();

        sprintf(hname,"pullMass_%.2f", mcSigTemplVal.at(i));
        toyPull->Clone(hname)->Write();

        sprintf(hname,"LL_%.2f",       mcSigTemplVal.at(i));
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

void RPoissonAnalysis::run(char* dataFileName) {
    setup();
    runCalibration(1000);
    fitAll();
    minimize(false);
}

void RPoissonAnalysis::save(char* outFileName) {

}
