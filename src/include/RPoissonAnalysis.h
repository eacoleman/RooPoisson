#ifndef RPOISSONANALYSIS
#define RPOISSONANALYSIS

// C Library inclusion
#include <cstdlib>
#include <vector>
#include <map>
#include <ctime>
#include <iostream>
#include <fstream>

// ROOT utility inclusion
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <TSystem.h>
#include <TLatex.h>
#include <TROOT.h>

// RooFit inclusion
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

// RooPoisson inclusion

class RPoissonAnalysis
{
    //friend class RPoissonFitHandler;

    public:
        /////////////
        // Methods //
        /////////////
        TH1F* getTemplHisto(string process, int iProp);
        TH1F* getBkgHisto(string process);

        void typeTag(char* nameToTag);

        int  generateToy(int templToUse);
        void doToys(int nExp, int iTemplate);

        void setOverflowBins(TH1F*);

        void assembleDatasets();

        // fitting methods
        double fitPoint(int index);
        int    fitAll();
        pair<double, double> minimize(bool fake, 
                                      bool all = true,
                                      int  pointsToUse = 2);

        void getCalibration(int numberOfExps = 1000);
        void calibrate(char tag[20]);

        ///////////////////////
        // Utility variables //
        ///////////////////////
        TCanvas *c_min;
        
        //
        vector<int> nTotSample;
        int totGenSig,
            totGenBkg;

        //
        vector<double> chiSquared;
        
        //
        int fittedTempl = -1;        

        //
        float fitVal,
              calVal;

        // fitting
        TRandom3 _random;

        //
        bool bkgsyst = false;

        /////////////////////////
        // Data and histograms //
        /////////////////////////

        map<string, TH1*>   datasets;
        map<string, int>    genSig, genBkg;

        map<string, TH1F*>  mcBkgHistosScaled,
                            mcTotalBkgHistosScaled,
                            mcTotalBkgHistosScaled_gen,
                            mcSigTemplHistosScaled,
                            mcSigTemplHistosScaled_gen;

        TH1F *dataHisto,
             *chi2Result,
             *calibHisto;

        TH1F *toyDataHisto,
             *toyMean,
             *toyBias,
             *toyError,
             *toyPull;
        TH2F *toyLL;

        TH1F *templMean,
             *templRMS;

        int  nFitFailed = 0,
             nFitTried  = 0;

        TGraph *gr;
        TGraphErrors *grc;
        

        //////////////////////
        // RooFit variables //
        //////////////////////
        
        RooWorkspace *workspace;

        RooRealVar   *bkgMean, 
                     *bkgWidth,
                     *propVal;

        RooDataHist  *data,
                     *bkgHisto;

        RooFitResult *fitResults;

        RooCategory  *sample;

        RooAbsPdf    *pdffit,
                     *bkgHistoPDF;

        map<string, RooRealVar *> nbrBkgEvts;

    public:
        /////////////
        // Methods //
        /////////////
        
        // construct/destruct
        RPoissonAnalysis();
        ~RPoissonAnalysis() {};

        // utils
        void setup();
        void run();
        void save(char* outFileName);


        ////////////
        // Labels //
        ////////////

        // property name, analysis acronym
        string sProp   = "width";
        string sAcro   = "mlbwa";
        string sUnit   = "GeV";

        // luminosity and 'b^{-1}' prefix
        float  lLumi   = 19.7;
        string lUnit   = "f";

        // labels for the signal process, subprocesses, background processes
        string         sigProcess  =   "TTbar";
        vector<string>  processes   = { "E", "EE", "EM", "MM", "M" };
        vector<TString> mcBkgLabels = { "SingleTop", 
                                        "WJets", 
                                        "DrellYan", 
                                        "Diboson", 
                                        "QCD" };


        ///////////////////////////////////
        // Variables to control settings //
        ///////////////////////////////////
       
        // location of your data files, output calib file
        string dataFileLoc = "./samples/2012_combined_EACMLB.root",
               systFileLoc = "./samples/calibration_19.700000762939453.root";
        
        // the signal template propVal values
        vector<float> mcSigTemplVal = { 1.50, 3.00, 4.50, 6.00, 7.50 };
        
        // whether to calibrate min/max of mcSigTemplVal[i]
        bool calibLastPts    = true;

        // cut information
        float minPropVal     = 0;
        float maxPropVal     = 10;
        
        // what is the index of the nominal template in mcSigTemplVal?
        int nomTemplIndex      = 0;

        // how many pseudoexperiments to perform in the calibration?
        int nPseudoexperiments = 100;

};

#endif
