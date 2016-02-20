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

class RPoissonAnalysis
{
    public:
        /////////////
        // Methods //
        /////////////

        // simple access functions. TODO: not used. keep?
        TH1F* getTemplHisto(string process, int iProp);
        TH1F* getBkgHisto(string process);
        void setOverflowBins(TH1F*);

        // methods to generate individual toys, get all toys
        int  generateToy(int templToUse);
        void doToys(int nExp, int iTemplate);

        // once we've input datasets, assemble them into useful arrays
        void assembleDatasets();

        // fitting methods: fit individual propVal points, get all fits
        double fitPoint(int index);
        int    fitAll();

        // minimize the -log(L/L_max)
        pair<double, double> minimize(bool fake, 
                                      bool all = true,
                                      int  pointsToUse = 2);

        // methods to retrieve and apply the bias correction
        void getCalibration(int numberOfExps = 1000);
        void calibrate();

        ///////////////////////
        // Utility variables //
        ///////////////////////

        // the canvas object to print to
        TCanvas *c_min;
        
        // statistics on sample sizes, sig/bkg-specific
        vector<int> nTotSample;
        int totGenSig,
            totGenBkg;

        // chisquared statistics generated in fitAll()
        vector<double> chiSquared;
        
        // the last template to be fitted
        int fittedTempl = -1;        

        // keeping track of the fit values/errors, both calibrated and uncalibrated
        float fitVal, fitErr,
              calVal, calErr;

        // randomization object for Poisson/Gaussian simulations
        TRandom3 _random;

        /////////////////////////
        // Data and histograms //
        /////////////////////////

        // maps of dataset histos and their signal/bkg integrals
        map<string, TH1*>   datasets;
        map<string, int>    genSig, genBkg;

        // histos which are stored nicely using assembleDatasets() (?)
        map<string, TH1F*>  mcBkgHistosScaled,
                            mcTotalBkgHistosScaled,
                            mcTotalBkgHistosScaled_gen,
                            mcSigTemplHistosScaled,
                            mcSigTemplHistosScaled_gen;

        // histogram for storage of toyDataHisto (when it needs to be modified)
        TH1F *dataHisto;

        // toy histograms: stores statistics across pseudoexperiments
        TH1F *toyDataHisto,
             *toyMean,          // mean of sample
             *toyBias,          // bias  = 
             *toyError,         // error = 
             *toyPull;          // pull  = 
        TH2F *toyLL;            // LL residuals = (TODO HERE)

        // keep track of how many fits fail and succeed
        int  nFitFailed = 0,
             nFitTried  = 0;

        // LL graph (TODO: not really used. remove?)
        TGraph *gr;
        TGraphErrors *grc;
        

        //////////////////////
        // RooFit variables //
        //////////////////////
        
        RooWorkspace *workspace;   // the workspace for our analysis

        RooRealVar   *bkgMean,     // the mean of the background sample we're looking at
                     *bkgWidth,    // the width of the gaussian fit to the bkg sample
                     *propVal;     // the actual sProp value we're optimizing

        RooDataHist  *data,        // histogram to keep track of data samples
                     *bkgHisto;    // histogram to keep track of bkg samples

        RooFitResult *fitResults;  // storage of LL (?) fit results (mostly for chi2)

        RooCategory  *sample;      // basically a constant for this analysis

        RooAbsPdf    *pdffit;      // how we get the fit of the data stored into fitResults

        map<string, RooRealVar *> nbrBkgEvts; // keep track of our background sample sizes

    public:
        /////////////
        // Methods //
        /////////////
        
        // construct/destruct
        RPoissonAnalysis();
        ~RPoissonAnalysis() {};

        // usage methods
        void setup();   // if the user wants to apply new settings
        void run();     // if the user wants to run the analysis


        ////////////
        // Labels //
        ////////////

        // property name, analysis acronym
        string sProp   = "width"; // the property we're trying to find (one word)
        string sAcro   = "mlbwa"; // the acronym describing this analysis
        string sUnit   = "GeV";   // the units of sProp

        // luminosity and 'b^{-1}' prefix
        float  lLumi   = 19.7;    // luminosity of samples
        string lUnit   = "f";     // SI prefix of inverse barns ("f" = "femto")

        // labels for the signal process, subprocesses, background processes
        string          sigProcess  =   "TTbar";
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
        string  dataFileLoc = "./samples/2012_combined_EACMLB.root",
               calibFileLoc = "./samples/calibration_19.700000762939453.root";
        
        // the signal template propVal values
        vector<float> mcSigTemplVal = { 1.50, 3.00, 4.50, 6.00, 7.50 };
        
        // what is the index of the nominal template in mcSigTemplVal?
        int nomTemplIndex      = 0;
        
        // whether to calibrate min/max of mcSigTemplVal[i]
        bool calibLastPts    = true;

        // cut information
        float minPropVal     = 0;       // minimum sProp to consider
        float maxPropVal     = 10;      // maximum sProp to consider
        int   minNumBkg      = 0;       // minimum number of background events
        int   maxNumBkg      = 10000;   // maximum number of background events

        // how many pseudoexperiments to perform in the calibration?
        int nPseudoexperiments = 100;

        // canvas sizes
        float canvWid = 600;
        float canvHei = 600;
        int   iPeriod =   4;     // for TStyleHandler().CMS_lumi()

        // plot binnings
        //  - bins = # bins to use in histogram/graph
        //  - delt = delta(x) to plot about the value of interest
        //  - min  = minimum value for axis
        //  - max  = maximum value for axis
        int   toyMeanBins  =  100;
        float toyMeanDelt  =  2.5;
        int   toyBiasBins  =  140;
        float toyBiasDelt  =  3.5;
        int   toyPullBins  =  200;
        float toyPullDelt  =  10.;
        
        int   toyErrorBins = 2000;
        float toyErrorMin  =    0;
        float toyErrorMax  =  0.4;

        int   toyLLResXBin =    9;
        float toyLLResXMin = -0.5;
        float toyLLResXMax =  8.5;
        int   toyLLResYBin =  200;
        float toyLLResYMin = -100;
        float toyLLResYMax =  100;

        // in toyDataHisto, the number of bins is 
        //      floor[(maxPropVal - minPropVal)/(this number)]
        float toyDataHistoDiv = 10;
};


#endif
