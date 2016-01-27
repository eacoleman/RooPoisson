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
//#include "RPoissonFitHandler.h"

class RPoissonAnalysis
{
    //friend class RPoissonFitHandler;

    private:
        /////////////
        // Methods //
        /////////////
        TH1F* getTemplHisto(string process, int iProp);
        TH1F* getBkgHisto(string process);

        void typeTag(char* nameToTag);

        void doToys(int nExp, int iTemplate);

        void setOverflowBins(TH1F*);

        void assembleDatasets();

        // Fitting (TODO: Migrate to RPoissonFitHandler)
        double fitPoint(int index);
        int    fitAll();
        pair<double, double> minimize(bool fake, 
                                      bool all = true,
                                      int  pointsToUse = 2);

        void runCalibration(int numberOfExps = 1000);

        ///////////////////////
        // Utility variables //
        ///////////////////////
        
        //
        //RPoissonFitHandler *fitHandler;
        TCanvas *c_min;
        
        //
        vector<int> nTotSample;
        int totGenSig,
            totGenBkg;

        //
        vector<double> chiSquared;
        
        //
        int fittedTempl;        

        //
        float fitVal,
              calVal;

        // fitting
        TRandom3 _random;

        /////////////////////////
        // Data and histograms //
        /////////////////////////
        
        map<string, TH1*>   datasets;
        map<string, int>    genSig, genBkg;

        map<string, TH1F *> mcBkgHistosScaled,
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

        TH1F *hFitFailed,
             *hFitFailedDiff;

        int  nFitFailed,
             nFitTried;

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
        void run(char* dataFileName);
        void save(char* outFileName);


        ////////////
        // Labels //
        ////////////

        // property name, analysis acronym
        string sProp   = "width";
        string sAcro   = "mlbwa";
        string sUnit   = "GeV";

        // luminosity and 'b^{-1}' prefix
        float lLumi   = 19.7;
        string lUnit   = "f";

        // labels for subprocesses, background processes
        vector<string> processes;
        vector<TString> mcBkgLabels;


        ///////////////////////////////////
        // Variables to control settings //
        ///////////////////////////////////
       
        string dataFileLoc,
               systFileLoc;
        
        vector<float> mcSigTemplVal;
        
        // whether to use systematics, PDF version
        //
        bool systematics     = true;
        bool systematicsPDF  = false;
        bool useRatio        = true;
        bool fixedSample     = true;
        bool bkgSyst         = false;

        // cut information
        float upperCut       = 400.;
        float lowerCut       = 0.;
        float upperWeightCut = 5000000.;
        float lowerWeightCut = 5.;

        float minPropVal     = 0;
        float maxPropVal     = 10;
        
        // what is the index of the nominal template
        // in mcSigTemplVal?
        int nomTemplIndex    = 0;

        // antiquated - consider removing
        int fixBkg           = 3;

};

#endif
