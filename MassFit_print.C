/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MassFit_print.C                                                             *
 * Author: Richard Nally et al.                                                *
 *         2015                                                                *
 *                                                                             *
 * Purpose:                                                                    *
 * Uses 2012_combined MC data to analytically determine the mass of the top    *
 * quark. Reads from ./2012_combined.root                                      *
 *                                                                             *
 * Usage:                                                                      *
 * Launch root via "root -l". Enter the following:                             *
 *                                                                             *
 *      .L MassFit_print.C                                                     *
 *      j = new MassFit()                                                      *
 *      j->fitPoint(int)  // for a single mass                                 *
 *      j->fitAll()       // for all masses                                    *
 *      j->findMin()      // to find the minimum after fitAll()                *
 *      j->printLL()      // print and be done!                                *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TSystem.h>
//#include <Cintex/Cintex.h>
#include <utility>
#include "TH1F.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TChain.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <map>
#include <ctime>
#include <TKey.h>
#include <TRegexp.h>
#include <TLatex.h>
#include <TROOT.h>
#include "TStyle.h"
#include <TMath.h>
#include "TASImage.h"
#include "CMS_lumi.C"

using namespace std;
using namespace RooFit;

static const bool systematics = true;
static const bool systematicsPDF = false;

static const bool useRatio = true;

static const int maxType = 5;
static const int maxSignalType = 2;
static const char type [maxType][7] =//{"EE_1W", "EE_2W"};
                                     //{"EE","EM","MM"};
                                     //{"EE_1B", "EM_1B", "MM_1B", "EE_2B", "EM_2B", "MM_2B"};
                                     //{"EE_2B"};
                                     //{"EE_0B", "EM_0B", "MM_0B"};
                                     //{"EE_1B", "EM_1B", "MM_1B"};
                                     //{"EE_2B", "EM_2B", "MM_2B"};
                                     //{"0B","1B","2B"};
                                     //{"2B"};
                                     //{"EE_0B", "EE_1B", "EE_2B",
                                     // "EM_0B", "EM_1B", "EM_2B",
                                     // "MM_0B", "MM_1B", "MM_2B"};
                                     //{"EM"};
                                     //{"MM;1", "EM;1", "EE;2", "EE;1"};
                                     //{"2B_All"};
                                       {"M", "MM", "EM", "EE", "E"};

static const bool fixedSample = true;
static int nTotSample[maxType];

static const float upperMassCut = 400.;
static const float lowerMassCut = 0;//100.;
static const float upperWeightCut = 5000000.;
static const float lowerWeightCut = 5.;

static const bool debug = true;

static char absolutePath[100] = "./";

class MassFit {
    public:
        MassFit();
        MassFit(string genFile);
        ~MassFit();
        double fitPoint(int i = 15);
        double fitData();
        int fitAll();
        TStyle* setTDRStyle();
        pair<double,double> findMin(bool all=true, int pointsToUse=2);
        pair<double,double> findMinFake(bool all=true, int pointsToUse=2);
        int fixBackground();
        void fixBackground(int fix);
        void fixBackground(float e, float ee, float em, float mm, float m);
        void quietFit(bool quiet) {quietFit_ = (quiet?-1:0);}
        void write_toy(TString name);
        void calib(char fn[100] = "calibration_19.700000762939453.root",char tag[20] = "");

        void do_toys(int n_exp, int templateToUse);
        int generate_toy(int templateToUse);
        void calibration(int number = 1000);
        void set_overflow_bins(TH1F * h);
        TH1F* templateHisto(const char* type, int i=0);
        TH1F* backgroundHisto(const char* type);
        void printFit(bool withFit = true, bool toFile = false, char* name = "");
        void printFit(int point, bool toFile = false, char* name = "");
        void checkMassRange(string histoName);
        void printMassRange();
        void typeTag(char *name);
        void printLL(bool toFile = false, char* name = "");

        void texEvents(int templateToUse, float lower = lowerMassCut, float upper = upperMassCut);

        int templateNr, nbrTemplates;
        int minTMass, maxTMass;
        float fitMass, calMass;

        string DataFileLocation, MCFileName, SystFileLocation;
        vector<TString> mcBackgroundLabels;
        vector<float> mcSignalTemplMass;
        TH1F* dataHisto;
        map<string, TH1F *> mcBackgroundHistosScaled;
        map<string, TH1F *> mcTotalBackgroundHistoScaled, mcTotalBackgroundHistoScaled_gen;
        map<string, TH1F *> mcSignalTemplHistosScaled, mcSignalTemplHistosScaled_gen;
        vector <double> chiSquared;
        TH1F* chi2Result;
        float ilumi;
        double signalXs;
        int minTemplate, maxTemplate, nominalTemplate, toyTemplate;
        TH1F *toyDataHisto;
        TH1F *toy_mean, *toy_bias;
        TH1F *toy_error;
        TH1F *toy_pull;
        TH2F *toy_LL;
        TH1F *calibrationH;
        TGraph *gr;
        TGraphErrors * grc;
        int nFitFailed;
        int nFitTried;
        TH1F *hFitFailed;
        TH1F *hFitFailedDiff;
        TFile *FitFailedFile;
        TH1F *templ_mean, *templ_rms ;
        int totalGeneratedSignal, totalGeneratedBkg;
        bool drawSwitch;
        int counter;
        RooRealVar * topMass;
        int fittedTemplate;

        RooDataHist *histo_bck;
        RooAbsPdf  *histo_bck_pdfff;
        RooFitResult* myFitResults_all;
        RooRealVar *bkgmean, *bkgwidth;
        RooAbsPdf* pdffit;
        RooExtendPdf *extBackgroundPdf;
        RooDataHist *data;
        bool bkgsyst;
        RooWorkspace *w;
        map<string, TH1*> datasets;
        map<string, RooRealVar *> nbrBackgroundEvents;
        map<string, int> generatedSignal, generatedBkg;
        RooCategory *sample;

        void pdfCalibration(int number);
        int generate_toy_pdf(int templateToUse);
        void cmsprelim();

    private:
        void assembleDatasets();
        void setup();
        void defaultData();
        void defaultMCsignal();
        void defaultMCbackground();
        void setMCsignal(string defaultFile, string genFile);
        void setMCsignal(string genFile);

        TRandom3 _random;
        int fixBckg;
        TCanvas *c_min;
        int quietFit_;

};

MassFit::MassFit()
{
    ilumi = 19.7;
    defaultData();
    defaultMCsignal();
    defaultMCbackground();
    setup();
}

MassFit::MassFit(string genFile)
{
    ilumi = 19.7;
    defaultData();
    defaultMCsignal();
    defaultMCbackground();
    SystFileLocation = genFile;
    setup();
}

void MassFit::defaultData()
{
    DataFileLocation = "./2012_combined_EACMLB.root";
}

TStyle* MassFit::setTDRStyle() {
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    // For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(600); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //Position on screen
    tdrStyle->SetCanvasDefY(0);

    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);

    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);

    // For the histo:
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    tdrStyle->SetEndErrorSize(2);
    tdrStyle->SetMarkerStyle(20);

    //For the fit/function:
    tdrStyle->SetOptFit(0);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);

    //For the date:
    tdrStyle->SetOptDate(0);

    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS: SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.5g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);

    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.16);
    tdrStyle->SetPadRightMargin(0.02);

    // For the Global title:
    tdrStyle->SetOptTitle(0);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);

    // For the axis titles:
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");
    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);

    // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XYZ");

    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);

    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);

    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
    tdrStyle->SetHatchesLineWidth(5);
    tdrStyle->SetHatchesSpacing(0.05);

    tdrStyle->cd();
    return tdrStyle;
//    gROOT->SetStyle("tdrStyle");
}

void MassFit::calib(char fn[100],char tag[20]){
    float massPoints[5] = { 1.5, 3.0, 4.5, 6.0, 7.5 };// { 169.5, 171.5, 172.5, 173.5, 175.5};
    float fitMass = 1.5;
    bool reduced = true;
    float fitUnc = 0.16;

    TCanvas* c_min = new TCanvas("c_min","", 600, 600);
    int W = 800;
    int H = 600;
    int H_ref = 600;
    int W_ref = 800;

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.16*W_ref;
    float R = 0.04*W_ref;

    c_min = new TCanvas("c_min","", 600, 600);

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

    writeExtraText = true;                   // if extra text
    extraText  = "Preliminary";               // default extra text is "Preliminary"
    //lumi_8TeV  = "19.7 fb^{-1}, #sqrt{s}=8"; // default is "19.7 fb^{-1}"
    //lumi_7TeV  = "4.9 fb^{-1}";              // default is "5.1 fb^{-1}"
    int iPeriod = 0;                         // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV

    char name[200], hname[50];

    TFile* theFile = new TFile (fn) ;
    int points = 0,minP,maxP;
    float minTMass = 0, maxTMass = 0;
    for (int i = 0;i<5;++i) {
        float mass = massPoints[i];
        sprintf(hname,"meanMass_%.2f", mass);
        TH1F*   toy_mean  = (TH1F*) gDirectory->Get(hname);
        cout << hname << endl;

        if (toy_mean == 0) continue;
        if (reduced && (mass<0.5||mass>8.5)) continue;

        ++points;
        if (minTMass==0.) {minTMass = mass; minP= i;}
        maxTMass = mass;
        maxP= i;
    }
    cout << points << " " << minTMass<< " " << maxTMass << endl;

    TVectorD massV(points), meanV(points), biasV(points), pullV(points), pullWV(points);
    TVectorD massErrV(points), meanErrV(points), biasErrV(points), pullErrV(points), pullWErrV(points);

    int pts = 0;
    for (int i = minP;i<=maxP;++i,++pts) {
        float mass = massPoints[i];
        sprintf(hname,"meanMass_%.2f", mass);
        TH1F* toy_mean  = (TH1F*) gDirectory->Get(hname);
        toy_mean->Fit("gaus","Q");

        sprintf(hname,"biasMass_%.2f", mass);
        TH1F* toy_bias= (TH1F*) gDirectory->Get(hname) ;
        sprintf(hname,"errMass_%.2f", mass);
        TH1F* toy_error= (TH1F*) gDirectory->Get(hname) ;
        sprintf(hname,"pullMass_%.2f", mass);
        TH1F* toy_pull= (TH1F*) gDirectory->Get(hname) ;

        toy_bias->Fit("gaus","Q");
        toy_pull->Rebin(5);
        toy_pull->Fit("gaus","Q");

        massV(pts)=mass;
        meanV(pts)=toy_mean->GetFunction("gaus")->GetParameter(1);
        biasV(pts)=toy_bias->GetFunction("gaus")->GetParameter(1);
        pullV(pts)=toy_pull->GetFunction("gaus")->GetParameter(1);
        pullWV(pts)=toy_pull->GetFunction("gaus")->GetParameter(2);
        massErrV(pts)=0.;
        meanErrV(pts)=toy_bias->GetFunction("gaus")->GetParameter(2)/sqrt(1);
        biasErrV(pts)=toy_bias->GetFunction("gaus")->GetParameter(2)/sqrt(1);
        pullErrV(pts)=toy_pull->GetFunction("gaus")->GetParError(1)/sqrt(1);
        pullWErrV(pts)=toy_pull->GetFunction("gaus")->GetParError(2)/sqrt(1);

        cout << "Template mass "<< mass<< " - Mass: " << meanV[pts] <<" +/- "<<meanErrV[pts];
        cout << "\t  Unc: " << toy_error->GetMean();
        cout << "\t  Pull: " << pullV[pts] <<" / "<<pullWV[pts]<<endl;
    }

    // Mass plot & fit
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TGraphErrors *massGraph = new TGraphErrors(massV, meanV, massErrV, meanErrV);
    massGraph->SetName("massGraph");
    TF1* meanFit = new TF1("meanFit","pol1",minTMass-10, maxTMass+10);
    meanFit->SetLineStyle(1);
    massGraph->Fit("meanFit");
    TF1* f2 = new TF1("f2","pol1",minTMass-10, maxTMass+10);
    f2->SetParameter(0,0.);
    f2->SetParameter(1,1.);
    f2->SetLineColor(4);
    massGraph->GetXaxis()->SetNdivisions(50205);
    massGraph->GetXaxis()->SetTitle("Generated mass [GeV]");
    massGraph->GetYaxis()->SetTitleOffset(1.22);
    massGraph->GetYaxis()->SetTitle("Estimated mass [GeV]");
    massGraph->Draw("apz");
    f2->Draw("same");
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_mass_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_mass_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_mass_%s.png", tag);
    c_min->Print(hname);

    // Pull plot & fit
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TGraphErrors *pullGraph = new TGraphErrors(massV, pullV, massErrV, pullErrV);
    pullGraph->SetName("pullGraph");
    TF1* pullFit = new TF1("pullFit","pol0",minTMass-10, maxTMass+10);
    pullGraph->Fit("pullFit");
    TF1* f3 = new TF1("f3","pol1",minTMass-10, maxTMass+10);
    f3->SetParameter(0,0.);
    f3->SetParameter(1,0.);
    f3->SetLineColor(4);
    pullGraph->SetMinimum(-1);
    pullGraph->SetMaximum(+1);
    pullGraph->GetXaxis()->SetNdivisions(50205);
    pullGraph->GetXaxis()->SetTitle("Generated mass [GeV]");
    pullGraph->GetYaxis()->SetTitleOffset(1.22);
    pullGraph->GetYaxis()->SetTitle("pull mean [GeV]");
    pullGraph->Draw("apz");
    f3->Draw("same");
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_pull_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pull_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pull_%s.png", tag);
    c_min->Print(hname);

    TGraphErrors *pullWGraph = new TGraphErrors(massV, pullWV, massErrV, pullWErrV);
    pullWGraph->SetName("pullWGraph");
    TF1* f4 = new TF1("f3","pol1",minTMass-10, maxTMass+10);
    f4->SetParameter(0,1.);
    f4->SetParameter(1,0.);
    f4->SetLineColor(4);
    f4->SetLineStyle(2);
    TF1* pullWFit = new TF1("pullWFit","pol0",minTMass-10, maxTMass+10);
    pullWGraph->GetXaxis()->SetNdivisions(50205);
    pullWGraph->GetXaxis()->SetTitle("Generated width [GeV]");
    pullWGraph->GetYaxis()->SetTitleOffset(1.22);
    pullWGraph->GetYaxis()->SetTitle("Pull width");
    pullWGraph->SetMinimum(0.75);
    pullWGraph->SetMaximum(1.25);
    pullWGraph->Draw("apz");
    f4->Draw("same");
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_pullW_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pullW_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_pullW_%s.png", tag);
    c_min->Print(hname);

    // bias plot & fit
    TGraphErrors *biasGraph = new TGraphErrors(massV, biasV, massErrV, biasErrV);
    biasGraph->SetName("biasGraph");
    TF1* biasFit = new TF1("biasFit","pol1",minTMass-10, maxTMass+10);
    biasGraph->Fit("biasFit");
    f3->SetLineStyle(2);

    biasGraph->GetXaxis()->SetNdivisions(50205);
    biasGraph->GetXaxis()->SetTitle("Generated width [GeV]");
    biasGraph->GetYaxis()->SetTitleOffset(1.22);
    biasGraph->GetYaxis()->SetTitle("Bias [GeV]");
    biasGraph->SetMinimum(-0.5);
    biasGraph->SetMaximum(+0.5);
    biasGraph->Draw("apz");
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"cal_bias_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_bias_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"cal_bias_%s.png", tag);
    c_min->Print(hname);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(kTRUE);

    cout << "Properties at 1.5:\n";
    TH1F*   toy_err  = (TH1F*) gDirectory->Get("errMass_1.50") ;
    cout << "Toy err is "<<toy_err<<endl;
    toy_err->GetXaxis()->SetTitle("Uncertainty [GeV]");
    toy_err->GetXaxis()->SetNdivisions(505);
    toy_err->GetYaxis()->SetTitleOffset(1.4);
    toy_err->Draw();
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"errMass_1.5_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"errMass_1.5_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"errMass_1.5_%s.png", tag);
    c_min->Print(hname);
    cout << "Mean uncertainty "<<toy_err->GetMean()<<endl;

    gStyle->SetOptFit(1111);
    TH1F*   toy_mean  = (TH1F*) gDirectory->Get("meanMass_1.50") ;
    toy_mean->Fit("gaus","Q");
    toy_mean->GetXaxis()->SetTitle("Mass [GeV]");
    toy_mean->GetYaxis()->SetTitleOffset(1.2);
    toy_mean->GetYaxis()->SetTitle("Events/bin");

    toy_mean->Draw();
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"meanMass_1.5_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"meanMass_1.5_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"meanMass_1.5_%s.png", tag);
    c_min->Print(hname);

    gStyle->SetOptFit(1111);
    TH1F*   toy_pull  = (TH1F*) gDirectory->Get("pullMass_1.50") ;
    toy_pull->Fit("gaus","Q");
    toy_pull->GetXaxis()->SetTitle("Pull");
    toy_pull->GetYaxis()->SetTitle("Events/bin");
    toy_pull->GetYaxis()->SetTitleOffset(1.4);
    toy_pull->Draw();
    CMS_lumi( c_min, iPeriod, 0 );
    sprintf(hname,"pullMass_1.5_%s.pdf", tag);
    c_min->Print(hname);
    sprintf(hname,"pullMass_1.5_%s.C", tag);
    c_min->Print(hname);
    sprintf(hname,"pullMass_1.5_%s.png", tag);
    c_min->Print(hname);

    cout << "Inversion\n";
    float a = meanFit->GetParameter(1);
    float b = meanFit->GetParameter(0);
    float ae = meanFit->GetParError(1);
    float be = meanFit->GetParError(0);
    cout <<a<<" "<<b<<endl;
    cout << "Template mass "<< fitMass << " - Fit: " << a <<" / "<<b<<endl;
    cout << "calib mass: " << (fitMass-b)/a<<endl;
    calMass = (fitMass-b)/a;
    float err = sqrt(ae*ae*((fitMass-b)/(a*a))*((fitMass-b)/(a*a)) + be*be/(a*a));
    cout << "calib err: " << err<<endl;

    cout << "Pull at 1.5: "<< pullWFit->Eval(1.5)<<endl;
    cout << "calbrated stat unc: "<< pullWFit->Eval(1.5)*fitUnc<<endl;

    theFile->Close();
}

void MassFit::defaultMCsignal()
{
    MCFileName = "./2012_combined_EACMLB.root";
    SystFileLocation = "./2012_combined_EACMLB.root";
}

void MassFit::setMCsignal(string defaultFile, string genFile)
{
    MCFileName = defaultFile;
    SystFileLocation = genFile;
}

//select type of mc, set xsec
void MassFit::defaultMCbackground()
{
    // Background templates
    mcBackgroundLabels.push_back("SingleTop");
    mcBackgroundLabels.push_back("WJets");
    mcBackgroundLabels.push_back("DrellYan");
    mcBackgroundLabels.push_back("Diboson");
    mcBackgroundLabels.push_back("QCD");
}

void MassFit::setup()
{
    double scale;
    minTMass = 151;maxTMass = 199; //intervalTMass = 3;

    setTDRStyle();

    int W = 800;
    int H = 600;

    // Simple example of macro: plot with CMS name and lumi text
    //  (this script does not pretend to work in all configurations)
    // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
    // For instance: 
    //               iPeriod = 3 means: 7 TeV + 8 TeV
    //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
    // Initiated by: Gautier Hamel de Monchenault (Saclay)
    int H_ref = 600; 
    int W_ref = 800; 

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref; 
    float L = 0.16*W_ref;
    float R = 0.04*W_ref;

    // setup a new canvas
    c_min = new TCanvas("c_min","", 600, 600);
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
    nFitFailed = 0;
    nFitTried = 0;
    hFitFailed = new TH1F("hFitFailed","hFitFailed", 11, 0, 9);
    hFitFailedDiff = new TH1F("hFitFailedDiff","hFitFailedDiff", 60, -30., 30.);
    fittedTemplate = -1;

    // initialize settings
    drawSwitch = false;
    fixBckg = 0;
    bkgsyst = false;
    quietFit_ = 0; //Default print level of the fit change with quiet()

    signalXs = 157.5;        //TTbar

    char dataHistoName[20] = "PeakMassTree_";
    char histoName[20] = "mlbwa_";
    char evHistoName[15] = "h_nEvents";

    topMass = new RooRealVar("mass", "Reconstructed top mass", lowerMassCut, upperMassCut);
    sample = new RooCategory("sample","sample") ;
    for (int itype = 0;itype!=maxType;++itype) sample->defineType(type[itype]); //set types of sample

    char name[300], hname[150], tag[50], sname[150];

    TFile* theFile;
    TH1::AddDirectory(kFALSE);

    // open data file
    theFile = new TFile (DataFileLocation.c_str());
    for (int ihisto =0 ; ihisto < maxType; ++ihisto) {
        // get the data histograms
        sprintf(sname, "%s", type[ihisto]);
        sprintf(hname, "%s_Data_%s", histoName, type[ihisto]);
        cout << hname << endl;
        datasets[sname] = (TH1F*) gDirectory->Get(hname);
        if (datasets[sname]==0) assert(false);
        cout << "Got dataset " << type[ihisto] << " " << datasets[sname]<<endl;
    }
    // cleanup
    theFile->Close();
    delete theFile;

    assembleDatasets();
    for (int i = 0; i!=maxType; ++i) {
        if (fixedSample) nTotSample[i] = datasets[type[i]]->GetEntries();
        cout << "in "<< type[i] << " " <<nTotSample[i]<<endl;
    }
    toyDataHisto = new TH1F("toyDataHisto", "toyDataHisto", 40, 0., 400);
    toyDataHisto->Reset();
    data->fillHistogram(toyDataHisto, *topMass);
    dataHisto = (TH1F*) toyDataHisto->Clone("dataHisto");

    // Get signal tt mass templates
    templ_mean   = new TH1F("tmean"  ,"Mean of the templates"  ,50, 0., 9.);
    templ_rms    = new TH1F("trms"   ,"RMS of the templates"   ,50, 0., 9.);

    // define template masses
    mcSignalTemplMass.push_back(1.50);
    mcSignalTemplMass.push_back(3.00);
    mcSignalTemplMass.push_back(4.50);
    mcSignalTemplMass.push_back(6.00);
    mcSignalTemplMass.push_back(7.50);

    minTemplate = 0;
    maxTemplate = mcSignalTemplMass.size();
    nominalTemplate = 0;
    toyTemplate = nominalTemplate;

    scale = 1;

    // open data file
    theFile = new TFile (DataFileLocation.c_str());
    // loop through the interaction types we want to check
    for (int itype =0; itype < maxType; ++itype) {
        // loop through the top masses we want to check
        for (unsigned int imass = 0; imass <  mcSignalTemplMass.size(); ++imass) {
            // get the histogram for the interaction and mass we want
            sprintf(hname, "mlbwa__TTbar_%.2f_%s",mcSignalTemplMass[imass], type[itype]);
            sprintf(tag, "%s%.2f", type[itype], mcSignalTemplMass[imass]);
            if (debug) cout << "Signal template " << hname << " "<< tag << endl;
            if (debug) printf("itype %d imass %d \n", itype, imass);
            TH1F* histo = (TH1F*) gDirectory->Get(hname);
            if (histo==0) { 
                cout << "Histo does not exist\n";
                histo = new TH1F(hname,hname,100,0,200); 
            }
            histo->Rebin(5);
            histo->SetLineColor(4);
            mcSignalTemplHistosScaled[tag] = histo;
        }

    }
    // cleanup
    theFile->Close();
    delete theFile;

    // plot the systematics GEN templates if they are available
    if (systematics) {
        cout << "Systematics Signal GEN template from "  <<SystFileLocation<<endl;
        theFile = new TFile (SystFileLocation.c_str()) ;

        unsigned int maxTemplates;
        if (systematicsPDF) maxTemplates = 41;
        else maxTemplates = mcSignalTemplMass.size();
        // loop through the top masses
        for (unsigned int i = 0; i!= maxTemplates; ++i){

            float mass = 0.;
            if (!systematicsPDF) mass = mcSignalTemplMass[i];

            scale = 1.;
            //loop through the interaction types
            for (int itype = 0;itype!=maxType;++itype) {
                // get the histogram name
                if (systematicsPDF) {
                    sprintf(tag,"%s_pdf%i",type[itype],i);
                    sprintf(hname,"%s_%s_pdf%i", histoName,type[itype],i);
                } else {
                    sprintf(tag,"%s%.2f",type[itype],mass);
                    sprintf(hname, "mlbwa__TTbar_%.2f_%s",mass, type[itype]);
                }

                // get the histogram and store in mcSignalHistos_gen
                if (debug) cout << "Signal GEN template "  << type[itype] << " " <<i<<" "<< tag <<" "<< hname <<endl;
                TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
                histo->Rebin(5);
                histo->SetLineColor(4);
                mcSignalTemplHistosScaled_gen[tag] = histo;
            }
        }
        // cleanup
        theFile->Close();
        delete theFile;
    }

    templateNr = mcSignalTemplHistosScaled.size();
    cout << "Have "<<templateNr << " signal template histos\n";

    // Get the background templates
    cout << "Will now retrieve the background templates\n";

    // open the histogram file
    theFile = new TFile (DataFileLocation.c_str());
    // iterate through types
    for (int itype = 0; itype < maxType; ++itype) {
        //iterate through background MC types
        for (unsigned int bkgType = 0;bkgType!=mcBackgroundLabels.size();++bkgType) {
            // format the name of the histo we want to get
            sprintf(hname, "mlbwa__%s_%s", mcBackgroundLabels[bkgType].Data(), type[itype]);
            if (debug) cout << "hname " << hname<<endl;
            TH1F* histo  = (TH1F*) gDirectory->Get(hname);
            if (histo==0) { 
                cout << "Histo does not exist\n";
                histo = new TH1F(hname,hname,100,0,200); 
            }

            // rebin the histogram and add to the total backgrounds histogram
            histo->Rebin(5);
            histo->SetLineColor(2);
            mcBackgroundHistosScaled[tag] = histo;
            if (mcTotalBackgroundHistoScaled.find(type[itype]) == mcTotalBackgroundHistoScaled.end()) {
                mcTotalBackgroundHistoScaled[type[itype]] = (TH1F*) histo->Clone("mcTotalBackgroundHistoScaled");
            } else {
                mcTotalBackgroundHistoScaled[type[itype]]->Add(histo);
            }
        }
        cout << "Total mean: " << mcTotalBackgroundHistoScaled[type[itype]]->GetMean()
             << " entries: " << mcTotalBackgroundHistoScaled[type[itype]]->Integral()
             << endl;
    }

    // clean up a bit, report to console
    theFile->Close();
    delete theFile;
    cout << "Got "<< mcBackgroundHistosScaled.size() << "background histos\n";

    if (systematics) {
        cout << "Systematics Background GEN template from "  <<SystFileLocation<<endl;
        theFile = new TFile (SystFileLocation.c_str());
        // iterate through the background MC types
        for (unsigned int i = 0; i!=mcBackgroundLabels.size();++i) {
            // iterate through types
            for (int itype = 0;itype!=maxType;++itype) {
                //format the name of the histo we want to get
                sprintf(hname, "mlbwa__%s_%s", mcBackgroundLabels[i].Data(), type[itype]);
                if (debug) cout << hname<<endl;
                TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
                if (histo==0) { 
                    cout << "Histo does not exist\n";
                    histo = new TH1F(hname,hname,100,0,200); 
                }

                //rebin the histogram and add to the total backgrounds_GEN histogram
                histo->Rebin(5);
                if (mcTotalBackgroundHistoScaled_gen.find(type[itype]) == mcTotalBackgroundHistoScaled_gen.end()) {
                    mcTotalBackgroundHistoScaled_gen[type[itype]] = (TH1F*) histo->Clone("mcTotalBackgroundHistoScaled_gen");
                } else {
                    mcTotalBackgroundHistoScaled_gen[type[itype]]->Add(histo);
                }
            }
        }

        // clean up
        theFile->Close();
        delete theFile;
    }

    // re-initialize some variables, establish the workspace
    gr = 0; toy_mean = 0; grc=0;
    myFitResults_all = 0;

    w = new RooWorkspace("w","workspace") ;
    w->import(*topMass);
    w->import(*sample);

    // define the datasets initially
    // loop through the interaction types
    for (int itype = 0;itype!=maxType;++itype) {

#ifdef GAUSS
        sprintf(hname,"Gaussian::background%s(mass,140,20)",type[itype]);
        w->factory(hname);
#else
        sprintf(hname,"histo_bck%s",type[itype]);
        if (debug) cout << hname<<endl;
        w->import(*(new RooDataHist(hname, hname, *topMass, mcTotalBackgroundHistoScaled[type[itype]])));
        sprintf(hname,"HistPdf::background%s(mass,histo_bck%s)", type[itype], type[itype]);
        if (debug) cout << hname<<endl;
        w->factory(hname);

        if (systematics) {
            sprintf(hname,"histo_bck_gen%s",type[itype]);
            w->import(*(new RooDataHist(hname, hname, *topMass, mcTotalBackgroundHistoScaled_gen[type[itype]])));
            sprintf(hname,"HistPdf::background_gen%s(mass,histo_bck%s)", type[itype], type[itype]);
            w->factory(hname);
        }
#endif
        //
        float totalBackground = (float) mcTotalBackgroundHistoScaled[type[itype]]->Integral();
        sprintf(hname,"nbrBackgroundEvents%s",type[itype]);
        if (debug) cout << hname<<endl;
        nbrBackgroundEvents[type[itype]] = new RooRealVar(hname,hname,totalBackground,0,10000);
        w->import(*nbrBackgroundEvents[type[itype]]);

        // loop through the top masses
        for (unsigned int i = 0; i!=mcSignalTemplMass.size();++i) {
            float mass = mcSignalTemplMass[i];
            if (debug) cout << "Build "  << type[itype] << " pdf for " <<mass<<endl;
            sprintf(tag,"%s%.2f",type[itype],mass);

#ifdef GAUSS
            sprintf(hname,"Gaussian::signal%s%.2f(mass,%.2f,20)",type[itype],mass,mass);
            w->factory(hname);
#else
            sprintf(tag,"%s%.2f",type[itype],mass);
            sprintf(hname,"histo_sgn%s%.2f",type[itype],mass);
            if (debug) cout << hname<< " " << tag<<endl;
            w->import( *(new RooDataHist(hname, hname, *topMass, mcSignalTemplHistosScaled[tag]) ));
            sprintf(hname,"HistPdf::signal%s%.2f(mass,histo_sgn%s%.2f)",type[itype],mass,type[itype],mass);
            if (debug) cout << hname<<endl;
            w->factory(hname);
            if (systematics && !systematicsPDF) {
                sprintf(hname,"histo_sgn_gen%s%.2f",type[itype],mass);
                w->import( *(new RooDataHist(hname, hname, *topMass, mcSignalTemplHistosScaled_gen[tag]) ));
                sprintf(hname,"HistPdf::signal_gen%s%.2f(mass,histo_sgn_gen%s%.2f)",type[itype],mass,type[itype],mass);
                w->factory(hname);
            }

#endif
            if (useRatio) {
                sprintf(hname,"SUM::model%s%.2f( ratio%s%.2f[0,1]*signal%s%.2f, background%s )",
                        type[itype], mass, type[itype], mass, type[itype], mass, type[itype]);
            } else {
                sprintf(hname,"SUM::model%s%.2f( Nsig%s%.2f[0,1000]*signal%s%.2f, nbrBackgroundEvents%s*background%s )",
                        type[itype], mass, type[itype], mass, type[itype], mass, type[itype], type[itype]);
            }
            if (debug) cout << hname<<endl;
            w->factory(hname);
        }

        if (systematics && systematicsPDF) {
            for (int i = 0; i!= 41; ++i){
                sprintf(hname,"histo_sgn_gen%s%i",type[itype],i);
                sprintf(tag,"%s_pdf%i",type[itype],i);
                w->import( *(new RooDataHist(hname, hname, *topMass, mcSignalTemplHistosScaled_gen[tag]) ));
                sprintf(hname,"HistPdf::signal_gen%s%i(mass,histo_sgn_gen%s%i)",type[itype],i,type[itype],i);
                cout << hname<<endl;
                w->factory(hname);
            }
        }

    }

    //This is only for the pdf systematics templates:
    char name2[300];
    for (unsigned int i = 0; i!=mcSignalTemplMass.size();++i) {
        float mass = mcSignalTemplMass[i];
        if (debug) cout << "Build simultaneous pdf for " <<mass<<endl;
        sprintf(name2,"SIMUL::model%.2f(sample",mass);
        for (int j =0; j<maxType;++j) {
            sprintf(name,"%s,%s=model%s%.2f", name2, type[j], type[j], mass);
            sprintf(name2,"%s", name);
        }
        sprintf(name,"%s)",name2);
        if (debug) cout << name<<endl;

        w->factory(name);
    }

    if (bkgsyst) {
        bkgmean  = new RooRealVar("bkgmean","bkg mass",180.) ;
        bkgwidth = new RooRealVar("bkgwidth","bkg width",20.) ;
        histo_bck_pdfff = new RooGaussian("background","background PDF",*topMass,*bkgmean,*bkgwidth) ;
    }

    toy_error=0;
    printMassRange();

    toy_LL     = new TH2F("LL"  ,"LL residuals",9, -0.5, 8.5,200,-100,100);

}

void MassFit::assembleDatasets()
{
    char sname[150], mname[150];

    // create a dictionary of TH1s to operate on
    if (debug) cout << "Combine "<< maxType<<" categories\n";
    map<string,TH1*> mapToImport;

    // iterate through all the types
    for (int itype = 0; itype < maxType; ++itype) {

        sprintf(sname, "%s", type[itype]);
        if (debug) cout << type[itype] << "\n";
        if (debug) printf("891 datasets %f \n", datasets[sname]->GetMean());
        mapToImport[sname] = datasets[sname];
    }

    // create dataset histogram and print it out
    data = new RooDataHist("data","combined data", *topMass, Index(*sample), Import(mapToImport)) ;
    data->Print();
}

TH1F* MassFit::templateHisto(const char* typeC, int i)
{
    char tag[50];
    sprintf(tag,"%s%.2f",typeC,mcSignalTemplMass[i]);
    cout << "Template mass "<< tag<<endl;
    return mcSignalTemplHistosScaled[tag];
}

TH1F* MassFit::backgroundHisto(const char* typeC)
{
    return mcTotalBackgroundHistoScaled[typeC];
}

double MassFit::fitPoint(int i)
{
    cout << "Fit with template mass "<< mcSignalTemplMass[i]<<endl;
    fittedTemplate = i;
    if (myFitResults_all)  {
        delete myFitResults_all;
        if (bkgsyst) {
            delete pdffit;
            delete extBackgroundPdf;
        }
    }

    char hname[50];
    sprintf(hname,"model%.2f", mcSignalTemplMass[i]);
    pdffit = w->pdf(hname) ;
    pdffit->Print();

    if (useRatio && fixBckg==3) {
        char tag[50];
        for (int itype = 0;itype!=maxType;++itype) {
            sprintf(tag,"%s%.2f",type[itype],mcSignalTemplMass[i]);
            sprintf(hname,"ratio%s%.2f", type[itype], mcSignalTemplMass[i]);
            w->var(hname)->setVal(mcSignalTemplHistosScaled[tag]->Integral()/
                    (mcSignalTemplHistosScaled[tag]->Integral()+mcTotalBackgroundHistoScaled[type[itype]]->Integral()));
            w->var(hname)->setConstant(1);
        }
    } else if (useRatio && fixBckg==1) {
        for (int itype = 0;itype!=maxType;++itype) {
            sprintf(hname,"ratio%s%.2f", type[itype], mcSignalTemplMass[i]);
            w->var(hname)->setVal(1.0);
            w->var(hname)->setConstant(1);
        }
    }

    myFitResults_all = pdffit->fitTo(*data, Save(), PrintLevel(quietFit_)) ;
    double fit_chi2 = myFitResults_all->minNll();
    cout << "LL " << fit_chi2 <<endl;
    return fit_chi2;
}


void MassFit::fixBackground(float e, float ee, float em, float mm, float m)
{
    if (maxType!=3){ cout << "Needs to be modified - Exiting\n";return;}

    cout << "Fix background to " << e << "for E," << ee <<" for EE, "
         << em <<" for EM, "<< mm <<" for MM " << m << "for M\n";

    w->var("nbrBackgroundEventsE")->setVal(e);
    w->var("nbrBackgroundEventsEE")->setVal(ee);
    w->var("nbrBackgroundEventsEM")->setVal(em);
    w->var("nbrBackgroundEventsMM")->setVal(mm);
    w->var("nbrBackgroundEventsM")->setVal(m);

    w->var("nbrBackgroundEventsE")->setConstant(kTRUE);
    w->var("nbrBackgroundEventsEE")->setConstant(kTRUE);
    w->var("nbrBackgroundEventsEM")->setConstant(kTRUE);
    w->var("nbrBackgroundEventsMM")->setConstant(kTRUE);
    w->var("nbrBackgroundEventsM")->setConstant(kTRUE);
    fixBckg = 2;
}

void MassFit::fixBackground(int fix)
{
    char hname[50];
    switch ( fix )
    {
        case 0:
            for (int itype = 0;itype!=maxType;++itype) {
                sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
                w->var(hname)->setConstant(0);
            }
            break;
        case 1:
            for (int itype = 0;itype!=maxType;++itype) {
                sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
                w->var(hname)->setConstant(1);
                w->var(hname)->setVal(0.);
            }
            break;
        case 2: 
            cout <<"Use method fixBackground(e, ee, em, mm, m) - Nothing done now\n";
            return;
            break;
        case 3:
            if (!useRatio) {
                for (int itype = 0;itype!=maxType;++itype) {
                    sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
                    w->var(hname)->setVal(mcTotalBackgroundHistoScaled[type[itype]]->Integral());
                    w->var(hname)->setConstant(1);
                }
            }
            break;
        default:
            cout << "input " << fix << " unknown!\n";
    }
    fixBckg = fix;
    fixBackground();
}

int MassFit::fixBackground()
{
    cout << "Choice can be: \n";
    cout << "   0: Not fixed\n";
    cout << "   1: Fixed to 0.\n";
    cout << "   2: Fixed to external input\n";
    cout << "   3: Fixed to MC\n";
    cout << "Present setting: " << fixBckg<<endl;
    if (!useRatio && fixBckg>=2) {
        char hname[50];
        cout << "   Fixed to: ";
        for (int itype = 0;itype!=maxType;++itype) {
            sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
            cout << w->var(hname)->getVal() <<" for "<<type[itype]<<", ";
        }
        cout <<endl;
    }
    return fixBckg;
}

int MassFit::fitAll() {
    chiSquared.clear();
    chiSquared.resize(maxTemplate+1);
    for (int i = minTemplate; i!=maxTemplate;++i) {
        chiSquared[i] = fitPoint(i);
        if (std::isinf(chiSquared[i]) || std::isnan(chiSquared[i])) return 1;
    }
    return 0;
}

pair<double,double> MassFit::findMin(bool all, int pointsToUse) {
    int pts=0;
    double minLL = 999999999;
    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i]>=0) {
            ++pts;
            if (chiSquared[i] < minLL) minLL = chiSquared[i];
        }
    }

    TVectorD x(pts), ex(pts);
    TVectorD y(pts), ey(pts);
    pts=0;
    int minPt=0;
    int failures=0;
    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i] >= 0){
            x[pts]=mcSignalTemplMass[i];
            y[pts]=chiSquared[i] - minLL;
            ey[pts]=9.4;
            cout << x[pts] << " "<<y[pts]<< " "<<y[minPt]<<endl;
            if (y[pts]<=y[minPt]) minPt=pts;
            ++pts;
        } else ++failures;
    }

    int points, min, max;
    if (all ) { points = pts; min=0;max=pts-1;}
    else {
        if (pointsToUse < 2) pointsToUse = 2;
        if (pointsToUse > pts) pointsToUse = 4;
        points = pointsToUse*2+1;
        if (minPt<pointsToUse) minPt=pointsToUse;
        if (minPt>maxTemplate-pointsToUse-1) minPt=maxTemplate-pointsToUse-1;
        min= minPt-pointsToUse;max=minPt+pointsToUse;
    }

    if(gr!=0) delete gr;

    gr = new TGraphErrors(x,y,ex,ey);
    gr->SetName("gr");
    gr->Fit("pol2","Q","",mcSignalTemplMass[min], mcSignalTemplMass[max]);
    gr->GetXaxis()->SetTitle("mlbwa  m(lb) [GeV]");
    gr->GetYaxis()->SetTitle("-log (L/L_{max})");
    gr->GetYaxis()->SetTitleOffset(1.25);

    Double_t a = gr->GetFunction("pol2")->GetParameter(2);
    fitMass = gr->GetFunction("pol2")->GetMinimumX(0,9);
    cout << "Min mass "<< gr->GetFunction("pol2")->GetMinimumX(0,9)<<" +/- "<< 1/sqrt(2*a)<<endl;
    if (failures) cout << "Fits failed: " <<failures <<endl;

    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i] >= 0){
            toy_LL->Fill(i,gr->GetFunction("pol2")->Eval(mcSignalTemplMass[i])-(chiSquared[i] - minLL));
            cout <<i<< " "<<chiSquared[i] << " "<<(chiSquared[i] - minLL)<< " "<<gr->GetFunction("pol2")->Eval(mcSignalTemplMass[i])-(chiSquared[i] - minLL)<<endl;
        }
    }

    if (a<0.) return pair<double,double>(0.,0.);

    return pair<double,double>(gr->GetFunction("pol2")->GetMinimumX(0,9),1/sqrt(2*a));
}

pair<double,double> MassFit::findMinFake(bool all, int pointsToUse) {
    int pts=0;
    double minLL = 999999999;
    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i]>=0) {
            ++pts;
            if (chiSquared[i] < minLL) minLL = chiSquared[i];
        }
    }

    TVectorD x(pts), ex(pts);
    TVectorD y(pts), ey(pts);
    pts=0;
    int minPt=0;
    int failures=0;
    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i] >= 0){
            x[pts]=mcSignalTemplMass[i];
            y[pts]=chiSquared[i] - minLL;
            ey[pts]=9.4;
            cout << x[pts] << " "<<y[pts]<< " "<<y[minPt]<<endl;
            if (y[pts]<=y[minPt]) minPt=pts;
            ++pts;
        } else ++failures;
    }

    int points, min, max;
    if (all ) { points = pts; min=0;max=pts-1;}
    else {
        if (pointsToUse < 2) pointsToUse = 2;
        if (pointsToUse > pts) pointsToUse = 4;
        points = pointsToUse*2+1;
        if (minPt<pointsToUse) minPt=pointsToUse;
        if (minPt>maxTemplate-pointsToUse-1) minPt=maxTemplate-pointsToUse-1;
        min= minPt-pointsToUse;max=minPt+pointsToUse;
    }

    if (gr!=0) {
        delete gr;
    }
    gr = new TGraphErrors(x,y,ex,ey);
    gr->SetName("gr");
    gr->Fit("pol2","Q","",mcSignalTemplMass[min], mcSignalTemplMass[max]);
    gr->GetXaxis()->SetTitle("MLBWA width [GeV]");
    gr->GetYaxis()->SetTitle("-log (L/L_{max})");
    gr->GetYaxis()->SetTitleOffset(1.25);

    Double_t a = gr->GetFunction("pol2")->GetParameter(2);
    cout << "Min mass "<< gr->GetFunction("pol2")->GetMinimumX(0,9)<<" +/- "<< 1/sqrt(2*a)<<endl;
    if (failures) cout << "Fits failed: " <<failures <<endl;

    for (int i = minTemplate; i!=maxTemplate;++i) {
        if (chiSquared[i] >= 0){
            toy_LL->Fill(i,gr->GetFunction("pol2")->Eval(mcSignalTemplMass[i])-(chiSquared[i] - minLL));
            cout <<i<< " "<<chiSquared[i] << " "<<(chiSquared[i] - minLL)<< " "<<gr->GetFunction("pol2")->Eval(mcSignalTemplMass[i])-(chiSquared[i] - minLL)<<endl;
        }
    }

    if (a<0.) return pair<double,double>(0.,0.);

    return pair<double,double>(gr->GetFunction("pol2")->GetMinimumX(0,9),1/sqrt(2*a));
}

void MassFit::printLL(bool toFile, char* name)
{

    writeExtraText = true;                // if extra text
    extraText  = "Preliminary";            // default extra text is "Preliminary"
    //lumi_8TeV  = "dilepton, 19.7 fb^{-1}"; // default is "19.7 fb^{-1}"
    //lumi_7TeV  = "4.9 fb^{-1}";            // default is "5.1 fb^{-1}"
    int iPeriod = 0;                       // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 

    c_min->cd();
    gr->Draw("apz");

    CMS_lumi( c_min, iPeriod, 0 );

    //cmsprelim();
    if (toFile) {
        char name2[200],name3[200];
        if (strlen(name)==0) {
            sprintf(name2,"%s/fit_result/ll", absolutePath);
            typeTag(name2);
        } else {
            sprintf(name2,"%s", name);
        }
        sprintf(name3,"%s.C",name2); c_min->Print(name3);
        sprintf(name3,"%s.png",name2); c_min->Print(name3);
        sprintf(name3,"%s.pdf",name2); c_min->Print(name3);
    }
    cout << "Chi2 : "<< gr->GetFunction("pol2")->GetChisquare()<<endl; // obtain chi^2
    cout << "NDOF : " <<gr->GetFunction("pol2")->GetNDF()<<endl;       // obtain ndf
    cout << "P(c2): "<< TMath::Prob(gr->GetFunction("pol2")->GetChisquare(),gr->GetFunction("pol2")->GetNDF())<<endl;
}

void MassFit::do_toys(int n_exp, int templateToUse)
{
    quietFit(true);
    if (!systematics || !systematicsPDF) cout << "Template mass "<< mcSignalTemplMass[templateToUse]<<endl;
    else cout << "PDF "<< templateToUse<<endl;
    if (toy_error!=0) {
        delete toy_mean; delete toy_error; delete toy_pull; delete toy_bias;
    }
    if (toy_LL!=0) delete toy_LL;
    float massPoint = 1.5;
    if (!systematics || !systematicsPDF) massPoint = mcSignalTemplMass[templateToUse];

    toy_mean   = new TH1F("mean"  ,"Top mass",100, massPoint-1.5, massPoint+7.5);
    toy_bias   = new TH1F("bias"  ,"Top mass bias",100, -3.5, 3.5);
    toy_error  = new TH1F("error" ,"top mass uncertainty",500, 0.1, 0.4);
    toy_pull   = new TH1F("pull"  ,"pull",200, -10, 10);
    toy_LL     = new TH2F("LL"  ,"LL residuals",9, -0.5, 8.5,200,-100,100);

    toy_mean->GetXaxis()->SetNdivisions(50205);
    toy_bias->GetXaxis()->SetNdivisions(50205);
    toy_mean->SetFillColor(44);
    toy_bias->SetFillColor(44);
    toy_error->SetFillColor(44);
    toy_pull->SetFillColor(44);
    int j, stat;
    for (int i=0;i<n_exp;i++){
        int countFailures = 0;
        do {
            do { j = generate_toy(templateToUse);} while (j==0);

            stat = fitAll();
            if (stat !=0) {
                ++countFailures;
                cout << "RooFit failure "<<countFailures<<endl;
            }
        } while ((stat !=0) && countFailures<10);
        if (stat !=0) {
            cout << "Too many consecutive failures\n";
            exit(1);
        }
        nFitFailed+=countFailures;

        pair<double,double> result = findMinFake();

        if (result.second>0.) {
            cout << result.first<<" "<<result.second<<endl;
            toy_mean->Fill(result.first);
            toy_bias->Fill(result.first-massPoint);
            toy_error->Fill(result.second);
            toy_pull->Fill((massPoint-result.first)/result.second);
        } else {
            ++nFitFailed;
        }
        ++nFitTried;
    }
    cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
}

void MassFit::texEvents(int templateToUse, float lower, float upper)
{
    int binLow, binHigh;
    char tag[50];
    cout <<"\\begin{table}[htp]\n";
    cout <<"\\begin{center}\n";

    cout <<"\\caption{Number of events, Mass range: "<<lower<<" - "<<upper
        << " ; weight cut="<<lowerWeightCut<<"}\n";

    cout <<"\\begin{tabular}{|l|";
    for (int itype = 0;itype!=maxType;++itype) cout <<"c|";
    cout <<"}\\hline\\hline\n";


    cout << "Sample ";
    for (int itype = 0;itype!=maxType;++itype) {
        cout << "\t& "<<type[itype];
    }
    cout <<"\\\\\n\\hline\n";

    cout << "Signal";
    for (int itype = 0;itype!=maxType;++itype) {
        sprintf(tag,"%s%.2f",type[itype],mcSignalTemplMass[templateToUse]);
        binLow = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(lower);
        binHigh = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(upper);
        printf ("\t& %7.1f", mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh));
    }
    cout <<"\\\\\n";

    for (unsigned int i = 0; i!=mcBackgroundLabels.size();++i) {
        cout << mcBackgroundLabels[i];
        for (int itype = 0;itype!=maxType;++itype) {
            sprintf(tag,"%s%s",mcBackgroundLabels[i].Data(),type[itype]);
            printf ("\t& %7.1f", mcBackgroundHistosScaled[tag]->Integral(binLow, binHigh));
        }
        cout <<"\\\\\n";
    }

    cout <<"\\hline\n";
    cout << "Total MC";
    for (int itype = 0;itype!=maxType;++itype) {
        sprintf(tag,"%s%.2f",type[itype],mcSignalTemplMass[templateToUse]);
        printf ("\t& %7.1f", mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh) +
                mcTotalBackgroundHistoScaled[type[itype]]->Integral(binLow, binHigh));
    }
    cout <<"\\\\\n\\hline\n";
    cout << "Data";
    for (int itype = 0;itype!=maxType;++itype) {
        cout << "\t& "
            << datasets[type[itype]]->GetEntries();
    }
    cout <<"\\\\\n\\hline\\hline\n";

    cout <<"\\end{tabular}\n";
    cout <<"\\end{center}\n";
    cout <<"\\end{table}\n";

}

int MassFit::generate_toy(int templateToUse)
{

    if (toyDataHisto!=0) {
        delete toyDataHisto;
    }
    toyDataHisto = (TH1F*)dataHisto->Clone("toyDataHisto");
    toyDataHisto->Reset();

    char hname[50], tag[50];
    int n;
    totalGeneratedSignal = totalGeneratedBkg = 0;
    for (int itype = 0;itype!=maxType;++itype) {
        delete datasets[type[itype]];
        if (!fixedSample) {
            if (!systematics || !systematicsPDF) sprintf(tag,"%s%.2f",type[itype],mcSignalTemplMass[templateToUse]);
            else sprintf(tag,"%s_pdf%i",type[itype],templateToUse);
            int binLow = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(lowerMassCut);
            int binHigh = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(upperMassCut);

            if (systematics) {
                n = _random.Poisson( (float) mcSignalTemplHistosScaled_gen[tag]->Integral(binLow, binHigh));
            } else {
                n = _random.Poisson( (float) mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh));
            }
            cout << "Generate "<< n << " " << type[itype]<<" signal events of ";
            if (!systematics || !systematicsPDF) cout << "mass "<<mcSignalTemplMass[templateToUse];
            else cout << "pdf "<<templateToUse;
            cout << ", with Poisson mean "<< (float) mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh) << endl;
            generatedSignal[type[itype]] = n; totalGeneratedSignal+=n;

            if (systematics) {
                n = _random.Poisson( (float) mcTotalBackgroundHistoScaled_gen[type[itype]]->Integral(binLow, binHigh));
            } else {
                n = _random.Poisson( (float) mcTotalBackgroundHistoScaled[type[itype]]->Integral(binLow, binHigh));
            }
            generatedBkg[type[itype]] = n; totalGeneratedBkg += n;
            cout << "Generate "<<n << " " << " background events";
            cout << ", with Poisson mean "<< (float) mcTotalBackgroundHistoScaled[type[itype]]->Integral(binLow, binHigh) << endl;
        } else {
            sprintf(tag,"%s%.2f",type[itype],mcSignalTemplMass[4]);
            int binLow = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(lowerMassCut);
            int binHigh = mcSignalTemplHistosScaled[tag]->GetXaxis()->FindBin(upperMassCut);
            cout << tag<<endl;
            double sigProb;
            if (systematics) {
                if (systematicsPDF) sprintf(tag,"%s_pdf%i",type[itype],templateToUse);
                sigProb = mcSignalTemplHistosScaled_gen[tag]->Integral(binLow, binHigh) /
                    (mcSignalTemplHistosScaled_gen[tag]->Integral(binLow, binHigh)+mcTotalBackgroundHistoScaled_gen[type[itype]]->Integral(binLow, binHigh));
            } else {
                sigProb = mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh) /
                    (mcSignalTemplHistosScaled[tag]->Integral(binLow, binHigh)+mcTotalBackgroundHistoScaled[type[itype]]->Integral(binLow, binHigh));
            }
            n = _random.Binomial(nTotSample[itype], sigProb);
            cout << "Generate "<< n << " " << type[itype]<<" signal events of ";
            if (!systematics || !systematicsPDF) cout << "mass "<<mcSignalTemplMass[templateToUse];
            else cout << "pdf "<<templateToUse;
            cout << ", with binomial prob "<< sigProb << endl;
            generatedSignal[type[itype]] = n; totalGeneratedSignal+=n;

            generatedBkg[type[itype]] = nTotSample[itype] - n;
            if (generatedBkg[type[itype]] < 0.){
                cout << "ERROR: Number of background events < 0.\n";
                assert(false);
            }
            totalGeneratedBkg += generatedBkg[type[itype]];
            cout << "Generate "<< generatedBkg[type[itype]] << " " << " background events\n";
        }

        if (systematics) {
            if (!systematicsPDF) sprintf(hname,"signal_gen%s%.2f", type[itype], mcSignalTemplMass[templateToUse]);
            else sprintf(hname,"signal_gen%s%i", type[itype], templateToUse);
        } else {
            sprintf(hname,"signal%s%.2f", type[itype], mcSignalTemplMass[templateToUse]);
        }
        cout << hname<<endl;
        datasets[type[itype]] = w->pdf(hname)->generateBinned(*topMass, generatedSignal[type[itype]])->createHistogram(hname,*topMass);
        cout << hname<<endl;

        if (systematics) {
            sprintf(hname,"background_gen%s", type[itype]);
        } else {
            sprintf(hname,"background%s", type[itype]);
        }
        datasets[type[itype]]->Add(w->pdf(hname)->generateBinned(*topMass, generatedBkg[type[itype]])->createHistogram(hname,*topMass));
    }

    if (data!=0) {
        delete data;
    }
    assembleDatasets();

    data->fillHistogram(toyDataHisto, *topMass);

    return data->numEntries();
}

int MassFit::generate_toy_pdf(int templateToUse)
{
    cout << "b\n";
    if (toyDataHisto!=0) {
        delete toyDataHisto;
    }
    toyDataHisto = (TH1F*)dataHisto->Clone("toyDataHisto");
    toyDataHisto->Reset();
    cout << "b\n";

    char hname[50], tag[50];
    int n;
    totalGeneratedSignal = totalGeneratedBkg = 0;
    for (int itype = 0;itype!=maxType;++itype) {

        delete datasets[type[itype]];

        if (!fixedSample) {

            sprintf(tag,"%s%i",type[itype],templateToUse);

            if (systematics) {
                n = _random.Poisson( (float) mcSignalTemplHistosScaled_gen[tag]->Integral());
            } else {
                n = _random.Poisson( (float) mcSignalTemplHistosScaled[tag]->Integral());
            }
            cout << "Generate "<< n << " " << type[itype]<<" signal events of mass "<<templateToUse;
            cout << ", with Poisson mean "<< (float) mcSignalTemplHistosScaled[tag]->Integral() << endl;
            generatedSignal[type[itype]] = n; totalGeneratedSignal+=n;

            if (systematics) {
                n = _random.Poisson( (float) mcTotalBackgroundHistoScaled_gen[type[itype]]->Integral());
            } else {
                n = _random.Poisson( (float) mcTotalBackgroundHistoScaled[type[itype]]->Integral());
            }
            generatedBkg[type[itype]] = n; totalGeneratedBkg += n;
            cout << "Generate "<<n << " " << " background events";
            cout << ", with Poisson mean "<< (float) mcTotalBackgroundHistoScaled[type[itype]]->Integral() << endl;
        } else {
            sprintf(tag,"%s_pdf%i",type[itype],templateToUse);
            cout << "c\n"<<tag<<endl;

            double sigProb;
            if (systematics) {
                sigProb = mcSignalTemplHistosScaled_gen[tag]->Integral() /
                    (mcSignalTemplHistosScaled_gen[tag]->Integral()+mcTotalBackgroundHistoScaled_gen[type[itype]]->Integral());
            } else {
                sigProb = mcSignalTemplHistosScaled[tag]->Integral() /
                    (mcSignalTemplHistosScaled[tag]->Integral()+mcTotalBackgroundHistoScaled[type[itype]]->Integral());
            }

            n = _random.Binomial(nTotSample[itype], sigProb);
            cout << "Generate "<< n << " " << type[itype]<<" signal events of mass "<<templateToUse;
            cout << ", with binomial prob "<< sigProb << endl;
            generatedSignal[type[itype]] = n; totalGeneratedSignal+=n;

            generatedBkg[type[itype]] = nTotSample[itype] - n;
            totalGeneratedBkg += generatedBkg[type[itype]];
            cout << "Generate "<< generatedBkg[type[itype]] << " " << " background events\n";
        }

        if (systematics) {
            sprintf(hname,"signal_gen%s%i", type[itype], templateToUse);
        } else {
            sprintf(hname,"signal%s%i", type[itype], templateToUse);
        }

        datasets[type[itype]] = w->pdf(hname)->generateBinned(*topMass, generatedSignal[type[itype]])->createHistogram(hname,*topMass);

        if (systematics) {
            sprintf(hname,"background_gen%s", type[itype]);
        } else {
            sprintf(hname,"background%s", type[itype]);
        }
        datasets[type[itype]]->Add(w->pdf(hname)->generateBinned(*topMass, generatedBkg[type[itype]])->createHistogram(hname,*topMass));
    }

    if (data!=0) {
        delete data;
    }
    assembleDatasets();

    data->fillHistogram(toyDataHisto, *topMass);

    return data->numEntries();
}

void MassFit::write_toy(TString name)
{
    name += ilumi;
    name.ReplaceAll ( " " , "" );
    cout <<name<<endl;
    TFile * out = new TFile(name+".root","RECREATE");
    toy_mean->Write();
    toy_bias->Write();
    toy_error->Write();
    toy_pull->Write();
    toy_LL->Write();
    out->Close();
}

void MassFit::calibration(int number)
{
    time_t start,end;
    time (&start);

    nFitFailed = 0;
    nFitTried = 0;
    int min=minTemplate, max=maxTemplate;
    TVectorD x(max-min+1), y(max-min+1), ex(max-min+1), ey(max-min+1);

    int pts=0;
    int i;
    char hname[50];
    TString name = TString("calibration_");
    name += ilumi;
    name.ReplaceAll ( " " , "" );
    cout <<name<<endl;
    quietFit(true);

    TFile * out = new TFile(name+".root","RECREATE");
    cout << "Opened TFile" << endl;

    for (i = min;i!=max;++i)
    {
        do_toys(number,i);
        toy_mean->Fit("gaus");

        x[pts]=mcSignalTemplMass[i]; ex[pts]=0.;
        y[pts]=toy_mean->GetFunction("gaus")->GetParameter(1);
        ey[pts]=toy_mean->GetFunction("gaus")->GetParameter(2)/sqrt(number);
        ++pts;

        sprintf(hname,"meanMass_%.2f", mcSignalTemplMass[i]);
        toy_mean->Clone(hname)->Write();
        sprintf(hname,"biasMass_%.2f", mcSignalTemplMass[i]);
        toy_bias->Clone(hname)->Write();
        sprintf(hname,"errMass_%.2f", mcSignalTemplMass[i]);
        toy_error->Clone(hname)->Write();
        sprintf(hname,"pullMass_%.2f", mcSignalTemplMass[i]);
        toy_pull->Clone(hname)->Write();
        sprintf(hname,"LL_%.2f", mcSignalTemplMass[i]);
        toy_LL->Clone(hname)->Write();
    }
    for (i = min, pts=0;i!=max;++i,++pts)
        cout << "Template mass "<< mcSignalTemplMass[i]<< " - Fit: " << y[pts] <<" / "<<ey[pts]<<endl;

    if (grc!=0) {
        delete grc;
    }
    grc = new TGraphErrors(x,y,ex,ey);
    grc->SetName("grc");
    TF1* f1 = new TF1("f1","pol1",minTMass, maxTMass);
    grc->Fit("f1");
    TF1* f2 = new TF1("f2","pol1",minTMass, maxTMass);
    f2->SetParameter(0,0.);
    f2->SetParameter(1,1.);
    f2->SetLineColor(4);

    grc->Draw("a*");
    f2->Draw("same");
    grc->Write();
    f1->Write();
    f2->Write();
    out->Close();
    cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
    time (&end);
    double dif = difftime (end,start);
    printf ("It took  %.2lf seconds to do the whole thing, %.2lf per loop.\n", dif , dif/number);
}


void MassFit::pdfCalibration(int number)
{
    time_t start,end;
    time (&start);

    nFitFailed = 0;
    nFitTried = 0;

    int pts=0;
    char hname[50];
    quietFit(true);
    TFile * out = new TFile("pdfSystematics.root","RECREATE");
    TVectorD x(41), y(41), ex(41), ey(41);

    for (int i = 0;i<41;++i)
    {
        cout <<"a\n";
        do_toys(number,i);
        toy_mean->Fit("gaus");

        x[pts]=mcSignalTemplMass[i]; ex[pts]=0.;
        y[pts]=toy_mean->GetFunction("gaus")->GetParameter(1);
        ey[pts]=toy_mean->GetFunction("gaus")->GetParameter(2)/sqrt(number);
        ++pts;

        sprintf(hname,"meanMass_%.2f", mcSignalTemplMass[i]);
        toy_mean->Clone(hname)->Write();
        sprintf(hname,"errMass_%.2f", mcSignalTemplMass[i]);
        toy_error->Clone(hname)->Write();
        sprintf(hname,"pullMass_%.2f", mcSignalTemplMass[i]);
        toy_pull->Clone(hname)->Write();
    }
    for (int i = 0;i<41;++i)
        cout << "Template pdf "<< i<< " - Fit: " << y[i] <<" / "<<ey[i]<<endl;

    out->Close();
    cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
    time (&end);
    double dif = difftime (end,start);
    printf ("It took  %.2lf seconds to do the whole thing, %.2lf per loop.\n", dif , dif/number);
}

void MassFit::set_overflow_bins(TH1F * h)
{
    Int_t _last_bin = h->GetNbinsX();
    Double_t _overflow = h->GetBinContent(_last_bin+1);
    Int_t n_entries = (Int_t)h->GetEntries();
    h->AddBinContent(1,h->GetBinContent(0));
    h->AddBinContent(_last_bin,_overflow);
    h->SetBinContent(0,0);
    h->SetBinContent(_last_bin+1,0);
    h->SetEntries(n_entries);
}

MassFit::~MassFit(){
}

void MassFit::printFit(int point, bool toFile, char* name)
{
    if (fittedTemplate!=point) fitPoint(point);
    printFit(true, toFile, name);
}

void MassFit::printFit(bool withFit, bool toFile, char* name)
{
    char hname[50];
    map<string, RooPlot*> frames;
    int bins = (int) (upperMassCut-lowerMassCut)/10;
    for (int itype = 0;itype!=maxType;++itype) {
        sprintf(hname,"Mass distribution %s", type[itype]);
        RooPlot* frame = topMass->frame(Bins(bins),Title(hname)) ;
        sprintf(hname,"sample==sample::%s", type[itype]);
        data->plotOn(frame,Cut(hname), DrawOption("ZpE") ) ;
        frames[type[itype]] = frame;
    }

    RooPlot* frame4 = topMass->frame(Bins(bins),Title("Mass distribution")) ;
    data->plotOn(frame4, DrawOption("Zp")) ;

    if (withFit && fittedTemplate==-1) {
        cout << "No fit performed yet\n";
    } else if (withFit) {

        for (int itype = 0;itype!=maxType;++itype) {
            sprintf(hname,"model%s%.2f", type[itype], mcSignalTemplMass[fittedTemplate]);
            pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data), LineColor(kBlack)) ;
            sprintf(hname,"signal%s%.2f", type[itype],mcSignalTemplMass[fittedTemplate]);
            pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDotted), LineColor(kBlue)) ;
            sprintf(hname,"background%s", type[itype]);
            pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDashed), LineColor(kRed) ) ;
        }

        sprintf(hname,"model%.2f", mcSignalTemplMass[fittedTemplate]);
        pdffit->plotOn(frame4,Components(hname),ProjWData(*sample,*data), LineColor(kBlack)) ;
    }

    TCanvas* c;
    if (maxType>3) {
        c = new TCanvas("mass","Top ",1200,1200) ;
        c->Divide(3,3) ;
    } else {
        c = new TCanvas("mass","Top ",800,800) ;
        c->Divide(2,2) ;
    }
    int i=1;
    for (int itype = 0;itype!=maxType;++itype) {
        frames[type[itype]]->SetMinimum(1e-5);
        frames[type[itype]]->GetYaxis()->SetTitle("Events / 10 GeV");
        c->cd(i);
        gPad->SetLeftMargin(0.15);
        sprintf(hname,"Reconstructed mass [GeV] %s", type[itype]);
        frames[type[itype]]->GetXaxis()->SetTitle(hname);
        frames[type[itype]]->GetYaxis()->SetTitleOffset(1.1) ;
        frames[type[itype]]->Draw() ;
        cout << frames[type[itype]]->chiSquare()<<endl;

        ++i;
        TLatex* text1 = new TLatex(3.570061,23.08044,"CMS Preliminary");
        text1->SetNDC();
        text1->SetTextAlign(13);
        text1->SetX(0.57);
        text1->SetY(0.928);
        text1->SetTextFont(42);
        text1->SetTextSizePixels(24);
        text1->Draw();

        sprintf(hname,"%g%s", ilumi, " pb^{-1} at #sqrt{s} = 8 TeV");
        TLatex* text2 = new TLatex(3.570061,23.08044,hname);
        text2->SetNDC();
        text2->SetTextAlign(13);
        text2->SetX(0.57);
        text2->SetY(0.88);
        text2->SetTextFont(42);
        text2->SetTextSizePixels(24);
        text2->Draw();
    }

    if (i<10){
        frame4->SetMinimum(1e-5);
        frame4->GetYaxis()->SetTitle("Events / 10 GeV");
        c->cd(i) ; gPad->SetLeftMargin(0.15) ; frame4->GetXaxis()->SetTitle("Reconstructed mass [GeV] - All");frame4->GetYaxis()->SetTitleOffset(1.1) ; frame4->Draw() ;
        cout << frame4->chiSquare()<<endl;
        TLatex* text1 = new TLatex(3.570061,23.08044,"CMS Preliminary");
        text1->SetNDC();
        text1->SetTextAlign(13);
        text1->SetX(0.57);
        text1->SetY(0.928);
        text1->SetTextFont(42);
        text1->SetTextSizePixels(24);
        text1->Draw();

        sprintf(hname,"%g%s", ilumi, " pb^{-1} at #sqrt{s} = 8 TeV");
        TLatex* text2 = new TLatex(3.570061,23.08044,hname);
        text2->SetNDC();
        text2->SetTextAlign(13);
        text2->SetX(0.57);
        text2->SetY(0.88);
        text2->SetTextFont(42);
        text2->SetTextSizePixels(24);
        text2->Draw();

    }

    if (toFile) {
        char name2[200],name3[200];
        if (strlen(name)==0) {
            if (withFit) sprintf(name2,"%s/fit_result/mass%.2f_", absolutePath,
                    mcSignalTemplMass[fittedTemplate]);
            else sprintf(name2,"%s/fit_result/mass_", absolutePath);
            typeTag(name2);
        } else {
            sprintf(name2,"%s", name);
        }
        sprintf(name3,"%s.C",name2); c->Print(name3);
        sprintf(name3,"%s.png",name2); c->Print(name3);
        sprintf(name3,"%s.pdf",name2); c->Print(name3);
    }
}

void MassFit::typeTag(char *name)
{
    char name2[200];
    for (int j =0; j<maxType;++j) {
        sprintf(name2,"%s%s", name, type[j]);
        if (j<(maxType-1)) sprintf(name,"%s_", name2);
        else sprintf(name,"%s", name2);
    }
    cout << name<<endl;
}

void MassFit::printMassRange()
{
    cout << "Templates: " <<mcSignalTemplMass.size()<<endl;
    for (vector<float>::iterator i = mcSignalTemplMass.begin(); i!= mcSignalTemplMass.end(); ++i)
        cout << i-mcSignalTemplMass.begin()<<" "<<*i<<endl;
    cout << "Nominal template (width="<<mcSignalTemplMass[nominalTemplate]<<"): " << nominalTemplate<<endl;
}

void MassFit::cmsprelim()
{
    TLatex* text1 = new TLatex(3.570061,25,"CMS Preliminary, 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
    text1->SetNDC();
    text1->SetTextAlign(13);
    text1->SetX(0.167);
    text1->SetY(0.995);
    text1->SetTextFont(42);
    text1->SetTextSizePixels(24);
    text1->Draw();
}
