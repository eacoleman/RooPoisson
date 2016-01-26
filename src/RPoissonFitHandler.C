#include "include/RPoissonFitHandler.h"

RPoissonFitHandler::RPoissonFitHandler() {

}

RPoissonFitHandler::~RPoissonFitHandler() {

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

}

pair<double,double> RPoissonFitHandler::minimize(bool all = true, 
                                                 int pointsToUse = 2) {

}

pair<double,double> RPoissonFitHandler::minFake(bool all = true, 
                                                int pointsToUse = 2) {

}

void RPoissonFitHandler::runCalibration(int numberOfExps = 1000) {

}

void RPoissonFitHandler::fit() {

}
