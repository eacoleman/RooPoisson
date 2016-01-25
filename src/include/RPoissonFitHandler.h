#ifndef RPOISSON_FIT_HANDLER
#define RPOISSON_FIT_HANDLER

// C Inclusions

// ROOT Inclusions
#include "TRandom3.h"

// RooFit Inclusions


// RooPoisson Inclusions
#include "RPoissonAnalysis.h"

class RPoissonFitHandler
{
    private:
        /////////////
        // Methods //
        /////////////
       
        //
        double fitPoint(int index);
        int    fitAll();

        //
        pair<double, double> minimize(bool all=true, int pointsToUse=2);
        pair<double, double>  minFake(bool all=true, int pointsToUse=2);

        //
        void runCalibration(int numberOfExps = 1000);

        ///////////////////////
        // Utility variables //
        ///////////////////////

        //RPoissonAnalysis *a;
        
        TRandom3 _random;



    public:
        /////////////
        // Methods //
        /////////////

        // constructor/destructor
        RPoissonFitHandler();
        ~RPoissonFitHandler();

        //
        void fit();

};

#endif
