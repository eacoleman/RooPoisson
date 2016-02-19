/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * rnallyDileptonAnalysis.c                                                    *
 * Richard Nally's 2015 analysis of the top quark mass, in cleaner code format.*
 *                                                                             *
 * Author: Evan Coleman                                                        *
 * Date:   February 2016                                                       *
 *                                                                             *
 * To run: root -l rnallyDileptonAnalysis.c                                    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "src/RPoissonAnalysis.C"

void rnallyDileptonAnalysis()
{
    
    // Create the analysis
    RPoissonAnalysis *r = new RPoissonAnalysis();

    //////////////////////////
    // Preliminary Settings //
    //////////////////////////

    // What property are we looking for?
    r->sProp = "mass";

    // What is the acronym for this analysis?
    r->sAcro = "amwt";

    // What are the units of the property we're searching for?
    r->sUnit = "GeV";

    // Where is the data? Where is the calibration to be stored/accessed?
    r->dataFileLoc = "./samples/2012_combined.root";
    r->systFileLoc = "./rnally_calibration_19.700000762939453.root";

    // Utility labels - our signal template values, process names, etc.
    r->mcSigTemplVal = { 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5 };
    r->processes = { "2B_All" };
    r->mcBkgLabels = { "WJets",
                       "SingleTop",
                       "DrellYan",
                       "Diboson" };

    // Set the cuts
    r->maxPropVal = 400;
    r->minPropVal = 100;

    // Do not calibrate the endpoints of mcSigTemplVal
    r->calibLastPts = false;

    // Set the nominal MC template index in mcSigTemplVal
    r->nomTemplIndex = 4;

    // Set the number of pseudoexperiments to run in the calibration
    r->nPseudoexperiments = 1000;
    
    // Plot binnings
    //  - bins = # bins to use in histogram/graph
    //  - delt = delta(x) to plot about the value of interest
    //  - min  = minimum value for axis
    //  - max  = maximum value for axis
    int   toyMeanBins  =  100;
    float toyMeanDelt  =  2.5;
    int   toyBiasBins  =  100;
    float toyBiasDelt  =  2.5;
    int   toyPullBins  =  100;
    float toyPullDelt  =   5.;
    
    int   toyErrorBins =  500;
    float toyErrorMin  =  0.1;
    float toyErrorMax  =  0.2;

    /////////////
    // Running //
    /////////////

    // Run the analysis!
    r->run();
}
