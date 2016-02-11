#include "src/RPoissonAnalysis.C"


void rnallyDileptonAnalysis()
{
    
    // Create the analysis
    ///gROOT->ProcessLine(".L src/RPoissonAnalysis.C");
    RPoissonAnalysis *r = new RPoissonAnalysis();


    // Set the preliminary settings

    r->sProp = "mass";
    r->sAcro = "amwt";
    r->sUnit = "GeV";

    r->dataFileLoc = "./samples/2012_combined.root";
    r->systFileLoc = "./rnally_calibration_19.700000762939453.root";

    r->mcSigTemplVal = { 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5 };
    r->processes = { "2B_All", "1B_All" };
    r->mcBkgLabels = { "WJets",
                       "SingleTop",
                       "DrellYan",
                       "Diboson" };

    r->nPseudoexperiments = 1000;

    
    
    r->run();
}
