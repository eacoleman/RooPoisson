{
    
    // Create the analysis
    gROOT->ProcessLine(".L src/RPoissonAnalysis.C");
    RPoissonAnalysis *r = new RPoissonAnalysis();


    // Set the preliminary settings

    r->systematics = false;

    r->dataFileLoc = "./samples/2012_combined.root";
    r->systFileLoc = "./rnally_calibration_19.700000762939453.root";

    r->mcSigTemplVal = { 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5 };

    r->sProp = "mass";
    r->sAcro = "amwt";
    r->sUnit = "GeV";

    cout << r->sAcro <<endl;

    r->processes = { "0B_EE", "1B_EE", "2B_EE",
                     "0B_EM", "1B_EM", "2B_EM",
                     "0B_MM", "1B_MM", "2B_MM" };

    r->mcBkgLabels = { "WJets",
                       "SingleTop",
                       "DrellYan",
                       "Diboson" };

    r->setup();
    r->runCalibration(100);
    r->systematics = true;
    r->fitAll();
    r->minimize(false);




}
