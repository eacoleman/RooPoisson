# RooPoisson
A set of ROOT utilities for Poisson likelihood fit analyses at CMS. It's nice 
not to have to reinvent the wheel every time we'd like to run one of these.

## Background
A Poisson likelihood analysis for a property P compares Monte Carlo (MC) samples 
generated at different values of P to data histograms. We extract a compatiblity,
or "likelihood," using chi-squared statistics for data-to-MC fits. We then fit
the -ln(likelihood) of the MC samples against their assumed values of P. 
Maximizing the compatibility means minimizing the -ln(likelihood) against P; we
minimize with a simple quadratic regression to determine an initial measurement
of P.

We calibrate for bias by performing dummy trials: let one of the MC samples
replace the data sample in our analysis. Theoretically, we should expect perfect
agreement between the output P and the P we assumed while generating the 
selected MC. If for some reason we get a value which does not agree with the
theoretical one, we can suppose there is some sort of linear bias affecting the
measurement. We perform these trials once for each MC sample, fit the plot of
output P to assumed P with a linear functional y = Ax + B, and then take our
measurement against data (call it m) and invert, so that our final measurement
(call it m') is just m' = (m - B)/A.

## Requirements
 - ROOT 6.06 and above
 - data and pre-generated MC samples

## Usage

The analysis process is fairly streamlined.

### Running

In its simplest form, an analysis looks like:

    // filename.c
    #include "src/RPoissonAnalysis.C"

    void filename() {
      RPoissonAnalysis *r = new RPoissonAnalysis();
      RPoissonAnalysis->run();
    }

You'll then get the default output, which runs on one of the given samples. 
A more in-depth example of the usage is given in the macro
`rnallyDileptonAnalysis.c`.

Options which you can modify as needed are as follows:

 o (string) sProp - the name of the property to measure (should be one word, no TeX)
 o (string) sUnit - the units of this property (ROOT TeX formatting acceptable)
 o (string) sAcro - the acronym of the analysis (one word, no TeX)
 o (float)  lLumi - the luminosity of your samples
 o (string) lUnit - the SI inverse barns prefix to include with the lumi (e.g. "f")
 o (string) dataFileLoc  - the location of the samples
 o (string) calibFileLoc - the place to store the calibration file
 o (string) sigProcess - the signal process in your samples
 o (vector<TString>) mcBkgLabels - the background processes in your samples
 o (vector<string>)  processes - the names of your subprocesses (additional final states)
 o (vector<float>) mcSigTemplVal - the template sProp values for your signal MC
 o (int) nomTemplIndex - the index of the nominal template in mcSigTemplVal
 o (float) minPropVal - the minimum sProp to consider
 o (float) maxPropVal - the maximum sProp to consider
 o (bool) calibLastPts - whether to include the highest and lowest sProp values when calibrating
 o (int) nPseudoexperiments - the number of toys to run in the analysis

You can additionally modify the following style options:

 o (float) canvWid, canvHei - canvas width and height
 o (int) iPeriod - CMS_lumi standard iPeriod configuration
 o (int) toy<PLOTNAME>Bins - number of bins to use in histogram toy<PLOTNAME>
   o PLOTNAME can be Mean, Bias, Pull
 o (float) toy<PLOTNAME>Delt - the delta(x) to plot around the value of interest in toy<PLOTNAME>
   o PLOTNAME can be Mean, Bias, Pull
 o toyErrorBins, toyErrorMax, toyErrorMin control the toyError histogram settings
 o toyLL is a TH2F, with corresponding settings toyLLRes<X/Y><Bin/Min/Max>

### Samples

The included .root files in `/samples/` show examples of what the input samples
need to look like.

The general guidelines are as follows:
 o Data histograms are titled: `<sAcro>__Data_<subprocess>`
 o Signal MC histograms: `<sAcro>__<signalProcess>_<templateValue>_<subprocess>`
 o Background MC histograms: `<sAcro>__<signalProcess>_<subprocess>`

Note that two underscores (`__`, not `_`) follow `sAcro` in both formats. 
Template values run strictly to the second decimal place, but if you need more 
or fewer digits, you can easily modify RPoissonAnalysis.C by replacing all
instances of `%.2f` with `%.<your precision>f`. (We are working on a way to
do this dynamically).