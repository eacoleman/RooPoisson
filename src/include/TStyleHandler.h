/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * TStyleHandler.h                                                             *
 * Author: Evan Coleman                                                        *
 *         2015                                                                *
 *                                                                             *
 * Purpose:                                                                    *
 *   Handles all style issues for this application. Includes standard          *
 *   definitions for CMS_lumi and tdrStyle.                                    *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TSTYLE_HANDLER 
#define TSTYLE_HANDLER 

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TImage.h"
#include "TASImage.h"
#include "TStyle.h"

class TStyleHandler 
{
    private:
        TStyleHandler() {}

        /////////////////////////////////////
        //             CMS_lumi            //
        /////////////////////////////////////

        static int iPeriod = 7;
        static int iPosX = 0;

        static TString cmsText     = "CMS";
        static float cmsTextFont   = 61;  // default is helvetic-bold

        static bool outOfFrame    = false;
        static bool writeExtraText = false;
        static TString extraText   = "Preliminary";
        static float extraTextFont = 52;  // default is helvetica-italics

        // text sizes and text offsets with respect to the top frame
        // in unit of the top margin size
        static float lumiTextSize     = 0.6;
        static float lumiTextOffset   = 0.2;
        static float cmsTextSize      = 0.75;
        static float cmsTextOffset    = 0.1;  // only used in outOfFrame version

        static float relPosX    = 0.045;
        static float relPosY    = 0.035;
        static float relExtraDY = 1.2;

        // ratio of "CMS" and extra text size
        static float extraOverCmsTextSize  = 0.76;

        static TString lumi_13TeV = "20.1 fb^{-1}";
        static TString lumi_8TeV  = "19.7 fb^{-1}";
        static TString lumi_7TeV  = "5.1 fb^{-1}";

        static bool drawLogo      = false;

        /////////////////////////////////////
        //             tdrStyle            //
        /////////////////////////////////////
        

    public:
        static void CMS_lumi();
        static void CMS_lumi(TPad* pad, int iPeriod, int iPosX);
        static TStyle* setTDRStyle();

        void setStyle(TPad *pad, int iPeriod, int iPosX);

};

#endif
