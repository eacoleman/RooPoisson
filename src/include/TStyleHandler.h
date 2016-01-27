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

        const int iPeriod = 7;
        const int iPosX = 0;

        const TString cmsText     = "CMS";
        const float cmsTextFont   = 61;  // default is helvetic-bold

        const bool outOfFrame    = false;
        const bool writeExtraText = false;
        const TString extraText   = "Preliminary";
        const float extraTextFont = 52;  // default is helvetica-italics

        // text sizes and text offsets with respect to the top frame
        // in unit of the top margin size
        const float lumiTextSize     = 0.6;
        const float lumiTextOffset   = 0.2;
        const float cmsTextSize      = 0.75;
        const float cmsTextOffset    = 0.1;  // only used in outOfFrame version

        const float relPosX    = 0.045;
        const float relPosY    = 0.035;
        const float relExtraDY = 1.2;

        // ratio of "CMS" and extra text size
        const float extraOverCmsTextSize  = 0.76;

        const TString lumi_13TeV = "20.1 fb^{-1}";
        const TString lumi_8TeV  = "19.7 fb^{-1}";
        const TString lumi_7TeV  = "5.1 fb^{-1}";

        const bool drawLogo      = false;

        /////////////////////////////////////
        //             tdrStyle            //
        /////////////////////////////////////
        

    public:
        static void CMS_lumi();
        static void CMS_lumi(TPad* pad, int iPeriod, int iPosX);
        static TStyle* setTDRStyle();

        static void setStyle(TPad *pad, int iPeriod, int iPosX);

};

#endif
