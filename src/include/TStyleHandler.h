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
        
        void CMS_lumi();
        void CMS_lumi(TPad* pad, int iPeriod, int iPosX);

        int iPeriod = 7;
        int iPosX = 0;

        TString cmsText     = "CMS";
        float cmsTextFont   = 61;  // default is helvetic-bold

        bool outOfFrame    = false;
        bool writeExtraText = false;
        TString extraText   = "Preliminary";
        float extraTextFont = 52;  // default is helvetica-italics

        // text sizes and text offsets with respect to the top frame
        // in unit of the top margin size
        float lumiTextSize     = 0.6;
        float lumiTextOffset   = 0.2;
        float cmsTextSize      = 0.75;
        float cmsTextOffset    = 0.1;  // only used in outOfFrame version

        float relPosX    = 0.045;
        float relPosY    = 0.035;
        float relExtraDY = 1.2;

        // ratio of "CMS" and extra text size
        float extraOverCmsTextSize  = 0.76;

        TString lumi_13TeV = "20.1 fb^{-1}";
        TString lumi_8TeV  = "19.7 fb^{-1}";
        TString lumi_7TeV  = "5.1 fb^{-1}";

        bool drawLogo      = false;

        /////////////////////////////////////
        //             tdrStyle            //
        /////////////////////////////////////
        
        TStyle* setTDRStyle();

    public:
        static void setStyle(TPad *pad, int iPeriod, int iPosX);

}

#endif
