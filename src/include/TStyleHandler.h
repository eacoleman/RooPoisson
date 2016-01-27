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

        static const int iPeriod;
        static const int iPosX;

        static const TString cmsText;
        static const float cmsTextFont; 

        static const bool outOfFrame;
        static const bool writeExtraText;
        static const TString extraText;
        static const float extraTextFont;

        // text sizes and text offsets with respect to the top frame
        // in unit of the top margin size
        static const float lumiTextSize;
        static const float lumiTextOffset;
        static const float cmsTextSize;
        static const float cmsTextOffset;

        static const float relPosY;
        static const float relExtraDY;

        // ratio of "CMS" and extra text size
        static const float extraOverCmsTextSize;

        static const TString lumi_13TeV;
        static const TString lumi_8TeV;
        static const TString lumi_7TeV;

        static const bool drawLogo;

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
