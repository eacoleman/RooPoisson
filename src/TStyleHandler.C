/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * TStyleHandler.C                                                             *
 * Author: Evan Coleman                                                        *
 *         2015                                                                *
 *                                                                             *
 * Purpose:                                                                    *
 *   Handles all style formatting. Includes standard                           *
 *   definitions for CMS_lumi and tdrStyle.                                    *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "include/TStyleHandler.h"
#include "TCanvas.h"


void TStyleHandler::CMS_lumi(TPad *pad, int iPeriod, int iPosX) {
    //if( iPosX/10==0 ) 
    //{
    //    outOfFrame = true;
    //}

    int alignY_=3;
    int alignX_=2;
    if( iPosX/10==0 ) alignX_=1;
    if( iPosX==0    ) alignX_=1;
    if( iPosX==0    ) alignY_=1;
    if( iPosX/10==1 ) alignX_=1;
    if( iPosX/10==2 ) alignX_=2;
    if( iPosX/10==3 ) alignX_=3;
    if( iPosX == 0  ) relPosX = 0.12;
    int align_ = 10*alignX_ + alignY_;

    float H = pad->GetWh();
    float W = pad->GetWw();
    float l = pad->GetLeftMargin();
    float t = pad->GetTopMargin();
    float r = pad->GetRightMargin();
    float b = pad->GetBottomMargin();
    float e = 0.025;

    pad->cd();

    TString lumiText;
    if( iPeriod==1 )
    {
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
    }
    else if ( iPeriod==2 )
    {
        lumiText += lumi_8TeV;
        lumiText += " (8 TeV)";
    }
    else if( iPeriod==3 ) 
    {
        lumiText = lumi_8TeV; 
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
    }
    else if ( iPeriod==4 )
    {
        lumiText += lumi_13TeV;
        lumiText += " (13 TeV)";
    }
    else if ( iPeriod==7 )
    { 
        if( outOfFrame ) lumiText += "#scale[0.85]{";
        lumiText += lumi_13TeV; 
        lumiText += " (13 TeV)";
        lumiText += " + ";
        lumiText += lumi_8TeV; 
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumi_7TeV;
        lumiText += " (7 TeV)";
        if( outOfFrame) lumiText += "}";
    }
    else if ( iPeriod==12 )
    {
        lumiText += "8 TeV";
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    

    float extraTextSize = extraOverCmsTextSize*cmsTextSize;

    latex.SetTextFont(42);
    latex.SetTextAlign(31); 
    latex.SetTextSize(lumiTextSize*t);    
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

    if( outOfFrame )
    {
        latex.SetTextFont(cmsTextFont);
        latex.SetTextAlign(11); 
        latex.SetTextSize(cmsTextSize*t);    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }

    pad->cd();

    float posX_=0;
    if( iPosX%10<=1 )
    {
        posX_ =   l + relPosX*(1-l-r);
    }
    else if( iPosX%10==2 )
    {
        posX_ =  l + 0.5*(1-l-r);
    }
    else if( iPosX%10==3 )
    {
        posX_ =  1-r - relPosX*(1-l-r);
    }
    float posY_ = 1-t - relPosY*(1-t-b);
    if( !outOfFrame )
    {
        if( drawLogo )
        {
            posX_ = l + 0.045*(1-l-r)*W/H;
            posY_ = 1 - t - 0.045*(1-t-b);
            float xl_0 = posX_;
            float yl_0 = posY_ - 0.15;
            float xl_1 = posX_ + 0.15*H/W;
            float yl_1 = posY_;
            TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
            TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
            pad_logo->Draw();
            pad_logo->cd();
            //CMS_logo->Draw("X");
            pad_logo->Modified();
            pad->cd();
        }
        else
        {
            latex.SetTextFont(cmsTextFont);
            latex.SetTextSize(cmsTextSize*t);
            latex.SetTextAlign(align_);
            latex.DrawLatex(posX_, posY_, cmsText);
            if( writeExtraText ) 
            {
                latex.SetTextFont(extraTextFont);
                latex.SetTextAlign(align_);
                latex.SetTextSize(extraTextSize*t);
                latex.DrawLatex(posX_,
                                posY_-
                                relExtraDY*cmsTextSize*t,
                                extraText);
            }
        }
    }
    else
        if( writeExtraText )
        {
            if( iPosX==0 ) 
            {
                posX_ = l + relPosX*(1-l-r);
                posY_ = 1 - t + lumiTextOffset*t;
            }
            latex.SetTextFont(extraTextFont);
            latex.SetTextSize(extraTextSize*t);
            latex.SetTextAlign(align_);
            latex.DrawLatex(posX_,
                            posY_,
                            extraText);      
        }
    return;
}


void TStyleHandler::CMS_lumi() {
    TCanvas *pad = new TCanvas("CMS_lumi_pad", "", 600, 600);
    CMS_lumi(pad, iPeriod, iPosX);
    return;
}


TStyle* TStyleHandler::setTDRStyle() {
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    // For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(600); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //Position on screen
    tdrStyle->SetCanvasDefY(0);

    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);

    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);

    // For the histo:
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    tdrStyle->SetEndErrorSize(2);
    tdrStyle->SetMarkerStyle(20);

    //For the fit/function:
    tdrStyle->SetOptFit(0);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);

    //For the date:
    tdrStyle->SetOptDate(0);

    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS: SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.5g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);

    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.16);
    tdrStyle->SetPadRightMargin(0.02);

    // For the Global title:
    tdrStyle->SetOptTitle(0);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);

    // For the axis titles:
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");
    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);

    // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XYZ");

    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on opposite side of frame
    tdrStyle->SetPadTickY(1);

    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);

    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
    tdrStyle->SetHatchesLineWidth(5);
    tdrStyle->SetHatchesSpacing(0.05);

    tdrStyle->cd();
    return tdrStyle;
}


void TStyleHandler::setStyle(TPad *tPad, int tPeriod, int tPosX) {
    CMS_lumi(tPad, tPeriod, tPosX);
    setTDRStyle();
}
