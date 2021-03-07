// TOF separation power = f(pT,y) for HL-LHC experiments

// Abstract :
//         The document addresses comparative performances expected for various configurations of Time-Of-Flight detectors,
//         for several CERN experiments (ALICE-1, ALICE-3, ATLAS, CMS).
//         The figures are based on analytical formulae derived for charged particles travelling in a solenoidal magnetic field.
//         While these performances correspond to ideal limits (with respect to real data or full simulations including particle transport),
//         the main point is to assess the theoretical separation among various identified species 
//         (e±, µ±, π±, K±, p±, d±, t±, 3He±, 4He±)
//         at the horizon of runs IV (~2027-2029) or V (> 2030) of the High-Luminosity LHC, HL-LHC.
//         The intent is thus to appreciate the experimental realm of possible physics cases.
//
// As of 2021-03, the corresponding document and figures can be found under arXiv: 
//          ...
// GitHub repository : 
//          https://github.com/maireiphc/TOFseparationPowerAsFuncPtY

// Origin : 
//          Wei Li, Rice University http://wl33.web.rice.edu/
//          i.e. Fig 1.5 p.8 in the Technical Design Report of MTD detector CERN-LHCC-2019-003 / CMS-TDR-020 = https://cds.cern.ch/record/2667167

// Versioning :
//          First implementation : 29 April 2020, Antonin Maire, antonin.maire_AT_cern.ch
//          Last modification : 07 March 2021, Antonin Maire


#include "Riostream.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphBentErrors.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TASImage.h"



static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);



enum gkColor{
    red     = 0,
    orange  = 1,
    yellow  = 2,
    kaki    = 3,
    green   = 4,
    cyan    = 5,
    azure   = 6,
    blue    = 7,
    violet  = 8,
    magenta = 9,
    black   = 10,
    nbColor
};



// Preferred colors and markers
const Int_t fillColors[] = {kRed-10,     kOrange-9,   kYellow-9,     kSpring+1,    kGreen-7,    kCyan-8,    kAzure-8,    kBlue-9,     kViolet-9,    kMagenta-9,   kGray+1     }; // for syst bands
const Int_t colors[]     = {kRed+1 ,     kOrange+1,   kYellow+0,     kSpring+9,    kGreen+2,    kCyan+2,    kAzure+2,    kBlue+1,     kViolet-5,    kMagenta+1,   kBlack      };
const Int_t markers[]    = {kFullCircle, kOpenSquare, kFullDotMedium,kOpenSquare,  kOpenSquare, kFullCross, kFullCircle, kFullCircle, kOpenCircle,  kOpenDiamond, kFullCircle };



enum gkParticle{
    kPDG_ePm,
    kPDG_muPm,
    kPDG_piPm,
    kPDG_Kpm,
    kPDG_K0s,
    kPDG_p,
    kPDG_2H,  
    kPDG_3H,  
    kPDG_3He, 
    kPDG_4He,   
    kPDG_phi1020,
    kPDG_Lambda,
    kPDG_XiPm,
    kPDG_Omega,
    kPDG_D0,
    kPDG_Dpm,
    kPDG_Dstar,
    kPDG_Ds,
    kPDG_LambdaC,
    kPDG_Jpsi,
    kPDG_Psi2S,
    kPDG_Upsilon1S,
    kPDG_Upsilon2S,
    kPDG_Upsilon3S,
    kPDG_Wpm,
    kPDG_Z0,
    nbPart
};



const TString Str_HadronLatexName[] = { "e^{#pm}", "#mu^{#pm}", "#pi^{#pm}","K^{#pm}", "K^{0}#kern[-1.5]{_{s}}", "p",  "d",  "t", "^{3}He", "^{4}He", "#phi(#scale[0.8]{1020})", "#Lambda", "#Xi^{#pm}", "#Omega^{#pm}", "D^{0}", "D^{#pm}", "D*(#scale[0.8]{2010})^{#pm}", "D^{+}#kern[-1.5]{_{s}}", "#Lambda^{+}#kern[-1.0]{_{c}}", "J/#psi", "#psi(2S)", "Y(1S)", "Y(2S)", "Y(3S)",  "W^{#pm}", "Z^{0}" };
const TString Str_HadronASCIIName[] = { "ePm", "muPm", "piPm", "Kpm", "K0s", "p", "2H", "3H", "3He", "4He",   "phi1020", "Lambda", "XiPm", "Omega", "D0", "Dpm", "Dstar", "Ds", "LambdaC", "Jpsi", "Psi2S", "Upsilon1S", "Upsilon2S", "Upsilon3S", "Wpm", "Z0"};
const Double_t mPDG[]          = {0.0005109989461, 0.1056583745, 0.13957018, 0.493677, 0.497611,                 0.938272081, 
    1.8756129, 2.8089211, 2.809413485, 3.728400129, 1.019461,                  1.115683,  1.31486,     1.67245,        1.86483, 1.86958,   2.01026,                       1.96827,                  2.28646,                        3.096900, 3.686097,   9.64030, 10.02326, 10.3552, 80.385,    91.1876 };
    // PDG 2016



void myLegendSetUp(TLegend *currentLegend = 0,
                    float currentTextSize = 0.07);

void myPadSetUp(   TPad *currentPad, 
                    Bool_t kLog         = 0,
                    float currentLeft   = 0.11, 
                    float currentTop    = 0.04, 
                    float currentRight  = 0.04, 
                    float currentBottom = 0.15);

void myHistoSetUp( TH1 *currentHisto             = nullptr, 
                    Size_t  currentMarkerSize    = 1.0,
                    Style_t currentMarkerStyle   = 20,
                    Color_t currentMarkerColor   = 0,
                    Width_t currentLineWidth     = 1, 
                    Style_t currentLineStyle     = 1, 
                    Color_t currentLineColor     = 0,
                    Style_t currentFillStyle     = 1001,
                    Color_t currentFillColor     = kGray);


void myGraphSetUp( TGraph  *currentGraph       = nullptr, 
                    Size_t  currentMarkerSize  = 1.0,
                    Style_t currentMarkerStyle = 20,  
                    Color_t currentMarkerColor = 0,
                    Width_t currentLineWidth   = 1,
                    Style_t currentLineStyle   = 1,
                    Color_t currentLineColor   = 0,
                    Style_t currentFillStyle   = 1001,
                    Color_t currentFillColor   = kGray);

void myFuncSetUp(  TF1* f1, 
                    Color_t lColor      = kGray+1, 
                    Style_t lLineStyle  = 1, 
                    Width_t lLineWidth  = 1);

void myOptions(     Bool_t kStat = 0);



enum gkCorner {
    kMidY_SW = 0,
    kMidY_NW = 1,
    kMidY_SEE,
    kMidY_NEE,
    kFwdY_NW,
    kFwdY_SEE,
    kFwdY_NEE,  
    kNall
};

TString Str_Corner[] = { 
    "kMidY_SW", 
    "kMidY_NW", 
    "kMidY_SEE", 
    "kMidY_NEE", 
    "kFwdY_NW", 
    "kFwdY_SEE",
    "kFwdY_NEE", 
    "kNall" };


    
Double_t funcEta2Rap(Double_t *x, Double_t *par){    
    // NOTE :  y = f(x=pT, param = m0,eta0), see https://en.wikipedia.org/wiki/Pseudorapidity
 
    Double_t pT  = x[0];
    Double_t m0  = par[0];
    Double_t eta = par[1]; 
    
    return TMath::Log(  ( TMath::Sqrt( m0*m0   +   pT*pT*TMath::CosH(eta)*TMath::CosH(eta) ) + pT * TMath::SinH(eta) ) / TMath::Sqrt( m0*m0 + pT*pT ));
}
    
    
    
Double_t funcPtMinBarrelDelta(Double_t *x, Double_t *par){// eq 45, PtMinBarrelDelta = f(eta) 
    Double_t eta = x[0];
    Double_t lPtMinBarrel      = par[0];
    Double_t lDeltaMinAngle    = par[1];
    return lPtMinBarrel / ( TMath::Sqrt( 1. - TMath::Power( TMath::Sin(lDeltaMinAngle), 2)* TMath::Power( TMath::CosH(eta), 2)) );    
}

    
void LocatePhaseSpaceCorner(    
                                std::vector<Int_t> &vCornerYBin ,
                                std::vector<Int_t> &vCornerPtBin,  
                                UInt_t kCase      = kNall,
                                Int_t lFixedYBin  = -1, 
                                Int_t lFixedPtBin = -1, 
                                TH2I *h2 = nullptr                          
){
    // Use cases :
    //  - 
    //  - 
       
    
    Int_t lNbc          = 0;
    Int_t lNbcPrev      = 0;
    Int_t lYrapBinReach = h2->FindLastBinAbove(0.5, 1, 1);  // NOTE : last bin above 0.5 (i.e. BinContent = 1) along the X axis (1), starting from bin 1 = spot abscissa of kFwdY_NW
    Int_t lPtBinReach   = h2->FindLastBinAbove(0.5, 2, 1);  // NOTE : last bin above 0.5 (i.e. BinContent = 1) along the Y axis (2), starting from bin 1 = spot ordinate of kFwdY_NW    

    
    Int_t lYMaxBin      = -1; 
    Int_t lNbPtBins     = -1;
    
    Int_t lInitYBin     = -1; 
    Int_t lInitPtBin    = -1;
    
    
    if(lFixedYBin > 0){
        lInitYBin = lFixedYBin;
        lYMaxBin  = lFixedYBin+1;  
    } 
    else{                        
        lInitYBin = 1;
        lYMaxBin  = h2->GetNbinsX();
    }
    
    if(lFixedPtBin > 0){
        lInitPtBin = lFixedPtBin;
        lNbPtBins  = lFixedPtBin+1;
    }
    else{
        lInitPtBin = 1;    
        lNbPtBins  = h2->GetNbinsY();
    }
        
        
    // Sanity checks among arguments    
    if(!h2)                                         { Printf("Missing valid histo h2... exit!");                                return;}
    if( lYMaxBin < 1        || lNbPtBins < 1 )      { Printf("Wrong values for nb of bins (<1)... exit!");                      return;}
    if( lFixedYBin == 0     || lFixedPtBin == 0 )   { Printf("Wrong values for fixed Y or Pt bin (=0 = underflow)... exit!");    return;}  
    
    if(kCase == kMidY_SW    && lFixedYBin != 1)  {Printf("Option case, mismatch in arguments (kMidY_SW  needs loop over column 1)... exit!"); return;}
    if(kCase == kMidY_NW    && lFixedYBin != 1)  {Printf("Option case, mismatch in arguments (kMidY_NW  needs loop over column 1)... exit!"); return;}
        
    
    Printf(" ");
    Printf("LocatePhaseSpaceCorner() invoked : Case (%d - %s)  with YBin = %+04d / PtBin = %+04d", kCase, Str_Corner[kCase].Data(), lFixedYBin, lFixedPtBin );
    Printf("        - Loop over Y  : from bin[%04d] to < bin[%04d]", lInitYBin , lYMaxBin );
    Printf("        - Loop over Pt : from bin[%04d] to < bin[%04d]", lInitPtBin, lNbPtBins);
    Printf("    . Y  last bin above in abscissa axis (%s) : yBin  = %04d", h2->GetName(), lYrapBinReach);    
    Printf("    . Pt last bin above in ordinate axis (%s) : PtBin = %04d", h2->GetName(), lPtBinReach);    
    
    // Printf("Check vector : bin[%03d,%03d] = (y = %6.4f, pT = %6.3f GeV/c)",
    //                             vCornerYBin  [kCase],
    //                             vCornerPtBin [kCase],
    //                             h2->GetXaxis()->GetBinCenter( vCornerYBin  [kCase] ), 
    //                             h2->GetYaxis()->GetBinCenter( vCornerPtBin [kCase] ) );    
    
    
    
    for(Int_t iYBin = lInitYBin; iYBin < lYMaxBin; ++iYBin){
        
        Short_t kLocationDone  = 0;
        
        
        for(Int_t iPtBin = lInitPtBin; iPtBin < lNbPtBins; ++iPtBin){ 
            lNbcPrev    = lNbc; // bin[n-1] content  in the loop
            lNbc        = h2->GetBinContent(iYBin, iPtBin); 
            
            Bool_t kHeavysideUp   = (lNbcPrev == 0 && lNbc  > 0); // step upwards   in two consecutive bins
            Bool_t kHeavysideDown = (lNbcPrev  > 0 && lNbc == 0); // step downwards in two consecutive bins
            
            Bool_t kConstantUp    = (lNbcPrev  > 0 && lNbc  > 0); // 1 in two consecutive bins
            Bool_t kConstantDown  = (lNbcPrev == 0 && lNbc == 0); // 0 in two consecutive bins 
            
            
            // Case kMidY_SW, kMidY_NW + (kMidY_NEE, part2)
            //      - kMidY_SW          : Loop over one column, the first column, i.e. y = 0
            //      - kMidY_NW          : idem
            //      - kMidY_NEE (part1) : spot the column with only zeros = the column next to the farthest yrap corner accessible
            //      - kMidY_NEE (part2) : Loop over the last filled column, mid y edge ~ 2.0, to flag 2 values 
            
                
                if(kCase == kMidY_SW    &&  kHeavysideUp){
                
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = 1 in fact
                        vCornerPtBin [kCase] = iPtBin;
                        kLocationDone = 1;
                }
                else
                if(kCase == kMidY_NW    &&  kHeavysideDown ){
                    
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = 1 in fact
                        vCornerPtBin [kCase] = iPtBin;
                        kLocationDone = 1;       
                }
                else                              
                if(kCase == kMidY_SEE   &&  kHeavysideUp    &&  iYBin == lYrapBinReach ){
                    
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = column spotted
                        vCornerPtBin [kCase] = iPtBin;
                        kLocationDone = 2;
                }      
                else    
                if(kCase == kMidY_NEE   &&  kHeavysideDown  &&  iYBin == lYrapBinReach ){
                    
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = column spotted
                        vCornerPtBin [kCase] = iPtBin-1;
                        kLocationDone = 2;
                }
                                        
                    
                    
                else    
                if(kCase == kFwdY_NW   &&  kHeavysideDown   &&  iPtBin == lPtBinReach+1 ){
            
                    
                    //     Printf(" . Detection (%d->%d) : bin[%03d,%03d] = (y = %6.4f, pT = %6.3f GeV/c)",
                    //                 lNbcPrev, lNbc,
                    //                 iYBin,
                    //                 iPtBin,
                    //                 h2->GetXaxis()->GetBinCenter( iYBin ), 
                    //                 h2->GetYaxis()->GetBinCenter( iPtBin ) );
                    

                        vCornerYBin  [kCase] = iYBin;
                        vCornerPtBin [kCase] = lPtBinReach;
                        kLocationDone = 3;                        
                }
                else                              
                if(kCase == kFwdY_SEE   &&  kHeavysideUp    &&  iYBin == lYrapBinReach ){
                    
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = column spotted
                        vCornerPtBin [kCase] = iPtBin;
                        kLocationDone = 2;
                }      
                else    
                if(kCase == kFwdY_NEE   &&  kHeavysideDown  &&  iYBin == lYrapBinReach ){
                    
                        vCornerYBin  [kCase] = iYBin; // = lFixedYBin here = column spotted
                        vCornerPtBin [kCase] = iPtBin;
                        kLocationDone = 2;
                }                
                
                
                

                if( kLocationDone == 1 ){
                    Printf(" . Detection (%d->%d) : bin[%03d,%03d] = (y = %6.4f, pT = %6.3f GeV/c)",
                                lNbcPrev, lNbc,
                                vCornerYBin  [kCase],
                                vCornerPtBin [kCase],
                                h2->GetXaxis()->GetBinCenter( iYBin ), 
                                h2->GetYaxis()->GetBinCenter( iPtBin ) );   // assume BinWidth = negligible wrt to targeted accuracy, i.e. %.2f, so taking Center, LowEdge, UpEdge does not matter
                    break; // save CPU no need to loop further
                }
                if( kLocationDone == 2 ) { // kCase == kMidY_NEE (part1)
                    Printf(" . Case[%d] - %9s : Column with max Y reach spotted : bin [%03d,%03d] -> column(i-1) %03d => y = %6.4f", 
                                kCase,
                                Str_Corner[kCase].Data(),
                                iYBin,
                                iPtBin,
                                lYrapBinReach,
                                h2->GetXaxis()->GetBinCenter( lYrapBinReach ) );
                    break; // break iPtBin loop, then need to break the iYBin loop (flag with lYrapBinReach)
                }
                if( kLocationDone == 3 ) { // kCase == kFwdY_NE (part1)
                    Printf(" . Case[%d] - %9s : Column with max Pt reach spotted : bin [%03d,%03d] -> (y = %6.4f, pT = %6.3f GeV/c)", 
                                kCase,
                                Str_Corner[kCase].Data(),
                                vCornerYBin  [kCase],
                                vCornerPtBin [kCase],
                                h2->GetXaxis()->GetBinCenter( iYBin-1 ),
                                h2->GetYaxis()->GetBinCenter( iPtBin  ) );
                    break; // break iPtBin loop, then need to break the iYBin loop (flag with lYrapBinReach)
                }

                
                
                
            
        }// end iPtBin
        
        
        lNbcPrev   = 0;
        lNbc       = 0;
     
        // Printf("Y [%03d], done", iYBin);
        
        if(kLocationDone) break;
        
        
    }// end iYBin
    
     Printf(" ");
    
}// end LocatePhaseSpaceCorner





TGraph* LocateLimitBetweenHalfTurns( TH2I* h2_ecap, TH2I* h2_loop, Int_t lNbHalfTurnUp){ 

    Int_t iPnt = 0;
    Int_t lNbc          = 0;
    Int_t lNbcPrev      = 0;    
    


    std::vector<Double_t> vYSpot ;     
    std::vector<Double_t> vPtSpot;

    
    
    // 1st loop : h2_ecap
    for(Int_t iPtBin = 1; iPtBin < h2_ecap->GetNbinsY() + 1; ++iPtBin){
        for(    Int_t iYBin  = 1; iYBin  < h2_ecap->GetNbinsX() + 1; ++iYBin){ 
        
   
            lNbcPrev    = lNbc; // bin[n-1] content  in the loop
            lNbc        = h2_ecap->GetBinContent(iYBin, iPtBin);
            
            Bool_t kLimitSpotted   = ( lNbcPrev == (10+lNbHalfTurnUp)   &&  lNbc  == (10+lNbHalfTurnUp)-1  ); 
                // NOTE : consecutive areas, 1 to 0, 2 to 1, etc
                //          Beware the conventions for histo filling "+10" : 0 = 10, 1 = 11, ... 1+ lAbsNbHalfTurns_Endcap1
                //         Possible cases : 1 half turn  -> 0 half turn, 
                //                          2 half turns -> 1 half turn  = naturally separated by acceptance hole, ~useless
                //                          3 half turns -> 2 half turns
                //                          4 half turns -> 3 half turns = naturally separated by acceptance hole, ~useless
                //                          5 half turns -> 4 half turns

            if( kLimitSpotted ) {
                    vYSpot .push_back( h2_ecap->GetXaxis()->GetBinCenter( iYBin  ) );
                    vPtSpot.push_back( h2_ecap->GetYaxis()->GetBinCenter( iPtBin ) );
                // Printf(" LocateLimitBetweenHalfTurns [%1d to %1d] - Ecap - y[%03d] = %5.3f , pt[%03d] = %5.3f GeV/c : lNbcPrev = %1d - lNbc = %1d", lNbHalfTurnUp, lNbHalfTurnUp-1, iPnt, vYSpot[ iPnt ], iPnt, vPtSpot[ iPnt ], lNbcPrev, lNbc);
                iPnt++;        
             }
        }// end iYBin
        
        lNbcPrev   = 0;
        lNbc       = 0;        
    }// end iPtBin
    
    
    // 2nd loop : h2_loop
        // NOTE :  the information is spread over the 2 histos, 
        //          but there should be a continuity between the two (ecap+loop), without overlap regions
        //          so that the final order of the joint search between ecap and loop should be natural in the std::vector
    for(Int_t iPtBin = 1; iPtBin < h2_loop->GetNbinsY() + 1; ++iPtBin){
        for(    Int_t iYBin  = 1; iYBin  < h2_loop->GetNbinsX() + 1; ++iYBin){ 
        
   
            lNbcPrev    = lNbc; // bin[n-1] content  in the loop
            lNbc        = h2_loop->GetBinContent(iYBin, iPtBin);
            
            
            Bool_t kLimitSpotted   = ( lNbcPrev == (10+lNbHalfTurnUp)   &&  lNbc  == (10+lNbHalfTurnUp)-1  ); 

            if( kLimitSpotted ) {
                    vYSpot .push_back( h2_loop->GetXaxis()->GetBinCenter( iYBin  ) );
                    vPtSpot.push_back( h2_loop->GetYaxis()->GetBinCenter( iPtBin ) );
                // Printf(" LocateLimitBetweenHalfTurns [%1d to %1d] - Loop - y[%03d] = %5.3f , pt[%03d] = %5.3f GeV/c : lNbcPrev = %1d - lNbc = %1d", lNbHalfTurnUp, lNbHalfTurnUp-1, iPnt, vYSpot[ iPnt ], iPnt, vPtSpot[ iPnt ], lNbcPrev, lNbc);
                iPnt++;        
             }
        }// end iYBin 
        
        lNbcPrev   = 0;
        lNbc       = 0;        
    }// end iPtBin
    
    Printf(" ");

    if( vYSpot.size() ) {
        TGraph *grHalfTurns = new TGraph( vYSpot.size(), vYSpot.data(), vPtSpot.data()  );
        return grHalfTurns;    
    }
    else return 0x0;
    
    
    
}// end LocateLimitBetweenHalfTurns





void Root_TOFseparationPowerAsFuncPtY(  TString Str_Exp = "CMS", 
                                        Double_t rInnerRadiusBarrel  = -1.0, // radial location in m
                                        Double_t rBmag               =  0.0, // B field in T
                                        Double_t rcDeltaT            = -1.0, // in ns
                                        TString rWrite = "0",
                                        Double_t rZendcapTOF         = -1.0, // z location in m
                                        Bool_t kLog = 1
){
    
// Str_Exp =    "CMS"    = CMS Phase 2 with ETL, BTL = MTD detector
//              "ATLAS"  = ATLAS phase 2 with HGTD
//              "ALICE-1"  = ALICE TOF run II
//              "ALICE-3" = TOF with SPAD Si sensors for Run V
        
// rWrite     = "eps", "pdf", "png", "tex"

    
    
// estimate BTL PID acceptance vs pT and eta
// L = pTcosh(eta)/(Bq) * acos(1-B^2*q^2*R^2/pT*2/2)
// 1/beta = sqrt(m^2/pT^2/cosh^2(eta)+1)
    
    TStopwatch stopWatch;
    stopWatch.Start();  
    
  Int_t lDebug = 0;    
    
  myOptions();
  gROOT->ForceStyle();

  TDatime now;
  int iDate = now.GetDate();
  int iYear=iDate/10000;
  int iMonth=(iDate%10000)/100;
  int iDay=iDate%100;
  TString cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"};


    TCanvas *myCanvas[30]     = {nullptr};
    TPad    *myPad   [30]     = {nullptr};
    TLegend *myLegend[30]     = {nullptr};
    TH1F    *myBlankHisto[30] = {nullptr};
    
      
    TH2I* h5sigm_ecap  [ 30 ] = {nullptr};
    TH2I* h4sigm_ecap  [ 30 ] = {nullptr};
    TH2I* h3sigm_ecap  [ 30 ] = {nullptr};
    TH2I* h2sigm_ecap  [ 30 ] = {nullptr};
    TH2I* h1sigm_ecap  [ 30 ] = {nullptr};                    
    
    TH2I* h5sigm_loop  [ 30 ] = {nullptr};
    TH2I* h4sigm_loop  [ 30 ] = {nullptr};
    TH2I* h3sigm_loop  [ 30 ] = {nullptr};
    TH2I* h2sigm_loop  [ 30 ] = {nullptr};
    TH2I* h1sigm_loop  [ 30 ] = {nullptr};
    
    
    TH2I* hdEdx_CMS    [ 30 ] = {nullptr};
    TH2I* h5sigm_barrel[ 30 ] = {nullptr};
    TH2I* h4sigm_barrel[ 30 ] = {nullptr};
    TH2I* h3sigm_barrel[ 30 ] = {nullptr};
    TH2I* h2sigm_barrel[ 30 ] = {nullptr};
    TH2I* h1sigm_barrel[ 30 ] = {nullptr};
    
    
    

    TGraph* grEta0pt5[ 30 ] = {nullptr};
    TGraph* grEta1pt0[ 30 ] = {nullptr};
    TGraph* grEta1pt5[ 30 ] = {nullptr};
    TGraph* grEta2pt0[ 30 ] = {nullptr};
    TGraph* grEta2pt5[ 30 ] = {nullptr};
    TGraph* grEta3pt0[ 30 ] = {nullptr};
    TGraph* grEta3pt5[ 30 ] = {nullptr};

    
    TGraph* grEtaMaxBarrel[ 30 ] = {nullptr};
    TGraph* grEtaMinEndcap[ 30 ] = {nullptr};
    TGraph* grEtaMaxEndcap[ 30 ] = {nullptr};
    
    TGraph* grPtMinBarrel[ 30 ] = {nullptr};
    TGraph* grPtMinEndcap[ 30 ] = {nullptr};
    TGraph* grPtMinCalo  [ 30 ] = {nullptr};
    
    TGraph* grEtaMinLambda[ 30 ] = {nullptr};
    TGraph* grEtaMinDelta [ 30 ] = {nullptr};
    
    
    TGraph* grHalfTurns_1to0[ 30 ] = {nullptr};
    TGraph* grHalfTurns_3to2[ 30 ] = {nullptr};    
    TGraph* grHalfTurns_5to4[ 30 ] = {nullptr};
    

    
    
    TPaveText*  lPaveTxt_Sep    [ 30 ] = {nullptr};
    TText*      lTxt_Sep        [ 30 ] = {nullptr};
    TLatex*     lLatex_Hyp      [ 30 ] = {nullptr};
    TPaveText*  lPaveTxt_Info   [ 30 ] = {nullptr};
    TLatex*     lLatexEta       [ 30 ] = {nullptr};
    
    TString Str_FigFilename     [ 30 ] = {nullptr};
  
  
    
    double lBz              = 0.;
    
    Bool_t kActivateBarrel  = kFALSE;
    double lRbarrelTOF      = 0.;
    double lEtaMaxBarrel    = 0.;
    double lZmaxBarrel      = 0.;
    
    Bool_t kActivateEndcap  = kFALSE;
    double lZendcapTOF      = 0.;
    double lEtaMinEndcap    = 0.;
    double lEtaMaxEndcap    = 0.;
    double lRmaxEndcap      = 0.;
    double lRminEndcap      = 0.;
    
    double lEtaMaxCalo      = 0.;
    double lRCalo           = 0;
    
    double cDeltaT          = 0.;  // time resolution in ns ! e.g. 0.03 ns  = 30 ps for CMS
    TString StrBarrelTOFname("");
    TString StrEndcapTOFname("");
    
    TString Str_LHCrun("");
        
    
    if( Str_Exp == "CMS" ){
        if( TMath::Abs(rBmag) > 1e-13 )  lBz = rBmag;     // Which field map exists for CMS  B = 1.9 T ?
        else                             lBz = -3.8  ;     // default CMS = 3.8
        kActivateBarrel = kTRUE;
        if(rInnerRadiusBarrel > 0)      lRbarrelTOF = rInnerRadiusBarrel;
        else                            lRbarrelTOF = 1.161;     // radial position on central barrel BTL, CERN-LHCC-2019-003, Section 1.4.1   
        lEtaMaxBarrel   = 1.48 ;     // 1.50 figure from Wei Li in the original code, ≠ consistent with CERN-LHCC-2019-003, Section 2.1 = 1.48
        kActivateEndcap = kTRUE;
        lZendcapTOF     = 3.04 ;     // figure from Wei Li in the original code, consistent with CERN-LHCC-2019-003, Section 1.4.2
        lEtaMinEndcap   = 1.6  ;
        lEtaMaxEndcap   = 3.0  ;
        lEtaMaxCalo     = 1.479;     // ECal TDR 1997, Section 1.6.1 + Fig 1.6 p.36 pdf, https://cms-docdb.cern.ch/cgi-bin/PublicDocDB/ShowDocument?docid=2713
        lRCalo          = 1.29 ;     // Ecal TDR 1997, p.6, Fig.3.1 p.74 pdf, https://cms-docdb.cern.ch/cgi-bin/PublicDocDB/ShowDocument?docid=2713
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.03 ;     // 30 ps Fig 1.5 TDR MTD CMS, CERN-LHCC-2019-003
        Str_LHCrun      = "HL-LHC run IV";
        StrBarrelTOFname= "BTL";
        StrEndcapTOFname= "ETL";
    }
    else if( Str_Exp == "ATLAS" ){
        if( TMath::Abs(rBmag) > 1e-13 )  lBz = rBmag;     // Which field map exists for ATLAS ?
        else                             lBz = -2.0  ;     // default ATLAS solenoid = 2.0
        kActivateBarrel = kTRUE;
        if(rInnerRadiusBarrel > 0)      lRbarrelTOF = rInnerRadiusBarrel;
        else                            lRbarrelTOF = 0.29;                 // r ~ 29 cm instead of pixel outer barrel or r = 100 cm instead of last strip ? To be tested
              // FIXME : hyp to be made for a potential SPAD layer instead of a tracker layer in ITk
              // FIXME : hyp, let's assume a SPAD layer (beware radiation tolerance as f(r)... See Seminar CERN detector, ITk, Kuehn, 2019/11/29, https://indico.cern.ch/event/865308/)   
        lEtaMaxBarrel   = 1.0  ;      // >1.1 at r = 1.0 m or 2.0 Pb cylinder Vs vertical planes vs tilted plane ?
        kActivateEndcap = kTRUE;
        
        if(TMath::Abs(rZendcapTOF - 1.5) < 1e-13) {
            lZendcapTOF     = rZendcapTOF;
            lEtaMinEndcap   = 1.22  ;     
            lEtaMaxEndcap   = 2.1  ;      
        }
        else {        
            lZendcapTOF     = 3.45 ;      // 3.45 m = HGTD LGad, https://cds.cern.ch/record/2623663 NB there are in fact 4 TOF disks per endcap... thickness, 4.0 cm fig 4.11 and tab 4.3
            lEtaMinEndcap   = 2.40 ;      // 2.40
            lEtaMaxEndcap   = 4.0  ;      // 4.0
        }
        lEtaMaxCalo     = 1.4  ;      // ECal TdR, ~ solenoid before ECal, Section 1.4.1 p.27, 1.4 value ≠ explicit but apparently recurring Fig. 1-4https://inspirehep.net/literature/432394
        lRCalo          = 1.15 ;      // ECal TdR, ~ solenoid before ECal, Section 1.4.1 p.27, https://inspirehep.net/literature/432394
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.032;      
        Str_LHCrun      = "HL-LHC run IV(+V?)";
        StrBarrelTOFname= "#kern[-0.2]{bTOF}^{#it{Hyp.}}";        
        if(TMath::Abs(rZendcapTOF - 1.5) < 1e-13)   StrEndcapTOFname= "#kern[-0.2]{eTOF}^{#it{Hyp.}}"; 
        else                                        StrEndcapTOFname= "HGTD"; 
        
    }
    else if( Str_Exp == "ALICE-1" ){
        if( TMath::Abs(rBmag) > 1e-13 )  lBz = rBmag; // existing field map for B(L3) = 0.2 T
        else                             lBz = -0.5  ; // default L3 0.5 T
        kActivateBarrel = kTRUE;
        lRbarrelTOF     = 3.8   ;     // Ravg(TOF) = 380 cm
        lEtaMaxBarrel   = 0.88  ;
        kActivateEndcap = kFALSE;
        lZendcapTOF     = NAN;        // TOF endcap inexistent... dummy, substitue by a double NaN C++11 (1.77 m = 3rd station of muon spectrometer)
        lEtaMinEndcap   = 2.5   ;
        lEtaMaxEndcap   = 4.0   ;
        lEtaMaxCalo     = 0.7   ;     // p.8 in EmCal TDR, https://cds.cern.ch/record/1121574
        lRCalo          = 4.25  ;     // Approximate value... from Fig 8.1 p.94 in EmCal TDR, https://cds.cern.ch/record/1121574
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.080 ;     // 56 ps in Pb-Pb run II to be checked with Bologna about pp MB ? See https://arxiv.org/abs/1806.03825        80 ps Run I, https://arxiv.org/pdf/1402.4476.pdf, p.51
        Str_LHCrun      = "LHC run II";
        StrBarrelTOFname= "TOF";
        StrEndcapTOFname= "None";        
    }
    else if( Str_Exp == "ALICE-3" ){
        if( TMath::Abs(rBmag) > 1e-13 )  lBz = rBmag; // existing field map for B(L3) = 0.5 or 0.2 T
        else                             lBz = -0.5  ; // default L3 0.5 T
        kActivateBarrel = kTRUE;
        if(rInnerRadiusBarrel > 0)      lRbarrelTOF = rInnerRadiusBarrel;    // Ravg(TOF) = 0.2 m, 0.3 m  inner SPAD layer = to be tested
        else                            lRbarrelTOF = 1.00;                  // Ravg(TOF) ~ 100 cm    outer SPAD layer
        lEtaMaxBarrel   = 1.4  ;   
        kActivateEndcap = kTRUE;
        lZendcapTOF     = 2.00 ;      // 2.0 m for ~std config, 
        lEtaMinEndcap   = 1.5  ;      // hermetic acceptance : connect to 1.4 barrel  or  µ spectro acceptance, 2.5 - 4.0 idea : dipole to be used ?
        lEtaMaxEndcap   = 4.0  ;
        lEtaMaxCalo     = 1.4  ;      // likely in line with the tracker extension
        lRCalo          = 1.15 ;      // approximate location of the Shower Pixel detector, https://arxiv.org/abs/1902.01211
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.02 ;      // hyp. SPAD resolution
        Str_LHCrun      = "HL-LHC run V";
        StrBarrelTOFname= "#kern[-0.2]{bTOF}^{#it{Hyp.}}";
        StrEndcapTOFname= "#kern[-0.2]{eTOF}^{#it{Hyp.}}";  
    }
    else{Printf("Exp case not foreseen... return without doing anything !"); return; }

    
    lZmaxBarrel     =  lRbarrelTOF * TMath::SinH(lEtaMaxBarrel) ;
    lRminEndcap     =  lZendcapTOF / TMath::SinH(lEtaMaxEndcap) ;   // CMS = 0.315;     //  CERN-LHCC-2019-003, p.14 Table 1.3 (Section 1.4.2)
    lRmaxEndcap     =  lZendcapTOF / TMath::SinH(lEtaMinEndcap) ;   // CMS = 1.27;      //  CERN-LHCC-2019-003, p.14 Table 1.3 (Section 1.4.2)

    
    
    if( lRbarrelTOF < lRminEndcap ){
            Printf("Weird detector configuration... Rtof[Barrel] below Rmin[Endcap] ?! to be re-checked... the code logic may fail under these circumstances. Exit for now ?!");
            // return;
    }
    
    if( lRbarrelTOF > lRCalo ){
        Printf("Fishy detector configuration... Rtof[Barrel] > RCalo ! unphysical. Exit...");
            return;
    }
    
    
    Printf("-- General settings ---------------------");        
    Printf("    Magnetic field B = %4.2f T", lBz            );
    Printf("    R[barrel](TOF)   = %4.2f m", lRbarrelTOF      );
    Printf("    Max eta[Barrel]  = %4.2f",   lEtaMaxBarrel    );
    Printf("    Max zBarrel      = %5.3f m", lZmaxBarrel      );
    Printf("    -");
    Printf("    zEndcap(TOF)     = %4.2f m", lZendcapTOF      );
    Printf("    Min eta[Endcap]  = %4.2f",   lEtaMinEndcap    );
    Printf("    Max eta[Endcap]  = %4.2f",   lEtaMaxEndcap    );
    Printf("    Min R[endcap]    = %5.3f m",   lRminEndcap      );
    Printf("    Max R[endcap]    = %5.3f m",   lRmaxEndcap      );
    Printf("    -");
    Printf("    TOF Timing Resol = %02f ps", cDeltaT*1e3      );
    Printf("    data period      : %s",      Str_LHCrun.Data()       );
    Printf("    Barrel TOF name  : %s",      StrBarrelTOFname.Data() );
    Printf("    Endcap TOF name  : %s",      StrEndcapTOFname.Data() );

  
    
        
   
    
    
    
      const Int_t lNbCanvas = 8;  // 8 canvas
      
        TString Str_SepType = "";
      
    
      // Display N lines of 3 columns of canvas
for(Int_t iCan = 0; iCan < lNbCanvas; iCan++){
    
    
        myCanvas[iCan] = new TCanvas( Form("myCanvas%02d", iCan ), Form("Canvas %02d - [%i %s %i]",iCan, iDay, cMonth[iMonth-1].Data(), iYear), (iCan%3)*650, TMath::Floor(iCan/3)*250, 700, 600 );
        myPadSetUp(myCanvas[iCan],kFALSE, 0.15,0.04, 0.04,0.16);
        myCanvas[iCan]->cd();
        myCanvas[iCan]->Draw();
        
        myPad[iCan] = new TPad( Form("myPad%02d", iCan), Form("pad - %02d", iCan),0,0,1,1);
        myPadSetUp(myPad[iCan],kFALSE, 0.135,0.04,0.145,0.15);  
        myPad[iCan]->Draw();      
        
        myPad[iCan]->Update();
        myCanvas[iCan]->ToggleEventStatus();
        myCanvas[iCan]->ToggleEditor();
        myCanvas[iCan]->ToggleToolBar();

    UInt_t lPart1 = -1;
    UInt_t lPart2 = -1; 
    
    double q1 = 0.;
    double q2 = 0.;
    
    switch(iCan){
        case  0 :   Str_SepType = "e-vs-pi"     ;   lPart1 = kPDG_piPm;      q1 = 1.;       lPart2 = kPDG_ePm;      q2 = 1.;    break;
        case  1 :   Str_SepType = "mu-vs-pi"    ;   lPart1 = kPDG_piPm;      q1 = 1.;       lPart2 = kPDG_muPm;     q2 = 1.;    break;
        case  2 :   Str_SepType = "pi-vs-K"     ;   lPart1 = kPDG_Kpm;       q1 = 1.;       lPart2 = kPDG_piPm;     q2 = 1.;    break;
        case  3 :   Str_SepType = "K-vs-p"      ;   lPart1 = kPDG_p;         q1 = 1.;       lPart2 = kPDG_Kpm;      q2 = 1.;    break;
        case  4 :   Str_SepType = "p-vs-2H"     ;   lPart1 = kPDG_2H;        q1 = 1.;       lPart2 = kPDG_p;        q2 = 1.;    break;
        case  5 :   Str_SepType = "2H-vs-3H"    ;   lPart1 = kPDG_3H;        q1 = 1.;       lPart2 = kPDG_2H;       q2 = 1.;    break;
        case  6 :   Str_SepType = "3H-vs-3He"   ;   lPart1 = kPDG_3He;       q1 = 2.;       lPart2 = kPDG_3H;       q2 = 1.;    break;
        case  7 :   Str_SepType = "3He-vs-4He"  ;   lPart1 = kPDG_4He;       q1 = 2.;       lPart2 = kPDG_3He;      q2 = 2.;    break;
                                    
        case 10 :   Str_SepType = "pi-vs-e"     ;   lPart2 = kPDG_piPm;      q2 = 1.;       lPart1 = kPDG_ePm;      q1 = 1.;   break;
        case 11 :   Str_SepType = "pi-vs-µ"     ;   lPart2 = kPDG_piPm;      q2 = 1.;       lPart1 = kPDG_muPm;     q1 = 1.;   break;
        case 12 :   Str_SepType = "K-vs-pi"     ;   lPart2 = kPDG_Kpm;       q2 = 1.;       lPart1 = kPDG_piPm;     q1 = 1.;   break;
        case 13 :   Str_SepType = "p-vs-K"      ;   lPart2 = kPDG_p;         q2 = 1.;       lPart1 = kPDG_Kpm;      q1 = 1.;   break;
        case 14 :   Str_SepType = "2H-vs-p"     ;   lPart2 = kPDG_2H;        q2 = 1.;       lPart1 = kPDG_p;        q1 = 1.;   break;
        case 15 :   Str_SepType = "3H-vs-2H"    ;   lPart2 = kPDG_3H;        q2 = 1.;       lPart1 = kPDG_2H;       q1 = 1.;   break;
        case 16 :   Str_SepType = "3He-vs-3H"   ;   lPart2 = kPDG_3He;       q2 = 2.;       lPart1 = kPDG_3H;       q1 = 1.;   break;
        case 17 :   Str_SepType = "4He-vs-3He"  ;   lPart2 = kPDG_4He;       q2 = 2.;       lPart1 = kPDG_3He;      q1 = 2.;   break;
        default :   Printf("Not foreseen PID separation case (mind the syntax ?) : %s... exit!", Str_SepType.Data() );  return;
    }
        
    
    
    if(iCan > lNbCanvas) { Printf("Not enough TCanvas are booked... exit !"); return; }
    myPad[ iCan ]->cd();
    
    Printf(" ");
    Printf("------- %02d ---", iCan);
    Printf("Chosen particle separation : m[%s] = %6.5f GeV/c² Vs m[%s] = %6.5f GeV/c²", Str_HadronASCIIName[lPart1].Data(), mPDG[lPart1], Str_HadronASCIIName[lPart2].Data(), mPDG[lPart2] );
  
    
    

    
    Double_t lMaxPtHist = -1.0;
    Double_t lMinPtHist = -1.0;
    Double_t lMaxPtEta  = -1.0;    
    
    if(kLog) {
        myPad[ iCan ]->SetLogy();
        lMaxPtHist = 100; // 50
        lMinPtHist = 0.009;
    }
    else{// linear
        lMaxPtHist = 11.;
        lMinPtHist = -0.4;    
    }
    
    myBlankHisto[ iCan ] = new TH1F( Form("myBlankHisto[%02d]", iCan),"Blank Histogram",   440,-0.4,4.0);
    myBlankHisto[ iCan ]->SetMaximum(lMaxPtHist);
    myBlankHisto[ iCan ]->SetMinimum(lMinPtHist);
    myBlankHisto[ iCan ]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
    myBlankHisto[ iCan ]->SetXTitle( Form("#it{y}_{%s} #scale[0.7]{#scale[0.7]{, particle rapidity}}",  Str_HadronLatexName[lPart1].Data() ) );
    myBlankHisto[ iCan ]->SetNdivisions(515,"y");
    myBlankHisto[ iCan ]->Draw();
    
    
    
    
    
    
    
    
    
    
    
    

        
        
    double clight = 0.299792458;      // = kB2C in AliRoot, [m.s^-1.V^-1]
    double qElem  = 1.602176634e-19;  // in Coulomb 

    
    double m1 = mPDG[lPart1];  // in GeV/c²
    double m2 = mPDG[lPart2];
    
    double m1InKg = m1 * 1.672649e-27 / mPDG[kPDG_p];  // in kg ! =  needed for the cyclotron frequency, use proton mass as reference in kg
    double m2InKg = m2 * 1.672649e-27 / mPDG[kPDG_p];
    
    
    
    double n5sigma = 5.;
    double n4sigma = 4.;
    double n3sigma = 3.;
    double n2sigma = 2.;
    double n1sigma = 1.;
    //  double p_dedx = 1.05; // pi/K separation in CMS with Si dE/dx
    double p_dedx = 1.7; // pT max in CMS for p/K separation in CMS with Si dE/dx = 1.7 GeV/c

    
    // gradual binning in pT
    double ptmax   = 26.0;
  
    std::vector<Double_t> vPtBinning;  
    
    for(Int_t iPt = -1; iPt <  90; iPt++){  vPtBinning.push_back( 0.01 + iPt * 0.001 );  } // from  0.009 to  0.1 GeV/c
    for(Int_t iPt =  0; iPt < 180; iPt++){  vPtBinning.push_back( 0.1  + iPt * 0.005 );  } // from  0.1   to  1   GeV/c
    for(Int_t iPt =  0; iPt < 450; iPt++){  vPtBinning.push_back( 1.   + iPt * 0.02  );  } // from  1.0   to 10   GeV/c
    for(Int_t iPt =  0; iPt <  81; iPt++){  vPtBinning.push_back(10.   + iPt * 0.2   );  } // from 10     to 26   GeV/c

    
    std::sort( vPtBinning.begin(), vPtBinning.end() );
    int    npt     = vPtBinning.size() - 1;  // 801
    // for(Int_t i = 0; i < npt; i++)  Printf("pT[%03d] = %06.3f GeV/c", i, vPtBinning[i]); 
    
    

    
        
        
    
    
    // Regular binning in eta
    Double_t etamax  = lEtaMaxEndcap;
    Double_t etastep = 0.005;      
    Int_t    neta    = etamax / etastep;
    
    Printf(" - ptstep  = from 0.001 to 0.2 GeV/c");
    Printf(" - npt     = %03d", npt    );   // 801 
    Printf(" - etastep = %6.4g", etastep);
    Printf(" - neta    = %03d", neta   );
    


    
    TH2I* hDummyMap  = new TH2I( Form("hDummyMap[%02d]", iCan),  "hDummyMap",    neta, 0,etamax,  npt, vPtBinning.data() );
    
    // TF1* funcEta2Rap = new TF1("funcEta2Rap","log((sqrt([0]*[0]+x*x*cosh([1])*cosh([1]))+x*sinh([1]))/sqrt([0]*[0]+x*x))",0,100);
    // NOTE :  y = f(x=pT, param = m0,eta0), see https://en.wikipedia.org/wiki/Pseudorapidity
    TF1* fnEta2Rap = new TF1("fnEta2Rap", funcEta2Rap,   0,10,  2);
   
    
    
  
  
    // Determination of pT thresholds, in barrel and endcap
    Double_t lPtMinBarrel = lRbarrelTOF/2                              * clight * TMath::Abs(lBz* q1);
    Double_t lPtMinEndcap = lZendcapTOF/(2*TMath::SinH(lEtaMaxEndcap)) * clight * TMath::Abs(lBz* q1);
    Double_t lPtMinCalo   = lRCalo/2                                   * clight * TMath::Abs(lBz* q1); 
    
    Bool_t kRelaxAngConstraints = kFALSE; // relax (lambda angle, delta angle) = turned down  / keep calo constraints

    
    
    
    
    // Minimal inclination angle for endcap      
    Int_t iBinUpperPtForEtaMinLambda = hDummyMap->GetYaxis()->FindBin ( lPtMinCalo );   // value which does a decent job to fix the upper idx point of grEtaMinLambda
    
    grEtaMinLambda[ iCan ] = new TGraph(iBinUpperPtForEtaMinLambda);    

    Double_t lLambdaMinAngle = 0;
    if(kRelaxAngConstraints)    lLambdaMinAngle =  5 * TMath::DegToRad() ; //  5° taken as minimal inclination angle for endcap = lambda_min loose
    else                        lLambdaMinAngle = 30 * TMath::DegToRad() ; // 30° taken as minimal inclination angle for endcap = lambda_min   FIXME
    
    
    Double_t EtaMinLambda = TMath::ASinH( TMath::Tan( lLambdaMinAngle )   );  // eq 48a
    
    
    fnEta2Rap->SetParameters(m1, EtaMinLambda);   
    for(int i=0;i<iBinUpperPtForEtaMinLambda;i++){
        
        grEtaMinLambda[ iCan ]->SetPoint(i, fnEta2Rap->Eval(hDummyMap->GetYaxis()->GetBinCenter(i+1)), 
                                            hDummyMap->GetYaxis()->GetBinCenter(i+1) );
        
    }
    
    
    Printf("Acceptance terms :");
    Printf(" - eta_min[lambda] for loopers = %5.3g", EtaMinLambda);
    
     
     
    // Minimal inclination angle for barrel     
    Double_t lDeltaMinAngle = 0;
    if(kRelaxAngConstraints)    lDeltaMinAngle = 10 * TMath::DegToRad(); //  5° taken as minimal inclination angle for barrel = delta_min loose 
    else                        lDeltaMinAngle = 30 * TMath::DegToRad(); // 30° taken as minimal inclination angle for barrel = delta_min  FIXME
    
    
    
    Double_t EtaMaxDelta    = TMath::ACosH(1./TMath::Sin( lDeltaMinAngle )  ); // corresponding max eta value, eq 46a

    
    Double_t lPtMinBarrelDelta  = 0.;   
    Int_t    iEta               = 0 ;
    Double_t iEtaVal            = 0.; 
    Double_t lEtaStepDelta1     = 0.01;
    
    
    TF1 *fnPtMinBarrelDelta = new TF1("fnPtMinBarrelDelta", funcPtMinBarrelDelta, 0, EtaMaxDelta, 2); // eq 45, PtMinBarrelDelta = f(eta)
    
    fnPtMinBarrelDelta->SetParameter( 0, lPtMinBarrel  );
    fnPtMinBarrelDelta->SetParameter( 1, lDeltaMinAngle);
    
    myFuncSetUp(fnPtMinBarrelDelta, kAzure+1, 7, 2);

    
    
    
    std::vector<Double_t>   vPtMinBarrelDelta;
    std::vector<Double_t>   vYBarrelDelta;

        // NOTE we do not know a priori what will be the nb of point needed to describe grEtaMinDelta...
        //          Issue = the steep rising when eta ~< EtaMaxDelta (vertical asymptote)
        //          1. Need a variating binning as a function of the steepness of the function 

    Double_t lDerivative = 0;
    while ( lDerivative <  15.0 ){            
        iEta++;        
                    iEtaVal     = iEta * lEtaStepDelta1;    
                    lDerivative = fnPtMinBarrelDelta->Derivative(iEtaVal);
        Double_t    lPt         = fnPtMinBarrelDelta->Eval( iEtaVal );
    
        fnEta2Rap->SetParameters(m1, iEtaVal);
        Double_t    lY          = fnEta2Rap      ->Eval(lPt);
        
        if(lPt > 9 || TMath::IsNaN(lPt) || TMath::IsNaN(lY) ) break;
        
        vPtMinBarrelDelta.push_back( lPt );
        vYBarrelDelta    .push_back( lY  );

        if(lDebug > 3) Printf("[%03d]  : lY = %8.6f, lPt = %8.6f  - (iEtaVal = %9.7f, lDerivative = %06.4g)", iEta, lY, lPt, iEtaVal, lDerivative);
    }
    
    Int_t     iEta_LastIn1stPart    = iEta;
    Double_t  iEtaVal_LastIn1stPart = iEtaVal;
    
    
    Int_t iEtaBis = 0;
    iEtaVal = EtaMaxDelta; // starts from the asymptotic value
    
    while(iEtaVal > iEtaVal_LastIn1stPart){// NOTE starts in reverse order, from the asymptote and get down with step that evolves quadratically, tiny at first then larger and larger step
        iEtaBis++;
                    iEtaVal     = EtaMaxDelta - iEtaBis*iEtaBis * 1e-5;    
                    lDerivative = fnPtMinBarrelDelta->Derivative(iEtaVal);
        Double_t    lPt         = fnPtMinBarrelDelta->Eval( iEtaVal );
        
        fnEta2Rap->SetParameters(m1, iEtaVal);
        Double_t    lY          = fnEta2Rap      ->Eval(lPt);
        
        if(lPt > 9 || TMath::IsNaN(lPt) || TMath::IsNaN(lY) ) continue;
        
        if(lDebug > 3) Printf("[%03d]' : lY = %8.6f, lPt = %8.6f  - (iEtaVal = %9.7f, EtaMaxDelta - iEtaVal = %9.7f, lDerivative = %06.4g)", iEta_LastIn1stPart+iEtaBis, lY, lPt, iEtaVal, EtaMaxDelta - iEtaVal, lDerivative);
        
        vPtMinBarrelDelta.push_back( lPt );
        vYBarrelDelta    .push_back( lY  );
    }
    
    
    std::sort( vPtMinBarrelDelta.begin(), vPtMinBarrelDelta.end() ); // makes sense because the function is strictly monotonous, default by increasing order
    std::sort( vYBarrelDelta    .begin(), vYBarrelDelta    .end() );
    
    
    Int_t lNbEtaDelta = vYBarrelDelta    .size();
    
    grEtaMinDelta[ iCan ] = new TGraph(lNbEtaDelta);

    for(Int_t iPnt = 0; iPnt < lNbEtaDelta; ++iPnt){       
        grEtaMinDelta[ iCan ]->SetPoint( iPnt, vYBarrelDelta[ iPnt ], vPtMinBarrelDelta[ iPnt ]);
    }

  
     
     Printf(" - eta_max[delta]      = %8.6f", EtaMaxDelta); 
     Printf("   NbEtaBins[delta]    = %03d",  lNbEtaDelta);
     Printf("   lEtaStepDelta1      = %6.4f", lEtaStepDelta1);
     Printf("   pTmin[Barrel]_delta = %5.3g to %5.3g GeV/c", vPtMinBarrelDelta[0], vPtMinBarrelDelta[ lNbEtaDelta-1 ] );
     

    
    
    // Determination of pT thresholds, in barrel and endcap (part2)

    fnEta2Rap->SetParameters(m1, lEtaMaxBarrel); // find the corresponding rapidity, ymax, for the upper edge of acceptance in eta along the horizontal line
    Double_t lYMaxForPtMinBarrel = fnEta2Rap->Eval(lPtMinBarrel);    
    
    fnEta2Rap->SetParameters(m1, lEtaMinEndcap);
    Double_t lYMinForPtMinEndcap = fnEta2Rap->Eval(lPtMinEndcap);

    fnEta2Rap->SetParameters(m1, lEtaMaxEndcap);
    Double_t lYMaxForPtMinEndcap = fnEta2Rap->Eval(lPtMinEndcap); 
     
     
    fnEta2Rap->SetParameters(m1, EtaMinLambda); 
    Float_t lYMinLooper = fnEta2Rap->Eval( lPtMinEndcap );
    // locate the crossing between horizontal straight line (lPtMinEndcap) and the curled vertical limit of grEtaMinLambda
     
     
     
    grPtMinBarrel[ iCan ] = new TGraph(2); 
    grPtMinEndcap[ iCan ] = new TGraph(2); 

    grPtMinBarrel[ iCan ]->SetPoint(0,  0,                   lPtMinBarrel);
    grPtMinBarrel[ iCan ]->SetPoint(1,  lYMaxForPtMinBarrel, lPtMinBarrel);

    grPtMinEndcap[ iCan ]->SetPoint(0,  lYMinLooper,         lPtMinEndcap);
    grPtMinEndcap[ iCan ]->SetPoint(1,  lYMaxForPtMinEndcap, lPtMinEndcap);  

    
    Printf(" - pT,min[Barrel]      = %5.3g GeV/c", lPtMinBarrel);
    Printf(" - pT,min[Endcap]      = %5.3g GeV/c", lPtMinEndcap);    

   
    // Barrel Calorimeter or solenoïd threshold
    // NOTE : the calorimeter, even if barely touched, can be seen as the end of looper track (Eloss ++)
    //         meaning that the redundancy between barrel TOF and endcap late hit of looper is forbidden 
    
    
    fnEta2Rap->SetParameters(m1, lEtaMaxCalo); // find the corresponding rapidity, ymax, for the upper edge of acceptance in eta along the horizontal line
    Double_t lYMaxForPtMinCalo = fnEta2Rap->Eval(lPtMinCalo);    
            // NOTE : lEtaMaxCalo does not necessarily coincide with lEtaMaxBarrel for TOF... TOF edge can be shorter, e.g. ATLAS
    
    fnEta2Rap->SetParameters(m1, EtaMinLambda); 
    Float_t lYMinForPtMinCalo = fnEta2Rap->Eval( lPtMinCalo );
    // locate the crossing between horizontal straight line (lPtMinCalo) and the curled vertical limit of grEtaMinLambda
    
    grPtMinCalo[ iCan ] = new TGraph(2); 
    grPtMinCalo[ iCan ]->SetPoint(0,  lYMinForPtMinCalo, lPtMinCalo);
    grPtMinCalo[ iCan ]->SetPoint(1,  lYMaxForPtMinCalo, lPtMinCalo);

    Printf(" - pT,min[Barrel Calo] = %5.3g GeV/c", lPtMinCalo); 
  
  
    
    Printf(" "); 
    Printf("- Histogram booking..."); 
    
  Int_t nEtaEndcap =  ( (lEtaMaxEndcap - lYMinLooper        ) / etastep);  // lEtaMaxEndcap - lYMinForPtMinEndcap
  Int_t nEtaLooper =  ( (lEtaMinEndcap - lYMinLooper        ) / etastep);
  Int_t nEtaBarrel =  ( (lEtaMaxBarrel - 0                  ) / etastep);
    //   Int_t nPtEndcap  =  ( (ptmax         - lPtMinEndcap       ) / ptstep );
    //   Int_t nPtLooper  =  ( (ptmax         - lPtMinEndcap       ) / ptstep ); 
    //   Int_t nPtBarrel  =  ( (ptmax         - lPtMinBarrel       ) / ptstep );

   
  
  Printf("   Reminder : etastep = %6.4g", etastep);   
  
  Printf("   EtaEndcap : Nb[bins] =  %04d / y(min)  : %4.2f - y(max)  : %05.2f", nEtaEndcap, lYMinLooper       , lEtaMaxEndcap  );
  Printf("   EtaLooper : Nb[bins] =  %04d / y(min)  : %4.2f - y(max)  : %05.2f", nEtaLooper, lYMinLooper       , lEtaMinEndcap  );
  Printf("   EtaBarrel : Nb[bins] =  %04d / y(min)  : %4.2f - y(max)  : %05.2f", nEtaBarrel, 0.                , lEtaMaxBarrel  );
    //   Printf("   PtEndcap  : Nb[bins] =  %04d / pT(min) : %4.2f - pT(max) : %05.2f", nPtEndcap , lPtMinEndcap      , ptmax          );
    //   Printf("   PtLooper  : Nb[bins] =  %04d / pT(min) : %4.2f - pT(max) : %05.2f", nPtLooper , lPtMinEndcap      , ptmax          );
    //   Printf("   PtBarrel  : Nb[bins] =  %04d / pT(min) : %4.2f - pT(max) : %05.2f", nPtBarrel , lPtMinBarrel      , ptmax          );
  Printf(" ");                             
  
   if ( lPtMinBarrel  - lPtMinEndcap < 0 ){
    Printf("Weird detector config... pT,min(barrel) < pT,min(endcap) !?  Exit for now...");
    Printf("Likely, TOF barrel at small radius (close to beam pipe), while endcap really far away in z. Better to be rechecked.");
    // return;       
    }
   
    // NOTE : beware we quote here eta values, but this is in fact a map in y
    //          So that this is "y" values that are read for the histo booking
    //          Issue : for any eta, y < eta
    //          To be dealt with properly, especially for lower bound in endcap
    //          We need to enlarge the rectangle window to make sure to encompass the necessary accetance (pT,y) for each case
    
        //   h5sigm_ecap[ iCan ]     = new TH2I( Form("h5sigm_ecap[%02d]",   iCan),  "h5sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  nPtEndcap, lPtMinEndcap,ptmax);  
        //   h4sigm_ecap[ iCan ]     = new TH2I( Form("h4sigm_ecap[%02d]",   iCan),  "h4sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  nPtEndcap, lPtMinEndcap,ptmax);
        //   h3sigm_ecap[ iCan ]     = new TH2I( Form("h3sigm_ecap[%02d]",   iCan),  "h3sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  nPtEndcap, lPtMinEndcap,ptmax);
        //   h2sigm_ecap[ iCan ]     = new TH2I( Form("h2sigm_ecap[%02d]",   iCan),  "h2sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  nPtEndcap, lPtMinEndcap,ptmax);
        //   h1sigm_ecap[ iCan ]     = new TH2I( Form("h1sigm_ecap[%02d]",   iCan),  "h1sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  nPtEndcap, lPtMinEndcap,ptmax);
        //   
        //   h5sigm_loop[ iCan ]     = new TH2I( Form("h5sigm_loop[%02d]",   iCan),  "h5sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,   nPtLooper, lPtMinEndcap,ptmax);  
        //   h4sigm_loop[ iCan ]     = new TH2I( Form("h4sigm_loop[%02d]",   iCan),  "h4sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,   nPtLooper, lPtMinEndcap,ptmax); 
        //   h3sigm_loop[ iCan ]     = new TH2I( Form("h3sigm_loop[%02d]",   iCan),  "h3sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,   nPtLooper, lPtMinEndcap,ptmax); 
        //   h2sigm_loop[ iCan ]     = new TH2I( Form("h2sigm_loop[%02d]",   iCan),  "h2sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,   nPtLooper, lPtMinEndcap,ptmax); 
        //   h1sigm_loop[ iCan ]     = new TH2I( Form("h1sigm_loop[%02d]",   iCan),  "h1sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,   nPtLooper, lPtMinEndcap,ptmax); 
        // 
        //   hdEdx_CMS[ iCan ]       = new TH2I( Form("hdEdx_CMS[%02d]",     iCan),  "hdEdx_CMS",        nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);
        //   h5sigm_barrel[ iCan ]   = new TH2I( Form("h5sigm_barrel[%02d]", iCan),  "h5sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);  
        //   h4sigm_barrel[ iCan ]   = new TH2I( Form("h4sigm_barrel[%02d]", iCan),  "h4sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);
        //   h3sigm_barrel[ iCan ]   = new TH2I( Form("h3sigm_barrel[%02d]", iCan),  "h3sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);
        //   h2sigm_barrel[ iCan ]   = new TH2I( Form("h2sigm_barrel[%02d]", iCan),  "h2sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);
        //   h1sigm_barrel[ iCan ]   = new TH2I( Form("h1sigm_barrel[%02d]", iCan),  "h1sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,   nPtBarrel, lPtMinBarrel,ptmax);
  
  h5sigm_ecap[ iCan ]     = new TH2I( Form("h5sigm_ecap[%02d]",   iCan),  "h5sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  npt, vPtBinning.data() );  
  h4sigm_ecap[ iCan ]     = new TH2I( Form("h4sigm_ecap[%02d]",   iCan),  "h4sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  npt, vPtBinning.data() );
  h3sigm_ecap[ iCan ]     = new TH2I( Form("h3sigm_ecap[%02d]",   iCan),  "h3sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  npt, vPtBinning.data() );
  h2sigm_ecap[ iCan ]     = new TH2I( Form("h2sigm_ecap[%02d]",   iCan),  "h2sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  npt, vPtBinning.data() );
  h1sigm_ecap[ iCan ]     = new TH2I( Form("h1sigm_ecap[%02d]",   iCan),  "h1sigm_ecap",      nEtaEndcap, lYMinLooper,lEtaMaxEndcap,  npt, vPtBinning.data() );
                                                                                                                                       
  h5sigm_loop[ iCan ]     = new TH2I( Form("h5sigm_loop[%02d]",   iCan),  "h5sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,  npt, vPtBinning.data() );  
  h4sigm_loop[ iCan ]     = new TH2I( Form("h4sigm_loop[%02d]",   iCan),  "h4sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,  npt, vPtBinning.data() ); 
  h3sigm_loop[ iCan ]     = new TH2I( Form("h3sigm_loop[%02d]",   iCan),  "h3sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,  npt, vPtBinning.data() ); 
  h2sigm_loop[ iCan ]     = new TH2I( Form("h2sigm_loop[%02d]",   iCan),  "h2sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,  npt, vPtBinning.data() ); 
  h1sigm_loop[ iCan ]     = new TH2I( Form("h1sigm_loop[%02d]",   iCan),  "h1sigm_loop",      nEtaLooper, lYMinLooper,lEtaMinEndcap,  npt, vPtBinning.data() ); 
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
  hdEdx_CMS[ iCan ]       = new TH2I( Form("hdEdx_CMS[%02d]",     iCan),  "hdEdx_CMS",        nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );
  h5sigm_barrel[ iCan ]   = new TH2I( Form("h5sigm_barrel[%02d]", iCan),  "h5sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );  
  h4sigm_barrel[ iCan ]   = new TH2I( Form("h4sigm_barrel[%02d]", iCan),  "h4sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );
  h3sigm_barrel[ iCan ]   = new TH2I( Form("h3sigm_barrel[%02d]", iCan),  "h3sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );
  h2sigm_barrel[ iCan ]   = new TH2I( Form("h2sigm_barrel[%02d]", iCan),  "h2sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );
  h1sigm_barrel[ iCan ]   = new TH2I( Form("h1sigm_barrel[%02d]", iCan),  "h1sigm_barrel",    nEtaBarrel, 0,          lEtaMaxBarrel,  npt, vPtBinning.data() );


    grEta0pt5[ iCan ] = new TGraph( (npt-0 )/10 );  // < 26 GeV/c
    grEta1pt0[ iCan ] = new TGraph( (npt-30)/10 );  // < 20 GeV/c
    grEta1pt5[ iCan ] = new TGraph( (npt-40)/10 );  // < 18 GeV/c
    grEta2pt0[ iCan ] = new TGraph( (npt-50)/10 );  // < 17 GeV/c
    grEta2pt5[ iCan ] = new TGraph( (npt-70)/10 );  // < 15 GeV/c
    grEta3pt0[ iCan ] = new TGraph( (npt-70)/10 );  // < 15 GeV/c
    grEta3pt5[ iCan ] = new TGraph( (npt-70)/10 );  // < 15 GeV/c = (801 -70)/20 = 731/20 = 36
    
        fnEta2Rap->SetParameters(m1,0.5);   for(int i=0; i < npt-0 ; i++){ if(i%10) continue; grEta0pt5[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }
        fnEta2Rap->SetParameters(m1,1.0);   for(int i=0; i < npt-30; i++){ if(i%10) continue; grEta1pt0[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }       
        fnEta2Rap->SetParameters(m1,1.5);   for(int i=0; i < npt-40; i++){ if(i%10) continue; grEta1pt5[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }
        fnEta2Rap->SetParameters(m1,2.0);   for(int i=0; i < npt-50; i++){ if(i%10) continue; grEta2pt0[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }
        fnEta2Rap->SetParameters(m1,2.5);   for(int i=0; i < npt-70; i++){ if(i%10) continue; grEta2pt5[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }
        fnEta2Rap->SetParameters(m1,3.0);   for(int i=0; i < npt-70; i++){ if(i%10) continue; grEta3pt0[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }
        fnEta2Rap->SetParameters(m1,3.5);   for(int i=0; i < npt-70; i++){ if(i%10) continue; grEta3pt5[ iCan ]->SetPoint((Int_t)i/10, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) ); }       


    grEtaMaxBarrel [ iCan ] = new TGraph(npt-70);       
    grEtaMinEndcap [ iCan ] = new TGraph(npt-70);       
    grEtaMaxEndcap [ iCan ] = new TGraph(npt-70);       
    
        fnEta2Rap->SetParameters(m1, lEtaMaxBarrel);   for(int i=0;i<npt-70;i++) grEtaMaxBarrel[ iCan ]->SetPoint(i, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) );
        fnEta2Rap->SetParameters(m1, lEtaMinEndcap);   for(int i=0;i<npt-70;i++) grEtaMinEndcap[ iCan ]->SetPoint(i, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) );    
        fnEta2Rap->SetParameters(m1, lEtaMaxEndcap);   for(int i=0;i<npt-70;i++) grEtaMaxEndcap[ iCan ]->SetPoint(i, fnEta2Rap->Eval( hDummyMap->GetYaxis()->GetBinCenter(i+1)), hDummyMap->GetYaxis()->GetBinCenter(i+1) );
        
  
    


     
    
    
  double L1=0;
  double L2=0;
  Double_t cT_diff = 0;
  Double_t lTransvRotAngle = 0;
  UShort_t lAbsNbHalfTurns_Barrel1 = 0;
  UShort_t lAbsNbHalfTurns_Endcap1 = 0;
  
  Int_t iDbleAcceptance = 0;
  
  
  
    Bool_t  kAbovePtMinEndcap             = kFALSE;        
    Bool_t  kAbovePtMinBarrel             = kFALSE; 
    Bool_t  kAbovePtMinCalo               = kFALSE;
    Bool_t  kAboveDeltaMinAngle           = kFALSE;
    Bool_t  kAboveLambdaMinAngle          = kFALSE;
    Bool_t  kLooper                       = kFALSE; 
    Bool_t  kWithinRZbarrelAccptce        = kFALSE;
    Bool_t  kWithinRZendcapAccptce        = kFALSE;
    Bool_t  IsPartForeseenInBarrel        = kFALSE;
    Bool_t  IsPartForeseenInEndcap        = kFALSE;
    Bool_t  IsPartForeseenBeyondEndcap    = kFALSE;
  
  
  
  for(Int_t iPt= 0; iPt < npt ; iPt++)
  {
    double pt = vPtBinning[ iPt ];
    
    for(Int_t iEta = 0; iEta < neta; iEta++)
    {
        
        // Reset booleans to start from clean ground at every sub-loop turn
        kAbovePtMinEndcap             = kFALSE;
        kAbovePtMinBarrel             = kFALSE;
        kAbovePtMinCalo               = kFALSE;
        kAboveDeltaMinAngle           = kFALSE;
        kAboveLambdaMinAngle          = kFALSE;
        kLooper                       = kFALSE;
        kWithinRZbarrelAccptce        = kFALSE;
        kWithinRZendcapAccptce        = kFALSE;
        IsPartForeseenInBarrel        = kFALSE;
        IsPartForeseenInEndcap        = kFALSE;
        IsPartForeseenBeyondEndcap    = kFALSE;
        
        lAbsNbHalfTurns_Barrel1 = 0;
        lAbsNbHalfTurns_Endcap1 = 0;
           
    
        if(pt < lPtMinEndcap) continue; // below the endcap minimum pT, so no detection at all allowed    
        else kAbovePtMinEndcap    = kTRUE;
    
       
        if(pt > lPtMinBarrel) kAbovePtMinBarrel    = kTRUE; // above the endcap minimum pT (loopers, regular endcap, regular barrel)
    
        if(pt > lPtMinCalo)   kAbovePtMinCalo      = kTRUE; // above the reach of barrel calorimeter -> forget about double acceptance (barrel +endcap) loopers

    
    
        Double_t eta            = iEta*etastep;
        Double_t pl             = pt * TMath::SinH(eta);
        Double_t ptot           = TMath::Sqrt(pl*pl + pt*pt);
    
      
        Double_t energy1        = TMath::Sqrt(pt*pt + pl*pl + m1*m1);
        Double_t energy2        = TMath::Sqrt(pt*pt + pl*pl + m2*m2);
        Double_t rapidity1      = 0.5*log((energy1+pl)/(energy1-pl));
        Double_t rapidity2      = 0.5*log((energy2+pl)/(energy2-pl));
        Double_t gamma1         = TMath::Sqrt(1. + ptot*ptot/m1/m1);
        Double_t gamma2         = TMath::Sqrt(1. + ptot*ptot/m2/m2);


        if( eta > EtaMinLambda)             kAboveLambdaMinAngle = kTRUE;
        
        if( eta < EtaMaxDelta-etastep  &&  pt > fnPtMinBarrelDelta->Eval(eta)  )  kAboveDeltaMinAngle  = kTRUE;
            // NOTE : due to discretisation of eta, one needs to remove one last etastep to avoid rounding effects between EtaMaxDelta and eta = i*etastep
        
        if(      kAboveLambdaMinAngle  
            &&   kAbovePtMinEndcap 
            &&   kAbovePtMinCalo == kFALSE
            &&   eta < lEtaMinEndcap
            &&   kActivateEndcap)           kLooper = kTRUE;       
    
        

        
        
        
        // Wei's Lee code      
        Bool_t kWeiCodeSwitch = kFALSE;
        
      
        if(!kWeiCodeSwitch){// Antonin Maire's math core of the code
                   
            Double_t AngFreq1       = (q1*qElem * lBz) / (gamma1 * m1InKg); // in rad.s-1, from var with international system unit q*e (C), m (kg), etc
            Double_t AngFreq2       = (q2*qElem * lBz) / (gamma2 * m2InKg);
            
            
            Double_t rho1           = pt / (clight * lBz * q1);
            Double_t rho2           = pt / (clight * lBz * q2);

            Double_t betaZ1         = pt * TMath::SinH( eta ) / ( gamma1 * m1 );
            Double_t betaZ2         = pt * TMath::SinH( eta ) / ( gamma2 * m2 );
            

            
            Double_t lPhit0 = 0.;  // Initial azimuthal angle, by default = 0,  could be anything in [0, 2pi] : TMath::Pi()/2. , ... 
            // NOTE : ok for the trajectory
            //             but still to be cured for the other displays : centre of curvature location, tan(t0)...
            //             The physics quantities : path legnth, TOF time, straight lengths, Rt, rotation angle ... should be left unchanged by this rotation
            // NOTE the helix will be rotated in block around the z axis, like a door on its hinges, by the side
            //             so that in terms of acceptance the initial phi(t0) has NO effect, if the detector is invariant by rotation
            //             the final TOF hit will end up always in the same circle, same radius, same z, just at a different location on it.
            //          -> it has been proved here with phi(t0) = 90°, 60°, 45° , all the output looks the same

            Double_t lTmpX = 0;
            Double_t lTmpY = 0;

            Double_t ltf1 = 0;
            Double_t ltf2 = 0;
            

            
            // Math for potential barrel case
            Double_t ltf_Barrel1    = TMath::Abs(1./AngFreq1) * 2 * TMath::ASin( lRbarrelTOF/ (2. * TMath::Abs(rho1)) );
            Double_t ltf_Barrel2    = TMath::Abs(1./AngFreq2) * 2 * TMath::ASin( lRbarrelTOF/ (2. * TMath::Abs(rho2)) );
            
            Double_t lXtf_Barrel1   = -rho1 *   TMath::Sin( AngFreq1 * ltf_Barrel1 )     ;
            Double_t lYtf_Barrel1   = -rho1 * ( TMath::Cos( AngFreq1 * ltf_Barrel1 ) - 1);
            Double_t lZtf_Barrel1   = betaZ1 * TMath::C()            * ltf_Barrel1 ;
              
                     lTmpX = lXtf_Barrel1;
                     lTmpY = lYtf_Barrel1;
                     
                     lXtf_Barrel1 = lTmpX * TMath::Cos( lPhit0 ) - lTmpY * TMath::Sin( lPhit0 );  // trigo formula cos(a+b) : x = r. cos (phit) -> x' = r. cos( phit + phi0 ) -> x' = x. cos(phi0) - y. sin(phi0)
                     lYtf_Barrel1 = lTmpY * TMath::Cos( lPhit0 ) + lTmpX * TMath::Sin( lPhit0 );  // trigo formula sin(a+b) : y = r. sin (phit) -> x' = r. sin( phit + phi0 ) -> y' = y. cos(phi0) + y. sin(phi0)       

            
            Double_t lXtf_Barrel2   = -rho2 *   TMath::Sin( AngFreq2 * ltf_Barrel2 )     ;       
            Double_t lYtf_Barrel2   = -rho2 * ( TMath::Cos( AngFreq2 * ltf_Barrel2 ) - 1);        
            Double_t lZtf_Barrel2   = betaZ2 * TMath::C()            * ltf_Barrel2 ;
            
                     lTmpX = lXtf_Barrel2;
                     lTmpY = lYtf_Barrel2;
                     
                     lXtf_Barrel2 = lTmpX * TMath::Cos( lPhit0 ) - lTmpY * TMath::Sin( lPhit0 );  // trigo formula cos(a+b) 
                     lYtf_Barrel2 = lTmpY * TMath::Cos( lPhit0 ) + lTmpX * TMath::Sin( lPhit0 );  // trigo formula sin(a+b)
            
            
            Double_t lRt_tfBarrel1  = TMath::Sqrt( lXtf_Barrel1 * lXtf_Barrel1  +  lYtf_Barrel1 * lYtf_Barrel1 );
            Double_t lRt_tfBarrel2  = TMath::Sqrt( lXtf_Barrel2 * lXtf_Barrel2  +  lYtf_Barrel2 * lYtf_Barrel2 );
            
            
            
            Double_t lTransvRotAngle_Barrel1 = 2 * TMath::ASin( lRt_tfBarrel1     /2. *1/(-rho1) );    // NOTE : beware = signed angle here, in accordance with right-angled frame (0;x,y,z)
                        // NOTE : should be valid for both TOF barrel and endcap... 
                        // NOTE : It is just that, in the barrel case, lRt_tf should be = lRbarrelTOF ...
                        // NOTE NOTE : but ... not valid if several half-turns ! This is only valid within [0, pi] = via visible final and concrete Rtf
                        //                  i.e. valid for a track that stays in the 1st quadrant, phi in [0;pi/2] = rotation around centre between [0;pi]
                        //                  That is typically true for that track meant to TOF barrel layer.
                        //                  But may not always be true for endcap (>180° or loopers, etc)
                        //              Hence the need to compute a specific angle for endcap...
                     lAbsNbHalfTurns_Barrel1 = TMath::Floor( TMath::Abs( lTransvRotAngle_Barrel1 *TMath::RadToDeg()/180.));   
                        // NOTE : number of completed half turn(s)
                        //          Should always be 0 in the physics case for barrel, no looper allowed
    

            
            
            
            
            
            // Math for potential endcap case
            Double_t ltf_Endcap1    = TMath::Abs( lZendcapTOF / ( TMath::C() * betaZ1  ) );
            Double_t ltf_Endcap2    = TMath::Abs( lZendcapTOF / ( TMath::C() * betaZ2  ) );  
    
            Double_t lXtf_Endcap1   = -rho1 *   TMath::Sin( AngFreq1 * ltf_Endcap1 )     ;
            Double_t lYtf_Endcap1   = -rho1 * ( TMath::Cos( AngFreq1 * ltf_Endcap1 ) - 1);
            Double_t lZtf_Endcap1   = betaZ1 * TMath::C()            * ltf_Endcap1 ; 
            
                     lTmpX = lXtf_Endcap1;
                     lTmpY = lYtf_Endcap1;
                     
                     lXtf_Endcap1 = lTmpX * TMath::Cos( lPhit0 ) - lTmpY * TMath::Sin( lPhit0 );  // trigo formula cos(a+b) 
                     lYtf_Endcap1 = lTmpY * TMath::Cos( lPhit0 ) + lTmpX * TMath::Sin( lPhit0 );  // trigo formula sin(a+b)
            
    
            
            Double_t lXtf_Endcap2   = -rho2 *   TMath::Sin( AngFreq2 * ltf_Endcap2 )     ;               
            Double_t lYtf_Endcap2   = -rho2 * ( TMath::Cos( AngFreq2 * ltf_Endcap2 ) - 1);                
            Double_t lZtf_Endcap2   = betaZ2 * TMath::C()            * ltf_Endcap2 ;
            
                     lTmpX = lXtf_Endcap2;
                     lTmpY = lYtf_Endcap2;
                     
                     lXtf_Endcap2 = lTmpX * TMath::Cos( lPhit0 ) - lTmpY * TMath::Sin( lPhit0 );  // trigo formula cos(a+b) 
                     lYtf_Endcap2 = lTmpY * TMath::Cos( lPhit0 ) + lTmpX * TMath::Sin( lPhit0 );  // trigo formula sin(a+b)
            
            
            Double_t lRt_tfEndcap1  = TMath::Sqrt( lXtf_Endcap1 * lXtf_Endcap1  +  lYtf_Endcap1 * lYtf_Endcap1 );        
            Double_t lRt_tfEndcap2  = TMath::Sqrt( lXtf_Endcap2 * lXtf_Endcap2  +  lYtf_Endcap2 * lYtf_Endcap2 );        
            

            Double_t lTransvRotAngle_Endcap1 = lZtf_Endcap1/(-rho1 * TMath::SinH(eta)); // e22 in South face p.9
                     lAbsNbHalfTurns_Endcap1 = TMath::Floor( TMath::Abs( lTransvRotAngle_Endcap1 *TMath::RadToDeg()/180.)); 
    
    


            
            

            
            if(     kLooper 
                &&  (lRt_tfEndcap1  < lRminEndcap ||  lRt_tfEndcap1  > lRmaxEndcap) 
                )                                                                       kLooper = kFALSE;  // acceptance check : track ends up below lRminEndcap, so undo the flag if necessary.
            
         // if(kLooper &&  lAbsNbHalfTurns_Endcap1 != 1) kLooper = kFALSE;  // For investigations on rotation for loopers : isolate a given number of half turn + apply l.1427 (~15 lines below)
            
            
            
            
            if(         TMath::Abs(lZtf_Barrel1)  < lZmaxBarrel 
                    &&  TMath::Abs(lRt_tfBarrel1 - lRbarrelTOF) < 2e-3 
                  //&&  lAbsNbHalfTurns_Barrel1 == 0
                ){ 
                                                                    kWithinRZbarrelAccptce  = kTRUE; }
                                                                    
                                                                    
            if(    lRt_tfEndcap1 < lRmaxEndcap                  // WARNING : not a else if, but an independent if to allow for double acceptance !
                    &&  lRt_tfEndcap1  > lRminEndcap 
                    &&  TMath::Abs(lZtf_Endcap1 - lZendcapTOF)  < 2e-3
                  //&&  lAbsNbHalfTurns_Endcap1 == 1     
                ){     
                                                                    kWithinRZendcapAccptce  = kTRUE; }
    
                
             
                
            if (  kWithinRZbarrelAccptce   && kWithinRZendcapAccptce && kActivateEndcap){
                ++iDbleAcceptance;
                if(lDebug > 1) Printf("Double acceptance (barrel + endcap) Point [pT = %6.4g GeV/c ; y1 = %6.4f / eta = %6.4f ] :", pt, rapidity1, eta );
                
                
                // 2nd thought needed : Reevaluate the assignment in terms of acceptance : only in barrel, only in endcap or in both, finally
                
                if( kWithinRZbarrelAccptce && kAbovePtMinCalo )                    kWithinRZendcapAccptce = kFALSE;
                        // NOTE : the track, first and foremost a barrel track, will stop in calorimeter, 
                        //          so for sure will not manage to further loop back towards the endcap plane
                        //          In practice, it can only be a barrel TOF measurement
                
                if( kWithinRZbarrelAccptce && kAboveLambdaMinAngle  == kFALSE )    kWithinRZendcapAccptce = kFALSE;
                        // NOTE : track, first and foremost a barrel track, is here a priori allowed to loop 
                        //              - above the TOF barrel radius
                        //              - but below the calorimeter one
                        //          and thus loop towards the endcap !
                        //          But then , the requirement about the minimal inclination angle is not met : lambda too small
                        //          In practice, it will remain a clean measurement only in barrel TOF.
                
                if(kRelaxAngConstraints)                                           kWithinRZendcapAccptce = kTRUE;
                        // NOTE : to visualise the full phase ~without constraints (calo, lambda,  delta angles)
                
                
                    
                if(lDebug > 2) Printf("Hyp. Barrel : tf1 = %8.6f s / x1 = %+8.6f m, y1 = %+8.6f m, z1 = %+8.6f m, r1 = %8.6f m / (NB : ZmaxBarrel = %5.3f m), Nb(1/2 completed turns) = %02d", 
                          ltf_Barrel1,
                         lXtf_Barrel1, 
                         lYtf_Barrel1, 
                         lZtf_Barrel1,
                        lRt_tfBarrel1,
                         lZmaxBarrel ,
                         lAbsNbHalfTurns_Barrel1
                    );              
                if(lDebug > 3) Printf("Hyp. Barrel : tf2 = %8.6f s / x2 = %+8.6f m, y2 = %+8.6f m, z2 = %+8.6f m, r2 = %8.6f m / (NB : ZmaxBarrel = %5.3f m)", 
                          ltf_Barrel2,
                         lXtf_Barrel2, 
                         lYtf_Barrel2, 
                         lZtf_Barrel2,
                        lRt_tfBarrel2,
                         lZmaxBarrel
                    );
                if(lDebug > 1) Printf("Hyp. Endcap : tf1 = %8.6f s / x1 = %+8.6f m, y1 = %+8.6f m, z1 = %+8.6f m, r1 = %8.6f m / (NB : lRmaxEndcap = %5.3f m), Nb(1/2 completed turns) = %02d", 
                          ltf_Endcap1,
                         lXtf_Endcap1, 
                         lYtf_Endcap1, 
                         lZtf_Endcap1,
                        lRt_tfEndcap1,
                         lRmaxEndcap ,
                         lAbsNbHalfTurns_Endcap1
                    );
                if(lDebug > 3) Printf("Hyp. Endcap : tf2 = %8.6f s / x2 = %+8.6f m, y2 = %+8.6f m, z2 = %+8.6f m, r2 = %8.6f m / (NB : lRmaxEndcap = %5.3f m)", 
                          ltf_Endcap2,
                         lXtf_Endcap2, 
                         lYtf_Endcap2, 
                         lZtf_Endcap2,
                        lRt_tfEndcap2,
                         lRmaxEndcap 
                    );
                if(lDebug > 2) Printf(" ");
            }
            
            
            
           
            if(kLooper  && kActivateEndcap){
                ltf1 = ltf_Endcap1;
                ltf2 = ltf_Endcap2;
            }            
            else if(kWithinRZendcapAccptce  && kActivateEndcap){  // simplistic case : eta > lEtaMinEndcap   or more advanced version : kWithinRZendcapAccptce
                IsPartForeseenInEndcap = kTRUE; 
                ltf1 = ltf_Endcap1;
                ltf2 = ltf_Endcap2;
            }
            else if(kWithinRZbarrelAccptce){  // simplistic case : eta < lEtaMaxBarrel  or more advanced version :  kWithinRZbarrelAccptce
                IsPartForeseenInBarrel = kTRUE;
                ltf1 = ltf_Barrel1;
                ltf2 = ltf_Barrel2;
            } 
            

            
           cT_diff = TMath::Abs( ltf1 - ltf2 ) * 1e9;  // time, from seconds to nanoseconds
           // if(IsPartForeseenInBarrel){
           // if( TMath::Abs(pt - 1.00) < 0.01 && TMath::Abs(eta - 0.65) < 1e-2   )
           //    Printf("Point [pT = %6.4g GeV/c ; y1 = %6.4f / y2 = %6.4f / eta = %6.4f ] : ctdiff = [tf1 = %8.6g] - [tf2 = %8.6g] = %8.6g ns", pt, rapidity1, rapidity2, eta, ltf1, ltf2, cT_diff);
           // }
         
            
        }    
        else{// Wei Li's core code
            if(eta<lEtaMaxBarrel){
                L1 = pt/clight*cosh(eta)/lBz/q1*acos(1-lBz*lBz*q1*q1*lRbarrelTOF*lRbarrelTOF/pt/pt*clight*clight/2.);
                L2 = pt/clight*cosh(eta)/lBz/q2*acos(1-lBz*lBz*q2*q2*lRbarrelTOF*lRbarrelTOF/pt/pt*clight*clight/2.);
            }
            else if(eta>lEtaMinEndcap && kActivateEndcap){
                L1 = lZendcapTOF / tanh(eta);
                L2 = lZendcapTOF / tanh(eta);
            }
            else{
                L1 = 0;
                L2 = 0;
            }

            double beta1_inv  = sqrt(m1*m1/pt/pt/cosh(eta)/cosh(eta)+1);
            double beta2_inv  = sqrt(m2*m2/pt/pt/cosh(eta)/cosh(eta)+1);
                   cT_diff    = fabs( beta1_inv*L1 - beta2_inv*L2 )/clight;
            
            if(eta<lEtaMaxBarrel)  IsPartForeseenInBarrel = kTRUE;
            if(eta>lEtaMinEndcap)  IsPartForeseenInEndcap = kTRUE;
            
            // NB : No looper anticipated in Wei's code
            
        }// end Wei's code    
            

            
        if(IsPartForeseenInBarrel && kAboveDeltaMinAngle){
            if( cT_diff > n5sigma*cDeltaT )     h5sigm_barrel[ iCan ]  ->SetBinContent( h5sigm_barrel[ iCan ]  ->FindBin(rapidity1,pt), 1);
            if( cT_diff > n4sigma*cDeltaT )     h4sigm_barrel[ iCan ]  ->SetBinContent( h4sigm_barrel[ iCan ]  ->FindBin(rapidity1,pt), 1);
            if( cT_diff > n3sigma*cDeltaT )     h3sigm_barrel[ iCan ]  ->SetBinContent( h3sigm_barrel[ iCan ]  ->FindBin(rapidity1,pt), 1);
            if( cT_diff > n2sigma*cDeltaT )     h2sigm_barrel[ iCan ]  ->SetBinContent( h2sigm_barrel[ iCan ]  ->FindBin(rapidity1,pt), 1);
            if( cT_diff > n1sigma*cDeltaT )     h1sigm_barrel[ iCan ]  ->SetBinContent( h1sigm_barrel[ iCan ]  ->FindBin(rapidity1,pt), 1);      
        }
            
        if(IsPartForeseenInEndcap && kAboveLambdaMinAngle){
            if( cT_diff > n5sigma*cDeltaT )     h5sigm_ecap  [ iCan ]  ->SetBinContent( h5sigm_ecap  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n4sigma*cDeltaT )     h4sigm_ecap  [ iCan ]  ->SetBinContent( h4sigm_ecap  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n3sigma*cDeltaT )     h3sigm_ecap  [ iCan ]  ->SetBinContent( h3sigm_ecap  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n2sigma*cDeltaT )     h2sigm_ecap  [ iCan ]  ->SetBinContent( h2sigm_ecap  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n1sigma*cDeltaT )     h1sigm_ecap  [ iCan ]  ->SetBinContent( h1sigm_ecap  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);    
        }
        
        if(kLooper){
            if( cT_diff > n5sigma*cDeltaT )     h5sigm_loop  [ iCan ]  ->SetBinContent( h5sigm_loop  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n4sigma*cDeltaT )     h4sigm_loop  [ iCan ]  ->SetBinContent( h4sigm_loop  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n3sigma*cDeltaT )     h3sigm_loop  [ iCan ]  ->SetBinContent( h3sigm_loop  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n2sigma*cDeltaT )     h2sigm_loop  [ iCan ]  ->SetBinContent( h2sigm_loop  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);
            if( cT_diff > n1sigma*cDeltaT )     h1sigm_loop  [ iCan ]  ->SetBinContent( h1sigm_loop  [ iCan ]  ->FindBin(rapidity1,pt), 10+lAbsNbHalfTurns_Endcap1);            
        }
        
        
        

        if(ptot<p_dedx) hdEdx_CMS[ iCan ]->SetBinContent(hdEdx_CMS[ iCan ]->FindBin(rapidity1,pt),1);
            
        
    }// end loop iEta
  }// end loop pt

  
    Printf(" Total double acceptance (barrel + endcap) counts = %03d", iDbleAcceptance);
  
  
  
    // kBlue-4        kGreen+2     kPink          kOrange+7
    // kAzure+7       kGreen+1     kPink-4        kOrange-3
    // kAzure+6       kSpring      kRed-9         kOrange-2
    myHistoSetUp( h5sigm_barrel[ iCan ], 1.0, 20, 0, 1, 1, 0, 1001, kPink     );
    myHistoSetUp( h4sigm_barrel[ iCan ], 1.0, 20, 0, 1, 1, 0, 1001, kRed-10   );   
    myHistoSetUp( h3sigm_barrel[ iCan ], 1.0, 20, 0, 1, 1, 0, 1001, kPink     );
    myHistoSetUp( h2sigm_barrel[ iCan ], 1.0, 20, 0, 1, 1, 0, 1001, kPink-4   );
    myHistoSetUp( h1sigm_barrel[ iCan ], 1.0, 20, 0, 1, 1, 0, 1001, kRed-9    );

    myHistoSetUp( h5sigm_ecap[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange+7 );    
    myHistoSetUp( h4sigm_ecap[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-4 );    
    myHistoSetUp( h3sigm_ecap[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange+7 );
    myHistoSetUp( h2sigm_ecap[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-3 );
    myHistoSetUp( h1sigm_ecap[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-2 );

    myHistoSetUp( h5sigm_loop[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange+7 ); 
    myHistoSetUp( h4sigm_loop[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-4 );  
    myHistoSetUp( h3sigm_loop[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange+7 ); 
    myHistoSetUp( h2sigm_loop[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-3 ); 
    myHistoSetUp( h1sigm_loop[ iCan ],   1.0, 20, 0, 1, 1, 0, 1001, kOrange-2 ); 
    
    
    //  hdEdx_CMS ->SetLineColor(kGray);  // kGreen-10
    //     
    //  latex2->SetTextSize(0.033);
    //  latex2->SetTextColor(kGreen+2);
    //  latex2->SetTextAngle(-30);
    //  latex2->DrawLatex(0.2,0.31,"dE/dx K/p limit");
    //  latex2->SetTextAngle(0);
    //  latex2->SetTextSize(0.045);
    //  latex2->DrawLatex(0.145,0.2,"dE/dx");   
    
    
    

    h1sigm_barrel[ iCan ]->Draw("BOXsame");
    h2sigm_barrel[ iCan ]->Draw("BOXsame");
    h3sigm_barrel[ iCan ]->Draw("BOXsame");
    h4sigm_barrel[ iCan ]->Draw("BOXsame");  
    h5sigm_barrel[ iCan ]->Draw("BOXsame");
    
    h1sigm_ecap  [ iCan ]->Draw("BOXsame");
    h2sigm_ecap  [ iCan ]->Draw("BOXsame");
    h3sigm_ecap  [ iCan ]->Draw("BOXsame");
    h4sigm_ecap  [ iCan ]->Draw("BOXsame");  
    h5sigm_ecap  [ iCan ]->Draw("BOXsame");
    
    h1sigm_loop  [ iCan ]->Draw("BOXsame");
    h2sigm_loop  [ iCan ]->Draw("BOXsame");
    h3sigm_loop  [ iCan ]->Draw("BOXsame");
    h4sigm_loop  [ iCan ]->Draw("BOXsame");  
    h5sigm_loop  [ iCan ]->Draw("BOXsame");
  
    
    
    //hdEdx_CMS->Draw("BOXsame");
    gPad->RedrawAxis("g");  // to have Grid on top of the previous curves/TH1
    
    
    
    myGraphSetUp( grEta0pt5[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta1pt0[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta1pt5[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta2pt0[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta2pt5[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta3pt0[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    myGraphSetUp( grEta3pt5[ iCan ], 1.0, kDot, kBlack, 2, 7, kGray+0, 0, kGray);
    

    grEta0pt5[ iCan ]->Draw("Csame");
    grEta1pt0[ iCan ]->Draw("Csame");
    grEta1pt5[ iCan ]->Draw("Csame");
    grEta2pt0[ iCan ]->Draw("Csame");
    grEta2pt5[ iCan ]->Draw("Csame");
    grEta3pt0[ iCan ]->Draw("Csame");
    grEta3pt5[ iCan ]->Draw("Csame");
    

    myGraphSetUp( grEtaMaxBarrel[ iCan ], 1.0, kDot,        kBlack,      2, 7, kPink-4,     0, kPink    );
    myGraphSetUp( grEtaMinEndcap[ iCan ], 1.0, kDot,        kBlack,      2, 7, kOrange-3,   0, kOrange+7);
    myGraphSetUp( grEtaMaxEndcap[ iCan ], 1.0, kDot,        kBlack,      2, 7, kOrange-3,   0, kOrange+7);
    
    myGraphSetUp( grEtaMinLambda[ iCan ], 0.8, kDot,        kAzure+1, +302, 7, kAzure+1, 3235, kAzure+1 );
    myGraphSetUp( grEtaMinDelta [ iCan ], 0.8, kDot,        kAzure+1, -302, 7, kAzure+1, 3253, kAzure+1 );
 // myGraphSetUp( grEtaMinDelta [ iCan ], 0.8, kOpenCircle, kAzure+1,    1, 7, kAzure+1,    0, kAzure+1 );
    myGraphSetUp( grPtMinBarrel [ iCan ], 0.8, kFullSquare, kTeal+3,     3, 7, kTeal+3,     0, kPink    );
    myGraphSetUp( grPtMinEndcap [ iCan ], 0.8, kFullSquare, kTeal+3,     3, 7, kTeal+3,     0, kOrange+7);
    myGraphSetUp( grPtMinCalo   [ iCan ], 0.6, kOpenSquare, kTeal+3,     3, 7, kTeal+3,     0, kOrange+7);

    
    

    
                        grEtaMaxBarrel[ iCan ]->Draw("Csame");
                        grEtaMinEndcap[ iCan ]->Draw("Csame");
                        grEtaMaxEndcap[ iCan ]->Draw("Csame");

                        grPtMinBarrel[ iCan ]->Draw("PCsame");
    if(kActivateEndcap) grPtMinEndcap[ iCan ]->Draw("PCsame");
                      //grPtMinCalo  [ iCan ]->Draw("PCsame");  // --> see below : point will be potentially updated, then Draw only at that stage

    
                        grEtaMinDelta [ iCan ]->Draw("Csame");
    if(kActivateEndcap) grEtaMinLambda[ iCan ]->Draw("Csame");
                      //fnPtMinBarrelDelta->Draw("LCsame"); // Beware the function is thought as f(eta) while here the abscissa are y ≠ eta!
    
    
    
    // Spot transition between consecutive nb of half turns in trajectory

    Printf(" LocateLimitBetweenHalfTurns ... 0-1, 2-3, 4-5");
    grHalfTurns_1to0[ iCan ] = LocateLimitBetweenHalfTurns( h5sigm_ecap[ iCan ], h5sigm_loop[ iCan ],  1);
    grHalfTurns_3to2[ iCan ] = LocateLimitBetweenHalfTurns( h5sigm_ecap[ iCan ], h5sigm_loop[ iCan ],  3);
    grHalfTurns_5to4[ iCan ] = LocateLimitBetweenHalfTurns( h5sigm_ecap[ iCan ], h5sigm_loop[ iCan ],  5);

    
    if( grHalfTurns_1to0[ iCan ] )  myGraphSetUp(   grHalfTurns_1to0[ iCan ], 0.15, kFullSquare, kWhite, 1, kSolid, kWhite, 0, kOrange+7);
    if( grHalfTurns_3to2[ iCan ] )  myGraphSetUp(   grHalfTurns_3to2[ iCan ], 0.15, kFullSquare, kWhite, 1, kSolid, kWhite, 0, kOrange+7);
    if( grHalfTurns_5to4[ iCan ] )  myGraphSetUp(   grHalfTurns_5to4[ iCan ], 0.15, kFullSquare, kWhite, 1, kSolid, kWhite, 0, kOrange+7);
    
    if(kActivateEndcap && grHalfTurns_1to0[ iCan ] ) grHalfTurns_1to0[ iCan ]->Draw("Psame");
    if(kActivateEndcap && grHalfTurns_3to2[ iCan ] ) grHalfTurns_3to2[ iCan ]->Draw("Psame");
    if(kActivateEndcap && grHalfTurns_5to4[ iCan ] ) grHalfTurns_5to4[ iCan ]->Draw("Psame");
    
  
    
    
    
    
    
    
    
    
    

    lPaveTxt_Sep[ iCan ] = new TPaveText(0.65, 0.89, 0.97, 0.97, "NDC ARC");
    lPaveTxt_Sep[ iCan ]->SetCornerRadius(0.22);
    lPaveTxt_Sep[ iCan ]->SetFillColor(kWhite);
    lPaveTxt_Sep[ iCan ]->SetShadowColor(kWhite);
    lPaveTxt_Sep[ iCan ]->SetBorderSize(0);  
        lTxt_Sep[ iCan ] = lPaveTxt_Sep[ iCan ]->AddText(Form("%s/%s #scale[0.7]{separation}", Str_HadronLatexName[lPart2].Data(), Str_HadronLatexName[lPart1].Data()) );
        lTxt_Sep[ iCan ]->SetTextSize(0.05);
        lTxt_Sep[ iCan ]->SetTextColor(kRed+1);    
    lPaveTxt_Sep[ iCan ]->Draw();
  

    lLatex_Hyp[ iCan ] = new TLatex();
    // lLatex_Hyp[ iCan ]->SetNDC();
    lLatex_Hyp[ iCan ]->SetTextSize(0.04);
    lLatex_Hyp[ iCan ]->SetTextColor(colors[red]);
    lLatex_Hyp[ iCan ]->DrawLatex(-0.2,55, Form("%s #scale[0.7]{config.} - #scale[0.7]{[ %s ]}", Str_Exp.Data(), Str_LHCrun.Data()) );
    
    
  
    lPaveTxt_Info[ iCan ] = new TPaveText(2.7, 0.012, 4.7, 0.09, "ARC");
    lPaveTxt_Info[ iCan ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ iCan ]->SetFillColor(kWhite);
    lPaveTxt_Info[ iCan ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ iCan ]->SetLineColor(kYellow-10);
    lPaveTxt_Info[ iCan ]->SetBorderSize(0);
    
    TText *lTxt_Info[4];
    lTxt_Info[0] = lPaveTxt_Info[ iCan ]->AddText( Form("Solenoidal B field : %3.1f T", lBz) );
    lTxt_Info[0]->SetTextSize(0.03);
    lTxt_Info[0]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[0]->SetTextColor(kGray+3);
    lTxt_Info[1] = lPaveTxt_Info[ iCan ]->AddText( Form("#sigma^{TOF}(#scale[0.7]{#color[%d]{EndCap} or #color[%d]{Barrel}}) : %2.0f ps", kOrange+7, kPink, cDeltaT*1000) ); // 807 = kOrange+7; 900 = kPink       
    lTxt_Info[1]->SetTextSize(0.03);
    lTxt_Info[1]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[1]->SetTextColor(kGray+3);
    
    lTxt_Info[2] = lPaveTxt_Info[ iCan ]->AddText( Form("<#it{r}>^{TOF}(Barrel) #kern[+4.3]{:} %4.2f m", lRbarrelTOF));       
    lTxt_Info[2]->SetTextSize(0.03);
    lTxt_Info[2]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[2]->SetTextColor(kPink);  
    
    if(kActivateEndcap)  lTxt_Info[3] = lPaveTxt_Info[ iCan ]->AddText( Form("<#it{z}>^{TOF}(EndCap) : %4.2f m", lZendcapTOF) );       
    else                 lTxt_Info[3] = lPaveTxt_Info[ iCan ]->AddText( "<#it{z}>^{TOF}(EndCap) : --");
    lTxt_Info[3]->SetTextSize(0.03);
    lTxt_Info[3]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[3]->SetTextColor(kOrange+7);

    lPaveTxt_Info[ iCan ]->Draw();
  


    lLatexEta[ iCan ] = new TLatex();
    // lLatexEta[ iCan ]->SetNDC();
    lLatexEta[ iCan ]->SetTextSize(0.02);
    lLatexEta[ iCan ]->SetTextAngle(65);
    lLatexEta[ iCan ]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lLatexEta[ iCan ]->SetTextColor(kGray+2);
    
    lLatexEta[ iCan ]->DrawLatex(0.56,  21, "#eta = 0.5");
    lLatexEta[ iCan ]->DrawLatex(1.05,  20," #eta = 1.0");
    lLatexEta[ iCan ]->DrawLatex(1.40,  25, "#eta = 1.5");
    lLatexEta[ iCan ]->DrawLatex(2.10,  15, "#eta = 2.0");
    lLatexEta[ iCan ]->DrawLatex(2.60,  15, "#eta = 2.5");
    lLatexEta[ iCan ]->DrawLatex(3.10,  12, "#eta = 3.0");
    lLatexEta[ iCan ]->DrawLatex(3.30,  12, "#eta = 3.5");
    
    
    lLatexEta[ iCan ]->SetTextSize(0.028);
    lLatexEta[ iCan ]->SetTextColor(kPink);    
    lLatexEta[ iCan ]->DrawLatex(lEtaMaxBarrel-0.10,  15.0, Form("#eta = %4.2f", lEtaMaxBarrel));
    
    lLatexEta[ iCan ]->SetTextColor(kOrange+7);
    if(kActivateEndcap) lLatexEta[ iCan ]->DrawLatex(lEtaMinEndcap+0.05,  15.0, Form("#eta = %4.2f", lEtaMinEndcap));
    if(kActivateEndcap) lLatexEta[ iCan ]->DrawLatex(lEtaMaxEndcap-0.25,  12.0, Form("#eta = %4.2f", lEtaMaxEndcap));
  

    
    lLatexEta[ iCan ]->SetTextAngle(0);
    lLatexEta[ iCan ]->SetTextSize(0.025);
    lLatexEta[ iCan ]->SetTextColor(kAzure+2);
    if(kActivateEndcap) lLatexEta[ iCan ]->DrawLatex(0.15, 0.012, Form("#eta^{min}_{#lambda} = %4.2f (#lambda>%02.f^{o})", EtaMinLambda, lLambdaMinAngle* TMath::RadToDeg() ));
    
    
    
    Double_t lPtCoord = 0., lYCoord = 0.;
    grEtaMinDelta[ iCan ] ->GetPoint(lNbEtaDelta-1, lYCoord, lPtCoord); // get the last point of the curve = the upper tip
    Printf(" ");
    Printf(" EtaDelta : Upper tip of the curve [y = %6.4f, pT = %6.4f GeV/c]", lYCoord, lPtCoord);
    
    lLatexEta[ iCan ]->SetTextAngle(0);
    lLatexEta[ iCan ]->DrawLatex(lYCoord * 1.10 , lPtCoord,      Form("p_{T,min}^{ barrel} x (k_{#delta})^{-1} (#delta>%02.f^{o})", lDeltaMinAngle* TMath::RadToDeg() ));
    lLatexEta[ iCan ]->DrawLatex(lYCoord * 1.27 , lPtCoord *0.6, Form("|#eta| < #eta^{max}_{#delta} = %4.2f", EtaMaxDelta));
    
   
    lLatexEta[ iCan ]->SetTextAngle(0);
    lLatexEta[ iCan ]->SetTextSize(0.025);
    lLatexEta[ iCan ]->SetTextColor(kTeal+3);
                        lLatexEta[ iCan ]->DrawLatex(lYMaxForPtMinBarrel-0.2, lPtMinBarrel*0.7,   Form("p_{T,min}^{barrel} = %5.3f GeV/c", lPtMinBarrel ));
    if(kActivateEndcap) lLatexEta[ iCan ]->DrawLatex(lYMaxForPtMinEndcap-0.1, lPtMinEndcap*0.7,   Form("p_{T,min}^{endcap} = %5.3f GeV/c", lPtMinEndcap ));
    
    
    Int_t    iLastPnt_HalfTurns_1to0 = -1;
    Double_t lYMaxECalEdge = -998;
    Double_t lPtDummy      = -997;
    if( grHalfTurns_1to0[ iCan ] ) {
        iLastPnt_HalfTurns_1to0 = grHalfTurns_1to0[ iCan ]->GetN() - 1;    
        grHalfTurns_1to0[ iCan ]->GetPoint( iLastPnt_HalfTurns_1to0, lYMaxECalEdge, lPtDummy);
    }
    
    Printf(" lYMaxECalEdge = %5.3f / lPtMinCalo = %5.3f GeV/c", lYMaxECalEdge, lPtDummy );
    
    //if(TMath::Abs(lYMaxECalEdge + 998) > 1e-13) grPtMinCalo  [ iCan ]->SetPoint(1,  lYMaxECalEdge, lPtMinCalo); 
        // NOTE : such a point update based on limit half-turn 1-0 does not make always sense, especially for TOF barrel in the middle of barrel tracker (ATLAS, ALICE-3)
        //  No unique solution...
        //  The point update is cool and suggestive for CMS while previous value (lEtaMaxCalo) is potentially confusing
        //  But at the end of the day, prefer to keep the point location based on lYMaxForPtMinCalo = lEtaMaxCalo
    
    grPtMinCalo  [ iCan ]->Draw("PCsame");    
    
    
    lLatexEta[ iCan ]->DrawLatex(lYMaxForPtMinCalo  +0.1,   lPtMinCalo*1.1,   Form("#scale[0.7]{ECal}"));  // Form("p_{T,min}^{ECal} = %5.3f GeV/c",   lPtMinCalo ));
    
    
    
    
    
    
  
    if(kActivateEndcap)   myLegend[ iCan ] = new TLegend(0.865,0.42,0.99,0.86);
    else                  myLegend[ iCan ] = new TLegend(0.865,0.62,0.99,0.86);
    myLegendSetUp(myLegend[ iCan ],0.03);
    myLegend[ iCan ]->AddEntry(myBlankHisto [ iCan ],  Form("#color[%d]{#font[%d]{#scale[1.2]{%s}}}", kPink-4, 62, StrBarrelTOFname.Data()), "");
    myLegend[ iCan ]->AddEntry(h1sigm_barrel[ iCan ],  "1-#sigma",     "F");
    myLegend[ iCan ]->AddEntry(h2sigm_barrel[ iCan ],  "2-#sigma",     "F");
    myLegend[ iCan ]->AddEntry(h3sigm_barrel[ iCan ],  "3-#sigma",     "F");
    myLegend[ iCan ]->AddEntry(h4sigm_barrel[ iCan ],  "4-#sigma",     "PF");
    myLegend[ iCan ]->AddEntry(h5sigm_barrel[ iCan ],  "5-#sigma",     "F");
    myLegend[ iCan ]->AddEntry(myBlankHisto [ iCan ],  "#color[0]{}",  "" ); // dummy entry, white spacer
    
    if(kActivateEndcap){
            myLegend[ iCan ]->AddEntry(myBlankHisto[ iCan ],    Form("#color[%d]{#font[%d]{#scale[1.2]{%s}}}", kOrange-3, 62, StrEndcapTOFname.Data()), "");
            myLegend[ iCan ]->AddEntry(h1sigm_ecap [ iCan ],    "1-#sigma", "F");
            myLegend[ iCan ]->AddEntry(h2sigm_ecap [ iCan ],    "2-#sigma", "F");
            myLegend[ iCan ]->AddEntry(h3sigm_ecap [ iCan ],    "3-#sigma", "F");
            myLegend[ iCan ]->AddEntry(h4sigm_ecap [ iCan ],    "4-#sigma", "PF");
            myLegend[ iCan ]->AddEntry(h5sigm_ecap [ iCan ],    "5-#sigma", "F");
    }
    myLegend[ iCan ]->Draw();

  
  
    

    

    // Spot specific values on the figure to get quantitative and explicit ideas
    // in some key corners of the TH2
  
    std::vector<Int_t> vCornerYBin;    vCornerYBin .resize(kNall);
    std::vector<Int_t> vCornerPtBin;   vCornerPtBin.resize(kNall);
      
    // kMidY_SW = 0,  kMidY_NW = 1,    kMidY_SEE,    kMidY_NEE,    kFwdY_NW,    kFwdY_NE,      kNall

    vCornerYBin  [0] = 999;
    
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kMidY_SW,   1,  -2, h4sigm_barrel[ iCan ] );
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kMidY_NW,   1,  -2, h4sigm_barrel[ iCan ] );
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kMidY_SEE, -2,  -2, h4sigm_barrel[ iCan ] );
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kMidY_NEE, -2,  -2, h4sigm_barrel[ iCan ] );

    if(kActivateEndcap){
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kFwdY_NW,  -2,  -2, h4sigm_ecap  [ iCan ] );
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kFwdY_SEE, -2,  -2, h4sigm_ecap  [ iCan ] );       
        LocatePhaseSpaceCorner(vCornerYBin, vCornerPtBin, kFwdY_NEE, -2,  -2, h4sigm_ecap  [ iCan ] );       
    }
    
   
    
    Printf(" ");  
    
    
    std::vector<TMarker*> vCornerMarker; vCornerMarker.resize(kNall);
    std::vector<TLatex*>  vLatexCoord;   vLatexCoord.resize(kNall);
    
    
    for(UInt_t iCorner = 0; iCorner < kNall; ++iCorner){
        
        Float_t lYrap = 0;
        Float_t lPt   = 0;
        
        switch(iCorner){
            case kMidY_SW   :
            case kMidY_NW   :
            case kMidY_SEE  :
            case kMidY_NEE  :   
                lYrap = h4sigm_barrel[ iCan ]->GetXaxis()->GetBinCenter( vCornerYBin  [iCorner] );
                lPt   = h4sigm_barrel[ iCan ]->GetYaxis()->GetBinCenter( vCornerPtBin [iCorner] );
                break;
                
                
            case kFwdY_NW   :
            case kFwdY_SEE  :
            case kFwdY_NEE  : 
                lYrap = h4sigm_ecap[ iCan ]->GetXaxis()->GetBinCenter( vCornerYBin  [iCorner] );
                lPt   = h4sigm_ecap[ iCan ]->GetYaxis()->GetBinCenter( vCornerPtBin [iCorner] );
                break;
        }
        
        
        if(iCorner == kFwdY_NW ) Printf(" ");
        
        Printf("Corner[%02d] = %9s - 4sigmaSep : bin[%03d,%03d] = (y = %+6.4f, pT = %+6.3f GeV/c)",
                                iCorner,
                                Str_Corner[iCorner].Data(),
                                vCornerYBin  [iCorner],
                                vCornerPtBin [iCorner],
                                lYrap, 
                                lPt); 
       
        
        if(iCorner == kMidY_SEE || iCorner == kFwdY_NEE ) continue;
        
        vCornerMarker[iCorner] = new TMarker(lYrap, lPt, kOpenCircle);
        vCornerMarker[iCorner]->SetMarkerSize(0.8);
        vCornerMarker[iCorner]->SetMarkerColor(kGray+3);
        if(iCorner == kFwdY_NW ) vCornerMarker[iCorner]->SetMarkerColor(kRed+3);
        
        
        vLatexCoord[iCorner] = new TLatex();
            vLatexCoord[iCorner]->SetTextSize(0.03);
            vLatexCoord[iCorner]->SetTextAngle(-18);
            vLatexCoord[iCorner]->SetTextColor(kGray+3);
            vLatexCoord[iCorner]->SetTextAlign(kHAlignCenter+kVAlignTop);
            Bool_t kUp = kFALSE;
            
            if( iCorner == kMidY_SW ) kUp = 1;
            if( iCorner == kMidY_NW ) kUp = 1;
            if( iCorner == kMidY_NEE) kUp = 1;
            
            if( iCorner == kFwdY_NW ) kUp = 0;  
            if( iCorner == kFwdY_SEE) kUp = 0;
            
            
            
        if( iCorner == kMidY_SW || iCorner == kMidY_NW || iCorner == kMidY_NEE){    
            vCornerMarker[iCorner]->Draw();
            vLatexCoord[iCorner]->DrawLatex(lYrap,  lPt*( kUp? 1.5:0.8),    Form("%4.2f #scale[0.7]{GeV/#it{c}}", lPt)); 
        }
        
        else if(kActivateEndcap && ( iCorner == kFwdY_NW || iCorner == kFwdY_SEE ) ){
            vCornerMarker[iCorner]->Draw();
            vLatexCoord[iCorner]->DrawLatex(lYrap,  lPt*( kUp? 1.5:0.8),    Form("%4.2f #scale[0.7]{GeV/#it{c}}", lPt)); 
        }
        
        
    }// end iCorner
    
    Printf(" ");
    
    
    Str_FigFilename[ iCan ] = Form("z%s_timing_%02d_%s_RbarrelTOF%4.2fm_ZendcapTOF%4.2fm_%2.0fps_WithBfield%+.1fT", 
                                   Str_Exp.Data(), 
                                   iCan, 
                                   Str_SepType.Data(), 
                                   kActivateBarrel? lRbarrelTOF:891,
                                   kActivateEndcap? lZendcapTOF:891,
                                   cDeltaT*1e3, 
                                   lBz);
    
    
    
    
    if (rWrite.Contains("eps") ) myCanvas[iCan]->SaveAs( Form("%s.eps", Str_FigFilename[ iCan ].Data() ) );   
    if (rWrite.Contains("pdf") ) myCanvas[iCan]->SaveAs( Form("%s.pdf", Str_FigFilename[ iCan ].Data() ) ); 
    if (rWrite.Contains("png") ) {gStyle->SetImageScaling(3.); myCanvas[iCan]->Update();  myCanvas[iCan]->SaveAs( Form("%s.png", Str_FigFilename[ iCan ].Data() ) ); }
    
    
    
    delete hDummyMap;
        
}// end iCan
    

    if (rWrite.Contains("tex") ){

    // 1 - figure placement
    std::ofstream TeXOutputFile;
    TString Str_TeXFileName(  Form("z%s_timing_AllSep_RbarrelTOF_%4.2fm_ZendcapTOF_%4.2fm_Timing_%2.0fps_WithBfield_%.1fT.tex", 
                                   Str_Exp.Data(), 
                                   kActivateBarrel? lRbarrelTOF:891,
                                   kActivateEndcap? lZendcapTOF:891,
                                   cDeltaT*1e3, 
                                   lBz) );
    
    TeXOutputFile.open ( Str_TeXFileName.Data() );

    std::ifstream TeXtemplate("Tex_Template_2x4Figures.tex");
    std::string readoutLine;
    std::string searchFigString[] = {"fake00.pdf", "fake01.pdf", "fake02.pdf", "fake03.pdf", "fake04.pdf", "fake05.pdf", "fake06.pdf", "fake07.pdf" };

    // 1 - Place Figures
    while( getline(TeXtemplate, readoutLine) ){
        size_t lposChar[8] = {0};
        for(Int_t iFig = 0; iFig< 8; ++iFig){ 
            lposChar[iFig] = readoutLine.find( searchFigString[iFig] );
            if ( lposChar[iFig] != std::string::npos ){                    
                readoutLine.replace( lposChar[iFig], 
                                searchFigString[iFig].length(),
                                Form("{%s}.png", Str_FigFilename[ iFig ].Data() )
                                );
            }// end if match
        }// end for iFig
        TeXOutputFile << readoutLine << std::endl;            
    }// end while
    TeXtemplate.clear();
    TeXtemplate.seekg(0, TeXtemplate.beg);
    
    TeXtemplate  .close();
    TeXOutputFile.close();   
    
    
    
    // 2 - Tune caption and labels
    gSystem->Exec( Form("sed -i 's|EXPMTid|%s|g' %s",           Str_Exp.Data()  , Str_TeXFileName.Data() ) );
    gSystem->Exec( Form("sed -i 's|RBARRELTOFid|%4.2f|g' %s",   lRbarrelTOF     , Str_TeXFileName.Data() ) );
    if(kActivateEndcap)
        gSystem->Exec( Form("sed -i 's|ZENDCAPTOFid|%4.2f|g' %s",   lZendcapTOF     , Str_TeXFileName.Data() ) );
    else 
        gSystem->Exec( Form("sed -i 's|ZENDCAPTOFid|--|g' %s",                        Str_TeXFileName.Data() ) );
    gSystem->Exec( Form("sed -i 's|TIMINGid|%2.0f|g' %s",       cDeltaT*1e3     , Str_TeXFileName.Data() ) );
    gSystem->Exec( Form("sed -i 's|BFIELDid|%0.1f|g' %s",       lBz             , Str_TeXFileName.Data() ) );   
    
    Printf(" ");
    Printf("-> %s now produced ...", Str_TeXFileName.Data() );

    
    }// end if LaTeX
 
 
    Printf(" ");
    stopWatch.Stop();
    stopWatch.Print(); 
  

}





void myLegendSetUp(TLegend *currentLegend,float currentTextSize){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  //currentLegend->SetFillStyle(0);
  currentLegend->SetFillStyle(1001);
  currentLegend->SetFillColor(0);
  currentLegend->SetLineWidth(1);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  // currentLegend->SetNColumns(2);  
  
  return;
}

void myPadSetUp(TPad *currentPad, Bool_t kLog, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  currentPad->SetFillColor(0);
  currentPad->SetBorderMode(0);
  currentPad->SetBorderSize(2);
  if(kLog) currentPad->SetLogy();
  // currentPad->SetTicks(); // in myOptions()
  currentPad->SetFrameBorderMode(0);
  currentPad->SetGridx();
  currentPad->SetGridy();
  
  return;
}


void myHistoSetUp( TH1 *currentHisto         , 
                   Size_t  currentMarkerSize ,
                   Style_t currentMarkerStyle,
                   Color_t currentMarkerColor,
                   Width_t currentLineWidth  , 
                   Style_t currentLineStyle  , 
                   Color_t currentLineColor  ,
                   Style_t currentFillStyle  ,
                   Color_t currentFillColor
                 ){
    
    currentHisto->SetStats(0);
    
    currentHisto->SetMarkerSize(currentMarkerSize);
    currentHisto->SetMarkerStyle(currentMarkerStyle);
    currentHisto->SetMarkerColor(currentMarkerColor);
    currentHisto->SetLineWidth(currentLineWidth);
    currentHisto->SetLineStyle(currentLineStyle);
    currentHisto->SetLineColor(currentLineColor);
    currentHisto->SetFillStyle(currentFillStyle);
    currentHisto->SetFillColor(currentFillColor);
    return;
}


void myGraphSetUp(TGraph *currentGraph, 
                  Size_t  currentMarkerSize  ,
                  Style_t currentMarkerStyle,  
                  Color_t currentMarkerColor,
                  Width_t currentLineWidth  ,           
                  Style_t currentLineStyle  ,     
                  Color_t currentLineColor  ,
                  Style_t currentFillStyle  ,
                  Color_t currentFillColor  
                 ){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineWidth(currentLineWidth);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  currentGraph->SetFillStyle(currentFillStyle);
  currentGraph->SetFillColor(currentFillColor);
  return;
}

void myFuncSetUp(TF1* f1, Color_t lColor, Style_t lLineStyle, Width_t lLineWidth){
    f1->SetLineColor(lColor);
    f1->SetLineStyle(lLineStyle);
    f1->SetLineWidth(lLineWidth);    
}

void myOptions(Bool_t kStat){
  // Set gStyle
  int font = 42;
  Bool_t kGrayPalette = 0;
  
  // From plain
  gStyle->Reset("Plain");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetGridColor(kGray);

  if(kGrayPalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleFillColor(kWhite);  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (kStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat("KSiouRMe");
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}





