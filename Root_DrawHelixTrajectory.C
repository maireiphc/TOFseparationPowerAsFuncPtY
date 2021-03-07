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
#include "TGraph2D.h"
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

void myGraph2DSetUp( TGraph2D  *currentGraph   = nullptr, 
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




void Root_DrawHelixTrajectory(          TString Str_Exp = "CMS", 
                                        Double_t    pT         = 100.0,  // [GeV/c]
                                        Double_t    eta        = 1.61,
                                        Int_t       qCharge    = 1,    // e units
                                        UInt_t      lPart1     = kPDG_piPm,
                                        Double_t rBmag               = -1.0,
                                        Double_t rInnerRadiusBarrel  = -1.0, 
                                        TString rWrite = "0"                                                   
){
    
      myOptions();
  gROOT->ForceStyle();

  TDatime now;
  int iDate = now.GetDate();
  int iYear=iDate/10000;
  int iMonth=(iDate%10000)/100;
  int iDay=iDate%100;
  TString cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"};


    TCanvas *myCanvas[30]       = {nullptr};
    TPad    *myPad   [30]       = {nullptr};
    TLegend *myLegend[30]       = {nullptr};
    TH3C    *myBlankHisto[30]   = {nullptr};
    
    TLatex*     lLatex_Hyp      [ 30 ] = {nullptr};
    TPaveText*  lPaveTxt_Info   [ 30 ] = {nullptr};
    
    TString      Str_FigFilename[ 30 ] = {nullptr};
  
      const Int_t lNbCanvas = 1;
      const Bool_t kLog = 0;
      
      // Display N lines of 3 columns of canvas
    for(Int_t iCan = 0; iCan < lNbCanvas; iCan++){
            myCanvas[iCan] = new TCanvas( Form("myCanvas%02d", iCan ), Form("Canvas %02d - [%i %s %i]",iCan, iDay, cMonth[iMonth-1].Data(), iYear), (iCan%3)*750, TMath::Floor(iCan/3)*440, 930, 700 );
            myPadSetUp(myCanvas[iCan],kFALSE, 0.0,0.04, 0.04,0.2);
            myCanvas[iCan]->cd();
            myCanvas[iCan]->Draw();
            
            myPad[iCan] = new TPad( Form("myPad%02d", iCan), Form("pad - %02d", iCan),0,0,1,1);
            myPadSetUp(myPad[iCan],kFALSE, 0.13,0.15,0.35,0.1);   
            myPad[iCan]->Draw();      
            
            myPad[iCan]->Update();
            myCanvas[iCan]->ToggleEventStatus();
            myCanvas[iCan]->ToggleEditor();
            myCanvas[iCan]->ToggleToolBar();
    }


  
  myPad[0]->cd();    
    
    
    
    double lBmag            = 0.;  // [Tesla]
    double lRbarrelTOF      = 0.;
    double lEtaMaxBarrel    = 0.;
    double lZendcapTOF      = 0.;
    double lEtaMinEndcap    = 0.;
    double lEtaMaxEndcap    = 0.;
    double lZmaxBarrel      = 0.;
    double lRmaxEndcap      = 0.;
    double lRminEndcap      = 0.;
    double cDeltaT          = 0.;  // time resolution in ns ! e.g. 0.03 ns  = 30 ps for CMS
    TString StrBarrelTOFname("");
    TString StrEndcapTOFname("");
    
    TString Str_LHCrun("");
    Double_t rcDeltaT       = -1.0;
        
    
    if( Str_Exp == "CMS" ){
        if( TMath::Abs(rBmag) > 0 )  lBmag = rBmag;     // Which field map exists for CMS  B = 1.9 T ?
        else                         lBmag = -3.8  ;     // default CMS = 3.8
        lRbarrelTOF     = 1.161;     // radial position on central barrel BTL, CERN-LHCC-2019-003, Section 1.4.1
        lEtaMaxBarrel   = 1.5  ;
        lZendcapTOF     = 3.04 ;     // figure from Wei Li in the original code, consistent with CERN-LHCC-2019-003, Section 1.4.2
        lEtaMinEndcap   = 1.6  ;
        lEtaMaxEndcap   = 3.0  ;
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.03 ;     // 30 ps Fig 1.5 TDR MTD CMS, CERN-LHCC-2019-003
        Str_LHCrun      = "HL-LHC run IV";
        StrBarrelTOFname= "BTL";
        StrEndcapTOFname= "ETL";
    }
    else if( Str_Exp == "ATLAS" ){
        if( TMath::Abs(rBmag) > 0 )  lBmag = rBmag;     // Which field map exists for ATLAS ?
        else                         lBmag = -2.0  ;     // default ATLAS solenoid = 2.0
        if(rInnerRadiusBarrel > 0)      lRbarrelTOF = rInnerRadiusBarrel;
        else                            lRbarrelTOF = 0.29;                 // r ~ 29 cm instead of pixel outer barrel or r = 100 cm instead of last strip ? To be tested
              // FIXME : hyp to be made for a potential SPAD layer instead of a tracker layer in ITk
              // FIXME : hyp, let's assume a SPAD layer (beware radiation tolerance as f(r)... See Seminar CERN detector, ITk, Kuehn, 2019/11/29, https://indico.cern.ch/event/865308/)   
        lEtaMaxBarrel   = 1.0  ;
        lZendcapTOF     = 3.45 ;      // HGTD LGad, https://cds.cern.ch/record/2623663
        lEtaMinEndcap   = 2.4  ;
        lEtaMaxEndcap   = 4.0  ;  
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.032;      
        Str_LHCrun      = "HL-LHC run IV(+V?)";
        StrBarrelTOFname= "#kern[-0.2]{bTOF}^{#it{Hyp.}}";
        StrEndcapTOFname= "HGTD"; 
    }
    else if( Str_Exp == "ALICE-1" ){
        if( TMath::Abs(rBmag) > 0 )  lBmag = rBmag; // existing field map for B(L3) = 0.2 T
        else                         lBmag = -0.5  ; // default L3 0.5 T
        lRbarrelTOF     = 3.8   ;     // Ravg(TOF) = 380 cm
        lEtaMaxBarrel   = 0.88  ;
        lZendcapTOF     = 3.7  ;     // inexistent... dummy
        lEtaMinEndcap   = 2.5   ;
        lEtaMaxEndcap   = 4.0   ;
        if(rcDeltaT > 0)  cDeltaT = rcDeltaT;
        else              cDeltaT = 0.080 ;     // 56 ps in Pb-Pb run II to be checked with Bologna about pp MB ? See https://arxiv.org/abs/1806.03825        80 ps Run I, https://arxiv.org/pdf/1402.4476.pdf, p.51
        Str_LHCrun      = "LHC run II";
        StrBarrelTOFname= "TOF";
        StrEndcapTOFname= "None";        
    }
    else if( Str_Exp == "ALICE-3" ){
        if( TMath::Abs(rBmag) > 0 )  lBmag = rBmag; // existing field map for B(L3) = 0.5 or 0.2 T
        else                         lBmag = -0.5  ; // default L3 0.5 T
        if(rInnerRadiusBarrel > 0)      lRbarrelTOF = rInnerRadiusBarrel;    // Ravg(TOF) = 0.2 m, 0.3 m  inner SPAD layer = to be tested
        else                            lRbarrelTOF = 1.00;                  // Ravg(TOF) ~ 100 cm    outer SPAD layer
        lEtaMaxBarrel   = 1.4  ;   
        lZendcapTOF     = 2.00 ;      // 2.0 m for ~std config, 
        lEtaMinEndcap   = 1.4  ;      // µ spectro acceptance, idea : dipole to be used ?
        lEtaMaxEndcap   = 4.0  ;
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
    
    
    
    
    
    
    
    
    Printf("-- General settings ---------------------");        
    Printf("    Magnetic field B = %4.2f T", lBmag            );
    Printf(" ");
    Printf("    Rbarrel(TOF)     = %4.2f m", lRbarrelTOF      );
    Printf("    Max eta_Barrel   = %4.2f",   lEtaMaxBarrel    );
    Printf("    Max Z_barrel     = %5.3f m", lZmaxBarrel      );
    Printf(" ");
    Printf("    zEndcap(TOF)     = %4.2f m", lZendcapTOF      );
    Printf("    Min eta_Endcap   = %4.2f",   lEtaMinEndcap    );
    Printf("    Max eta_Endcap   = %4.2f",   lEtaMaxEndcap    );
    Printf("    Max R_Endcap     = %5.3f m", lRmaxEndcap      );    
    Printf("    Min R_Endcap     = %5.3f m", lRminEndcap      );
    Printf(" ");
    Printf("    TOF Timing Resol = %02f ps", cDeltaT*1e3      );
    Printf("    data period      : %s",      Str_LHCrun.Data()       );
    Printf("    Barrel TOF name  : %s",      StrBarrelTOFname.Data() );
    Printf("    Endcap TOF name  : %s",      StrEndcapTOFname.Data() );

  
  myPad[0]->cd();

  myBlankHisto[0] = new TH3C("myBlankHisto[0]","Blank Histogram",   60,-lRbarrelTOF,lRbarrelTOF,   60,-lRbarrelTOF,lRbarrelTOF,   80, -0.1,lZendcapTOF);
  myBlankHisto[0]->SetXTitle("X coordinate (m)");
  myBlankHisto[0]->SetYTitle("Y coordinate (m)");
  myBlankHisto[0]->SetZTitle("Z coordinate (m)");
  myBlankHisto[0]->SetNdivisions(510,"xyz");
  
  myBlankHisto[0]->GetZaxis()->SetRangeUser(0, 10);
  myBlankHisto[0]->Draw();
  

     
    
    
    

    Double_t m0InGeV        = mPDG[lPart1]; // [GeV/c²]
    
    Double_t pTot           = pT * TMath::CosH( eta );
    Double_t pZ             = pT * TMath::SinH( eta );
    Double_t Etot           = TMath::Sqrt( pTot*pTot + m0InGeV*m0InGeV );
    Double_t betaZ          = pZ/Etot;      // vZ with c=1, then in fact, fraction of c, so betaZ
    Double_t betaT          = pT/Etot;      // vT with c=1, then in fact, fraction of c, so betaT
    
    
    Double_t gammaTot       = TMath::Sqrt( 1 + pTot/m0InGeV * pTot/m0InGeV );
    Double_t gammaT         = TMath::Sqrt( 1 + pT/m0InGeV * pT/m0InGeV );
    Double_t gammaZ         = TMath::Sqrt( 1 + pZ/m0InGeV * pZ/m0InGeV );
    Double_t betaTot        = pTot/Etot;
    
    Double_t m0InKg         = m0InGeV * 1.672649e-27 / mPDG[kPDG_p];  // needed for the cyclotron frequency, use proton mass as reference
    Double_t eCharge        = 1.602176565e-19;  // [Coulomb]
    Double_t lB2c_GeV       = TMath::C() / 1e9 ;// ~0.3 = c/ (1eV) = c/(1e-9) GeV
    
    
    Double_t lRho           = pT/(lB2c_GeV * qCharge *lBmag);  // curvature radius of the track, [m]   
//     Double_t AngFreq_Tot    = (TMath::Abs(qCharge)*eCharge) * TMath::Abs(lBmag) / (gammaTot  * m0InKg); // in rad.s-1, from var with international system unit q*e (C), m (kg), etc
//     Double_t AngFreq_T      = (TMath::Abs(qCharge)*eCharge) * TMath::Abs(lBmag) / (gammaT    * m0InKg); // in rad.s-1, from var with international system unit q*e (C), m (kg), etc
    Double_t AngFreq_Tot    = (qCharge*eCharge) * lBmag / (gammaTot  * m0InKg); // in rad.s-1, from var with international system unit q*e (C), m (kg), etc
    Double_t AngFreq_T      = (qCharge*eCharge) * lBmag / (gammaT    * m0InKg); // in rad.s-1, from var with international system unit q*e (C), m (kg), etc
    Double_t AngFreq        = AngFreq_Tot;
    
    Double_t pathLength      =  0.;
    UShort_t lAbsNbTurns     =  0 ;
    UShort_t lAbsNbHalfTurns =  0 ;
    
    
    Int_t lNbHelixPnts = 500;
    
    
    Double_t l_tf       = -1; // final timing, when the trajectory should cross a TOF layer or endcap (given in seconds)
    Double_t l_tf_WeiLi = -1; // final timing, when the trajectory should cross a TOF layer or endcap (given in seconds)
    Double_t l_t0 = 0.;
    
    
    // Check minimal acceptance in pT and pZ
    
    // Will the particle have the necessary pT to reach the TOF barrel radius
    // Is the sinus angle defined, i.e. sin in [-1;1] i.e. angle between -pi/2 and pi/2 ?
    //   With pT(t0) along x here, we lean towards the edge +1 (rotation >0 ...)
    //   so TMath::ASin will simply silently return pi/2...
    
    

    
    Bool_t IsPartForeseenInBarrel        = kFALSE;
    Bool_t IsPartForeseenInEndcap        = kFALSE;
    Bool_t IsPartForeseenBeyondEndcap    = kFALSE;
    TString Str_KineCase("");
    

        // NOTE : beware the order of the cases tested matter to be exhaustive and avoid wrong case assignment (test on pZ Vs transverse reach of the track) !
    if( (TMath::Abs(eta) > lEtaMaxEndcap) ){
        l_tf       = lZendcapTOF / (betaZ   * TMath::C() ) ;    // pZ ≠ necessarily different from 0
        l_tf_WeiLi = lZendcapTOF / (betaTot * TMath::C() * TMath::TanH(eta) );
        IsPartForeseenBeyondEndcap = kTRUE;
        Str_KineCase = "Case 1 - Very forward particle near beam pipe (eta < min endcap at t0)";
    }
    else if( (TMath::Abs(eta) > lEtaMinEndcap) && (TMath::Abs(eta) < lEtaMaxEndcap)  ){
        l_tf       = lZendcapTOF / (betaZ   * TMath::C() ) ;    // pZ ≠ necessarily different from 0
        l_tf_WeiLi = lZendcapTOF / (betaTot * TMath::C() * TMath::TanH(eta) );
        IsPartForeseenInEndcap = kTRUE;
        Str_KineCase = "Case 2 - Forward particle in endcap (pT, eta in acceptance at t0)";
    }
    
    
    else if(    TMath::Abs(lRbarrelTOF/2. * 1./lRho) > 1.0  &&  pZ > 0)    {// 1. looper track at y ≠ 0, Rho too short, pT too small 
        l_tf       = TMath::Abs( lZendcapTOF / (betaZ   * TMath::C() ) );                //      -> let's consider the trajectory up to endcap
        l_tf_WeiLi = TMath::Abs( lZendcapTOF / (betaTot * TMath::C() * TMath::TanH(eta) ) );
        IsPartForeseenInEndcap = kTRUE;
        Str_KineCase = "Case 3 - looper at y ≠ 0 (small pT + eta ≠ 0)";
    }                                                                        
    else if(    TMath::Abs(lRbarrelTOF/2. * 1./lRho) > 1.0  &&  pZ < 1e-13){// 2. looper track at y = 0, Rho too short, pT too small 
        l_tf       = TMath::Abs( 1./AngFreq * TMath::TwoPi() );                           //      -> let's consider the description of the 1st entire loop only, one period...
        l_tf_WeiLi = TMath::Abs( TMath::TwoPi() * lRho * TMath::CosH(eta)/ (betaTot * TMath::C() ) );
        IsPartForeseenInBarrel = kTRUE;
        Str_KineCase = "Case 4 - looper at y = 0 (small pT + eta = 0)";
    }
    else if(  TMath::Abs(lRbarrelTOF/2. * 1./lRho) < 1.0  &&  pZ > 0 )     {// 3. Normal case y ≠ 0, we reach TOF barrel
        l_tf       = TMath::Abs( 1./AngFreq * 2. * TMath::ASin( lRbarrelTOF/2. * 1./lRho  ) );
        l_tf_WeiLi = TMath::Abs( TMath::ACos( 1. - lRbarrelTOF*lRbarrelTOF/(2*lRho*lRho)  ) * lRho * TMath::CosH(eta)/ ( betaTot * TMath::C() ) );
        IsPartForeseenInBarrel = kTRUE;        
        Str_KineCase = "Case 5 - Std Barrel at y ≠ 0 (pT in acceptance / eta ≠ 0)";
    }
    else if(  TMath::Abs(lRbarrelTOF/2. * 1./lRho) < 1.0  &&  pZ < 1e-13 ) {// 4. Normal case y = 0, we reach TOF barrel
        l_tf       = TMath::Abs( 1./AngFreq * 2. * TMath::ASin( lRbarrelTOF/2. * 1./lRho  ) );
        l_tf_WeiLi = TMath::Abs( TMath::ACos( 1. - lRbarrelTOF*lRbarrelTOF/(2*lRho*lRho)  ) * lRho * 1.              / ( betaTot * TMath::C() ) );
        IsPartForeseenInBarrel = kTRUE;        
        Str_KineCase = "Case 6 - Specific Barrel at y = 0 (pT in acceptance / eta = 0)";
    }        

         
    // Hack : to force propagation up to TOF endcap, if need be for investigations on loopers    
    Bool_t kForcePropagationToEndcap = kFALSE;
    if(kForcePropagationToEndcap)  l_tf       = lZendcapTOF / (betaZ   * TMath::C() ) ;
    
    
    if( !IsPartForeseenInBarrel && !IsPartForeseenInEndcap && !IsPartForeseenBeyondEndcap ){
        Printf("Problem : particle foreseen nowhere in the detector ! ... return.");
        return;
    }
        
    
    std::vector<Double_t> xHelix;    xHelix .resize(lNbHelixPnts+1); 
    std::vector<Double_t> yHelix;    yHelix .resize(lNbHelixPnts+1);
    std::vector<Double_t> zHelix;    zHelix .resize(lNbHelixPnts+1);
    

    Printf(" "); 
    Printf("-- Specific quantities ---------------------");     
    Printf("    pT          = %4.2f GeV/c",     pT);
    Printf("    qCharge     = %d e units",      qCharge);
//  Printf("           ");  
    Printf("    lBmag       = %4.2f T",         lBmag     );
    Printf("    m0InGeV     = %6.4g GeV/c",     m0InGeV    );
    Printf("    eta         = %4.2f",           eta        );
    Printf("    y           = %4.2f",           1./2*TMath::Log( (Etot + pZ) / (Etot - pZ) ) );
    Printf(" ");            
    Printf("    pTot        = %6.4g GeV/c",     pTot       );
    Printf("    Etot (sqrt(pTot²+m²) = %6.4g GeV",       Etot       );
    Printf("    Etot (gamma.m)       = %6.4g GeV",       gammaTot*m0InGeV );
    Printf("    betaTot     = %10.8g (/c)",     betaTot    );
    Printf("    pZ          = %6.4g GeV/c (via sinh(eta) )",  pZ         );
    Printf("    betaZ       = %10.8g (/c)",     betaZ         );
    Printf("    betaT       = %10.8g (/c)",     betaT         );
    Printf(" ");            
    Printf("    vz          : betaZ.c          = %12.10e m.s-1", betaZ * TMath::C()               );
    Printf("    vz          : pZ.c/gammaTot.m0 = %12.10e m.s-1", pZ* TMath::C() / gammaTot/m0InGeV);
    Printf("    vz          : pZ.c/gammaZ.m0   = %12.10e m.s-1", pZ* TMath::C() / gammaZ  /m0InGeV);
    Printf(" ");  
    Printf("    vT          : betaT.c          = %12.10e m.s-1", betaT * TMath::C()               );
    Printf("    vT          : pT.c/gammaTot.m0 = %12.10e m.s-1", pT* TMath::C() / gammaTot/m0InGeV);
    Printf("    vT          : AngFreq_Tot.rho  = %12.10e m.s-1", AngFreq_Tot * lRho               );
    Printf("    vT          : pT.c/gammaT.m0   = %12.10e m.s-1", pT* TMath::C() / gammaT  /m0InGeV);
    Printf("    vT          : AngFreq_T  .rho  = %12.10e m.s-1", AngFreq_T   * lRho               );
    
    
    
    Printf(" ");            
    Printf("    gammaTot            = %10.8g",           gammaTot   );
    Printf("    gammaT              = %10.8g",           gammaT     );
    Printf("    gammaZ              = %10.8g",           gammaZ     );
    Printf("    gammaT*gammaZ       = %10.8g",           gammaT*gammaZ   );
    Printf(" ");            
    Printf("    1/gammaTot          = %10.8g",           1./gammaTot);
    Printf("    1/gammaT            = %10.8g",           1./gammaT  );
    Printf("    1/gammaZ            = %10.8g",           1./gammaZ  );    
    Printf("    1/gammaZ*1/gammaT   = %10.8g",           1./gammaZ * 1./gammaT );    
    Printf(" ");            
    Printf("    1/betaTot   = %14.12g (via betaTot)",        1./betaTot);
    Printf("    1/betaTot   = %14.12g (via m0, pTot)",       TMath::Sqrt( 1 + m0InGeV / pTot * m0InGeV / pTot) ); 
    
    
    Printf(" ");
    Printf("    m0InKg          = %6.4g kg",        m0InKg     );
    Printf("    eCharge         = %10.8g C",        eCharge    );
    Printf("    lB2c_GeV        = %11.9g",          lB2c_GeV   );
    Printf(" ");
    
    Printf("    lRho            = %6.4g m",         lRho);
    Printf("    AngFreq_Tot     = %6.4g rad.s-1",   AngFreq_Tot); 
    Printf("    AngFreq_T       = %6.4g rad.s-1",   AngFreq_T  ); 
    Printf("    -> AngFreq      = %6.4g rad.s-1",   AngFreq    );      
    Printf(" ");
    Printf("    Particle foreseen in : Barrel = %d / Endcap = %d / Beyond Endcap = %d", IsPartForeseenInBarrel, IsPartForeseenInEndcap, IsPartForeseenBeyondEndcap);
    Printf("    Kinematic case  = %s",               Str_KineCase.Data() );
    Printf("    -> tf           = %14.12g s",       l_tf              );
    Printf("    -> tf(Wei)      = %14.12g s",       l_tf_WeiLi        );
    Printf("    -> tf / tf(Wei) = %14.12f  ",       l_tf / l_tf_WeiLi );
    
    
    Printf(" ");
    
    
       

    
    
    Double_t lPhit0 = 0.;  // Initial azimuthal angle, by default = 0,  could be anything in [0, 2pi] : TMath::Pi()/2. , ...
            // FIXME TBD : ok for the trajectory
            //             but still to be cured for the other displays : centre of curvature location, tan(t0)...
            //             The physics quantities : path legnth, TOF time, straight lengths, Rt, rotation angle ... should be left unchanged by this rotation
            // NOTE the helix will be rotated in block around the z axis, like a door on its hinges, by the side
            //             so that in terms of acceptance the initial phi(t0) has NO effect, if the detector is invariant by rotation
            //             the final TOF hit will end up always in the same circle, same radius, same z, just at a different location on it.
    
    
    TGraph2D *gHelix = new TGraph2D(lNbHelixPnts+1);    gHelix->SetName("gHelix");
    
    for(Int_t iStep = 0; iStep < lNbHelixPnts+1; ++iStep){
       Double_t iTime = (l_tf - l_t0)/lNbHelixPnts * iStep;
       
       if( (iStep >0) && TMath::Floor( TMath::Abs( AngFreq * iTime *TMath::RadToDeg()/360.) ) > lAbsNbTurns)                    lAbsNbTurns++; 
            // everytime the angle change the multiple of 360° (= modulus+1), i.e. enter a new turn, increase the absolute turn counter
            // Beware the sign of AngFreq and the direction of the rotation , + or -...

       if( (iStep >0) && TMath::Floor( TMath::Abs( AngFreq * iTime *TMath::RadToDeg()/180.) ) > lAbsNbHalfTurns)                lAbsNbHalfTurns++;
            // everytime the angle change the multiple of 180°, i.e. enter a new hemisphere, increase the absolute half-turn counter
            // Beware the sign of AngFreq and the direction of the rotation , + or -...
                   


        xHelix[iStep] = lRho *   TMath::Sin( AngFreq * iTime )     ;
        yHelix[iStep] = lRho * ( TMath::Cos( AngFreq * iTime ) - 1);
        zHelix[iStep] = betaZ* TMath::C() * iTime ;
        
        
        Double_t lTmpX = xHelix[iStep];
        Double_t lTmpY = yHelix[iStep];
        
        
        xHelix[iStep] = lTmpX * TMath::Cos( lPhit0 ) - lTmpY * TMath::Sin( lPhit0 );  // trigo formula cos(a+b) : x = r. cos (phit) -> x' = r. cos( phit + phi0 ) -> x' = x. cos(phi0) - y. sin(phi0)
        yHelix[iStep] = lTmpY * TMath::Cos( lPhit0 ) + lTmpX * TMath::Sin( lPhit0 );  // trigo formula sin(a+b) : y = r. sin (phit) -> x' = r. sin( phit + phi0 ) -> y' = y. cos(phi0) + y. sin(phi0)       
            
            
        
        
        
        if(iStep >0) pathLength += TMath::Sqrt( TMath::Power(xHelix[iStep] - xHelix[iStep-1], 2) +  
                                                TMath::Power(yHelix[iStep] - yHelix[iStep-1], 2) +
                                                TMath::Power(zHelix[iStep] - zHelix[iStep-1], 2)  );
       
     gHelix->SetPoint(iStep, xHelix[iStep], yHelix[iStep], zHelix[iStep] );
     
     if(!(iStep%100) || iStep == lNbHelixPnts) 
            Printf("Helix trajectory at [t = %10.8e s] : x[%04d] = %10.8f, y[%04d] = %10.8f, z[%04d] = %10.8f", iTime, iStep, xHelix[iStep], iStep, yHelix[iStep], iStep, zHelix[iStep]);
              
    }// for gHelix
    

    
    
    TGraph2D *gArcCircleRPhi = new TGraph2D(lNbHelixPnts+1);    gHelix->SetName("gArcCircleRPhi");
    
    for(Int_t iStep = 0; iStep < lNbHelixPnts+1; ++iStep){
        gArcCircleRPhi->SetPoint(iStep, xHelix[iStep], yHelix[iStep], 0. );
    }
    
    
    
    TGraph2D *gCurvatureCentre = new TGraph2D(8);
        gCurvatureCentre->SetName("gCurvatureCentre");
        
        gCurvatureCentre->SetPoint(0,  0,                    0,                    0);
        gCurvatureCentre->SetPoint(1,  0,                    -lRho,                0);      // NOTE : beware the opposite sign of Rho, needed here
        gCurvatureCentre->SetPoint(2,  0,                    -lRho,                zHelix[lNbHelixPnts]);
        gCurvatureCentre->SetPoint(3,  xHelix[lNbHelixPnts], yHelix[lNbHelixPnts], zHelix[lNbHelixPnts]);
        gCurvatureCentre->SetPoint(4,  xHelix[lNbHelixPnts], yHelix[lNbHelixPnts], 0);
        gCurvatureCentre->SetPoint(5,  0,                    0,                    0);
        gCurvatureCentre->SetPoint(6,  0,                    0,                    zHelix[lNbHelixPnts]);
        gCurvatureCentre->SetPoint(7,  xHelix[lNbHelixPnts], yHelix[lNbHelixPnts], zHelix[lNbHelixPnts]);



    Double_t lRt_tf =  TMath::Sqrt( xHelix[lNbHelixPnts] * xHelix[lNbHelixPnts] + yHelix[lNbHelixPnts] * yHelix[lNbHelixPnts]);
    Double_t lZ_tf  =  zHelix[lNbHelixPnts];


        
    TGraph2D *gStraightLineTanInO = new TGraph2D(3);        // _OptA : straight line with ptot(t0) prolongated to TOF (barrel or endcap)
        gStraightLineTanInO->SetName("gStraightLineTanInO");
        
        
        Double_t lRtofFromInitialEta = 0;
        Double_t lZtofFromInitialEta = 0;
        if(IsPartForeseenInBarrel){
            lRtofFromInitialEta = lRt_tf;                           // Rt_tf = here lRbarrelTOF       
            lZtofFromInitialEta = lRt_tf * TMath::SinH( eta ) ;     // eq. e4 App. D or F            
        }
        else if(IsPartForeseenInEndcap || IsPartForeseenBeyondEndcap){
            lRtofFromInitialEta = lZ_tf * 1./TMath::SinH( eta );    // eq e3 App. D;
            lZtofFromInitialEta = lZ_tf;                            // lZ_tf = here lZendcapTOF
        }
        
        gStraightLineTanInO->SetPoint(0,  0,                    0,  0);
        gStraightLineTanInO->SetPoint(1,  lRtofFromInitialEta,  0,  lZtofFromInitialEta ); // eq. e4 App. D or F; Rt_tf = here lRbarrelTOF       
        gStraightLineTanInO->SetPoint(2,  lRtofFromInitialEta,  0,  0);
        

        
    TGraph2D *gStraightLineOtoTOFhit = new TGraph2D(2);      // _OptB : staight line between (x,y,z at tf = last point = TOF hit) - origin (0,0,0)
        gStraightLineOtoTOFhit->SetName("gStraightLineOtoTOFhit");
        
        gStraightLineOtoTOFhit->SetPoint(0,  0,                    0,                    0);
        gStraightLineOtoTOFhit->SetPoint(1,  xHelix[lNbHelixPnts], yHelix[lNbHelixPnts], zHelix[lNbHelixPnts]);        
        
        
                                            
    Double_t lStraightLength_OptA = 0;  // straight line with ptot(t0) prolongated to TOF (barrel or endcap)
        if(IsPartForeseenInBarrel)                                      lStraightLength_OptA = lRt_tf * TMath::CosH( eta ); //  eq. e5 App. F
        else if(IsPartForeseenInEndcap || IsPartForeseenBeyondEndcap)   lStraightLength_OptA = lZ_tf  / TMath::TanH( eta ); //  eq. e8 App. F = curled path length !! = pathLength1 here                                      
    
    
    Double_t lStraightLength_OptB = TMath::Sqrt( xHelix[lNbHelixPnts] * xHelix[lNbHelixPnts] + 
                                            yHelix[lNbHelixPnts] * yHelix[lNbHelixPnts] +
                                            zHelix[lNbHelixPnts] * zHelix[lNbHelixPnts] ); // length : staight line between (x,y,z at tf = last point = TOF hit) - origin (0,0,0)

    


        
    Double_t lTransvRotAngle_Barrel = 2 * TMath::ASin( lRt_tf     /2. *1/(-lRho) );    // NOTE : beware = signed angle here, in accordance with right-angled frame (0;x,y,z)
        // NOTE : should be valid for both TOF barrel and endcap... 
        // NOTE : It is just that, in the barrel case, lRt_tf should be = lRbarrelTOF ...
        // NOTE NOTE : but ... not valid if several half-turns ! This is only valid within [0, pi] = via visible final and concrete Rtf
        //                  i.e. valid for a track that stays in the 1st quadrant, phi in [0;pi/2] = rotation around centre between [0;pi]
        //                  That is typically true for that track meant to TOF barrel layer.
        //                  But may not always be true for endcap (>180° or loopers, etc)
        //              Hence the need to compute a specific angle for endcap...
    
    
    Double_t lTransvRotAngle_Endcap = lZ_tf/(-lRho * TMath::SinH(eta)); // e22 in South face p.9
    
    Double_t lTransvRotAngle = 0;
    
    if(IsPartForeseenInBarrel)                                      lTransvRotAngle = lTransvRotAngle_Barrel;
    else if(IsPartForeseenInEndcap || IsPartForeseenBeyondEndcap)   lTransvRotAngle = lTransvRotAngle_Endcap;
    
    Float_t lSignedRot       = -AngFreq/ TMath::Abs( AngFreq ); // ±1, positive or negative rotation ? 
    Int_t lSignedNbTurns     = lSignedRot * lAbsNbTurns;
    Int_t lSignedNbHalfTurns = lSignedRot * lAbsNbHalfTurns;
    
                                          
    Printf(" ");                                      
    Printf(" -- TOF hit");
    Printf("         Rt(tf)                     = %10.8g m (last point in helix)", lRt_tf );
    Printf("         z(tf)                      = %10.8g m (last point in helix)", lZ_tf );
    Printf("         z(tf)                      = %10.8g m (equation)",            lTransvRotAngle*(-lRho)* TMath::SinH(eta) );
    Printf("         1. Angle (~barrel)"           );
    Printf("            TransvRotAngle          = %8.6g rad",           lTransvRotAngle_Barrel);
    Printf("            TransvRotAngle          = %8.6g °",             lTransvRotAngle_Barrel * TMath::RadToDeg() );
    Printf("         2. Angle (~endcap)"           );
    Printf("            TransvRotAngle          = %8.6g rad",           lTransvRotAngle_Endcap );
    Printf("            TransvRotAngle          = %8.6g °",             lTransvRotAngle_Endcap * TMath::RadToDeg() );
    Printf("         -> Angle final choice"        );
    Printf("            TransvRotAngle          = %8.6g °",             lTransvRotAngle        * TMath::RadToDeg()  );
    Printf("         Signed Rotation            = %+1f",                lSignedRot);
    Printf("         Nb of full turns           = %02d turn(s)",        lSignedNbTurns);
    Printf("         Nb of half-turns           = %02d half-turn(s)",   lSignedNbHalfTurns);
    Printf("         TransvRotAngle [0;2pi]     = %8.6g °",             (lTransvRotAngle - lSignedNbTurns * TMath::TwoPi())* TMath::RadToDeg() );
    Printf("         'final' eta                = %8.6f",               TMath::ASinH( lZ_tf / (lTransvRotAngle * -lRho) ) );  // e22 in South Face
    Printf("         Apparent final 'eta'       = %8.6f",               TMath::ASinH( lZ_tf / lRt_tf )                    );  // e4  in Appendix D
    Printf(" ");    
    
    Double_t lCheck_Recompute_tf    = pathLength/(betaTot * TMath::C()); 
    Double_t lCheck_StraightA_tf     = lStraightLength_OptA / (betaTot * TMath::C());
    Double_t lCheck_StraightB_tf     = lStraightLength_OptB / (betaTot * TMath::C());
    
    
    
    Double_t pathLength1 = lZ_tf /  TMath::TanH(eta);
    Double_t pathLength2 = betaTot * TMath::C() * l_tf;
    Double_t pathLength3 = lTransvRotAngle * -lRho * TMath::CosH(eta); // e20 south face p.7
    
    Printf(" -- Path lengths");
    Printf("        0. Integrated pathLength (int dl)                 = %14.12f m",   pathLength);
    Printf("        1. Integrated pathLength [z(tf) / tanh(eta)]      = %14.12f m",   pathLength1 );    
    Printf("        2. Integrated pathLength [vTot * tf]              = %14.12f m",   pathLength2 );   
    Printf("        3. Integrated pathLength [TrRotAng.rho.cosh(eta)] = %14.12f m",   pathLength3 );
    Printf("        -> ratio option2 / option0                        = %14.12f",     pathLength2 / pathLength);
    Printf("        -> ratio option3 / option0                        = %14.12f",     pathLength3 / pathLength);
    Printf(" ");    
    Printf("        Straight line (Origin to TOF hit)                 = %14.12f m",   lStraightLength_OptB);
    Printf("        -> Ratio of Length (lineB Vs Curled)              = %14.12f  ",     lStraightLength_OptB / pathLength);
    Printf("        -> Delta in Length (lineB  - Curled)              = %14.12g m",     lStraightLength_OptB - pathLength);
    Printf(" ");    
    Printf("        Straight line (tan(t0) prolong. to TOF)           = %14.12f m",   lStraightLength_OptA);
    Printf("        -> Ratio of Length (lineA Vs Curled)              = %14.12f  ",     lStraightLength_OptA / pathLength);
    Printf("        -> Delta in Length (lineA  - Curled)              = %14.12g m",     lStraightLength_OptA - pathLength);    
    Printf("        -> Ratio of Length (lineA Vs lineB)               = %14.12f  ",     lStraightLength_OptA / lStraightLength_OptB);
    Printf("        -> Delta in Length (lineA  - lineB)               = %14.12g m",     lStraightLength_OptA - lStraightLength_OptB);
    
    Printf(" ");    
    Printf(" -- Final time");
    Printf("        TOF time tf    (reminder, from equation)          = %14.12g s",   l_tf);
    Printf("        TOF time tf    (curled [int dl]/betaTot)          = %14.12g s",   lCheck_Recompute_tf );
    Printf("        -> Ratio of tf (l_tf/lCheck_Recompute_tf)         = %14.12f",     l_tf/ lCheck_Recompute_tf);
    Printf(" ");    
    Printf("        TOF time tf    (reminder, from equation)          = %14.12g s",   l_tf);
    Printf("        TOF time tf(B) (straight L1/betaTot)              = %14.12g s",   lCheck_StraightB_tf);
    Printf("        -> Ratio of tf (StraightB Vs Curled)              = %14.12f  ",   lCheck_StraightB_tf / l_tf);
    Printf("        -> Delta in tf (StraightB  - Curled)              = %14.12g s",   lCheck_StraightB_tf - l_tf);
    Printf(" ");  
    Printf("        TOF time tf(A) (straight L2/betaTot)              = %14.12g s",   lCheck_StraightA_tf);
    Printf("        -> Ratio of tf (StraightA Vs Curled)              = %14.12f  ",   lCheck_StraightA_tf / l_tf);
    Printf("        -> Delta in tf (StraightA  - Curled)              = %14.12g s",   lCheck_StraightA_tf - l_tf);
    Printf("        -> Ratio of tf (StraightA Vs StraightB)           = %14.12f  ",   lCheck_StraightA_tf / lCheck_StraightB_tf);
    Printf("        -> Delta in tf (StraightA  - StraightB)           = %14.12g s",   lCheck_StraightA_tf - lCheck_StraightB_tf);


    Int_t  lPhiSamplingDense  = 360;
    Int_t  lPhiSamplingMild   = 72;
    UInt_t lZSampling         = 15; 


    std::vector<TGraph2D *> gTOFbarrel;           gTOFbarrel .resize( lZSampling +1 );
    
    for(UInt_t iZ = 0; iZ < lZSampling +1; ++iZ){
        gTOFbarrel[iZ] = new TGraph2D( lPhiSamplingDense );
        gTOFbarrel[iZ]->SetName( Form("gTOFbarrel[%d]", iZ) );
    }

    TGraph2D *gTOFendcapRmax = new TGraph2D( lPhiSamplingMild );    gTOFendcapRmax->SetName("gTOFendcapRmax");
    TGraph2D *gTOFendcapRmin = new TGraph2D( lPhiSamplingMild );    gTOFendcapRmax->SetName("gTOFendcapRmin");
       
  
    std::vector<Double_t> xTOFbarrel;    xTOFbarrel .resize(lPhiSamplingDense); 
    std::vector<Double_t> yTOFbarrel;    yTOFbarrel .resize(lPhiSamplingDense);
    std::vector<Double_t> zTOFbarrel;    zTOFbarrel .resize(lPhiSamplingDense);
    
    
    std::vector<Double_t> xTOFendcap;    xTOFendcap .resize(lPhiSamplingMild); 
    std::vector<Double_t> yTOFendcap;    yTOFendcap .resize(lPhiSamplingMild);
    std::vector<Double_t> zTOFendcap;    zTOFendcap .resize(lPhiSamplingMild);    

    
   
    
    // Barrel geometry
    for(Int_t iPnt = 0; iPnt < lPhiSamplingDense; ++iPnt){
        
        Double_t lAngle = iPnt * 360./lPhiSamplingDense - 50; // start angle at -50° to have subsequent points filled in TGraph2D
        
        if(lAngle > 160 && lAngle < 320) continue; // FIXME
        
        xTOFbarrel[iPnt] = lRbarrelTOF * TMath::Cos(lAngle * TMath::DegToRad() );
        yTOFbarrel[iPnt] = lRbarrelTOF * TMath::Sin(lAngle * TMath::DegToRad() );
        
        // Double_t lCurrentTurn = TMath::Floor( (lAngle) / (2*TMath::Pi()*TMath::RadToDeg()) ); // angle/360° = current nb of completed turn : 0, 1, 3
        // zTOFbarrel[iPnt] = lCurrentTurn * lZmaxBarrel/lZSampling ; 

            
       for(UInt_t iZ = 0; iZ < lZSampling +1; ++iZ){
           
           if(iZ == 0)  zTOFbarrel[iPnt] =  1e-6;
           else         zTOFbarrel[iPnt] = iZ * lZmaxBarrel/lZSampling ;
           // Printf("[%02d] : x=%+6.4f, y=%+6.4f, z=%+6.4f", iPnt, xTOFbarrel[iPnt], yTOFbarrel[iPnt], zTOFbarrel[iPnt]);    
           gTOFbarrel[iZ]->SetPoint(iPnt, xTOFbarrel[iPnt], yTOFbarrel[iPnt], zTOFbarrel[iPnt] );
       }
    }// end phi sampling
    
    
    
    // Endcap geometry
    for(Int_t iPnt = 0; iPnt < lPhiSamplingMild; ++iPnt){ 
        
        xTOFendcap[iPnt] = TMath::Cos(iPnt * 360./lPhiSamplingMild * TMath::DegToRad() );
        yTOFendcap[iPnt] = TMath::Sin(iPnt * 360./lPhiSamplingMild * TMath::DegToRad() );
        zTOFendcap[iPnt] = lZendcapTOF;
                
        gTOFendcapRmax->SetPoint(iPnt, lRmaxEndcap * xTOFendcap[iPnt], lRmaxEndcap * yTOFendcap[iPnt], zTOFendcap[iPnt] );
        gTOFendcapRmin->SetPoint(iPnt, lRminEndcap * xTOFendcap[iPnt], lRminEndcap * yTOFendcap[iPnt], zTOFendcap[iPnt] );     
    }    
    
    
    // Eta highlights        
    std::vector<Double_t>   vEtaValInTOF {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, lEtaMaxBarrel, lEtaMinEndcap, 2.0, 2.5, lEtaMaxEndcap};
    
    std::vector<TGraph2D *> gEtaMapInTOF;           gEtaMapInTOF        .resize( vEtaValInTOF.size() );
    
    for(UInt_t iEta = 0; iEta < vEtaValInTOF.size(); ++iEta){
        gEtaMapInTOF[iEta] = new TGraph2D( lPhiSamplingDense );
        gEtaMapInTOF[iEta]->SetName( Form("gEtaMapInTOF[%d]", iEta) );
    }

    
    
    
    // Eta highlights
    for(Int_t iPnt = 0; iPnt < lPhiSamplingMild; ++iPnt){// Eta highlights in barrel 
        
        Double_t lAngle = iPnt * 360./lPhiSamplingMild;
        
        if(lAngle > 90 && lAngle < 330) continue;
        
        for(UInt_t iEta = 0; iEta < 8; ++iEta){
            
            Double_t xTOFmap = lRbarrelTOF * TMath::Cos(lAngle * TMath::DegToRad() );
            Double_t yTOFmap = lRbarrelTOF * TMath::Sin(lAngle * TMath::DegToRad() );
            Double_t zTOFmap = lRbarrelTOF * TMath::SinH(  vEtaValInTOF[iEta] );
                        
            gEtaMapInTOF[iEta]->SetPoint(iPnt, xTOFmap, yTOFmap, zTOFmap);            
        }
    }// end : Eta highlights in barrel

    
    
    for(Int_t iPnt = 0; iPnt < lPhiSamplingMild; ++iPnt){// Eta highlights in endcap
        
        //if(iPnt * 360./lPhiSamplingMild > 140 && iPnt * 360./lPhiSamplingMild < 360) continue;
        
        for(UInt_t iEta = 8; iEta < vEtaValInTOF.size(); ++iEta){
                        
            Double_t yTOFmap = lZendcapTOF / TMath::SinH(  vEtaValInTOF[iEta] ) * TMath::Sin(iPnt * 360./lPhiSamplingMild * TMath::DegToRad() );
            Double_t xTOFmap = lZendcapTOF / TMath::SinH(  vEtaValInTOF[iEta] ) * TMath::Cos(iPnt * 360./lPhiSamplingMild * TMath::DegToRad() );
            Double_t zTOFmap = lZendcapTOF;
                        
            gEtaMapInTOF[iEta]->SetPoint(iPnt, xTOFmap, yTOFmap, zTOFmap);            
        }  
    }// end : Eta highlights in endcap
    
    
    
    
    TGraph2D *gAxisXbottom = new TGraph2D(2);   gAxisXbottom->SetName("gAxisXbottom");      gAxisXbottom->SetPoint(0,   -lRbarrelTOF, 0  , 0);        gAxisXbottom->SetPoint(1,  lRbarrelTOF, 0  , 0);   
    TGraph2D *gAxisXup     = new TGraph2D(2);   gAxisXup    ->SetName("gAxisXup"    );      gAxisXup    ->SetPoint(0,   -lRbarrelTOF, 0  , lZ_tf);    gAxisXup    ->SetPoint(1,  lRbarrelTOF, 0  , lZ_tf);
    TGraph2D *gAxisYbottom = new TGraph2D(2);   gAxisYbottom->SetName("gAxisYbottom");      gAxisYbottom->SetPoint(0,    0,  -lRbarrelTOF, 0);        gAxisYbottom->SetPoint(1,  0,   lRbarrelTOF, 0); 
    TGraph2D *gAxisYup     = new TGraph2D(2);   gAxisYup    ->SetName("gAxisYup");          gAxisYup    ->SetPoint(0,    0,  -lRbarrelTOF, lZ_tf);    gAxisYup    ->SetPoint(1,  0,   lRbarrelTOF, lZ_tf); 
    TGraph2D *gAxisZ       = new TGraph2D(2);   gAxisZ      ->SetName("gAxisZ");            gAxisZ      ->SetPoint(0,    0,   0  ,         -0.2);     gAxisZ      ->SetPoint(1,  0,   0  , lZendcapTOF);   



    
    
    myGraph2DSetUp( gAxisXbottom,     0.7, kFullCircle,    kGreen+2,   1, kDashed, kGreen+2,   0, kWhite);
    myGraph2DSetUp( gAxisXup,         0.7, kFullCircle,    kGreen+2,   1, kDashed, kGreen+2,   0, kWhite);
    myGraph2DSetUp( gAxisYbottom,     0.7, kFullCircle,    kAzure+2,   1, kDashed, kAzure+2,   0, kWhite);
    myGraph2DSetUp( gAxisYup,         0.7, kFullCircle,    kAzure+2,   1, kDashed, kAzure+2,   0, kWhite);
    myGraph2DSetUp( gAxisZ,           0.7, kFullCircle,    kRed+1,     2, kDashed, kRed+1,     0, kWhite);
    
    myGraph2DSetUp( gHelix,             0.6, kFullCircle,    kGreen+2  ,   1, kSolid , kGray+2,    0, kWhite );
    myGraph2DSetUp( gArcCircleRPhi,     0.6, kDot,           kGreen+2  ,   2, kDashed, kGreen+2,   0, kWhite ); 
    
    myGraph2DSetUp( gCurvatureCentre,       0.6, kFullSquare,    kGreen+2  ,   1, kDashed, kGreen+2 ,   0, kWhite );
    myGraph2DSetUp( gStraightLineOtoTOFhit, 0.6, kFullSquare,    kYellow+1 ,   2, 9      , kYellow+1,   0, kWhite );
    myGraph2DSetUp( gStraightLineTanInO,    0.6, kFullCircle,    kAzure+1  ,   2, kDotted, kAzure+1 ,   0, kWhite );
    
    
    
    
    
    
    
    for(UInt_t iZ = 0; iZ < lZSampling +1; ++iZ){    
        myGraph2DSetUp( gTOFbarrel[iZ],     1.0, kDot,  kGray ,    1, kDotted, kGray,      0, kWhite);    
    }
    myGraph2DSetUp( gTOFendcapRmax, 1.0, kDot,  kOrange+7, 1, kDotted, kOrange+7,  0, kWhite); 
    myGraph2DSetUp( gTOFendcapRmin, 1.0, kDot,  kOrange+7, 1, kDotted, kOrange+7,  0, kWhite);
    
    
    
    gAxisXbottom->Draw("P0 line same");
    gAxisYbottom->Draw("P0 line same");
    //if(IsPartForeseenInEndcap) {
        gAxisXup    ->Draw("line same");
        gAxisYup    ->Draw("line same");
    //}
    gAxisZ      ->Draw("P0 line same");
    
    for(UInt_t iZ = 0; iZ < lZSampling +1; ++iZ){
        gTOFbarrel[iZ]->Draw("P same");
    }
    gTOFendcapRmax->Draw("P line same");
    gTOFendcapRmin->Draw("P line same");
 
    for(UInt_t iEta = 0; iEta < vEtaValInTOF.size(); ++iEta){
        myGraph2DSetUp( gEtaMapInTOF[iEta], 1.3, kFullDotSmall,  kPink, 1, kSolid, kPink,  0, kWhite); 
        

        // Display only maps in concerned areas, barrel or (exclusive) endcap (for visibility)
        if(         IsPartForeseenInBarrel 
            &&  TMath::Abs(vEtaValInTOF[iEta]) < lEtaMaxBarrel )    gEtaMapInTOF[iEta]->Draw("P same");                                 

        else if(    IsPartForeseenInEndcap
            && (TMath::Abs(vEtaValInTOF[iEta]) > lEtaMinEndcap )
            && (TMath::Abs(vEtaValInTOF[iEta]) < lEtaMaxEndcap ))       gEtaMapInTOF[iEta]->Draw("P same");
        
        else if(   IsPartForeseenBeyondEndcap
            && (TMath::Abs(vEtaValInTOF[iEta]) > lEtaMaxEndcap ))       gEtaMapInTOF[iEta]->Draw("P same");
        
    }
    
    
    
    gCurvatureCentre        ->Draw("P line same");
    
    gArcCircleRPhi          ->Draw("line same");
    gHelix                  ->Draw("P line same");
    gStraightLineOtoTOFhit  ->Draw("P line same");
    gStraightLineTanInO     ->Draw("P line same");

    
    lLatex_Hyp[ 0 ] = new TLatex();
    lLatex_Hyp[ 0 ]->SetNDC();
    lLatex_Hyp[ 0 ]->SetTextSize(0.035);
    lLatex_Hyp[ 0 ]->SetTextColor(colors[red]);
    lLatex_Hyp[ 0 ]->DrawLatex(0.65,0.95, Form("%s #scale[0.7]{configuration [ %s ]}", Str_Exp.Data(), Str_LHCrun.Data()) );

 

    
    TText *lTxt_Info[32];
    


    lPaveTxt_Info[ 0 ] = new TPaveText(0.0, 0.80, 0.29, 0.98, "NDC ARC");
    lPaveTxt_Info[ 0 ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ 0 ]->SetFillColor(kWhite);
    lPaveTxt_Info[ 0 ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ 0 ]->SetLineColor(kWhite);
    lPaveTxt_Info[ 0 ]->SetBorderSize(0);
    

    lTxt_Info[0] = lPaveTxt_Info[ 0 ]->AddText( Form("#color[%d]{Particle :} %s (#scale[0.7]{#it{m}_{PDG} = %7.5f GeV/c^{2}})", colors[red], Str_HadronLatexName[lPart1].Data(), m0InGeV ) );
    lTxt_Info[0]->SetTextSize(0.03);
    lTxt_Info[0]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[0]->SetTextColor(kGray+3);

        
    lTxt_Info[1] = lPaveTxt_Info[ 0 ]->AddText( Form("  #it{p}_{T} = %4.2f GeV/#it{c}", pT ) );
    lTxt_Info[1]->SetTextSize(0.03);
    lTxt_Info[1]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[1]->SetTextColor(kGray+3);
    
    lTxt_Info[2] = lPaveTxt_Info[ 0 ]->AddText( Form("  #eta = %4.2f", eta ) );
    lTxt_Info[2]->SetTextSize(0.03);
    lTxt_Info[2]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[2]->SetTextColor(kGray+3);
    
    lTxt_Info[3] = lPaveTxt_Info[ 0 ]->AddText( Form("  charge = %+01d #it{e} unit", qCharge ) );
    lTxt_Info[3]->SetTextSize(0.03);
    lTxt_Info[3]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[3]->SetTextColor(kGray+3);

    lTxt_Info[4] = lPaveTxt_Info[ 0 ]->AddText( Form("  #rightarrow #it{#beta}_{tot} #approx %10.8f", betaTot) );
    lTxt_Info[4]->SetTextSize(0.025);
    lTxt_Info[4]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[4]->SetTextColor(kGray+3);  

    
    
    lPaveTxt_Info[ 1 ] = new TPaveText(0.53, 0.85, 0.75, 0.92, "NDC ARC");
    lPaveTxt_Info[ 1 ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ 1 ]->SetFillColor(kWhite);
    lPaveTxt_Info[ 1 ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ 1 ]->SetLineColor(kWhite);
    lPaveTxt_Info[ 1 ]->SetBorderSize(0);

    lTxt_Info[5] = lPaveTxt_Info[ 1 ]->AddText( Form("Solenoidal B_{z} field : %+4.2g T", lBmag) );
    lTxt_Info[5]->SetTextSize(0.025);
    lTxt_Info[5]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[5]->SetTextColor(kGray+3);

    lTxt_Info[6] = lPaveTxt_Info[ 1 ]->AddText( Form("#sigma^{TOF}(#scale[0.7]{#color[%d]{EndCap} or #color[%d]{Barrel}}) : %2.0f ps", kOrange+7, kPink, cDeltaT*1000) ); // 807 = kOrange+7; 900 = kPink       
    lTxt_Info[6]->SetTextSize(0.025);
    lTxt_Info[6]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[6]->SetTextColor(kGray+3); 
    
    
    lPaveTxt_Info[ 2 ] = new TPaveText(0.77, 0.85, 0.98, 0.92, "NDC ARC");
    lPaveTxt_Info[ 2 ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ 2 ]->SetFillColor(kWhite);
    lPaveTxt_Info[ 2 ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ 2 ]->SetLineColor(kWhite);
    lPaveTxt_Info[ 2 ]->SetBorderSize(0);
    
    
    
    lTxt_Info[7] = lPaveTxt_Info[ 2 ]->AddText( Form("<#it{r}>^{TOF}(Barrel) #kern[+4.3]{:} %4.2f m", lRbarrelTOF));       
    lTxt_Info[7]->SetTextSize(0.025);
    lTxt_Info[7]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[7]->SetTextColor(kPink);  

    if(lZendcapTOF > 0)  lTxt_Info[8] = lPaveTxt_Info[ 2 ]->AddText( Form("<#it{z}>^{TOF}(EndCap) : %4.2f m", lZendcapTOF) );       
    else                 lTxt_Info[8] = lPaveTxt_Info[ 2 ]->AddText( "<#it{z}>^{TOF}(EndCap) : --");
    lTxt_Info[8]->SetTextSize(0.025);
    lTxt_Info[8]->SetTextAlign(kHAlignLeft+kVAlignCenter);
    lTxt_Info[8]->SetTextColor(kOrange+7);

    
    TLine *lSepLine1 = new TLine();
    lSepLine1->SetLineColor(kGray);
    lSepLine1->SetLineStyle(kDashed);
    lSepLine1->SetLineWidth(1);
    lSepLine1->DrawLineNDC(0.55,0.84,0.98,0.84);
    
    
 
    
    
    lPaveTxt_Info[ 3 ] = new TPaveText(0.70, 0.55, 0.98, 0.82, "NDC");
    lPaveTxt_Info[ 3 ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ 3 ]->SetFillColor(kWhite);
    lPaveTxt_Info[ 3 ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ 3 ]->SetLineColor(kGray);
    lPaveTxt_Info[ 3 ]->SetBorderSize(0);
    
    
    Color_t colorTOFhit = 0;
    
    if(       IsPartForeseenInBarrel )       colorTOFhit = kPink;
    else if ( IsPartForeseenInEndcap )       colorTOFhit = kOrange+7;
    else if ( IsPartForeseenBeyondEndcap )   colorTOFhit = kGray+2;
    
    
    lTxt_Info[ 9] = lPaveTxt_Info[ 3 ]->AddText("TOF hit at #it{t}_{f} :");
    lTxt_Info[10] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.2]{#it{x} #approx %+6.4f m}",      xHelix[lNbHelixPnts]) );
    lTxt_Info[11] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.2]{#it{y} #approx %+6.4f m}",      yHelix[lNbHelixPnts]) );
    lTxt_Info[12] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.2]{#it{z} #approx %+6.4f m}",      zHelix[lNbHelixPnts]) );
    lTxt_Info[13] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.18]{#it{R}_{T} #approx %6.4f m}",  lRt_tf              ) );
    lTxt_Info[14] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.21]{Full turn(s) = %d}",     lSignedNbTurns) );
    lTxt_Info[15] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.09]{Rotation #approx %+6.3f #circ (#scale[0.8]{%+6.3f #circ})}", lTransvRotAngle * TMath::RadToDeg(), (lTransvRotAngle-lSignedNbTurns*TMath::TwoPi())*TMath::RadToDeg() ) );
    lTxt_Info[16] = lPaveTxt_Info[ 3 ]->AddText( Form("#kern[+0.21]{#rho #approx %+6.4f m}",     lRho) );

    for(Int_t iTxt = 9; iTxt < 17; ++iTxt){
        lTxt_Info[iTxt]->SetTextSize(0.025);
        lTxt_Info[iTxt]->SetTextAlign(kHAlignLeft+kVAlignCenter);
        if(iTxt == 9) lTxt_Info[iTxt]->SetTextColor(kGray+3);
        else lTxt_Info[iTxt]->SetTextColor(colorTOFhit); 
    }


    
    TLine *lSepLine2 = new TLine();
    lSepLine2->SetLineColor(kGray);
    lSepLine2->SetLineStyle(kSolid);
    lSepLine2->SetLineWidth(1);
    lSepLine2->DrawLineNDC(0.70,0.54,0.98,0.54);
    


    
    
        
    lPaveTxt_Info[ 4 ] = new TPaveText(0.70, 0.02, 0.98, 0.53, "NDC");
    lPaveTxt_Info[ 4 ]->SetCornerRadius(0.1);
    lPaveTxt_Info[ 4 ]->SetFillColor(kWhite);
    lPaveTxt_Info[ 4 ]->SetShadowColor(kWhite);
    lPaveTxt_Info[ 4 ]->SetLineColor(kGray);
    lPaveTxt_Info[ 4 ]->SetBorderSize(0);
    
    lTxt_Info[17] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{curled path, #it{#font[12]{L}} #kern[+0.14]{#approx %6.4f m}}", kGreen+2, betaTot * TMath::C() * l_tf) );
    lTxt_Info[18] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{straight path, #it{a}} #approx %6.4f m",   kAzure+1,  lStraightLength_OptA     ) );
    lTxt_Info[19] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{#font[12]{L} / a} #kern[+0.2]{#approx %6.4f}", (betaTot * TMath::C() * l_tf) / lStraightLength_OptA) );    
    lTxt_Info[20] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{#font[12]{L}} - #it{a} #kern[+0.16]{#approx %6.4f m}", (betaTot * TMath::C() * l_tf) - lStraightLength_OptA)); 

    lTxt_Info[21] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{straight path, #it{b}} #approx %6.4f m",   kYellow+1, lStraightLength_OptB     ) );    
    lTxt_Info[22] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{#font[12]{L} / b} #kern[+0.2]{#approx %6.4f}", (betaTot * TMath::C() * l_tf) / lStraightLength_OptB) );    
    lTxt_Info[23] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{#font[12]{L}} - #it{b} #kern[+0.16]{#approx %6.4f m}", (betaTot * TMath::C() * l_tf) - lStraightLength_OptB));        
    
    lTxt_Info[24] = lPaveTxt_Info[ 4 ]->AddText( " " );    
    lTxt_Info[25] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{#it{t}_{f} along #it{#font[12]{L}} #approx %6.4g ps}", kGreen+2, l_tf *1e12) );
    lTxt_Info[26] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{#it{t}_{f} along #it{a}} #approx %6.4g ps",  kAzure+1, lCheck_StraightA_tf *1e12 ) );
    lTxt_Info[27] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{t}_{f}(#it{#font[12]{L}}) / #it{t}_{f}(#it{a}) #approx %6.4g", l_tf / lCheck_StraightA_tf) );
    lTxt_Info[28] = lPaveTxt_Info[ 4 ]->AddText( Form("   #scale[0.9]{#it{t}_{f}(#it{#font[12]{L}}) - #it{t}_{f}(#it{a}) #approx %5.3g ps #approx %4.2f #sigma^{TOF}}", (l_tf - lCheck_StraightA_tf)*1e12, (l_tf -lCheck_StraightA_tf)*1e12/ (cDeltaT*1000) ));
    
    lTxt_Info[29] = lPaveTxt_Info[ 4 ]->AddText( Form("#color[%d]{#it{t}_{f} along #it{b}} #approx %6.4g ps", kYellow+1, lCheck_StraightB_tf *1e12 ) );
    lTxt_Info[30] = lPaveTxt_Info[ 4 ]->AddText( Form("   #it{t}_{f}(#it{#font[12]{L}}) / #it{t}_{f}(#it{b}) #approx %6.4g", l_tf / lCheck_StraightB_tf) );
    lTxt_Info[31] = lPaveTxt_Info[ 4 ]->AddText( Form("   #scale[0.9]{#it{t}_{f}(#it{#font[12]{L}}) - #it{t}_{f}(#it{b}) #approx %5.3g ps #approx %4.2f #sigma^{TOF}}", (l_tf - lCheck_StraightB_tf)*1e12, (l_tf -lCheck_StraightB_tf)*1e12/ (cDeltaT*1000) ));
    
    
    
    for(Int_t iTxt = 17; iTxt < 32; ++iTxt){
        lTxt_Info[iTxt]->SetTextSize(0.025);
        lTxt_Info[iTxt]->SetTextAlign(kHAlignLeft+kVAlignCenter);
        lTxt_Info[iTxt]->SetTextColor(kGray+3);
    }
        
   

    
    
    
    lPaveTxt_Info[ 0 ]->Draw();
    lPaveTxt_Info[ 1 ]->Draw();
    lPaveTxt_Info[ 2 ]->Draw();
    lPaveTxt_Info[ 3 ]->Draw();
    lPaveTxt_Info[ 4 ]->Draw();


    
    Str_FigFilename[0] = Form("zTrajectory_%s_Particle[%s%s]_pT[%4.2fGeVc]_eta[%4.2f]_RbarrelTOF[%4.2fm]_TimingRes[%2.0fps]_WithBfield[%.1fT]", 
                                    Str_Exp.Data(), 
                                    Str_HadronASCIIName[lPart1].Data(), 
                                    qCharge>0?"Pos":"Neg",
                                    pT,
                                    eta,
                                    lRbarrelTOF, 
                                    cDeltaT*1e3, 
                                    lBmag);
    
    
    if (rWrite.Contains("eps") ) myCanvas[0]->SaveAs( Form("%s.eps", Str_FigFilename[0].Data() ) );   
    if (rWrite.Contains("pdf") ) myCanvas[0]->SaveAs( Form("%s.pdf", Str_FigFilename[0].Data() ) ); 
    if (rWrite.Contains("png") ) myCanvas[0]->SaveAs( Form("%s.png", Str_FigFilename[0].Data() ) ); 
    
            
    
    
    
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
  //currentLegend->SetNColumns(2);  
  
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


void myGraphSetUp(TGraphErrors *currentGraph, 
                  Size_t currentMarkerSize  ,
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


void myGraph2DSetUp(TGraph2D *currentGraph, 
                  Size_t currentMarkerSize  ,
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
  gStyle->SetLabelSize(0.035,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.8,"xyz");  
  gStyle->SetTitleSize(0.04,"xyz");
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

