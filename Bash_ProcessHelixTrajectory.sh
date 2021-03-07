


# void Root_DrawHelixTrajectory(          TString Str_Exp = "CMS", 
#                                         Double_t    pT         = 100.0,  // [GeV/c]
#                                         Double_t    eta        = 1.61,
#                                         Int_t       qCharge    = 1,    // e units
#                                         UInt_t      lPart1     = kPDG_piPm,
#                                         Double_t rBmag               = -1.0,
#                                         Double_t rInnerRadiusBarrel  = -1.0, 
#                                         TString rWrite = "0"   

root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 4.20, 1.2, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[2.45GeVc]_eta[1.20]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 2.45, 1.2, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[2.45GeVc]_eta[1.20]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 1.00, 1.2, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[1.00GeVc]_eta[1.20]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 1.00, 1.2, -1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pNeg]_pT[1.00GeVc]_eta[1.20]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 0.65, 2.0, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[0.65GeVc]_eta[2.00]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 0.50, 0.6, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[0.50GeVc]_eta[0.60]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
root -l -b -q Root_DrawHelixTrajectory.C++g'("CMS", 0.50, 0.9, +1, kPDG_p, -3.8, 1.16, "png+pdf")' | tee zTrajectory_CMS_Particle[pPos]_pT[0.50GeVc]_eta[0.90]_RbarrelTOF[1.16m]_TimingRes[30ps]_WithBfield[-3.8T].log
