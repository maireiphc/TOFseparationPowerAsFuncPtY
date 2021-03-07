
# void Root_TOFseparationPowerAsFuncPtY(  TString Str_Exp = "CMS", 
#                                         Double_t rInnerRadiusBarrel  = -1.0, 
#                                         Double_t rBmag               =  0.0,
#                                         Double_t rcDeltaT            = -1.0,
#                                         TString rWrite = "0",
#                                         Double_t rZendcapTOF         = -1.0,
#                                         Bool_t kLog = 1

echo "Launching time :"
date 

root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("CMS",      -1 , -3.8, 0.03,   "png+tex")' | tee zCMS_timing_Log01_RbarrelTOF_1.16m_ZendcapTOF_3.04m_Timing_30ps_WithBfield_3.8T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("CMS",      -1 , -1.9, 0.03,   "png+tex")' | tee zCMS_timing_Log02_RbarrelTOF_1.16m_ZendcapTOF_3.04m_Timing_30ps_WithBfield_1.9T.log 

root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   1.00, -2.0, 0.032,  "png+tex")' | tee zATLAS_timing_Log01_RbarrelTOF_1.00m_ZendcapTOF_3.45m_Timing_32ps_WithBfield_2.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   1.00, -1.0, 0.032,  "png+tex")' | tee zATLAS_timing_Log02_RbarrelTOF_1.00m_ZendcapTOF_3.45m_Timing_32ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   0.29, -1.0, 0.032,  "png+tex")' | tee zATLAS_timing_Log03_RbarrelTOF_0.29m_ZendcapTOF_3.45m_Timing_32ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   0.29, -2.0, 0.032,  "png+tex")' | tee zATLAS_timing_Log04_RbarrelTOF_0.29m_ZendcapTOF_3.45m_Timing_32ps_WithBfield_2.0T.log


root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   1.00, -2.0, 0.020,  "png+tex", 1.5)' | tee zATLAS_timing_Log05_RbarrelTOF_1.00m_ZendcapTOF_1.50m_Timing_20ps_WithBfield_2.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   1.00, -1.0, 0.020,  "png+tex", 1.5)' | tee zATLAS_timing_Log06_RbarrelTOF_1.00m_ZendcapTOF_1.50m_Timing_20ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   0.29, -1.0, 0.020,  "png+tex", 1.5)' | tee zATLAS_timing_Log07_RbarrelTOF_0.29m_ZendcapTOF_1.50m_Timing_20ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ATLAS",   0.29, -2.0, 0.020,  "png+tex", 1.5)' | tee zATLAS_timing_Log08_RbarrelTOF_0.29m_ZendcapTOF_1.50m_Timing_20ps_WithBfield_2.0T.log




root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-1",  -1 , -0.5, 0.080,  "png+tex")' | tee zALICE-1_timing_Log01_RbarrelTOF_3.80m_ZendcapTOF_NaNm_Timing_80ps_WithBfield_0.5T.log # run I
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-1",  -1 , -0.5, 0.056,  "png+tex")' | tee zALICE-1_timing_Log02_RbarrelTOF_3.80m_ZendcapTOF_NaNm_Timing_56ps_WithBfield_0.5T.log # run II Pb-Pb, pp
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-1",  -1 , -0.2, 0.056,  "png+tex")' | tee zALICE-1_timing_Log03_RbarrelTOF_3.80m_ZendcapTOF_NaNm_Timing_56ps_WithBfield_0.2T.log # run II Xe-Xe 

root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.5, 0.030,  "png+tex")' | tee zALICE-3_timing_Log01_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.5, 0.020,  "png+tex")' | tee zALICE-3_timing_Log02_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.5, 0.010,  "png+tex")' | tee zALICE-3_timing_Log03_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.2, 0.030,  "png+tex")' | tee zALICE-3_timing_Log04_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.2, 0.020,  "png+tex")' | tee zALICE-3_timing_Log05_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -0.2, 0.010,  "png+tex")' | tee zALICE-3_timing_Log06_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -1.0, 0.030,  "png+tex")' | tee zALICE-3_timing_Log07_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -1.0, 0.020,  "png+tex")' | tee zALICE-3_timing_Log08_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 1.00, -1.0, 0.010,  "png+tex")' | tee zALICE-3_timing_Log09_RbarrelTOF_1.00m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.5, 0.030,  "png+tex")' | tee zALICE-3_timing_Log10_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.5, 0.020,  "png+tex")' | tee zALICE-3_timing_Log11_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.5, 0.010,  "png+tex")' | tee zALICE-3_timing_Log12_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_0.5T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.2, 0.030,  "png+tex")' | tee zALICE-3_timing_Log13_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.2, 0.020,  "png+tex")' | tee zALICE-3_timing_Log14_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -0.2, 0.010,  "png+tex")' | tee zALICE-3_timing_Log15_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_0.2T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -1.0, 0.030,  "png+tex")' | tee zALICE-3_timing_Log16_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_30ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -1.0, 0.020,  "png+tex")' | tee zALICE-3_timing_Log17_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_20ps_WithBfield_1.0T.log
root -l -b -q Root_TOFseparationPowerAsFuncPtY.C++g'("ALICE-3", 0.20, -1.0, 0.010,  "png+tex")' | tee zALICE-3_timing_Log18_RbarrelTOF_0.20m_ZendcapTOF_2.00m_Timing_10ps_WithBfield_1.0T.log

echo "Finishing time :"
date 
