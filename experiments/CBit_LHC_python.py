from .ColliderBit_analysis import Analysis

# # All
# "ATLAS_13TeV_4LEP_36invfb",
# "ATLAS_13TeV_RJ3L_2Lep2Jets_36invfb",
# "ATLAS_13TeV_RJ3L_3Lep_36invfb",
# "ATLAS_13TeV_MultiLEP_2Lep0Jets_36invfb",
# "ATLAS_13TeV_MultiLEP_2LepPlusJets_36invfb",
# "ATLAS_13TeV_MultiLEP_3Lep_36invfb",
# "ATLAS_13TeV_3b_24invfb",
# "ATLAS_13TeV_3b_discoverySR_24invfb",
# "ATLAS_13TeV_3b_36invfb",
# "ATLAS_13TeV_3b_discoverySR_36invfb"
# "CMS_13TeV_1LEPbb_36invfb",
# "CMS_13TeV_MONOJET_36invfb"
# "CMS_13TeV_2LEPsoft_36invfb",
# "CMS_13TeV_2LEPsoft_36invfb_nocovar",
# "CMS_13TeV_MultiLEP_2SSLep_36invfb",
# "CMS_13TeV_MultiLEP_3Lep_36invfb",
# "CMS_13TeV_2OSLEP_36invfb",
# "CMS_13TeV_2OSLEP_36invfb_nocovar"
# 
# # Skip
# "CMS_13TeV_2LEPsoft_36invfb_nocovar",
# "CMS_13TeV_2OSLEP_36invfb_nocovar",
# "ATLAS_13TeV_3b_discoverySR_24invfb",
# "ATLAS_13TeV_3b_discoverySR_36invfb",
# "ATLAS_13TeV_3b_36invfb",
# "CMS_13TeV_MONOJET_36invfb"
# 
# # Final:
# "ATLAS_13TeV_4LEP_36invfb",
# "ATLAS_13TeV_RJ3L_2Lep2Jets_36invfb",
# "ATLAS_13TeV_RJ3L_3Lep_36invfb",
# "ATLAS_13TeV_3b_24invfb",
# "ATLAS_13TeV_MultiLEP_2Lep0Jets_36invfb",
# "ATLAS_13TeV_MultiLEP_2LepPlusJets_36invfb",
# "ATLAS_13TeV_MultiLEP_3Lep_36invfb",
# "CMS_13TeV_1LEPbb_36invfb",
# "CMS_13TeV_2LEPsoft_36invfb",
# "CMS_13TeV_MultiLEP_2SSLep_36invfb",
# "CMS_13TeV_MultiLEP_3Lep_36invfb",
# "CMS_13TeV_2OSLEP_36invfb",

allow_corr = False # Use the analyses with correlation matrices (much faster to leave them out for testing)
 
analyses = []

a = Analysis("ATLAS_13TeV_3b_24invfb")
a.SR_names = ["meff160_ETmiss0__i0", "meff160_ETmiss20__i1", "meff200_ETmiss0__i2", "meff200_ETmiss20__i3", "meff200_ETmiss45__i4", "meff200_ETmiss70__i5", "meff260_ETmiss0__i6", "meff260_ETmiss20__i7", "meff260_ETmiss45__i8", "meff260_ETmiss70__i9", "meff260_ETmiss100__i10", "meff340_ETmiss0__i11", "meff340_ETmiss20__i12", "meff340_ETmiss45__i13", "meff340_ETmiss70__i14", "meff340_ETmiss100__i15", "meff340_ETmiss150__i16", "meff340_ETmiss200__i17", "meff440_ETmiss0__i18", "meff440_ETmiss20__i19", "meff440_ETmiss45__i20", "meff440_ETmiss70__i21", "meff440_ETmiss100__i22", "meff440_ETmiss150__i23", "meff440_ETmiss200__i24", "meff560_ETmiss0__i25", "meff560_ETmiss20__i26", "meff560_ETmiss45__i27", "meff560_ETmiss70__i28", "meff560_ETmiss100__i29", "meff560_ETmiss150__i30", "meff560_ETmiss200__i31", "meff700_ETmiss0__i32", "meff700_ETmiss20__i33", "meff700_ETmiss45__i34", "meff700_ETmiss70__i35", "meff700_ETmiss100__i36", "meff700_ETmiss150__i37", "meff700_ETmiss200__i38", "meff860_ETmiss0__i39", "meff860_ETmiss20__i40", "meff860_ETmiss45__i41", "meff860_ETmiss70__i42", "meff860_ETmiss100__i43", "meff860_ETmiss150__i44", "meff860_ETmiss200__i45", ]
a.SR_n     = [20, 3, 1503, 1137, 65, 0, 1329, 2877, 951, 150, 2, 373, 873, 444, 164, 40, 3, 0, 121, 304, 170, 62, 31, 3, 1, 40, 95, 75, 20, 15, 2, 2, 17, 30, 22, 12, 6, 2, 2, 2, 7, 10, 5, 2, 4, 1, ]
a.SR_b     = [16.21, 0.6503, 1480, 1088, 58.05, 0.2691, 1297, 2860, 991, 149.4, 2.024, 390.1, 884.6, 472.6, 171.1, 36.24, 1.457, 0.006531, 130.3, 310.8, 176.6, 65.1, 22.16, 3.895, 0.4816, 43.46, 102.6, 68.03, 30.72, 14.13, 2.358, 1.08, 13.56, 32.67, 23.78, 12.47, 5.549, 1.728, 0.8551, 2.816, 7.766, 8.968, 4.297, 2.785, 0.9345, 0.4297, ]
a.SR_b_sys = [0.11, 0.0747, 26, 7, 0.39, 0.0547, 8, 36, 6.5, 1, 1.426, 2.6, 13.1, 3, 1.1, 0.24, 0.111, 0.004409, 0.8, 9.5, 1.2, 1.1, 6.03, 0.14, 0.0551, 0.29, 6.6, 0.45, 0.2, 3.19, 1.02, 0.23, 0.09, 3.39, 0.15, 0.08, 0.873, 0.879, 0.1211, 0.246, 2.114, 2.332, 0.335, 0.29, 0.2345, 0.0719, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_3b_36invfb")
a.SR_names = ["SR-3b-meff1-A__i0", "SR-3b-meff2-A__i1", "SR-3b-meff3-A__i2", "SR-4b-meff1-A__i3", "SR-4b-meff1-B__i4", "SR-4b-meff2-A__i5", "SR-4b-meff2-B__i6", ]
a.SR_n     = [4, 3, 0, 1, 2, 1, 0, ]
a.SR_b     = [2.5, 2, 0.8, 0.43, 2.6, 0.43, 1.3, ]
a.SR_b_sys = [1, 0.5, 0.5, 0.31, 0.9, 0.27, 0.6, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("ATLAS_13TeV_3b_discoverySR_24invfb")
a.SR_names = ["low-SR-MET0meff440__i0", "low-SR-MET150meff440__i1", ]
a.SR_n     = [1063, 17, ]
a.SR_b     = [1100, 12, ]
a.SR_b_sys = [25, 8, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("ATLAS_13TeV_3b_discoverySR_36invfb")
a.SR_names = ["SR-4b-meff1-A-disc__i0", ]
a.SR_n     = [2, ]
a.SR_b     = [0.7, ]
a.SR_b_sys = [0.5, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("ATLAS_13TeV_4LEP_36invfb")
a.SR_names = ["SR0A__i0", "SR0B__i1", "SR0C__i2", "SR0D__i3", ]
a.SR_n     = [13, 2, 47, 10, ]
a.SR_b     = [10.2, 1.31, 37, 4.1, ]
a.SR_b_sys = [2.1, 0.24, 9, 0.7, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_MultiLEP_2Lep0Jets_36invfb")
a.SR_names = ["SR2_SF_loose__i0", "SR2_SF_tight__i1", "SR2_DF_100__i2", "SR2_DF_150__i3", "SR2_DF_200__i4", "SR2_DF_300__i5", ]
a.SR_n     = [153, 9, 78, 11, 6, 2, ]
a.SR_b     = [133, 9.8, 68, 11.5, 2.1, 0.6, ]
a.SR_b_sys = [22, 2.9, 7, 3.1, 1.9, 0.6, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_MultiLEP_2LepPlusJets_36invfb")
a.SR_names = ["SR2_int__i0", "SR2_high__i1", "SR2_low__i2", ]
a.SR_n     = [2, 0, 11, ]
a.SR_b     = [4.1, 1.6, 4.2, ]
a.SR_b_sys = [2.6, 1.6, 3.4, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_MultiLEP_3Lep_36invfb")
a.SR_names = ["SR3_slep_a__i0", "SR3_slep_b__i1", "SR3_slep_c__i2", "SR3_slep_d__i3", "SR3_slep_e__i4", "SR3_WZ_0Ja__i5", "SR3_WZ_0Jb__i6", "SR3_WZ_0Jc__i7", "SR3_WZ_1Ja__i8", "SR3_WZ_1Jb__i9", "SR3_WZ_1Jc__i10", ]
a.SR_n     = [4, 3, 9, 0, 0, 21, 1, 2, 1, 3, 4, ]
a.SR_b     = [2.2, 2.8, 5.4, 1.4, 1.1, 21.7, 2.7, 1.6, 2.2, 1.8, 1.3, ]
a.SR_b_sys = [0.8, 0.4, 0.9, 0.4, 0.2, 2.9, 0.5, 0.3, 0.5, 0.3, 0.3, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_RJ3L_2Lep2Jets_36invfb")
a.SR_names = ["2L2JHIGH__i0", "2L2JINT__i1", "2L2JLOW__i2", "2L2JCOMP__i3", ]
a.SR_n     = [0, 1, 19, 11, ]
a.SR_b     = [1.9, 2.4, 8.4, 2.7, ]
a.SR_b_sys = [0.8, 0.9, 5.8, 2.7, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("ATLAS_13TeV_RJ3L_3Lep_36invfb")
a.SR_names = ["3LHIGH__i0", "3LINT__i1", "3LLOW__i2", "3LCOMP__i3", ]
a.SR_n     = [2, 1, 20, 12, ]
a.SR_b     = [1.1, 2.3, 10, 3.9, ]
a.SR_b_sys = [0.5, 0.5, 2, 1, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("CMS_13TeV_1LEPbb_36invfb")
a.SR_names = ["SRA__i0", "SRB__i1", ]
a.SR_n     = [11, 7, ]
a.SR_b     = [7.5, 8.7, ]
a.SR_b_sys = [2.5, 2.2, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("CMS_13TeV_MONOJET_36invfb")
a.SR_names = ["sr-0__i0", "sr-1__i1", "sr-2__i2", "sr-3__i3", "sr-4__i4", "sr-5__i5", "sr-6__i6", "sr-7__i7", "sr-8__i8", "sr-9__i9", "sr-10__i10", "sr-11__i11", "sr-12__i12", "sr-13__i13", "sr-14__i14", "sr-15__i15", "sr-16__i16", "sr-17__i17", "sr-18__i18", "sr-19__i19", "sr-20__i20", "sr-21__i21", ]
a.SR_n     = [136865, 74340, 42540, 25316, 15653, 10092, 8298, 4906, 2987, 2032, 1514, 926, 557, 316, 233, 172, 101, 65, 46, 26, 31, 29, ]
a.SR_b     = [134500, 73400, 42320, 25490, 15430, 10160, 8480, 4865, 2970, 1915, 1506, 844, 526, 325, 223, 169, 107, 88.1, 52.8, 25, 25.5, 26.9, ]
a.SR_b_sys = [3700, 2000, 810, 490, 310, 170, 140, 95, 49, 33, 32, 18, 14, 12, 9, 8, 6, 5.3, 3.9, 2.5, 2.6, 2.8, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("CMS_13TeV_2LEPsoft_36invfb")
a.SR_names = ["SR1__i0", "SR2__i1", "SR3__i2", "SR4__i3", "SR5__i4", "SR6__i5", "SR7__i6", "SR8__i7", "SR9__i8", "SR10__i9", "SR11__i10", "SR12__i11", ]
a.SR_n     = [2, 15, 19, 18, 1, 0, 3, 1, 2, 1, 2, 0, ]
a.SR_b     = [3.5, 12, 17, 11, 1.6, 3.5, 2, 0.51, 1.4, 1.5, 1.5, 1.2, ]
a.SR_b_sys = [1, 2.3, 2.4, 2, 0.7, 0.9, 0.7, 0.52, 0.7, 0.6, 0.8, 0.6, ]
a.cov = [[1.29, 0.33, 0.45, 0.49, 0.06, 0.09, 0.12, 0.08, 0.12, 0.09, 0.07, 0.12],
 [0.33, 5.09, 1.01, 0.62, 0.12, 0.13,  0.2, 0.12, 0.12, 0.11, 0.15, 0.13],
 [0.45, 1.01, 6.44, 0.78, 0.21, 0.19, 0.18,  0.1, 0.18, 0.18, 0.15, 0.19],
 [0.49, 0.62, 0.78,  3.6, 0.09, 0.07, 0.12, 0.19, 0.19, 0.13, 0.17, 0.32],
 [0.06, 0.12, 0.21, 0.09, 0.59, 0.03, 0.06, 0.03, 0.02, 0.03, 0.03, 0.03],
 [0.09, 0.13, 0.19, 0.07, 0.03, 0.72, 0.03, 0.03, 0.03, 0.04, 0.03, 0.01],
 [0.12,  0.2, 0.18, 0.12, 0.06, 0.03,  0.6, 0.05, 0.04, 0.05, 0.04, 0.05],
 [0.08, 0.12,  0.1, 0.19, 0.03, 0.03, 0.05, 0.17, 0.05, 0.03, 0.04, 0.06],
 [0.12, 0.12, 0.18, 0.19, 0.02, 0.03, 0.04, 0.05, 0.26, 0.05, 0.07, 0.07],
 [0.09, 0.11, 0.18, 0.13, 0.03, 0.04, 0.05, 0.03, 0.05, 0.32, 0.05, 0.04],
 [0.07, 0.15, 0.15, 0.17, 0.03, 0.03, 0.04, 0.04, 0.07, 0.05,  0.2, 0.06],
 [0.12, 0.13, 0.19, 0.32, 0.03, 0.01, 0.05, 0.06, 0.07, 0.04, 0.06, 0.28]]
a.N_SR = len(a.SR_names)
if allow_corr:
    analyses += [a]

a = Analysis("CMS_13TeV_2LEPsoft_36invfb_nocovar")
a.SR_names = ["SR1__i0", "SR2__i1", "SR3__i2", "SR4__i3", "SR5__i4", "SR6__i5", "SR7__i6", "SR8__i7", "SR9__i8", "SR10__i9", "SR11__i10", "SR12__i11", ]
a.SR_n     = [2, 15, 19, 18, 1, 0, 3, 1, 2, 1, 2, 0, ]
a.SR_b     = [3.5, 12, 17, 11, 1.6, 3.5, 2, 0.51, 1.4, 1.5, 1.5, 1.2, ]
a.SR_b_sys = [1, 2.3, 2.4, 2, 0.7, 0.9, 0.7, 0.52, 0.7, 0.6, 0.8, 0.6, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("CMS_13TeV_2OSLEP_36invfb")
a.SR_names = ["SR1__i0", "SR2__i1", "SR3__i2", "SR4__i3", "SR5__i4", "SR6__i5", "SR7__i6", ]
a.SR_n     = [57, 29, 2, 0, 9, 5, 1, ]
a.SR_b     = [54.9, 21.6, 6, 2.5, 7.6, 5.6, 1.3, ]
a.SR_b_sys = [7, 5.6, 1.9, 0.9, 2.8, 1.6, 0.4, ]
a.cov = [[52.8, 12.7,    3,  1.2,  4.5,  5.1,  1.2],
 [12.7, 41.4,  3.6,    2,  2.5,    2,  0.7],
 [   3,  3.6,  1.6,  0.6,  0.4,  0.3,  0.1],
 [ 1.2,    2,  0.6,  1.1,  0.3,  0.1,  0.1],
 [ 4.5,  2.5,  0.4,  0.3,  6.5,  1.8,  0.4],
 [ 5.1,    2,  0.3,  0.1,  1.8,  2.4,  0.4],
 [ 1.2,  0.7,  0.1,  0.1,  0.4,  0.4,  0.2]]
a.N_SR = len(a.SR_names)
if allow_corr:
    analyses += [a]

a = Analysis("CMS_13TeV_2OSLEP_36invfb_nocovar")
a.SR_names = ["SR1__i0", "SR2__i1", "SR3__i2", "SR4__i3", "SR5__i4", "SR6__i5", "SR7__i6", ]
a.SR_n     = [57, 29, 2, 0, 9, 5, 1, ]
a.SR_b     = [54.9, 21.6, 6, 2.5, 7.6, 5.6, 1.3, ]
a.SR_b_sys = [7, 5.6, 1.9, 0.9, 2.8, 1.6, 0.4, ]
a.N_SR = len(a.SR_names)
#analyses += [a]

a = Analysis("CMS_13TeV_MultiLEP_2SSLep_36invfb")
a.SR_names = ["SR1__i0", "SR2__i1", ]
a.SR_n     = [13, 18, ]
a.SR_b     = [12, 18, ]
a.SR_b_sys = [3, 4, ]
a.N_SR = len(a.SR_names)
analyses += [a]

a = Analysis("CMS_13TeV_MultiLEP_3Lep_36invfb")
a.SR_names = ["SR3__i0", "SR4__i1", "SR5__i2", "SR6__i3", "SR7__i4", "SR8__i5", ]
a.SR_n     = [19, 128, 18, 2, 82, 166, ]
a.SR_b     = [19, 142, 22, 1, 109, 197, ]
a.SR_b_sys = [4, 34, 5, 0.6, 28, 42, ]
a.N_SR = len(a.SR_names)
analyses += [a]

# Extra hacking! Set background uncertainties in each analysis to be highly correlated!
#for a in analyses:
#    if a.cov is None:
#        # No covariance matrix: so add one! 0.9 correlation coefficient (1 gives singular matrix)
#        a.cov = [[0.9*a.SR_b_sys[i]*a.SR_b_sys[j] for i in range(a.N_SR)] for j in range(a.N_SR)]
#        for k in range(a.N_SR):
#            a.cov[k][k] = a.SR_b_sys[k]**2 # Need diagonal to be correct variance though
#
#    print(a.cov)
