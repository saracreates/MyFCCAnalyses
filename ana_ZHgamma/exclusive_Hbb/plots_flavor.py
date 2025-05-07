import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow #gamma H'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

#outdir         = './outputs/plots/flavor/' 
#inputDir       = './outputs/histmaker/flavor/' 

outdir         = '/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/plots/ZHgamma_btag' 
inputDir       = '/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma_btag/' 

plotStatUnc    = True

colors = {}
colors['AH'] = ROOT.kRed
colors['Acc'] = ROOT.kBlue+1
colors['Aqq'] = ROOT.kGreen+2
colors['Abb'] = ROOT.kYellow+3
colors['WW'] = ROOT.kCyan
colors['ZZ'] = ROOT.kAzure-9
colors['Aee'] = ROOT.kViolet+3
colors['Atautau'] = ROOT.kOrange
colors['Amumu'] = ROOT.kMagenta
colors['ZH'] = ROOT.kGray+2

#procs = {}
#procs['signal'] = {'ZH':['wzp6_ee_mumuH_ecm240']}
#procs['backgrounds'] =  {'WW':['p8_ee_WW_ecm240'], 'ZZ':['p8_ee_ZZ_ecm240']}
procs = {}
procs['signal'] = {'AH':['p8_ee_Hgamma_ecm240']}
procs['backgrounds'] =  {'Aqq':['p8_ee_qqgamma_ecm240'], 
                        'Acc':['p8_ee_ccgamma_ecm240'], 
                        'Abb':['p8_ee_bbgamma_ecm240'], 
                        'Atautau':['p8_ee_tautaugamma_ecm240'], 
                        'Amumu':['p8_ee_mumugamma_ecm240'], 
                        'Aee':['p8_ee_eegamma_ecm240'], 
                        # 'WW':['p8_ee_WW_ecm240'], 
                        # 'ZZ':['p8_ee_ZZ_ecm240'], 
                        'ZH':['p8_ee_ZH_ecm240']}


legend = {}
legend['AH'] = '#gamma H'
legend['Aqq'] = '#gamma q#bar{q}'
legend['Acc'] = '#gamma c#bar{c}'
legend['Abb'] = '#gamma b#bar{b}'
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'
legend['Aee'] = '#gamma e^{+} e^{-}'
legend['Atautau'] = '#gamma #tau^{+} #tau^{-}'
legend['Amumu'] = '#gamma #mu^{+} #mu^{-}'
legend['ZH'] = 'ZH'


hists = {}
hists2D = {}





hists["cutFlow"] = {
    "input":   "cutFlow",
    "output":   "cutFlow",
    "logy":     True,
    "stack":   True,
    "xmin":     0,
    "xmax":     8,
    "ymin":     1e4,
    "ymax":     1e11,
    #"xtitle":   ["All events", "iso < 0.2", "60  < p_{#gamma} < 100 ", "|cos(#theta)_{#gamma}|<0.9", "n particles > 5"],
    "xtitle":   ["All events", "iso < 0.2", "60  < p_{#gamma} < 100 ", "|cos(#theta)_{#gamma}|<0.9", "n particles > 9", "110 < m_{recoil} < 140 ", "b score sum > 1",  "120 < m_{recoil} < 132 "],
    "ytitle":   "Events ",
}

# hists["gamma_recoil_m"] = {
#     "input":   "gamma_recoil_m",
#     "output":   "gamma_recoil_m",
#     "logy":     False,
#     "stack":    True,
#     "xmin":     110,
#     "xmax":     150,
#     "xtitle":   "Recoil (GeV)",
#     "ytitle":   "Events ",
#     "density": False,
#     "scaleSig": 1000,

# }


hists["b_tags_sum"] = {
    "input":   "b_tags_sum",
    "output":   "b_tags_sum",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "b_tags_sume",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["gamma_recoil_m_LHF"] = {
    "input":   "gamma_recoil_m_LHF",
    "output":   "gamma_recoil_m_LHF",
    "logy":     False,
    "stack":    True,
    "xmin":     110,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "scaleSig": 1000,
    "density": False
}

hists["gamma_recoil_m_cut_3"] = {
    "input":   "gamma_recoil_m_cut_3",
    "output":   "gamma_recoil_m_cut_3",
    "logy":     False,
    "stack":    True,
    "xmin":     100,
    "xmax":     170,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "scaleSig": 1000,
    "density": False
}



hists["gamma_recoil_m_cut_4"] = {
    "input":   "gamma_recoil_m_cut_4",
    "output":   "gamma_recoil_m_cut_4",
    "logy":     False,
    "stack":    True,
    "xmin":     80,
    "xmax":     200,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
   # "scaleSig": 1000,
    "density": True
}



hists["gamma_recoil_m_norm_cut_4"] = {
    "input":   "gamma_recoil_m_cut_4",
    "output":   "gamma_recoil_m_norm_cut_4",
    "logy":     False,
    "stack":    False,
    "xmin":     80,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
   # "scaleSig": 1000,
    "density": True

}

hists["gamma_recoil_m_signal_cut"] = {
    "input":   "gamma_recoil_m_signal_cut",
    "output":   "gamma_recoil_m_signal_cut",
    "logy":     False,
    "stack":    True,
    "xmin":     110,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,
}

hists["electrons_p_baseline"] = {
    "input":   "electrons_p_baseline",
    "output":   "electrons_p_baseline",
    "logy":     True,
    "stack":    False,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "electrons_p_baseline",
    "ytitle":   "Events ",
    "density": False,
    "density": True
}


hists["photons_p_cut_0"] = {
    "input":   "photons_p_cut_0",
    "output":   "photons_p_cut_0",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "photons_p_cut_0",
    "ytitle":   "Events ",
    "density": False,
    "density": True
}


hists["photons_p_cut_1"] = {
    "input":   "photons_p_cut_1",
    "output":   "photons_p_cut_1",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "p(\gamma) [GeV]",
    "ytitle":   "Normalized Events",
    "density": True,
}

hists["photons_p_cut_2"] = {
    "input":   "photons_p_cut_2",
    "output":   "photons_p_cut_2",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "photons_p_cut_2",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,
    "density": False
}

hists["photons_p_cut_3"] = {
    "input":   "photons_p_cut_3",
    "output":   "photons_p_cut_3",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "photons_p_cut_3",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,
    "density": False
}

hists["electrons_n_baseline"] = {
    "input":   "electrons_n_baseline",
    "output":   "electrons_n_baseline",
    "logy":     False,
    "stack":    False,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "electrons_n_baseline",
    "ytitle":   "Events ",
    "density": True
}

hists["photons_n_cut_0"] = {
    "input":   "photons_n_cut_0",
    "output":   "photons_n_cut_0",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_n_cut_0",
    "ytitle":   "Events ",
    "density": True
}

hists["photons_n_cut_1"] = {
    "input":   "photons_n_cut_1",
    "output":   "photons_n_cut_1",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_n_cut_1",
    "ytitle":   "Events ",
    "density": True
}

hists["photons_n_cut_2"] = {
    "input":   "photons_n_cut_2",
    "output":   "photons_n_cut_2",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_n_cut_2",
    "ytitle":   "Events ",
    "density": True
}

hists["photons_n_cut_3"] = {
    "input":   "photons_n_cut_3",
    "output":   "photons_n_cut_3",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_n_cut_3",
    "ytitle":   "Events ",
    "density": True
}

hists["recopart_no_gamma_n_cut_0"] = {
    "input":   "recopart_no_gamma_n_cut_0",
    "output":   "recopart_no_gamma_n_cut_0",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     60,
    "xtitle":   "# reco particles",
    "ytitle":   "Normalized Events",
    "density": True
}

hists["recopart_no_gamma_n_cut_4"] = {
    "input":   "recopart_no_gamma_n_cut_4",
    "output":   "recopart_no_gamma_n_cut_4",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     60,
    "xtitle":   "recopart_no_gamma_n_cut_4",
    "ytitle":   "Events ",
    "density": False
}

hists["photon_isolation"] = {
    "input":   "photon_isolation",
    "output":   "photon_isolation",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "iso(\gamma)_{\Delta R< 0.5}",
    "ytitle":   "Normalized Events ",
    "density": True
}


hists["photons_cos_theta_cut_0"] = {
    "input":   "photons_cos_theta_cut_0",
    "output":   "photons_cos_theta_cut_0",
    "logy":     False,
    "stack":    False,
    "xmin":     -1,
    "xmax":     1,
    "xtitle":   "photons_cos_theta_cut_0",
    "ytitle":   "Events ",
     "scaleSig": 10000,
     "density": True
}



hists["photons_cos_theta_cut_1"] = {
    "input":   "photons_cos_theta_cut_1",
    "output":   "photons_cos_theta_cut_1",
    "logy":     False,
    "stack":    False,
    "xmin":     -1,
    "xmax":     1,
    "xtitle":   "photons_cos_theta_cut_1",
    "ytitle":   "Events ",
     "scaleSig": 10000,
     "density": True
}

hists["photons_cos_theta_cut_2"] = {
    "input":   "photons_cos_theta_cut_2",
    "output":   "photons_cos_theta_cut_2",
    "logy":     False,
    "stack":    True,
    "xmin":     -1,
    "xmax":     1,
    "xtitle":   "cos(\Theta)_{\gamma}",
    "ytitle":   "Normalized Events",
   #  "scaleSig": 10000,
     "density": True
}


hists["photons_cos_theta_cut_3"] = {
    "input":   "photons_cos_theta_cut_3",
    "output":   "photons_cos_theta_cut_3",
    "logy":     False,
    "stack":    False,
    "xmin":     -1,
    "xmax":     1,
    "xtitle":   "photons_cos_theta_cut_3",
    "ytitle":   "Events ",
     "scaleSig": 10000,
     "density": True
}

# hists["recojet_isB0"] = {
#     "input":   "recojet_isB0",
#     "output":   "recojet_isB0",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     1,
#     "xtitle":   "recojet_isB0",
#     "ytitle":   "Events ",
#     "density": True,
#     #"scaleSig": 1000,

# }

# hists["recojet_isB1"] = {
#     "input":   "recojet_isB1",
#     "output":   "recojet_isB1",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     1,
#     "xtitle":   "recojet_isB1",
#     "ytitle":   "Events ",
#     "density": True,
#     #"scaleSig": 1000,

# }

# hists["scoresum_B"] = {
#     "input":   "scoresum_B",
#     "output":   "scoresum_B",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_B",
#     "ytitle":   "Normalized Events",
#     "density": True,
#     #"scaleSig": 1000,

# }

# hists["scorediv_B"] = {
#     "input":   "scorediv_B",
#     "output":   "scorediv_B",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     3,
#     "xtitle":   "scorediv_B",
#     "ytitle":   "Events ",
#     "density": False,
#     #"scaleSig": 1000,

# }

# hists["scoresum_C"] = {
#     "input":   "scoresum_C",
#     "output":   "scoresum_C",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_C",
#     "ytitle":   "Events ",
#     "density": False,
#     #"scaleSig": 1000,

# }

# hists["scoresum_U"] = {
#     "input":   "scoresum_U",
#     "output":   "scoresum_U",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_U",
#     "ytitle":   "Events ",
#     "density": False,
#     #"scaleSig": 1000,

# }

# hists["scoresum_D"] = {
#     "input":   "scoresum_D",
#     "output":   "scoresum_D",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_D",
#     "ytitle":   "Events ",
#     "density": False,
#     #"scaleSig": 1000,

# }

# hists["scoresum_S"] = {
#     "input":   "scoresum_S",
#     "output":   "scoresum_S",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_S",
#     "ytitle":   "Events ",
#     "density": False,
#     #"scaleSig": 1000,

# }

# hists["scoresum_G"] = {
#     "input":   "scoresum_G",
#     "output":   "scoresum_G",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_G",
#     "ytitle":   "Events ",
#     "density": True,
#     #"scaleSig": 1000,

# }

# hists["scoresum_Tau"] = {
#     "input":   "scoresum_Tau",
#     "output":   "scoresum_Tau",
#     "logy":     True,
#     "stack":    True,
#     "xmin":     0,
#     "xmax":     2,
#     "xtitle":   "scoresum_Tau",
#     "ytitle":   "Events ",
#     "density": False,
#     "scaleSig": 1000,

# }






