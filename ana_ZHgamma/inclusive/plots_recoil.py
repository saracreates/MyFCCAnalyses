# from Lena, 28. April 2025: https://github.com/herrmannlena/FCCAnalyses/blob/higgsgamma/myanalysis/plots_recoil.py 

import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow #gamma H'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = '/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/plots/ZHgamma/' 
inputDir       = '/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma' 

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
# procs['signal'] = {'AH':['reco_higgsgamma_test_REC.edm4hep']} # for testing purpose

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
    "xmax":     7,
    "ymin":     1e4,
    "ymax":     1e11,
    #"xtitle":   ["All events", "iso < 0.2", "60  < p_{#gamma} < 100 ", "|cos(#theta)_{#gamma}|<0.9", "n particles > 5"],
    "xtitle":   ["All events", "iso < 0.2", "60  < p_{#gamma} < 100 ", "|cos(#theta)_{#gamma}|<0.9", "n particles > 9", "110 < m_{recoil} < 140 ", "114 < m_{recoil} < 128 "],
    "ytitle":   "Events ",
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
"""
hists["gamma_recoil_m_signal_cut"] = {
    "input":   "gamma_recoil_m_signal_cut",
    "output":   "gamma_recoil_m_signal_cut",
    "logy":     False,
    "stack":    True,
    #"xmin":     110,
    "xmin":     116,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "density": False

}
"""
hists["gamma_recoil_m_tight_cut"] = {
    "input":   "gamma_recoil_m_tight_cut",
    "output":   "gamma_recoil_m_tight_cut",
    "logy":     False,
    "stack":    True,
    "xmin":     115,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "scaleSig": 1000,
    "density": False

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
"""
hists["photons_p_cut_4"] = {
    "input":   "photons_p_cut_4",
    "output":   "photons_p_cut_4",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     130,
    "xtitle":   "photons_p_cut_4",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,
    "density": True
}
"""
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
"""
hists["photons_n_cut_4"] = {
    "input":   "photons_n_cut_4",
    "output":   "photons_n_cut_4",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_n_cut_4",
    "ytitle":   "Events ",
    "density": True
}
"""
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
"""
hists["photons_cos_theta_cut_4"] = {
    "input":   "photons_cos_theta_cut_4",
    "output":   "photons_cos_theta_cut_4",
    "logy":     False,
    "stack":    False,
    "xmin":     -1,
    "xmax":     1,
    "xtitle":   "photons_cos_theta_cut_4",
    "ytitle":   "Events ",
     "scaleSig": 10000,
     "density": True
}
"""
"""
hists["photons_boosted_p"] = {
    "input":   "photons_boosted_p",
    "output":   "photons_boosted_p",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     100,
    "xtitle":   "photons_boosted_p",
    "ytitle":   "Events ",
    "density": True
}



hists["photons_boosted_n"] = {
    "input":   "photons_boosted_n",
    "output":   "photons_boosted_n",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     10,
    "xtitle":   "photons_boosted_n",
    "ytitle":   "Events ",
}
"""


"""



hists["zmumu_recoil_m"] = {
    "output":   "zmumu_recoil_m",
    "logy":     False,
    "stack":    True,
    "rebin":    100,
    "xmin":     120,
    "xmax":     140,
    "ymin":     0,
    "ymax":     2500,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events / 100 MeV",
}

hists["zmumu_p"] = {
    "output":   "zmumu_p",
    "logy":     False,
    "stack":    True,
    "rebin":    2,
    "xmin":     0,
    "xmax":     80,
    "ymin":     0,
    "ymax":     2000,
    "xtitle":   "p(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events ",
}

hists["zmumu_m"] = {
    "output":   "zmumu_m",
    "logy":     False,
    "stack":    True,
    "rebin":    2,
    "xmin":     86,
    "xmax":     96,
    "ymin":     0,
    "ymax":     3000,
    "xtitle":   "m(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events ",
}"""