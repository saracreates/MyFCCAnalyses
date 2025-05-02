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

outdir         = './outputs/plots/flavor/extended' 
inputDir       = './outputs/histmaker/flavor/extended' 

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
procs['backgrounds'] =  {'Aqq':['p8_ee_qqgamma_ecm240'], 'Acc':['p8_ee_ccgamma_ecm240'], 'Abb':['p8_ee_bbgamma_ecm240'], 'Atautau':['p8_ee_tautaugamma_ecm240'], 'Amumu':['p8_ee_mumugamma_ecm240'], 'Aee':['p8_ee_eegamma_ecm240'], 'WW':['p8_ee_WW_ecm240'], 'ZZ':['p8_ee_ZZ_ecm240'], 'ZH':['p8_ee_ZH_ecm240']}
#procs['backgrounds'] =  {'Abb':['p8_ee_bbgamma_ecm240']}


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
    "xmax":     5,
    "ymin":     1e4,
    "ymax":     1e11,
    #"xtitle":   ["All events", "iso < 0.2", "60  < p_{#gamma} < 100 ", "|cos(#theta)_{#gamma}|<0.9", "n particles > 5"],
    "xtitle":   ["cut 1 to 5", "validate", "B score sum > 1 ", "recoil mass tight", ""],
    "ytitle":   "Events ",
}


hists["recojet_isB_size"] = {
    "input":   "recojet_isB_size",
    "output":   "recojet_isB_size",
    "logy":     False,
    "stack":    True,
    "xmin":     0,
    "xmax":     5,
    "xtitle":   "recojet_isB_size",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["recojet_isB0"] = {
    "input":   "recojet_isB0",
    "output":   "recojet_isB0",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     1,
    "xtitle":   "recojet_isB0",
    "ytitle":   "Events ",
    "density": True,
    #"scaleSig": 1000,

}

hists["recojet_isB1"] = {
    "input":   "recojet_isB1",
    "output":   "recojet_isB1",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     1,
    "xtitle":   "recojet_isB1",
    "ytitle":   "Events ",
    "density": True,
    #"scaleSig": 1000,

}

hists["scoresum_B"] = {
    "input":   "scoresum_B",
    "output":   "scoresum_B",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_B",
    "ytitle":   "Normalized Events",
    "density": True,
    #"scaleSig": 1000,

}

hists["scorediv_B"] = {
    "input":   "scorediv_B",
    "output":   "scorediv_B",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     3,
    "xtitle":   "scorediv_B",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["scoresum_C"] = {
    "input":   "scoresum_C",
    "output":   "scoresum_C",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_C",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["scoresum_U"] = {
    "input":   "scoresum_U",
    "output":   "scoresum_U",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_U",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["scoresum_D"] = {
    "input":   "scoresum_D",
    "output":   "scoresum_D",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_D",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["scoresum_S"] = {
    "input":   "scoresum_S",
    "output":   "scoresum_S",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_S",
    "ytitle":   "Events ",
    "density": False,
    #"scaleSig": 1000,

}

hists["scoresum_G"] = {
    "input":   "scoresum_G",
    "output":   "scoresum_G",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_G",
    "ytitle":   "Events ",
    "density": True,
    #"scaleSig": 1000,

}

hists["scoresum_Tau"] = {
    "input":   "scoresum_Tau",
    "output":   "scoresum_Tau",
    "logy":     True,
    "stack":    True,
    "xmin":     0,
    "xmax":     2,
    "xtitle":   "scoresum_Tau",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,

}




hists["gamma_recoil_m"] = {
    "input":   "gamma_recoil_m",
    "output":   "gamma_recoil_m",
    "logy":     False,
    "stack":    True,
    "xmin":     110,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,

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
    #"scaleSig": 1000,

}

hists["gamma_recoil_m_cut_3"] = {
    "input":   "gamma_recoil_m_cut_3",
    "output":   "gamma_recoil_m_cut_3",
    "logy":     False,
    "stack":    True,
    "xmin":     110,
    "xmax":     150,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
    "density": False,
    "scaleSig": 1000,

}



