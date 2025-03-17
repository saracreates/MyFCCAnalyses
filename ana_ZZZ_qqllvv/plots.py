import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.6 ab^{-1}"
ana_tex        = 'e^{+}e^{-}#rightarrow ZH; H#rightarrow ZZ; Z#rightarrow ll; l=e,#mu'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = './outputs/plots/ZZZqqllvv/' 
inputDir       = './outputs/histmaker/ZZZqqllvv/'

plotStatUnc    = True

colors = {}
colors['qqH-HZZ'] = ROOT.kRed
colors['qqH-HWW'] = ROOT.kBlue+1
colors['ee-ZZ'] = ROOT.kGreen+2
colors['ee-WW'] = ROOT.kOrange+1
colors['qqH-Hbb'] = ROOT.kMagenta+1
colors['qqH-Htautau'] = ROOT.kCyan+1


procs = {}
procs['signal'] = {'qqH-HZZ':['wzp6_ee_qqH_HZZ_ecm240']} 
procs['backgrounds'] =  {'qqH-HWW':['wzp6_ee_qqH_HWW_ecm240'], 
                        'ee-ZZ':['p8_ee_ZZ_ecm240'], 
                        'ee-WW':['p8_ee_WW_ecm240'], 
                        'qqH-Hbb':['wzp6_ee_qqH_Hbb_ecm240'], 
                        'qqH-Htautau':['wzp6_ee_qqH_Htautau_ecm240'] 
}


legend = {}
legend['qqH-HZZ'] = 'qqH-HZZ'
legend['qqH-HWW'] = 'qqH-HWW'
legend['ee-ZZ'] = 'ee-ZZ'
legend['ee-WW'] = 'ee-WW'
legend['qqH-Hbb'] = 'qqH-Hbb'
legend['qqH-Htautau'] = 'qqH-Htautau'




hists = {}

hists["m_ll"] = {
    "output":   "m_ll",
    "input":    "m_ll",
    "logy":     False,
    "stack":    True,
    "logy":     True,
    # "rebin":    100,
    "xmin":     0,
    "xmax":     200,
    # "ymin":     0,
    # "ymax":     2500,
    "xtitle":   "Mass (GeV)",
    "ytitle":   "Events",
}


hists["cutFlow"] = {
    "output":   "cutFlow",
    "input":    "cutFlow",
    "logy":     True,
    "stack":    False,
    "xmin":     0,
    "xmax":     9,
    "ymin":     1e4,
    "ymax":     1e11,
    "xtitle":   ["All events", "1 lepton pair", "122 GeV < recoil_{jj} < 130 GeV", "80 GeV < m_{ll} < 100 GeV", "6 GeV < p_{miss} < 60 GeV", "85 < m_{jj} < 105 GeV", "40 < p_{jj} < 55 GeV", "10 < recoil_{jjll} < 50 GeV", "cut on p_{miss} * p_{jets}"],
    "ytitle":   "Events ",
    # "scaleSig": 10
}