import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-}#rightarrow Z(qq)H; H#rightarrow Z(ll)Z*(#nu#nu)' #; l=e,#mu'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = './outputs/plots/ZZZqqllvv_LHF/' 
inputDir       = './outputs/histmaker_LHF/ZZZqqllvv'

plotStatUnc    = True

colors = {}
colors['qqH-HZZ'] = ROOT.kRed
colors['qqH-HWW'] = ROOT.kBlue+1
colors['ee-ZZ'] = ROOT.kGreen+2
colors['ee-WW'] = ROOT.kOrange+1
#colors['qqH-Hbb'] = ROOT.kBlack
colors['qqH-Htautau'] = ROOT.kCyan+1
colors['Zqq'] = ROOT.kMagenta+1
colors['llH-HZZ'] = ROOT.kBlack
# colors['eeH-HZZ'] = ROOT.kViolet+1
# colors['mumuH-HZZ'] = ROOT.kViolet+2
# colors['nunuH-HZZ'] = ROOT.kViolet+3


procs = {}
#procs['signal'] = {'qqH-HZZ':['wzp6_ee_qqH_HZZ_ecm240']} 
procs['signal'] = {'qqH-HZZ':['wzp6_ee_qqH_HZZ_llvv_ecm240']}
procs['backgrounds'] =  {'qqH-HWW':['wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240', 'wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240' ], 
                        'ee-ZZ':['p8_ee_ZZ_ecm240'], 
                        'ee-WW':['p8_ee_WW_ecm240'], 
                        # 'qqH-Hbb':['wzp6_ee_qqH_Hbb_ecm240', 'wzp6_ee_ssH_Hbb_ecm240', 'wzp6_ee_ccH_Hbb_ecm240', 'wzp6_ee_bbH_Hbb_ecm240'], 
                        'qqH-Htautau':['wzp6_ee_qqH_Htautau_ecm240', 'wzp6_ee_ssH_Htautau_ecm240', 'wzp6_ee_ccH_Htautau_ecm240', 'wzp6_ee_bbH_Htautau_ecm240'],
                        'Zqq':['p8_ee_Zqq_ecm240'], 
                        'llH-HZZ': ['wzp6_ee_eeH_HZZ_ecm240', 'wzp6_ee_mumuH_HZZ_ecm240'],
                        # 'eeH-HZZ':['wzp6_ee_eeH_HZZ_ecm240'],
                        # 'mumuH-HZZ':['wzp6_ee_mumuH_HZZ_ecm240'],
                        # 'nunuH-HZZ':['wzp6_ee_nunuH_HZZ_ecm240'], # none left
}


legend = {}
legend['qqH-HZZ'] = 'qqH-HZZ'
legend['qqH-HWW'] = 'qqH-HWW'
legend['ee-ZZ'] = 'ee-ZZ'
legend['ee-WW'] = 'ee-WW'
# legend['qqH-Hbb'] = 'qqH-Hbb'
legend['qqH-Htautau'] = 'qqH-Htautau'
legend['Zqq'] = 'Zqq'
legend['llH-HZZ'] = 'llH-HZZ'
# legend['eeH-HZZ'] = 'eeH-HZZ'
# legend['mumuH-HZZ'] = 'mumu-HZZ'
# legend['nunuH-HZZ'] = 'nunuH-HZZ'




hists = {}

hists["recoil_mass_LHF"] = {
    "output":   "recoil_mass_LHF",
    "input":    "recoil_mass_LHF",
    "logy":     False,
    "stack":    True,
    # "rebin":    100,
    "xmin":     100,
    "xmax":     170,
    # "ymin":     0,
    # "ymax":     2500,
    "xtitle":   "recoil mass m_{qq} (GeV)",
    "ytitle":   "Events",
}

hists["mva_score_signal"] = {
    "output":   "mva_score_signal",
    "input":    "mva_score_signal",
    "logy":     True,
    "stack":    True,
    "xtitles":  ["mva score for signal"],
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
    "xtitle":   ["All events", "1 lepton pair", "100 GeV < recoil_{jj} < 170 GeV", "80 GeV < m_{ll} < 100 GeV", "85 < m_{jj} < 105 GeV", "40 < p_{jj} < 55 GeV", "10 < recoil_{jjll} < 50 GeV", "5 < p_{T, miss} < 50 GeV", "cos(#theta)_{p_{miss}, p_{had}} < 0.3", "120 GeV < recoil_{jj} < 140 GeV", "score(ZZ) < 0.5"],
    "ytitle":   "Events ",
    # "scaleSig": 10
}