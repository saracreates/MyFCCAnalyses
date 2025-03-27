import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-}#rightarrow Z(qq)H; H#rightarrow Z(ll)Z*(#nu#nu)' #; l=e,#mu'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = './outputs/mva/ZZZ_qqllvv/plots' 
inputDir       = './outputs/mva/ZZZ_qqllvv/final_selection'

plotStatUnc    = True

colors = {}
colors['qqH-HZZ'] = ROOT.kRed
colors['qqH-HWW'] = ROOT.kBlue+1
colors['ee-ZZ'] = ROOT.kGreen+2
colors['ee-WW'] = ROOT.kOrange+1
# colors['qqH-Hbb'] = ROOT.kBlack
colors['qqH-Htautau'] = ROOT.kCyan+1
colors['Zqq'] = ROOT.kMagenta+1
colors['llH-HZZ'] = ROOT.kBlack
colors['eeH-HZZ'] = ROOT.kViolet+1
colors['mumuH-HZZ'] = ROOT.kViolet+2
colors['nunuH-HZZ'] = ROOT.kViolet+3

i = 0


procs = {}
#procs['signal'] = {'qqH-HZZ':['wzp6_ee_qqH_HZZ_ecm240']} 
procs['signal'] = {'qqH-HZZ':[f'wzp6_ee_qqH_HZZ_llvv_ecm240_sel{i}_histo']}
procs['backgrounds'] =  {'qqH-HWW':[f'wzp6_ee_qqH_HWW_ecm240_sel{i}_histo', f'wzp6_ee_ssH_HWW_ecm240_sel{i}_histo', f'wzp6_ee_ccH_HWW_ecm240_sel{i}_histo', f'wzp6_ee_bbH_HWW_ecm240_sel{i}_histo' ], 
                        'ee-ZZ':[f'p8_ee_ZZ_ecm240_sel{i}_histo'], 
                        'ee-WW':[f'p8_ee_WW_ecm240_sel{i}_histo'], 
                        #'qqH-Hbb':[f'wzp6_ee_qqH_Hbb_ecm240_sel{i}_histo', f'wzp6_ee_ssH_Hbb_ecm240_sel{i}_histo', f'wzp6_ee_ccH_Hbb_ecm240_sel{i}_histo', f'wzp6_ee_bbH_Hbb_ecm240_sel{i}_histo'], 
                        'qqH-Htautau':[f'wzp6_ee_qqH_Htautau_ecm240_sel{i}_histo', f'wzp6_ee_ssH_Htautau_ecm240_sel{i}_histo', f'wzp6_ee_ccH_Htautau_ecm240_sel{i}_histo', f'wzp6_ee_bbH_Htautau_ecm240_sel{i}_histo'],
                        'Zqq':[f'p8_ee_Zqq_ecm240_sel{i}_histo'], 
                        'llH-HZZ': [f'wzp6_ee_eeH_HZZ_ecm240_sel{i}_histo', f'wzp6_ee_mumuH_HZZ_ecm240_sel{i}_histo'],
                        'eeH-HZZ':[f'wzp6_ee_eeH_HZZ_ecm240_sel{i}_histo'],
                        'mumuH-HZZ':[f'wzp6_ee_mumuH_HZZ_ecm240_sel{i}_histo'],
                        'nunuH-HZZ':[f'wzp6_ee_nunuH_HZZ_ecm240_sel{i}_histo'], #none left
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
legend['eeH-HZZ'] = 'eeH-HZZ'
legend['mumuH-HZZ'] = 'mumu-HZZ'
legend['nunuH-HZZ'] = 'nunuH-HZZ'




hists = {}

hists["mva_score"] = {
    "output":   "mva_score",
    "input":    "mva_score",
    "logy":     True,
    "stack":    False,
    # "rebin":    100,
    "xmin":     0,
    "xmax":     1,
    # "ymin":     0,
    # "ymax":     2500,
    "xtitle":   "mva score",
    "ytitle":   "Events",
    # "scaleSig": 100
}


# hists["cutFlow"] = {
#     "output":   "cutFlow",
#     "input":    "cutFlow",
#     "logy":     True,
#     "stack":    False,
#     "xmin":     0,
#     "xmax":     9,
#     "ymin":     1e4,
#     "ymax":     1e11,
#     "xtitle":   ["All events", "1 lepton pair", "100 GeV < recoil_{jj} < 170 GeV", "80 GeV < m_{ll} < 100 GeV", "85 < m_{jj} < 105 GeV", "40 < p_{jj} < 55 GeV", "10 < recoil_{jjll} < 50 GeV", "5 < p_{T, miss} < 50 GeV", "cos(#theta)_{p_{miss}, p_{had}} < 0.3", "120 GeV < recoil_{jj} < 140 GeV"],
#     "ytitle":   "Events ",
#     # "scaleSig": 10
# }