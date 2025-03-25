import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-}#rightarrow ZH; H#rightarrow ZZ; Z#rightarrow ll; l=e,#mu'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = './outputs/plots/ZZZ6l/' 
inputDir       = './outputs/histmaker/ZZZ6l/' 

plotStatUnc    = True

colors = {}
colors['ZH-6l'] = ROOT.kRed
colors['llH-qqll'] = ROOT.kBlue+1
colors['qqH-4l'] = ROOT.kGreen+2
#colors['ZZll-6l'] = ROOT.kOrange+1
colors['llZ(Z)-4l'] = ROOT.kMagenta+1
colors['qqZ(Z)-4l'] = ROOT.kCyan+1

procs = {}
procs['signal'] = {'ZH-6l':['wzp6_ee_llH_HZZ_llll_ecm240']} # 60 events
procs['backgrounds'] =  {'llH-qqll':['wzp6_ee_llH_HZZ_qqll_ecm240'], # 5k?
                        'qqH-4l':['wzp6_ee_qqH_HZZ_llll_ecm240'], # 380 events 
                        #'ZZll-6l':['p8_ee_llZZ_Zll_ecm240'], # 1k events
                        'llZ(Z)-4l': ['p8_ee_llZZ_Zll_single_os_Z_ecm240'], # 9 events
                        'qqZ(Z)-4l': ['p8_ee_jjZZ_Zll_ecm240', 'p8_ee_jjllZ_Zll_ecm240'], #
                        }


legend = {}
legend['ZH-6l'] = 'ZH-6l'
legend['llH-qqll'] = 'llH-qqll'
legend['qqH-4l'] = 'qqH-4l'
#legend['ZZll-6l'] = 'ZZll-6l'
legend['llZ(Z)-4l'] = 'llZ(Z)-4l'
legend['qqZ(Z)-4l'] = 'qqZ(Z)-4l'



hists = {}

hists["higgs_mass"] = {
    "output":   "higgs_mass",
    "input":    "higgs_mass",
    "logy":     False,
    "stack":    True,
    # "rebin":    100,
    "xmin":     115,
    "xmax":     130,
    "ymin":     0,
    "ymax":     2500,
    "xtitle":   "Mass (GeV)",
    "ytitle":   "Events / 100 MeV",
}

hists["recoil_mass_cut4"] = {
    "output":   "recoil_mass_cut4",
    "input":    "recoil_mass_cut4",
    "logy":     False,
    "stack":    True,
    # "rebin":    100,
    # "xmin":     110,
    # "xmax":     160,
    # "ymin":     0,
    # "ymax":     2500,
    "xtitle":   "recoil mass m_{ll} (GeV)",
    "ytitle":   "Events / 0.5 GeV",
}


hists["cutFlow"] = {
    "output":   "cutFlow",
    "input":    "cutFlow",
    "logy":     True,
    "stack":    False,
    "xmin":     0,
    "xmax":     6,
    "ymin":     1e4,
    "ymax":     1e11,
    "xtitle":   ["All events", "N(e) + N(#mu) >=6", "N(e_{iso}) + N(#mu_{iso}) = 6", "l pairs", "115 < m_{H} < 127", "122 < m_{recoil} < 160"], #"40 < p_{Z, prod} < 60"],
    "ytitle":   "Events",
    # "scaleSig": 10
}