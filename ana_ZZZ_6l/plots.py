import ROOT

# global parameters
intLumi        = 1.
intLumiLabel   = "L = 10.6 ab^{-1}"
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
colors['Z(Z)ll-6l'] = ROOT.kMagenta+1
colors['qqZZ-4l'] = ROOT.kCyan+1
colors['qqllZ-2l'] = ROOT.kBlack

procs = {}
procs['signal'] = {'ZH-6l':['wzp6_ee_llH_HZZ_llll_ecm240']} # 60 events
procs['backgrounds'] =  {'llH-qqll':['wzp6_ee_llH_HZZ_qqll_ecm240'], # 5k?
                        'qqH-4l':['wzp6_ee_qqH_HZZ_llll_ecm240'], # 380 events 
                        #'ZZll-6l':['p8_ee_llZZ_Zll_ecm240'], # 1k events
                        'Z(Z)ll-6l': ['p8_ee_llZZ_Zll_single_os_Z_ecm240'], # 9 events
                        'qqZZ-4l': ['p8_ee_jjZZ_Zll_ecm240'], #
                        'qqllZ-2l': ['p8_ee_jjllZ_Zll_ecm240'], #
                        }


legend = {}
legend['ZH-6l'] = 'ZH-6l'
legend['llH-qqll'] = 'llH-qqll'
legend['qqH-4l'] = 'qqH-4l'
#legend['ZZll-6l'] = 'ZZll-6l'
legend['Z(Z)ll-6l'] = 'Z(Z)ll-6l'
legend['qqZZ-4l'] = 'jjZZ-4l'
legend['qqllZ-2l'] = 'jjllZ-2l'



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

hists["recoil_mass"] = {
    "input":    "recoil_mass",
    "output":   "recoil_mass",
    "logy":     False,
    "stack":    True,
    # "rebin":    100,
    "xmin":     120,
    "xmax":     135,
    "ymin":     0,
    "ymax":     2500,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events ",
}

#bins_recoil = (2000, 110, 180) # 1 MeV bins 
#bins_higgs = (2000, 100, 135) # 100 MeV bins



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
    "ytitle":   "Events ",
    # "scaleSig": 10
}