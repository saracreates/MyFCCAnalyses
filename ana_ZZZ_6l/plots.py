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
colors['ZZZ-6l'] = ROOT.kRed
colors['llH-qqll'] = ROOT.kBlue+1
colors['qqH-4l'] = ROOT.kGreen+2

procs = {}
procs['signal'] = {'ZZZ-6l':['wzp6_ee_llH_HZZ_llll_ecm240']}
procs['backgrounds'] =  {'llH-qqll':['wzp6_ee_llH_HZZ_qqll_ecm240'], 
                        'qqH-4l':['wzp6_ee_llH_HZZ_llll_ecm240'],
                        }


legend = {}
legend['ZZZ-6l'] = 'ZZZ-6l'
legend['llH-qqll'] = 'llH-qqll'
legend['qqH-4l'] = 'qqH-4l'



hists = {}

hists["higgs_mass"] = {
    "output":   "higgs_mass",
    "input":    "higgs_mass",
    "logy":     False,
    "stack":    True,
    "rebin":    100,
    "xmin":     120,
    "xmax":     140,
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
    "rebin":    100,
    "xmin":     120,
    "xmax":     140,
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
    "stack":    True,
    "xmin":     0,
    "xmax":     6,
    "ymin":     1e4,
    "ymax":     1e11,
    "xtitle":   ["All events", "N(e) + N(#mu) >=6", "N(e_{iso}) + N(#mu_{iso}) = 6", "l pairs", "124 < m_H < 125.5"],
    "ytitle":   "Events ",
    "scaleSig": 10
}