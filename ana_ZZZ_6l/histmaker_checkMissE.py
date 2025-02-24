# list of processes (mandatory)
processList = {
    'wzp6_ee_llH_HZZ_llll_ecm240':    {'fraction':1},
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["ZHfunctions.h"]

# Define the input dir (optional)
inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker_missE/ZZZ6l/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10600000 # 10.6 /ab


# define some binning for various histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (2000, 0, 200) # 100 MeV bins
bins_recoil = (200000, 0, 200) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, -5, 5)
bins_eta = (600, -3, 3)
bins_phi = (500, -5, 5)

bins_count = (10, 0, 10)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 3)

bins_missE = (2000, 0, 200) # 100 MeV bins
bins_missP = (2400, 0,80)



# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):
    """
    Goal of this exercise is to check if applying a cut on missing momentum/energy cuts away all the tau events. We expect this, because taus generate neutrinos.
    """

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    # define some aliases to be used later on
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")


    # get all the leptons from the collection
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")

    
    # select leptons with momentum > 5 GeV

    # muons
    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    # electrons
    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")
    
    # compute the muon isolation and store muons with an isolation cut of 0.25 in a separate column muons_sel_iso
    df = df.Define("muons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(muons, ReconstructedParticles)")
    df = df.Define("muons_sel_iso", "FCCAnalyses::ZHfunctions::sel_iso(0.25)(muons, muons_iso)")
    df = df.Define("muons_sel_iso_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)")
    # same for electrons
    df = df.Define("electrons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(electrons, ReconstructedParticles)")
    df = df.Define("electrons_sel_iso", "FCCAnalyses::ZHfunctions::sel_iso(0.25)(electrons, electrons_iso)")
    df = df.Define("electrons_sel_iso_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_sel_iso)")
    
        
    # baseline histograms, before any selection cuts (store with _cut0)
    results.append(df.Histo1D(("muons_p_cut0", "", *bins_p_mu), "muons_p"))
    results.append(df.Histo1D(("muons_theta_cut0", "", *bins_theta), "muons_theta"))
    results.append(df.Histo1D(("muons_phi_cut0", "", *bins_phi), "muons_phi"))
    results.append(df.Histo1D(("muons_q_cut0", "", *bins_charge), "muons_q"))
    results.append(df.Histo1D(("muons_no_cut0", "", *bins_count), "muons_no"))
    results.append(df.Histo1D(("muons_iso_cut0", "", *bins_iso), "muons_iso"))

    # same for electrons with same bins too
    results.append(df.Histo1D(("electrons_p_cut0", "", *bins_p_mu), "electrons_p"))
    results.append(df.Histo1D(("electrons_theta_cut0", "", *bins_theta), "electrons_theta"))
    results.append(df.Histo1D(("electrons_phi_cut0", "", *bins_phi), "electrons_phi"))
    results.append(df.Histo1D(("electrons_q_cut0", "", *bins_charge), "electrons_q"))
    results.append(df.Histo1D(("electrons_no_cut0", "", *bins_count), "electrons_no"))
    results.append(df.Histo1D(("electrons_iso_cut0", "", *bins_iso), "electrons_iso"))
    

    #########
    ### CUT 0: all events (lepton momentum > 5 GeV)
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))

    #########
    ### CUT 1: cut on missing energy
    #########
    df = df.Define("part_missingEnergy", "FCCAnalyses::ZHfunctions::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "FCCAnalyses::ReconstructedParticle::get_e(part_missingEnergy)")
    df = df.Define("missingMomentum", "FCCAnalyses::ReconstructedParticle::get_p(part_missingEnergy)")

    results.append(df.Histo1D(("missingEnergy", "", *bins_missE), "missingEnergy"))
    results.append(df.Histo1D(("missingMomentum", "", *bins_missP), "missingMomentum"))

    df = df.Filter("missingMomentum[0] < 5")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    # cut again on lower missing energy
    df = df.Filter("missingMomentum[0] < 2")
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    # again
    df = df.Filter("missingMomentum[0] < 1")
    df = df.Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))



    

    return results, weightsum