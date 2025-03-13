# list of processes (mandatory)
processList = {
    # cross sections given on the webpage: https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/ 
    'wzp6_ee_qqH_HZZ_ecm240':    {'fraction':1}, # 0.001409 pb -> 15200 events
    # 'wzp6_ee_qqH_HWW_ecm240':   {'fraction':1}, # 0.01148 pb  -> 186000 events
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["MCfunctions.h"]

# Define the input dir (optional)
# inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmakerMC/ZZZqqllvv/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1

# scale the histograms with the cross-section and integrated luminosity
# doScale = True
intLumi = 10800000 # 10.8 /ab


# define some binning for various histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (2000, 0, 200) # 100 MeV bins
bins_recoil = (2000, 110, 180) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, -5, 5)
bins_eta = (600, -3, 3)
bins_phi = (500, -5, 5)

bins_count = (10, 0, 10)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 3)

bins_higgs = (2000, 100, 135) # 100 MeV bins
bins_Z = (2000, 0, 200) # 100 MeV bins



# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    # define some aliases to be used later on
    df = df.Alias("Particle0", "Particle#0.index") # index of the parents particles
    df = df.Alias("Particle1", "Particle#1.index") # index of the daughters particles 
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")
    

    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))

    # want to plot the momentum distribution of the leptons of the on shell Z from Higgs
    df = df.Define("lep_from_os_Z", "FCCAnalyses::MCfunctions::get_leptons_from_onshell_Z_decay(Particle, Particle0, Particle1)") # Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers
    df = df.Define("p_lep_from_os_Z", "FCCAnalyses::MCParticle::get_p(lep_from_os_Z)")

    results.append(df.Histo1D(("p_lep_from_os_Z", "", *bins_p_mu), "p_lep_from_os_Z"))

    """
    THE PLAN:
    Find events with
    - two jets with p around 53 GeV
    - two leptons with inv m ~ Z mass
    - missing energy
    """


    return results, weightsum