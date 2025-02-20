# list of processes (mandatory)
processList = {
    'wzp6_ee_llH_HZZ_llll_ecm240':    {'fraction':1},
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["ZHMCfunctions.h"]

# Define the input dir (optional)
inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker_MC/ZZZ6l/"


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

bins_pdg = (2000, -1000, 1000)



# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):
    """
    Goal of this exercise is to check the MC truth information for the ZZZ->6l decay.
    - print out all particles produced and their attributes
    - make p_Z histograms for Z from H and Z from production
    - make m_ll histograms for Z from H and Z from production
    """

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    # define some aliases to be used later on
    df = df.Alias("Particle0", "Particle#0.index") # Vec_i parents
    df = df.Alias("Particle1", "Particle#1.index") # Vec_i daugthers
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")


    # get all the leptons from the collection
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    # print info on MC particles
    df = df.Define("pdg_mc", "FCCAnalyses::ZHMCfunctions::print_MC_info(Particle)")

    # add histogram
    results.append(df.Histo1D(("pdg_mc", "", *bins_pdg), "pdg_mc"))



    

    return results, weightsum