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
bins_recoil = (3000, 120, 150)#(200000, 0, 200) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, -5, 5)
bins_eta = (600, -3, 3)
bins_phi = (500, -5, 5)

bins_count = (10, 0, 10)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 3)

bins_pdg = (2000, -1000, 1000)
bins_missP = (2400, 0,80)



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
    # df = df.Define("pdg_mc", "FCCAnalyses::ZHMCfunctions::print_MC_info(Particle, Particle0, Particle1)")

    df = df.Define("leptons_from_Z", "FCCAnalyses::ZHMCfunctions::get_leptons_from_Z(Particle, Particle0, Particle1)")
    df = df.Define("recoil_from_Z", "FCCAnalyses::ZHMCfunctions::get_recoil_from_Z(leptons_from_Z, 240)")
    df = df.Define("m_recoil_from_Z", "FCCAnalyses::MCParticle::get_mass(recoil_from_Z)")
    
    # get the 3 Z bosons
    df = df.Define("Z", "FCCAnalyses::ZHMCfunctions::get_Z_from_leptons(leptons_from_Z)")
    df = df.Define("Z_onshell_from_H", "FCCAnalyses::ZHMCfunctions::get_Z_from_H(Particle, Particle0, Particle1, 0)")
    df = df.Define("Z_offshell_from_H", "FCCAnalyses::ZHMCfunctions::get_Z_from_H(Particle, Particle0, Particle1, 1)")

    # get the mass of the Z bosons
    df = df.Define("m_Z", "FCCAnalyses::MCParticle::get_mass(Z)")
    df = df.Define("m_Z_onshell_from_H", "FCCAnalyses::MCParticle::get_mass(Z_onshell_from_H)")
    df = df.Define("m_Z_offshell_from_H", "FCCAnalyses::MCParticle::get_mass(Z_offshell_from_H)")

    # get the momentum of the Z bosons
    df = df.Define("p_Z", "FCCAnalyses::MCParticle::get_p(Z)")
    df = df.Define("p_Z_onshell_from_H", "FCCAnalyses::MCParticle::get_p(Z_onshell_from_H)")
    df = df.Define("p_Z_offshell_from_H", "FCCAnalyses::MCParticle::get_p(Z_offshell_from_H)")


    # add histogram
    # results.append(df.Histo1D(("pdg_mc", "", *bins_pdg), "pdg_mc")) # for printing!!!

    # masses
    results.append(df.Histo1D(("m_Z", "", *bins_m_ll), "m_Z"))
    results.append(df.Histo1D(("m_Z_onshell", "", *bins_m_ll), "m_Z_onshell_from_H"))
    results.append(df.Histo1D(("m_Z_offshell", "", *bins_m_ll), "m_Z_offshell_from_H"))
    results.append(df.Histo1D(("m_recoil_from_Z", "", *bins_recoil), "m_recoil_from_Z"))

    # momentum
    results.append(df.Histo1D(("p_Z", "", *bins_p_mu), "p_Z"))
    results.append(df.Histo1D(("p_Z_onshell", "", *bins_p_mu), "p_Z_onshell_from_H"))
    results.append(df.Histo1D(("p_Z_offshell", "", *bins_p_mu), "p_Z_offshell_from_H"))

    # caluclate the missing momentum in these events 

    df = df.Define("ll_Z_onshell", "FCCAnalyses::ZHMCfunctions::get_leptons(Z_onshell_from_H, Particle1, Particle)")
    df = df.Define("ll_Z_offshell", "FCCAnalyses::ZHMCfunctions::get_leptons(Z_offshell_from_H, Particle1, Particle)")

    df = df.Define("part_missingE", "FCCAnalyses::ZHMCfunctions::create_lepton_missingE(leptons_from_Z, ll_Z_onshell, ll_Z_offshell, Particle1, Particle)")
    df = df.Define("part_missingE_emu_only", "FCCAnalyses::ZHMCfunctions::create_lepton_missingE_emu_only(leptons_from_Z, ll_Z_onshell, ll_Z_offshell, Particle1, Particle)")
    df = df.Define("miss_p", "FCCAnalyses::MCParticle::get_p(part_missingE)")
    df = df.Define("miss_p_emu_only", "FCCAnalyses::MCParticle::get_p(part_missingE_emu_only)")

    # plot histogram

    results.append(df.Histo1D(("missing_momentum", "", *bins_missP), "miss_p"))
    results.append(df.Histo1D(("missing_momentum_emu_only", "", *bins_missP), "miss_p_emu_only"))




    

    return results, weightsum