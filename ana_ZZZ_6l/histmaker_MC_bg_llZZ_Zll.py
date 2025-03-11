"""
Comment:
This file should check if the Z has a natural width when creating on shell Z with madgraph and then pythia8. 
Turns out 
- IT DOES NOT if the decay is done in pythia8
- IT DOES if the decay is done in madgraph

"""



# list of processes (mandatory)
processList = {
    'p8_ee_llZZ_ecm240': {'fraction':1}, # Z decay to ll done in pythia -> NO natural width
    'p8_ee_llZZ_Zll_ecm240': {'fraction':1}, # Z decay to ll done in madgraph -> natural width of Z!! 
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["ZHMCfunctions.h"]

# Define the input dir (optional)
# inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
inputDir = "/afs/cern.ch/work/s/saaumill/public/tmp_madgraph_output/edm4hep_data"

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker_bg_MC/ZZZ6l/"


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
    # df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    # df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")

    # require to have mc particles
    df = df.Filter("Particle.size() > 0")

    # get the two Z from the event  ee-> llZZ; Z->ll
    df = df.Define("Zs", "FCCAnalyses::ZHMCfunctions::get_all_Zs_bg(Particle, Particle0, Particle1)")
    df = df.Define("num_Z_daughters", "FCCAnalyses::ZHMCfunctions::n_daughters(Zs)")
    # cut whole event if number of Z daughters is not 



    df = df.Define("Z1", "Vec_mc{Zs[0]}")
    df = df.Define("Z2", "Vec_mc{Zs[1]}")
    df = df.Define("m_Z1", "FCCAnalyses::MCParticle::get_mass(Z1)")
    df = df.Define("m_Z2", "FCCAnalyses::MCParticle::get_mass(Z2)")

    # histograms
    results.append(df.Histo1D(("m_Z1", "m_Z1", 50, 80, 100), "m_Z1"))
    results.append(df.Histo1D(("m_Z2", "m_Z2", 50, 80, 100), "m_Z2"))


    df = df.Define("ll_Z1", "FCCAnalyses::ZHMCfunctions::get_leptons(Z1, Particle1, Particle)")
    df = df.Define("ll_Z2", "FCCAnalyses::ZHMCfunctions::get_leptons(Z2, Particle1, Particle)")
    df = df.Define("inv_mass_ll_Z1", "FCCAnalyses::ZHMCfunctions::calculate_inv_mass_of_two_leptons(ll_Z1)")
    df = df.Define("inv_mass_ll_Z2", "FCCAnalyses::ZHMCfunctions::calculate_inv_mass_of_two_leptons(ll_Z2)")

    # plot the histograms

    results.append(df.Histo1D(("inv_mass_ll_Z1", "inv_mass_ll_Z1", 200, 0, 200), "inv_mass_ll_Z1"))
    results.append(df.Histo1D(("inv_mass_ll_Z2", "inv_mass_ll_Z2", 200, 0, 200), "inv_mass_ll_Z2"))






    

    return results, weightsum