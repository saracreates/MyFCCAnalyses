# list of processes (mandatory)
processList = {
    # change x sections!! 
    'wzp6_ee_llH_HZZ_llll_ecm365':    {'fraction':1, 'crossSection': 0.00000517}, # 5.17 ab (from Louis) 
    'wzp6_ee_llH_HZZ_qqll_ecm365':    {'fraction':1, 'crossSection': 0.000025}, # 25 ab 
    'wzp6_ee_qqH_HZZ_llll_ecm365':    {'fraction':1, 'crossSection': 0.000038}, # 38 ab 
    # 'p8_ee_llZZ_ecm240_edm4hep': {'fraction':1, 'crossSection': 0.00008746} # 8.746e-05 +- 1.176e-07 pb
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["ZHfunctions.h"]

# Define the input dir (optional)
# inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
#inputDir = "/afs/cern.ch/work/s/saaumill/public/madgraph_to_edm4hep_pipeline/pythia_to_delphes_to_edm4hep/"

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker/ZZZ6l/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1

# scale the histograms with the cross-section and integrated luminosity
# doScale = True
intLumi = 10600000 # 10.6 /ab


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
    df = df.Alias("Particle0", "Particle#0.index") # index of the daughter particles
    df = df.Alias("Particle1", "Particle#1.index") # index of the mother particles 
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")


    # get all the leptons from the collection
    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")

    
    # select leptons with momentum > 5 GeV
    # df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(muons_all)")
    # df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(electrons_all)")

    # muons
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    # electrons
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")
    
    # compute the muon isolation and store muons with an isolation cut of 0.25 in a separate column muons_sel_iso
    
    df = df.Define("muons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(muons, ReconstructedParticles)") # (0.01, 0.5) too tight?
    df = df.Define("muons_sel_iso", "FCCAnalyses::ZHfunctions::sel_iso(0.5)(muons, muons_iso)") # 0.25 too tight

    df = df.Define("muons_sel_iso_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)")
    # same for electrons
    df = df.Define("electrons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(electrons, ReconstructedParticles)")
    df = df.Define("electrons_sel_iso", "FCCAnalyses::ZHfunctions::sel_iso(0.5)(electrons, electrons_iso)") # 0.25 too tight
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
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))

    #########
    ### CUT 1: There must be equal or more than 6 electrons or muons in the event
    #########

    df = df.Filter("muons.size() + electrons.size() >= 6")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    #########
    ### CUT 2: There must be 6 isolated leptons (electrons or muons) in the event
    #########
    df = df.Filter("muons_sel_iso.size() + electrons_sel_iso.size() == 6")
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    # histograms after the first cut
    results.append(df.Histo1D(("muons_p_cut2", "", *bins_p_mu), "muons_p"))
    results.append(df.Histo1D(("muons_q_cut2", "", *bins_charge), "muons_q"))
    results.append(df.Histo1D(("muons_no_cut2", "", *bins_count), "muons_no"))
    results.append(df.Histo1D(("muons_iso_cut2", "", *bins_iso), "muons_iso"))
    # and electrons
    results.append(df.Histo1D(("electrons_p_cut2", "", *bins_p_mu), "electrons_p"))
    results.append(df.Histo1D(("electrons_q_cut2", "", *bins_charge), "electrons_q"))
    results.append(df.Histo1D(("electrons_no_cut2", "", *bins_count), "electrons_no"))
    results.append(df.Histo1D(("electrons_iso_cut2", "", *bins_iso), "electrons_iso"))

    #########
    ### CUT 3: 3 lepton pairs with opposite charge and right flavour
    #########
    df = df.Filter("Sum(muons_sel_iso_q) == 0 && Sum(electrons_sel_iso_q) == 0")
    df = df.Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))

    # histograms after the third cut
    results.append(df.Histo1D(("muons_p_cut3", "", *bins_p_mu), "muons_p"))
    results.append(df.Histo1D(("muons_q_cut3", "", *bins_charge), "muons_q"))
    results.append(df.Histo1D(("muons_no_cut3", "", *bins_count), "muons_no"))

    results.append(df.Histo1D(("electrons_p_cut3", "", *bins_p_mu), "electrons_p"))
    results.append(df.Histo1D(("electrons_q_cut3", "", *bins_charge), "electrons_q"))
    results.append(df.Histo1D(("electrons_no_cut3", "", *bins_count), "electrons_no"))

    # now we build the Z resonance based on the isolated leptons. 
    df = df.Define("zbuilder_result", "FCCAnalyses::ZHfunctions::resonanceBuilder(91.2)(electrons_sel_iso, muons_sel_iso)")
    df = df.Define("res1", "Vec_rp{zbuilder_result[0]}") # the first resonance
    df = df.Define("res2", "Vec_rp{zbuilder_result[1]}") # the second resonance
    df = df.Define("res3", "Vec_rp{zbuilder_result[2]}") # the third resonance

    # caluclate the Higgs mass from H->ZZ->4l
    df = df.Define("Higgs_results", "FCCAnalyses::ZHfunctions::higgsmassBuilder(125.0)(res1, res2, res3)")
    df = df.Define("Higgs", "Vec_rp{Higgs_results[0]}") # the Higgs candidate
    df = df.Define("Z_onshell_from_H", "Vec_rp{Higgs_results[1]}") 
    df = df.Define("Z_offshell_from_H", "Vec_rp{Higgs_results[2]}") 
    df = df.Define("Z_onshell", "Vec_rp{Higgs_results[3]}") 
    df = df.Define("m_higgs", "FCCAnalyses::ReconstructedParticle::get_mass(Higgs)[0]") # Higgs mass

    # Z mass
    df = df.Define("m_Z_onshell_from_H", "FCCAnalyses::ReconstructedParticle::get_mass(Z_onshell_from_H)") # Z mass
    df = df.Define("m_Z_offshell_from_H", "FCCAnalyses::ReconstructedParticle::get_mass(Z_offshell_from_H)") # Z mass
    df = df.Define("m_Z_onshell", "FCCAnalyses::ReconstructedParticle::get_mass(Z_onshell)") # Z mass

    # Z momentum 
    df = df.Define("p_Z_onshell_from_H", "FCCAnalyses::ReconstructedParticle::get_p(Z_onshell_from_H)") # Z momentum
    df = df.Define("p_Z_offshell_from_H", "FCCAnalyses::ReconstructedParticle::get_p(Z_offshell_from_H)") # Z momentum
    df = df.Define("p_Z_onshell", "FCCAnalyses::ReconstructedParticle::get_p(Z_onshell)") # Z momentum - peaks at 50 GeV like expected!! 

    # check the recoil mass
    df = df.Define("recoil", "FCCAnalyses::ZHfunctions::recoilBuilder(125, 240)(Z_onshell)")
    df = df.Define("m_recoil", "FCCAnalyses::ReconstructedParticle::get_mass(recoil)") # recoil mass

    # save in histogram 
    results.append(df.Histo1D(("higgs_mass", "", *bins_higgs), "m_higgs"))
    results.append(df.Histo1D(("recoil_mass", "", *bins_recoil), "m_recoil"))
    results.append(df.Histo1D(("m_Z_onshell_from_H", "", *bins_Z), "m_Z_onshell_from_H"))
    results.append(df.Histo1D(("m_Z_offshell_from_H", "", *bins_Z), "m_Z_offshell_from_H"))
    results.append(df.Histo1D(("m_Z_onshell", "", *bins_Z), "m_Z_onshell"))
    results.append(df.Histo1D(("p_Z_onshell_from_H", "", *bins_p_mu), "p_Z_onshell_from_H"))
    results.append(df.Histo1D(("p_Z_offshell_from_H", "", *bins_p_mu), "p_Z_offshell_from_H"))
    results.append(df.Histo1D(("p_Z_onshell", "", *bins_p_mu), "p_Z_onshell"))

    #########
    ### CUT 4: higgs_mass 124 < m < 125.5
    #########
    df = df.Filter("m_higgs > 124 && m_higgs < 125.5")
    df = df.Define("cut4", "4")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))

    # histograms after the fourth cut - hist on higgs mass 
    results.append(df.Histo1D(("higgs_mass_cut4", "", *bins_higgs), "m_higgs"))
    results.append(df.Histo1D(("recoil_mass_cut4", "", *bins_recoil), "m_recoil"))
    results.append(df.Histo1D(("m_Z_onshell_from_H_cut4", "", *bins_Z), "m_Z_onshell_from_H"))
    results.append(df.Histo1D(("m_Z_offshell_from_H_cut4", "", *bins_Z), "m_Z_offshell_from_H"))
    results.append(df.Histo1D(("m_Z_onshell_cut4", "", *bins_Z), "m_Z_onshell"))
    results.append(df.Histo1D(("p_Z_onshell_from_H_cut4", "", *bins_p_mu), "p_Z_onshell_from_H"))
    results.append(df.Histo1D(("p_Z_offshell_from_H_cut4", "", *bins_p_mu), "p_Z_offshell_from_H"))
    results.append(df.Histo1D(("p_Z_onshell_cut4", "", *bins_p_mu), "p_Z_onshell"))


    return results, weightsum