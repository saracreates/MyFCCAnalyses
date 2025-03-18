import os, copy

# list of processes (mandatory)
processList = {
    # cross sections given on the webpage: https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/ 
    'wzp6_ee_qqH_HZZ_ecm240':    {'fraction':1}, # 0.001409 pb -> 15200 events
    'wzp6_ee_qqH_HWW_ecm240':   {'fraction':1}, # 0.01148 pb  -> 186000 events
    'p8_ee_ZZ_ecm240':          {'fraction':0.1},
    'p8_ee_WW_ecm240':          {'fraction':0.01},
    # 'wzp6_ee_qqH_Hbb_ecm240':  {'fraction':1},
    'wzp6_ee_qqH_Htautau_ecm240':  {'fraction':1},
    # 'p8_ee_Zqq_ecm240':         {'fraction':0.01},
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker/ZZZqqllvv/"

# additional/costom C++ functions, defined in header files (optional)
includePaths = ["functions.h"]

## latest particle transformer model, trained on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc_v1"

## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

## model files locally stored on /eos
model_dir = (
    "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
)
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

## get local file, else download from url
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)


weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import (
    ExclusiveJetClusteringHelper,
)

jetFlavourHelper = None
jetClusteringHelper = None

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
bins_cosThetaMiss = (1000, 0.0, 1)

bins_theta = (500, -5, 5)
bins_eta = (600, -3, 3)
bins_phi = (500, -5, 5)

bins_count = (10, 0, 10)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 3)

bins_higgs = (2000, 100, 135) # 100 MeV bins
bins_Z = (2000, 0, 200) # 100 MeV bins

bin_njets = (200, 0, 200)
bin_dotprod = (100, -1, 1)
bin_dotprod_cut = (2, 0, 2)

bin_miss_pt = (1000, 0, 100) # 100 MeV bins



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
    df = df.Define(
        "muons",
        "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)",
    )
    df = df.Define(
        "electrons",
        "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)",
    )

    # compute the muon isolation and store muons with an isolation cut of 0df = df.25 in a separate column muons_sel_iso
    df = df.Define(
        "muons_iso",
        "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(muons, ReconstructedParticles)",
    )
    df = df.Define(
        "muons_sel_iso",
        "FCCAnalyses::ZHfunctions::sel_iso(0.5)(muons, muons_iso)",
    )
    df = df.Define(
        "muons_sel_q",
        "FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)",
    )
    df = df.Define(
        "electrons_iso",
        "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(electrons, ReconstructedParticles)",
    )
    df = df.Define(
        "electrons_sel_iso",
        "FCCAnalyses::ZHfunctions::sel_iso(0.5)(electrons, electrons_iso)",
    )
    df = df.Define(
        "electrons_sel_q",
        "FCCAnalyses::ReconstructedParticle::get_charge(electrons_sel_iso)",
    )

    # plot iso values
    results.append(df.Histo1D(("muons_iso", "", *bins_iso), "muons_iso"))
    results.append(df.Histo1D(("electrons_iso", "", *bins_iso), "electrons_iso"))

    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))
    

    #########
    ### CUT 1: exactely two isolated leptons of opposite charge and same flavor
    #########

    ## FALSCH, koennten auch 4 iso leptons sein 
    df = df.Filter("((muons_sel_iso.size() == 2 && Sum(muons_sel_q) == 0) || (electrons_sel_iso.size() == 2 && Sum(electrons_sel_q) == 0)) && (muons_sel_iso.size() + electrons_sel_iso.size() == 2)")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    # get the two leptons
    df = df.Define("l1", "muons_sel_iso.size() == 2 ? muons_sel_iso[0] : electrons_sel_iso[0]")
    df = df.Define("l2", "muons_sel_iso.size() == 2 ? muons_sel_iso[1] : electrons_sel_iso[1]")

    df = df.Define("res_ll", "FCCAnalyses::ZHfunctions::get_two_lep_res(l1, l2)")
    df = df.Define("m_ll", "FCCAnalyses::ReconstructedParticle::get_mass(res_ll)[0]")

    # plot m_ll
    results.append(df.Histo1D(("m_ll", "", *bins_m_ll), "m_ll"))


    ## here cluster jets in the events but first remove muons from the list of
    ## reconstructed particles

    ## create a new collection of reconstructed particles removing muons with p>20
    df = df.Define(
        "ReconstructedParticlesNoMuons",
        "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, muons_sel_iso)",
    )
    df = df.Define(
        "ReconstructedParticlesNoLeptons",
        "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticlesNoMuons, electrons_sel_iso)",
    )

    ## perform N=2 jet clustering
    global jetClusteringHelper
    global jetFlavourHelper

    ## define jet and run clustering parameters
    ## name of collections in EDM root files
    collections = {
        "GenParticles": "Particle",
        "PFParticles": "ReconstructedParticles",
        "PFTracks": "EFlowTrack",
        "PFPhotons": "EFlowPhoton",
        "PFNeutralHadrons": "EFlowNeutralHadron",
        "TrackState": "EFlowTrack_1",
        "TrackerHits": "TrackerHits",
        "CalorimeterHits": "CalorimeterHits",
        "dNdx": "EFlowTrack_2",
        "PathLength": "EFlowTrack_L",
        "Bz": "magFieldBz",
    }

    collections_noleptons = copy.deepcopy(collections)
    collections_noleptons["PFParticles"] = "ReconstructedParticlesNoLeptons"


    # https://github.com/HEP-FCC/FCCAnalyses/blob/51488171d3492afd407b331083b2bd0bb231417a/addons/FastJet/python/jetClusteringHelper.py#L6 
    jetClusteringHelper = ExclusiveJetClusteringHelper(collections_noleptons["PFParticles"], 2, "N2")
    df = jetClusteringHelper.define(df)

    ## define jet flavour tagging parameters

    jetFlavourHelper = JetFlavourHelper(
        collections_noleptons,
        jetClusteringHelper.jets,
        jetClusteringHelper.constituents,
    )

    ## define observables for tagger
    df = jetFlavourHelper.define(df)

    ## tagger inference
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)


    # plot number of jets - not possible! Only the value of how good the two jet clustering was! 

    df = df.Define("y23", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet_N2, 2))")  # dmerge from 3 to 2
    df = df.Define("y34", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet_N2, 3))")  # dmerge from 4 to 3

    i = 2
    for j in range(1, 3):
        df = df.Define(f"jet{j}_nconst_N{i}", f"jet_nconst_N{i}[{j-1}]")


    df = df.Define(
        "jets_p4",
        "JetConstituentsUtils::compute_tlv_jets({})".format(
            jetClusteringHelper.jets
        ),
    )

    df = df.Define(
        "m_jj",
        "JetConstituentsUtils::InvariantMass(jets_p4[0], jets_p4[1])",
    )

    # add two lorentz vectors for the two jets to get Z resonance?
    df = df.Define("jet1", "jets_p4[0]")
    df = df.Define("jet2", "jets_p4[1]")
    df = df.Define("res_jj", "FCCAnalyses::ZHfunctions::get_two_jets_res(jet1, jet2)")
    df = df.Define("p_res_jj", "FCCAnalyses::ReconstructedParticle::get_p(res_jj)[0]")

    # plot the invariant mass of the two jets
    results.append(df.Histo1D(("m_jj", "", *bins_Z), "m_jj"))
    results.append(df.Histo1D(("p_res_jj", "", *bins_p_mu), "p_res_jj"))
    # calculate recoil mass of the two jets
    df = df.Define("recoil", "FCCAnalyses::ZHfunctions::get_recoil_jets(240.0, jet1, jet2)")
    df = df.Define("recoil_mass", "FCCAnalyses::ReconstructedParticle::get_mass(recoil)[0]")

    results.append(df.Histo1D(("recoil_mass", "", *bins_m_ll), "recoil_mass"))


    results.append(df.Histo1D(("y23", "", *bin_njets), "y23"))
    results.append(df.Histo1D(("y34", "", *bin_njets), "y34"))

    results.append(df.Histo1D(("jet1_nconst_N2", "", *bin_njets), "jet1_nconst_N2"))
    results.append(df.Histo1D(("jet2_nconst_N2", "", *bin_njets), "jet2_nconst_N2"))


    df = df.Define("missP", "FCCAnalyses::ZHfunctions::missingParticle(240.0, ReconstructedParticles)")
    df = df.Define("miss_p", "FCCAnalyses::ReconstructedParticle::get_p(missP)[0]")
    df = df.Define("miss_e", "FCCAnalyses::ReconstructedParticle::get_e(missP)[0]")
    df = df.Define("miss_pz", "FCCAnalyses::ReconstructedParticle::get_pz(missP)[0]")
    df = df.Define("miss_theta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(missP)")
    df = df.Define("miss_pT", "FCCAnalyses::ZHfunctions::miss_pT(missP)")


    results.append(df.Histo1D(("miss_p_cut1", "", *bins_p_mu), "miss_p"))
    results.append(df.Histo1D(("miss_pT_cut1", "", *bins_p_mu), "miss_pT"))
    results.append(df.Histo1D(("miss_pz_cut1", "", *bins_p_mu), "miss_pz"))





    #########
    ### CUT 2: recoil mass of the two jets must match Higgs mass
    #########

    df = df.Filter("recoil_mass > 120 && recoil_mass < 140")
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))



    #########
    ### CUT 3: lepton invariant mass around Z mass
    #########
    df = df.Filter("m_ll > 80 && m_ll < 100") # WW has background smaller 80/90 GeV
    df = df.Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    #########
    ### CUT 4: m_jj between 85 and 105 GeV - should cut away WW background
    #########
    df = df.Filter("m_jj > 85 && m_jj < 105")
    df = df.Define("cut4", "4")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #########
    ### CUT 5: p_jj > 43 GeV and < 55 GeV
    #########

    df = df.Filter("p_res_jj > 40 && p_res_jj < 55")
    df = df.Define("cut5", "5")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5")) 


    #########
    ### CUT 6: inv mass of ll jj system must be the offshell Z
    #########

    # calculate invariant mass of llqq system which gives the offshell Z and apply prop cut on recoil mass ~ 10-50 GeV 
    df = df.Define("recoil_jjll", "FCCAnalyses::ZHfunctions::get_recoil_from_lep_and_jets(240.0, jet1, jet2, l1, l2)")
    df = df.Define("recoil_mass_jjll", "FCCAnalyses::ReconstructedParticle::get_mass(recoil_jjll)[0]")

    results.append(df.Histo1D(("recoil_mass_jjll", "", *bins_m_ll), "recoil_mass_jjll"))


    df = df.Filter("recoil_mass_jjll > 10 && recoil_mass_jjll < 50")
    df = df.Define("cut6", "6")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))



    # look at miss system
    results.append(df.Histo1D(("miss_e_cut6", "", *bins_p_mu), "miss_e"))
    results.append(df.Histo1D(("miss_theta_cut6", "", *bins_cosThetaMiss), "miss_theta"))
    results.append(df.Histo1D(("miss_p_cut6", "", *bins_p_mu), "miss_p"))
    results.append(df.Histo1D(("miss_pT_cut6", "", *bins_p_mu), "miss_pT"))
    results.append(df.Histo1D(("miss_pz_cut6", "", *bins_p_mu), "miss_pz"))


    # look at ll system
    df = df.Define("Zll_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_ll)")
    df = df.Define("Zll_p", "FCCAnalyses::ReconstructedParticle::get_p(res_ll)[0]")
    df = df.Define("Zll_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_ll)[0]")

    results.append(df.Histo1D(("Zll_costheta_cut6", "", *bins_cosThetaMiss), "Zll_costheta"))
    results.append(df.Histo1D(("Zll_p_cut6", "", *bins_p_mu), "Zll_p"))
    results.append(df.Histo1D(("Zll_pT_cut6", "", *bins_p_mu), "Zll_pT"))

    # look at jj system
    df = df.Define("Zjj_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_jj)")
    df = df.Define("Zjj_p", "FCCAnalyses::ReconstructedParticle::get_p(res_jj)[0]")
    df = df.Define("Zjj_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_jj)[0]")
    results.append(df.Histo1D(("Zjj_costheta_cut6", "", *bins_cosThetaMiss), "Zjj_costheta"))
    results.append(df.Histo1D(("Zjj_p_cut6", "", *bins_p_mu), "Zjj_p"))
    results.append(df.Histo1D(("Zjj_pT_cut6", "", *bins_p_mu), "Zjj_pT"))



    #########
    ### CUT 7: miss pT > 5 GeV and pT < 50 GeV
    #########

    df = df.Filter("miss_pT > 5 && miss_pT < 50")
    df = df.Define("cut7", "7")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))



    df = df.Define("dot_prod_had", "FCCAnalyses::ZHfunctions::dot_prod_had(missP, jet1, jet2)")
    df = df.Define("dot_prod_lep", "FCCAnalyses::ZHfunctions::dot_prod_lep(missP, l1, l2)")
    results.append(df.Histo1D(("dot_prod_had", "", *bin_dotprod), "dot_prod_had"))
    results.append(df.Histo1D(("dot_prod_lep", "", *bin_dotprod), "dot_prod_lep"))



    #########
    ### CUT 8: dot product of hadronic system and missing momentum
    #########
    # in signal, the dot product is smaller than 0.3 aka the neutrinos go into the opposite direction of the jets
    # in the background, the missing momentum comes from decays inside the jets, so it's more likely to be larger than 0.3 because it's aligned with the jets

    df = df.Filter("dot_prod_had < 0.3")
    df = df.Define("cut8", "8")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))

    results.append(df.Histo1D(("dot_prod_lep_cut8", "", *bin_dotprod), "dot_prod_lep"))


    #########
    ### CUT 9: dot product of leptonic system and missing momentum
    #########

    df = df.Filter("dot_prod_lep > -0.8")
    df = df.Define("cut9", "9")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))


    ## check miss p after miss pT cut
    results.append(df.Histo1D(("miss_p_cut9", "", *bins_p_mu), "miss_p"))


    #########
    ### CUT 10: 5 < p_miss < 50 GeV
    #########

    df = df.Filter("miss_p > 5 && miss_p < 50")
    df = df.Define("cut10", "10")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut10"))









    return results, weightsum