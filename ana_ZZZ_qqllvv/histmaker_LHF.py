import os, copy
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

# list of processes (mandatory)
processList = {
    # cross sections given on the webpage: https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/ 
    'wzp6_ee_qqH_HZZ_llvv_ecm240': {'fraction':1, 'crossSection': 0.00015, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"}, # 
    'p8_ee_ZZ_ecm240':          {'fraction':1},
    'wzp6_ee_qqH_Htautau_ecm240': {'fraction': 1}, 
    'wzp6_ee_ssH_Htautau_ecm240': {'fraction': 1},
    'wzp6_ee_ccH_Htautau_ecm240': {'fraction': 1},
    'wzp6_ee_bbH_Htautau_ecm240': {'fraction': 1},
    # signal permutation
    'wzp6_ee_eeH_HZZ_ecm240': {'fraction': 1},
    'wzp6_ee_mumuH_HZZ_ecm240': {'fraction': 1},
    'wzp6_ee_nunuH_HZZ_ecm240': {'fraction': 1},
    # no BDT trained on 
    'wzp6_ee_qqH_HWW_ecm240':   {'fraction':1}, # q = u, d
    'wzp6_ee_ssH_HWW_ecm240':   {'fraction':1}, # s
    'wzp6_ee_ccH_HWW_ecm240':   {'fraction':1}, # c
    'wzp6_ee_bbH_HWW_ecm240':   {'fraction':1}, # b
    'p8_ee_WW_ecm240':          {'fraction':1},
    'p8_ee_Zqq_ecm240':         {'fraction':1}, # q = u,d,s,c,b,t 
    'wzp6_ee_qqH_Hbb_ecm240':  {'fraction':1}, # q = u, d
    'wzp6_ee_ssH_Hbb_ecm240':  {'fraction':1}, # s
    'wzp6_ee_ccH_Hbb_ecm240':  {'fraction':1}, # c
    'wzp6_ee_bbH_Hbb_ecm240':  {'fraction':1}, # b
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json" # QUESTION: is this correct?

#Optional: output directory, default is local running directory
outputDir   = "./outputs/histmaker_LHF/ZZZqqllvv/"

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

bins_count = (15, 0, 15)
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

    df = df.Filter("((muons_sel_iso.size() == 2 && Sum(muons_sel_q) == 0) || (electrons_sel_iso.size() == 2 && Sum(electrons_sel_q) == 0)) && (muons_sel_iso.size() + electrons_sel_iso.size() == 2)")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    # get the two leptons
    df = df.Define("l1", "muons_sel_iso.size() == 2 ? muons_sel_iso[0] : electrons_sel_iso[0]")
    df = df.Define("l2", "muons_sel_iso.size() == 2 ? muons_sel_iso[1] : electrons_sel_iso[1]")

    # save these params
    df = df.Define("l1_p", "FCCAnalyses::ReconstructedParticle::get_p(Vec_rp{l1})[0]")
    df = df.Define("l2_p", "FCCAnalyses::ReconstructedParticle::get_p(Vec_rp{l2})[0]")
    df = df.Define("l1_theta", "FCCAnalyses::ReconstructedParticle::get_theta(Vec_rp{l1})[0]")
    df = df.Define("l2_theta", "FCCAnalyses::ReconstructedParticle::get_theta(Vec_rp{l2})[0]")

    df = df.Define("res_ll", "FCCAnalyses::ZHfunctions::get_two_lep_res(l1, l2)")
    df = df.Define("m_ll", "FCCAnalyses::ReconstructedParticle::get_mass(res_ll)[0]")
    df = df.Define("recoil_ll", "FCCAnalyses::ZHfunctions::get_recoil_lep(240.0, l1, l2)")
    df = df.Define("m_recoil_ll", "FCCAnalyses::ReconstructedParticle::get_mass(recoil_ll)[0]")

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

    # calculate recoil mass of the two jets
    df = df.Define("recoil", "FCCAnalyses::ZHfunctions::get_recoil_jets(240.0, jet1, jet2)")
    df = df.Define("recoil_mass", "FCCAnalyses::ReconstructedParticle::get_mass(recoil)[0]")


    df = df.Define("missP", "FCCAnalyses::ZHfunctions::missingParticle(240.0, ReconstructedParticles)")
    df = df.Define("miss_p", "FCCAnalyses::ReconstructedParticle::get_p(missP)[0]")
    df = df.Define("miss_e", "FCCAnalyses::ReconstructedParticle::get_e(missP)[0]")
    df = df.Define("miss_pz", "FCCAnalyses::ReconstructedParticle::get_pz(missP)[0]")
    df = df.Define("miss_theta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(missP)")
    df = df.Define("miss_pT", "FCCAnalyses::ZHfunctions::miss_pT(missP)")







    #########
    ### CUT 2: rough cut on recoil mass of the two jets must match Higgs mass
    #########

    df = df.Filter("recoil_mass > 100 && recoil_mass < 170")
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



    # look at ll system
    df = df.Define("Zll_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_ll)")
    df = df.Define("Zll_p", "FCCAnalyses::ReconstructedParticle::get_p(res_ll)[0]")
    df = df.Define("Zll_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_ll)[0]")

    # look at jj system
    df = df.Define("Zjj_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_jj)")
    df = df.Define("Zjj_p", "FCCAnalyses::ReconstructedParticle::get_p(res_jj)[0]")
    df = df.Define("Zjj_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_jj)[0]")




    #########
    ### CUT 7: miss pT > 5 GeV and pT < 50 GeV
    #########

    df = df.Filter("miss_pT > 5 && miss_pT < 50")
    df = df.Define("cut7", "7")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))



    df = df.Define("dot_prod_had", "FCCAnalyses::ZHfunctions::dot_prod_had(missP, jet1, jet2)")
    df = df.Define("dot_prod_lep", "FCCAnalyses::ZHfunctions::dot_prod_lep(missP, l1, l2)")
    df = df.Define("dot_prod_ll", "FCCAnalyses::ZHfunctions::dot_prod_ll(l1, l2)")



    #########
    ### CUT 8: dot product of hadronic system and missing momentum
    #########
    # in signal, the dot product is smaller than 0.3 aka the neutrinos go into the opposite direction of the jets
    # in the background, the missing momentum comes from decays inside the jets, so it's more likely to be larger than 0.3 because it's aligned with the jets

    df = df.Filter("dot_prod_had < 0.3")
    df = df.Define("cut8", "8")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))


    # recoil mass for LHF
    bin_lhf = (70, 100, 170)
    results.append(df.Histo1D(("recoil_mass_LHF", "", *bin_lhf), "recoil_mass"))

    #########
    ### CUT 9: tighter cut on the recoil mass
    #########

    df = df.Filter("recoil_mass > 120 && recoil_mass < 140")
    df = df.Define("cut9", "9")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))



    doInference = False
    if doInference:
        # tmva_helper = TMVAHelperXGB("outputs/mva_multi/ZZZ_qqllvv/bdt_model_multi.root", "bdt_model") # read the XGBoost training
        # df = tmva_helper.run_inference(df, col_name="mva_score") # by default, makes a new column mva_score
        # df = df.Define("mva_score_ZZ", "mva_score[0]")
        # df = df.Define("mva_score_Htautau", "mva_score[1]")
        # df = df.Define("mva_score_signal", "mva_score[2]")

        tmva_helper = TMVAHelperXGB("outputs/mva_multi/ZZZ_qqllvv/bdt_model_multi_all_bkg.root", "bdt_model") # read the XGBoost training
        df = tmva_helper.run_inference(df, col_name="mva_score") # by default, makes a new column mva_score
        df = df.Define("mva_score_ZZ", "mva_score[3]")
        df = df.Define("mva_score_perm_sig", "mva_score[2]")
        df = df.Define("mva_score_Htautau", "mva_score[1]")
        df = df.Define("mva_score_signal", "mva_score[0]")

        #########
        ### CUT 9: cut on the mva score ZZ
        #########
        bins_mva = (100, 0, 1)
        results.append(df.Histo1D(("mva_score_ZZ", "", *bins_mva), "mva_score_ZZ"))
        results.append(df.Histo1D(("mva_score_Htautau", "", *bins_mva), "mva_score_Htautau"))
        results.append(df.Histo1D(("mva_score_signal", "", *bins_mva), "mva_score_signal"))
        results.append(df.Histo1D(("mva_score_perm_sig", "", *bins_mva), "mva_score_perm_sig"))


        df = df.Filter("mva_score_ZZ < 0.5")
        df = df.Define("cut9", "9")
        results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))


        #########
        ### CUT 10: cut on the mva score perm sig
        #########

        df = df.Filter("mva_score_perm_sig < 0.5")
        df = df.Define("cut10", "10")
        results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut10"))


        ### save some hists for LHF

        # signal score
        results.append(df.Histo1D(("mva_score_signal_LHF", "", *bins_mva), "mva_score_signal"))


    



    return results, weightsum