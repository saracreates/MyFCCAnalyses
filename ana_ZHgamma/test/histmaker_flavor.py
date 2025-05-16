# from Lena 11.April 2025, see: https://github.com/herrmannlena/FCCAnalyses/blob/higgsgamma/myanalysis/histmaker_recoil.py 

import os, copy


"""
from Juraj 03. Mai 2025: 
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/ZH/CLD_o2_v05/rec/00016944
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/eegamma/CLD_o2_v05/rec/00016945
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/mumugamma/CLD_o2_v05/rec/00016946
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/tautaugamma/CLD_o2_v05/rec/00016947
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/qqgamma/CLD_o2_v05/rec/00016948
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/ccgamma/CLD_o2_v05/rec/00016949
/eos/experiment/fcc/prod/fcc/ee/test_spring2024/240gev/bbgamma/CLD_o2_v05/rec/00016950

The tagger was included in CLDConfig when creating this data.

As these are not nice to link, I've created simlinks here: /afs/cern.ch/work/s/saaumill/public/analyses/Hgamma_fullsim_simlink/exlusive_bb
This folder has folders with the names of the processes, and inside are the simlinks to the files.
"""

input_base = "/afs/cern.ch/work/s/saaumill/public/analyses/Hgamma_fullsim_simlink/exclusive_bb"

# list of processes (mandatory)
processList = {
    'p8_ee_Hgamma_ecm240':    {'fraction':1, 'crossSection': 8.20481e-05, 'inputDir': input_base},  #what are the exact values here?
    #'reco_higgsgamma_test_REC.edm4hep': {'fraction':1, 'crossSection': 8.20481e-05, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/tmp_fullsim_output/reco_higgsgamma"},
    # 'p8_ee_qqgamma_ecm240':    {'fraction':1, 'crossSection': 6.9, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_ccgamma_ecm240':    {'fraction':1, 'crossSection': 2.15, 'inputDir': input_base},  #what are the exact values here?
    'p8_ee_bbgamma_ecm240':    {'fraction':1, 'crossSection': 2.35, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_ZH_ecm240':              {'fraction':1, 'crossSection': 0.2, 'inputDir': input_base},  #what are the exact values here?
    # # 'p8_ee_WW_ecm240':    {'fraction':1},  
    # # 'p8_ee_ZZ_ecm240':    {'fraction':1}, 
    # 'p8_ee_eegamma_ecm240':    {'fraction':1, 'crossSection': 190, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_tautaugamma_ecm240':    {'fraction':1, 'crossSection': 0.77, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_mumugamma_ecm240':    {'fraction':1, 'crossSection': 0.8, 'inputDir': input_base},  #what are the exact values here?
}

ecm= 240
# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

# additional/custom C++ functions, defined in header files (optional)
includePaths = ["../functions.h"]

# Define the input dir (optional)
#inputDir    = "outputs/FCCee/higgs/mH-recoil/mumu/stage1"
#inputDir    = "/afs/cern.ch/work/l/lherrman/private/HiggsGamma/data"

#Optional: output directory, default is local running directory

# outputDir   = "/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma_btag/"
outputDir   = "/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma_btag_test/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000  # 10.8 /ab

## latest particle transformer model, trained on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc"

## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)


## model files locally stored on /eos
# model_dir = (
#     "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"
# )
# local_preproc = "{}/{}.json".format(model_dir, model_name)
# local_model = "{}/{}.onnx".format(model_dir, model_name)

### TEST 

local_model = "/eos/experiment/fcc/ee/jet_flavour_tagging/fullsim_test_spring2024/fullsimCLD240_2mio.onnx"
local_preproc = "/eos/experiment/fcc/ee/jet_flavour_tagging/fullsim_test_spring2024/preprocess_fullsimCLD240_2mio.json"



## get local file, else download from url
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        #urllib.request.urlretrieve(url, os.path.basename(url))
        #return os.path.basename(url)
        raise ValueError("Model not available locally")


weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import (
    ExclusiveJetClusteringHelper,
)

jetFlavourHelper = None
jetClusteringHelper = None


# define some binning for various histograms
bins_a_p = (100, 0, 500) # 100 MeV bins
bins_a_n = (10, 0, 10) # 100 MeV bins

bins_count = (10, 0, 10)

# helper dictionary 
collections = {
    "GenParticles": "Particle",
    "PFParticles": "ReconstructedParticles",
    "PFTracks": "EFlowTrack",
    "PFPhotons": "EFlowPhoton",
    "PFNeutralHadrons": "EFlowNeutralHadron",
    "Tracks": "SiTracks_Refitted",
    "TrackStates": "_SiTracks_Refitted_trackStates",
    "TrackerHits": "TrackerHits",
    "CalorimeterHits": "CalorimeterHits",
    "dNdx": "EFlowTrack_2",
    "PathLength": "EFlowTrack_L",
    "Bz": "magFieldBz",
}



# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):

    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    # full sim defines

    # in FullSim, both the reco and gen particles are produced with a crossing angle
    df = df.Define("ReconstructedParticles", "FCCAnalyses::unBoostCrossingAngle(PandoraPFOs, -0.015)")
    df = df.Define("Particle", "FCCAnalyses::unBoostCrossingAngle(MCParticles, -0.015)")

    df = df.Define("photons_all", "FCCAnalyses::sel_type(22, ReconstructedParticles)")
    df = df.Define("electrons_all", f"FCCAnalyses::sel_type(11, ReconstructedParticles)")


    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)") 
    df = df.Define("photons_n","FCCAnalyses::ReconstructedParticle::get_n(photons_all)")  #number of photons per event
    df = df.Define("photons_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_all))")
    

    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)") 
    df = df.Define("electrons_n","FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")  #number of photons per event
    df = df.Define("electrons_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(electrons_all))")



    # order the cos theta values, and return arrays, when filter require length 2 for cut!

    # get cos theta from electrons
    df = df.Define("electrons_ordered_cos_theta","FCCAnalyses::ZHfunctions::ee_costheta_max(electrons_cos_theta)")
   
    #print and check
    #df = df.Define("photons_print", "FCCAnalyses::ZHfunctions::print_momentum(electrons_all)")
    #results.append(df.Histo1D(("photons_print", "", 100, 0, 100), "photons_print"))




    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    #Baseline selection
    results.append(df.Histo1D(("photons_p_cut_0", "", 130, 0, 130), "photons_p"))
    results.append(df.Histo1D(("photons_n_cut_0", "", *bins_a_n), "photons_n"))
    results.append(df.Histo1D(("photons_cos_theta_cut_0", "", 50, -1, 1), "photons_cos_theta"))

    results.append(df.Histo1D(("electrons_p_baseline", "", 130, 0, 130), "electrons_p"))
    results.append(df.Histo1D(("electrons_n_baseline", "", *bins_a_n), "electrons_n"))
    results.append(df.Histo1D(("electrons_cos_theta", "", 50, -1, 1), "electrons_ordered_cos_theta"))

    results.append(df.Histo1D(("photons_p_all", "", 40, 60, 100), "photons_p"))

    # fix the energy of the photons by shifting it by 2.6 GeV

    df = df.Define("photons_all_shifted", "FCCAnalyses::shift_E_photons(photons_all)")

    df = df.Define("photons_all_shifted_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all_shifted)")
    results.append(df.Histo1D(("photons_p_shifted", "", 40, 60, 100), "photons_all_shifted_p"))

   

    #isolation cut
    df = df.Define("photons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(photons_all_shifted, ReconstructedParticles)")  # is this correct?
    df = df.Define("photons_sel_iso","FCCAnalyses::ZHfunctions::sel_iso(0.2)(photons_all_shifted, photons_iso)",) # and this??
   
    df = df.Define("photons_iso_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_sel_iso)") 
    df = df.Define("photons_iso_n","FCCAnalyses::ReconstructedParticle::get_n(photons_sel_iso)")  #number of photons per event
    df = df.Define("photons_iso_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_sel_iso))")

    results.append(df.Histo1D(("photon_isolation", "", 50, 0, 10), "photons_iso"))

     
    #########
    ### CUT 1: Photons must be isolated
    #########
    
    df = df.Filter("photons_sel_iso.size()>0 ")  
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))
    
    results.append(df.Histo1D(("photons_p_cut_1", "",  130, 0, 130), "photons_iso_p"))
    results.append(df.Histo1D(("photons_n_cut_1", "", *bins_a_n), "photons_iso_n"))
    results.append(df.Histo1D(("photons_cos_theta_cut_1", "", 50, -1, 1), "photons_iso_cos_theta"))
    results.append(df.Histo1D(("electrons_cos_theta_cut_1", "", 50, -1, 1), "electrons_ordered_cos_theta"))
 
    

    #sort in p  and select highest energetic one
    df = df.Define("iso_highest_p","FCCAnalyses::ZHfunctions::sort_by_energy(photons_sel_iso)")

    #print and check
    #df = df.Define("photons_print", "FCCAnalyses::ZHfunctions::print_momentum(iso_highest_p)")
    #results.append(df.Histo1D(("photons_print", "", 100, 0, 100), "photons_print"))



    #energy cut
    df = df.Define("photons_boosted", "FCCAnalyses::ReconstructedParticle::sel_p(60,100)(iso_highest_p)") # looked okay from photons all
    #df = df.Define("photons_boosted", "FCCAnalyses::ReconstructedParticle::sel_p(60,100)(iso_highest_p)")

    df = df.Define("photons_boosted_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_boosted)") # is this correct?
    df = df.Define("photons_boosted_n","FCCAnalyses::ReconstructedParticle::get_n(photons_boosted)") 
    df = df.Define("photons_boosted_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_boosted))")

    
    #########
    ### CUT 2: Photons energy > 50
    #########
    
    df = df.Filter("photons_boosted.size()>0 ")  
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))
    
    results.append(df.Histo1D(("photons_p_cut_2", "",  130, 0, 130), "photons_boosted_p"))
    results.append(df.Histo1D(("photons_n_cut_2", "", *bins_a_n), "photons_boosted_n"))
    results.append(df.Histo1D(("photons_cos_theta_cut_2", "", 50, -1, 1), "photons_boosted_cos_theta"))
    results.append(df.Histo1D(("electrons_cos_theta_cut_2", "", 50, -1, 1), "electrons_ordered_cos_theta"))
 

    
    #########
    ### CUT 3: Cos Theta cut
    #########
    df = df.Filter("ROOT::VecOps::All(abs(photons_boosted_cos_theta) < 0.9) ") 
   
    df = df.Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))
    
    results.append(df.Histo1D(("photons_p_cut_3", "", 130, 0, 130), "photons_boosted_p"))
    results.append(df.Histo1D(("photons_n_cut_3", "", *bins_a_n), "photons_boosted_n"))
    results.append(df.Histo1D(("photons_cos_theta_cut_3", "", 50, -1, 1), "photons_boosted_cos_theta"))
    results.append(df.Histo1D(("electrons_cos_theta_cut_3", "", 50, -1, 1), "electrons_ordered_cos_theta"))


    df = df.Define("recopart_no_gamma", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons_boosted)",)
    df = df.Define("recopart_no_gamma_n","FCCAnalyses::ReconstructedParticle::get_n(recopart_no_gamma)") 
   
 
    results.append(df.Histo1D(("recopart_no_gamma_n_cut_0", "", 60, 0, 60), "recopart_no_gamma_n"))

    """
    #########
    ### CUT 4: Cos Theta cut on ee to reduce bhabhar
    #########
    df = df.Filter("electrons_ordered_cos_theta.size()<2 || (abs(electrons_ordered_cos_theta)[0]<0.8 && abs(electrons_ordered_cos_theta)[1]<0.8)") 
    
    df = df.Define("cut4", "4")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))
 
    results.append(df.Histo1D(("photons_p_cut_4", "", 130, 0, 130), "photons_boosted_p"))
    results.append(df.Histo1D(("photons_n_cut_4", "", *bins_a_n), "photons_boosted_n"))
    results.append(df.Histo1D(("photons_cos_theta_cut_4", "", 50, -1, 1), "photons_boosted_cos_theta"))
    results.append(df.Histo1D(("electrons_cos_theta_cut_4", "", 50, -1, 1), "electrons_ordered_cos_theta"))
    """

    # checks before recoil plot
    results.append(df.Histo1D(("num_boosted_gamma_for_recoil", "", 5, 0, 5), "photons_boosted_n"))
    results.append(df.Histo1D(("p_boosted_gamma_for_recoil", "", 40, 60, 100), "photons_boosted_p"))
    # check gen particles
    df = df.Define("gen_photons", "FCCAnalyses::sel_type(22, Particle)")
    df = df.Define("gen_photons_n","FCCAnalyses::MCParticle::get_n(gen_photons)")  #number of photons per event
    df = df.Define("gen_photons_p", "FCCAnalyses::MCParticle::get_p(gen_photons)")
    results.append(df.Histo1D(("gen_photons_p", "", 40, 60, 100), "gen_photons_p"))
    results.append(df.Histo1D(("gen_photons_p_fullrange", "", 100, 0, 100), "gen_photons_p"))
    results.append(df.Histo1D(("gen_photons_n", "", 5, 0, 5), "gen_photons_n"))


    # recoil mass
    df = df.Define("gamma_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(photons_boosted)") 
    df = df.Define("gamma_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(gamma_recoil)[0]") # recoil mass
    results.append(df.Histo1D(("gamma_recoil_m_cut_3", "", 170, 80, 250), "gamma_recoil_m"))
    
    #########
    ### CUT 4: require at least 10 (fast sim: 6) reconstructed particles (except gamma)
    #########
    df = df.Filter(" recopart_no_gamma_n > 9") # 5 in fast sim
    
    df = df.Define("cut4", "4")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))
 
    results.append(df.Histo1D(("recopart_no_gamma_n_cut_4", "", 60, 0, 60), "recopart_no_gamma_n"))
    

    
    results.append(df.Histo1D(("gamma_recoil_m_cut_4", "", 170, 80, 250), "gamma_recoil_m"))
   


    #########
    ### CUT 5: gamma recoil cut
    #########
    df = df.Filter("110 < gamma_recoil_m && gamma_recoil_m < 150") 
    # df = df.Filter("120 < gamma_recoil_m && gamma_recoil_m < 132")
    #df = df.Filter("115 < gamma_recoil_m && gamma_recoil_m < 170") 

    df = df.Define("cut5", "5")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    results.append(df.Histo1D(("gamma_recoil_m_signal_cut", "", 40, 110, 150), "gamma_recoil_m"))
    #results.append(df.Histo1D(("gamma_recoil_m_signal_cut", "", 64, 116, 170), "gamma_recoil_m"))
   
    #########
    ### CUT 6: gamma recoil cut tight
    #########
    #df = df.Filter("123.5 < gamma_recoil_m && gamma_recoil_m < 126.5") 
    # df = df.Filter("123.5 < gamma_recoil_m && gamma_recoil_m < 126.5") 

    # df = df.Define("cut6", "6")
    # results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))

    # results.append(df.Histo1D(("gamma_recoil_m_tight_cut", "", 70, 80, 150), "gamma_recoil_m"))


    #########
    ### Cut 6: On b-tagging
    #########

    ### THE CLUSTERING & TAGGING 

    ## perform N=2 jet clustering
    global jetClusteringHelper
    global jetFlavourHelper

    collections_nogamma = copy.deepcopy(collections)
    collections_nogamma["PFParticles"] = "recopart_no_gamma"

    jetClusteringHelper = ExclusiveJetClusteringHelper(collections_nogamma["PFParticles"], 2, "N2")
    df = jetClusteringHelper.define(df)

    ## define jet flavour tagging parameters

    jetFlavourHelper = JetFlavourHelper(
        collections_nogamma,
        jetClusteringHelper.jets,
        jetClusteringHelper.constituents,
        sim_type="full",
    )

    ## define observables for tagger
    df = jetFlavourHelper.define(df)

    ## tagger inference
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

    #df = df.Define("y23", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet_N2, 2))")

    df = df.Define(
        "jets_p4",
        "JetConstituentsUtils::compute_tlv_jets({})".format(
            jetClusteringHelper.jets
        ),
    )
    df = df.Define(
        "jj_m",
        "JetConstituentsUtils::InvariantMass(jets_p4[0], jets_p4[1])",
    )



    ### THE TAGGGING CUTS

    df = df.Define("b_tags_sum", "recojet_isB[0] + recojet_isB[1]")
    results.append(df.Histo1D(("b_tags_sum", "", 100, 0, 2), "b_tags_sum"))

    # cut 
    # df = df.Filter("b_tags_sum > 1") # cut on btag score 
    df = df.Filter("b_tags_sum > 0.6") 
    df = df.Define("cut6", "6")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))

    results.append(df.Histo1D(("gamma_recoil_m_LHF", "", 40, 110, 150), "gamma_recoil_m"))

    # look at inv mass of two b system?

    df = df.Define("jets", "FCCAnalyses::get_jet_res(RefinedVertexJets)") 

    df = df.Define("jets_p", "FCCAnalyses::ReconstructedParticle::get_p(jets)[0]")
    results.append(df.Histo1D(("jets_p", "", 100, 0, 240), "jets_p"))

    df = df.Define("jets_m", "FCCAnalyses::ReconstructedParticle::get_mass(jets)[0]") # recoil mass
    results.append(df.Histo1D(("jets_m", "", 100, 0, 240), "jets_m"))

    # df = df.Define("jets_p4", "FCCAnalyses::return_p4_jets(RefinedVertexJets)") # recoil mass
    # df = df.Define("m_jj", "JetConstituentsUtils::InvariantMass(jets_p4[0], jets_p4[1])")
    # results.append(df.Histo1D(("m_jj", "", 100, 0, 240), "m_jj"))


    ########
    ### Cut 7: tighter cut around recoil mass
    ########
    # df = df.Filter("114 < gamma_recoil_m && gamma_recoil_m < 128") # photon E not calibrated
    df = df.Filter("120 < gamma_recoil_m && gamma_recoil_m < 132") # phototn E shifted to Higgs mass
    df = df.Define("cut7", "7")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))


   
    # define further variables for plotting
    # df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    # df = df.Define("photons_boosted_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_boosted)")
    # df = df.Define("photons_boosted_n","FCCAnalyses::ReconstructedParticle::get_n(photons_boosted)")  #number of photons per event
    
   
    #select highest energetic photon

    ########################
    # Final histograms
    ########################
    #results.append(df.Histo1D(("photons_all_p", "", 100, 0, 100), "photons_all_p"))
   # results.append(df.Histo1D(("photons_boosted_p", "", *bins_a_p), "photons_boosted_p"))
    #results.append(df.Histo1D(("photons_n", "", *bins_a_n), "photons_n"))
    #results.append(df.Histo1D(("photons_boosted_n", "", *bins_a_n), "photons_boosted_n"))
    
    #results.append(df.Histo1D(("zmumu_recoil_m", "", *bins_recoil), "zmumu_recoil_m"))   # see how recoil determined
   
    #need to select the highest energetic photon

    return results, weightsum

