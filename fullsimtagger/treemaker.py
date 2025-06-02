import os, copy


# list of processes
input_base = '/afs/cern.ch/work/s/saaumill/public/analyses/symlink_Hxx_fullsim/Hxx/'
processList = {
    # u,d ,s ,c ,b, g, tau
    'Hbb':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Hbb data
    'Hcc':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Hcc data
    'Hgg':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Hgg data
    'Htautau':{'fraction':0.01, 'inputDir': input_base},  # test tagger on Htautau data
    'Huu':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Huu data
    'Hdd':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Hdd data
    'Hss':    {'fraction':0.01, 'inputDir': input_base},  # test tagger on Hss data
    
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/treemaker/fullsimtagger/"

# Define the input dir (optional)
#inputDir    = "./localSamples/"

# additional/costom C++ functions, defined in header files (optional)
includePaths = ["functions.h"]

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

bins_count = (10, 0, 10)
# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:

    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
       
        # __________________________________________________________
        # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
        results = []
        # in FullSim, both the reco and gen particles are produced with a crossing angle
        df = df.Define("ReconstructedParticles", "FCCAnalyses::unBoostCrossingAngle(PandoraPFOs, -0.015)")
        df = df.Define("Particle", "FCCAnalyses::unBoostCrossingAngle(MCParticles, -0.015)")
        df = df.Alias("Particle0", "_MCParticles_parents.index")
        df = df.Alias("Particle1", "_MCParticles_daughters.index")

        df = df.Define("photons_all", "FCCAnalyses::sel_type(22, ReconstructedParticles)")      
        
        
        ### THE CLUSTERING & TAGGING 

        ## perform N=2 jet clustering
        global jetClusteringHelper
        global jetFlavourHelper

        jetClusteringHelper = ExclusiveJetClusteringHelper(collections["PFParticles"], 2, "N2")
        df = jetClusteringHelper.define(df)

        ## define jet flavour tagging parameters

        jetFlavourHelper = JetFlavourHelper(
            collections,
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

        # check the true jet flavour
        df = df.Define("MC_pdg_flavour", "FCCAnalyses::fullsimtagger::get_higgs_daughters_MC_pdg(Particle, Particle1)")
        # df = df.Define("dummy", "FCCAnalyses::fullsimtagger::print_scores(recojet_isB)")


        # Rename columns appropriately
        flavors = {
            "U": 1, "D": 2, "S": 3, "B": 5,
            "C": 4, "G": 21, "TAU": 15,
        }

        # Rename predictions
        for flav in flavors:
            old_col = f"recojet_is{flav}"
            new_col = f"score_recojet_is{flav}"
            df = df.Define(new_col, old_col)

        # Define truth labels
        for flav, pdg in flavors.items():
            df = df.Redefine(f"recojet_is{flav}", f"FCCAnalyses::fullsimtagger::is_of_flavor(MC_pdg_flavour, {pdg})")


        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = []
        #     "ReconstructedParticles",
        #     "photons_all",
        #     "MC_pdg_flavour",
        # ]

        # save input variables to the network
        branchList += ["pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", "pfcand_dptdpt", "pfcand_detadeta", "pfcand_dphidphi", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dxydz", "pfcand_dphidxy", "pfcand_dlambdadz", "pfcand_dxyc", "pfcand_dxyctgtheta", "pfcand_phic", "pfcand_phidz", "pfcand_phictgtheta", "pfcand_cdz", "pfcand_cctgtheta", "pfcand_mtof", "pfcand_dndx", "pfcand_charge", "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad", "pfcand_dxy", "pfcand_dz", "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal", "pfcand_btagSip3dSig", "pfcand_btagJetDistVal", "pfcand_btagJetDistSig", "pfcand_type",
        ]

        # save some jet related variables
        branchList += [
            "jet_nmu",
            "jet_nel",
            "jet_nchad", 
            "jet_nnhad",
        ]

        ## extras

        branchList += [
            "pfcand_e",
            "pfcand_p",
        ]

        ##  outputs jet tag properties

        flavors = ["U", "D", "S", "B", "C", "G", "TAU"]
        for flav in flavors:
            branchList += [
                f"recojet_is{flav}",
                f"score_recojet_is{flav}",
            ]

        print("RDataFrame is calling following columns:", branchList)

        return branchList