import os, copy

input_base = "/afs/cern.ch/work/s/saaumill/public/analyses/Hgamma_fullsim_simlink/exclusive_bb"
# list of processes
processList = {
    'p8_ee_Hgamma_ecm240':    {'fraction':1, 'crossSection': 8.20481e-05, 'inputDir': input_base},  #what are the exact values here?
    # #'reco_higgsgamma_test_REC.edm4hep': {'fraction':1, 'crossSection': 8.20481e-05, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/tmp_fullsim_output/reco_higgsgamma"},
    # 'p8_ee_qqgamma_ecm240':    {'fraction':1, 'crossSection': 6.9, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_ccgamma_ecm240':    {'fraction':1, 'crossSection': 2.15, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_bbgamma_ecm240':    {'fraction':1, 'crossSection': 2.35, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_ZH_ecm240':              {'fraction':1, 'crossSection': 0.2, 'inputDir': input_base},  #what are the exact values here?
    # # 'p8_ee_WW_ecm240':    {'fraction':1},  
    # # 'p8_ee_ZZ_ecm240':    {'fraction':1}, 
    # 'p8_ee_eegamma_ecm240':    {'fraction':1, 'crossSection': 190, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_tautaugamma_ecm240':    {'fraction':1, 'crossSection': 0.77, 'inputDir': input_base},  #what are the exact values here?
    # 'p8_ee_mumugamma_ecm240':    {'fraction':1, 'crossSection': 0.8, 'inputDir': input_base},  #what are the exact values here?
}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma_test/"

# Define the input dir (optional)
#inputDir    = "./localSamples/"

# additional/costom C++ functions, defined in header files (optional)
includePaths = ["./../functions.h"]

## latest particle transformer model, trained on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc"

## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)


## model files locally stored on /eos
model_dir = (
    "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"
)
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

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

        df = df.Define("photons_all", "FCCAnalyses::sel_type(22, ReconstructedParticles)")
        df = df.Define("electrons_all", f"FCCAnalyses::sel_type(11, ReconstructedParticles)")
        
        df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)") 
        df = df.Define("photons_n","FCCAnalyses::ReconstructedParticle::get_n(photons_all)")  #number of photons per event
        df = df.Define("photons_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_all))")

    
    
        #isolation cut
        df = df.Define("photons_iso", "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(photons_all, ReconstructedParticles)")  # is this correct?
        df = df.Define("photons_sel_iso","FCCAnalyses::ZHfunctions::sel_iso(0.2)(photons_all, photons_iso)",) # and this??
   
        df = df.Define("photons_iso_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_sel_iso)") 
        df = df.Define("photons_iso_n","FCCAnalyses::ReconstructedParticle::get_n(photons_sel_iso)")  #number of photons per event
        df = df.Define("photons_iso_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_sel_iso))")

        

        #########
        ### CUT 1: Photons must be isolated
        #########
    
        df = df.Filter("photons_sel_iso.size()>0 ")
      
        #sort in p  and select highest energetic one
        df = df.Define("iso_highest_p","FCCAnalyses::ZHfunctions::sort_by_energy(photons_sel_iso)")

        #energy cut
        df = df.Define("photons_boosted", "FCCAnalyses::ReconstructedParticle::sel_p(60,100)(iso_highest_p)") # looked okay from photons all
   
        df = df.Define("photons_boosted_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_boosted)") # is this correct?
        df = df.Define("photons_boosted_n","FCCAnalyses::ReconstructedParticle::get_n(photons_boosted)") 
        df = df.Define("photons_boosted_cos_theta","cos(FCCAnalyses::ReconstructedParticle::get_theta(photons_boosted))")
       
        #########
        ### CUT 2: Photons energy > 50
        #########
    
        df = df.Filter("photons_boosted.size()>0 ") 
        
        #########
        ### CUT 3: Cos Theta cut
        #########
        df = df.Filter("ROOT::VecOps::All(abs(photons_boosted_cos_theta) < 0.9) ") 
      
        ## create a new collection of reconstructed particles removing targeted photons
        df = df.Define("recopart_no_gamma", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons_boosted)",)
        df = df.Define("recopart_no_gamma_n","FCCAnalyses::ReconstructedParticle::get_n(recopart_no_gamma)") 

        df = df.Define("gamma_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(photons_boosted)") 
        df = df.Define("gamma_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(gamma_recoil)[0]") # recoil mass
        
       
        #########
        ### CUT 4: require at least 6 reconstructed particles (except gamma)
        #########
        df = df.Filter(" recopart_no_gamma_n > 5") 
        
      
        
        
        ## perform N=2 jet clustering
        global jetClusteringHelper
        global jetFlavourHelper

        ## define jet and run clustering parameters
        ## name of collections in EDM root files
        

        # fast sim 
        '''
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
        '''
        # full sim
        collections = {
            "GenParticles": "Particle",
            "PFParticles": "ReconstructedParticles",
            "PFTracks": "EFlowTrack",
            "PFPhotons": "EFlowPhoton",
            "PFNeutralHadrons": "EFlowNeutralHadron",
            "TrackState": "SiTracks_Refitted",

            # SiTracks: input_line_193:2:230: error: no viable conversion from 'RVec<edm4hep::TrackData>' to 'const RVec<edm4hep::TrackState>'

            "TrackerHits": "TrackerHits",
            "CalorimeterHits": "CalorimeterHits",
            "dNdx": "EFlowTrack_2",
            "PathLength": "EFlowTrack_L",
            "Bz": "magFieldBz",
        }

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

      #  df = df.Define(
      #      "btag",
      #      "JetTaggingUtils::get_btag({},90,90,90,90)".format(
      #          jetClusteringHelper.jets
      #      ),)
       
       
       
        return df

        #how do I get flavor

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            "ReconstructedParticles",
            "photons_all",
            "jj_m", #also require Higgs mass?
            #"photons_sel_iso",
            #"photons_boosted",
            #"recopart_no_gamma_n",
            #"photons_p"
            #"y23",
        ]

        ##  outputs jet properties
        # branchList += jetClusteringHelper.outputBranches()

        ## outputs jet scores and constituent breakdown
        branchList += jetFlavourHelper.outputBranches()

        print(branchList)

        return branchList