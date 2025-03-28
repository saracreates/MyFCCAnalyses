
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB
import os, copy
# import argparse

# # Create the argument parser
# parser = argparse.ArgumentParser(description="Run the analysis with configurable doInference parameter")
# parser.add_argument(
#     '--doInference',
#     type=bool,
#     default=False,  # Default value
#     help="Set to True to perform inference, False to skip it"
# )

# # Parse the arguments
# args = parser.parse_args()

# # Now use args.doInference to set the parameter
# doInference = args.doInference

doInference = False
# Output directory
outputDir   = f"outputs/mva_multi/ZZZ_qqvvll/preselection/"


processList = {
    # cross sections given on the webpage: https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/ 
    'wzp6_ee_qqH_HZZ_llvv_ecm240': {'fraction':0.0008, 'crossSection': 0.00015, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"}, # 3000 events .....
    # # load two files for training
    'wzp6_ee_qqH_HWW_ecm240':   {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # q = u, d
    'wzp6_ee_ssH_HWW_ecm240':   {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # s
    'wzp6_ee_ccH_HWW_ecm240':   {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # c
    'wzp6_ee_bbH_HWW_ecm240':   {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # b
    'p8_ee_ZZ_ecm240':          {'fraction':1, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"}, # 680 events
    'p8_ee_WW_ecm240':          {'fraction':1, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"},
    # no Hbb left after selection cuts
    # 'wzp6_ee_qqH_Hbb_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # q = u, d
    # 'wzp6_ee_ssH_Hbb_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # s
    # 'wzp6_ee_ccH_Hbb_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # c
    # 'wzp6_ee_bbH_Hbb_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # b
    'wzp6_ee_qqH_Htautau_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # q = u, d
    'wzp6_ee_ssH_Htautau_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # s
    'wzp6_ee_ccH_Htautau_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # c
    'wzp6_ee_bbH_Htautau_ecm240':  {'fraction':1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"}, # b
    'p8_ee_Zqq_ecm240':         {'fraction':1, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"}, # q = u,d,s,c,b,t 
    # add other signal as bkg
    # 'wzp6_ee_eeH_HZZ_ecm240': {'fraction': 1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"},
    # 'wzp6_ee_mumuH_HZZ_ecm240': {'fraction': 1, 'inputDir': "/afs/cern.ch/work/s/saaumill/public/symlink_qqllvv_data"},
    'wzp6_ee_nunuH_HZZ_ecm240': {'fraction': 0.1, 'inputDir': "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"},
}


# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

# Additional/custom C++ functions, defined in header files
includePaths = ["./../functions.h"]


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

# Multithreading: -1 means using all cores
nCPUS       = -1

# Batch settings
#runBatch    = False
#batchQueue  = "longlunch"
#compGroup = "group_u_FCC.local_gen"


class RDFanalysis():

    # encapsulate analysis logic, definitions and filters in the dataframe
    def analysers(df):

        # define some aliases to be used later on
        df = df.Alias("Particle0", "Particle#0.index") # index of the daughter particles
        df = df.Alias("Particle1", "Particle#1.index") # index of the mother particles 
        df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
        df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Electron0", "Electron#0.index")


        # get all the leptons from the collection
        df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)",)
        df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)",)

        # compute the muon isolation and store muons with an isolation cut of 0df = df.25 in a separate column muons_sel_iso
        df = df.Define("muons_iso","FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(muons, ReconstructedParticles)",)
        df = df.Define("muons_sel_iso","FCCAnalyses::ZHfunctions::sel_iso(0.25)(muons, muons_iso)",)
        df = df.Define("muons_sel_q","FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)",)
        df = df.Define("electrons_iso","FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.1)(electrons, ReconstructedParticles)",)
        df = df.Define("electrons_sel_iso","FCCAnalyses::ZHfunctions::sel_iso(0.25)(electrons, electrons_iso)",)
        df = df.Define("electrons_sel_q","FCCAnalyses::ReconstructedParticle::get_charge(electrons_sel_iso)",)

        # VETO: two leptons that form pair
        df = df.Filter("((muons_sel_iso.size() == 2 && Sum(muons_sel_q) == 0) || (electrons_sel_iso.size() == 2 && Sum(electrons_sel_q) == 0)) && (muons_sel_iso.size() + electrons_sel_iso.size() == 2)")

        # get the two leptons from Z decay
        df = df.Define("l1", "muons_sel_iso.size() == 2 ? muons_sel_iso[0] : electrons_sel_iso[0]")
        df = df.Define("l2", "muons_sel_iso.size() == 2 ? muons_sel_iso[1] : electrons_sel_iso[1]")

        # save these params
        df = df.Define("l1_p", "FCCAnalyses::ReconstructedParticle::get_p(Vec_rp{l1})[0]")
        df = df.Define("l2_p", "FCCAnalyses::ReconstructedParticle::get_p(Vec_rp{l2})[0]")
        df = df.Define("l1_theta", "FCCAnalyses::ReconstructedParticle::get_theta(Vec_rp{l1})[0]")
        df = df.Define("l2_theta", "FCCAnalyses::ReconstructedParticle::get_theta(Vec_rp{l2})[0]")
        df = df.Define("p_lep_sorted", "FCCAnalyses::ZHfunctions::sort_by_momentum(l1_p, l2_p)")
        df = df.Define("l_p_max", "p_lep_sorted[0]")
        # isolation value?


        df = df.Define("res_ll", "FCCAnalyses::ZHfunctions::get_two_lep_res(l1, l2)")
        df = df.Define("m_ll", "FCCAnalyses::ReconstructedParticle::get_mass(res_ll)[0]")
        df = df.Define("recoil_ll", "FCCAnalyses::ZHfunctions::get_recoil_lep(240.0, l1, l2)")
        df = df.Define("m_recoil_ll", "FCCAnalyses::ReconstructedParticle::get_mass(recoil_ll)[0]")


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


        df = df.Define("jets_p4", "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets),)

        df = df.Define("m_jj", "JetConstituentsUtils::InvariantMass(jets_p4[0], jets_p4[1])",)

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


        # calculate invariant mass of llqq system which gives the offshell Z and apply prop cut on recoil mass ~ 10-50 GeV 
        df = df.Define("recoil_jjll", "FCCAnalyses::ZHfunctions::get_recoil_from_lep_and_jets(240.0, jet1, jet2, l1, l2)")
        df = df.Define("recoil_mass_jjll", "FCCAnalyses::ReconstructedParticle::get_mass(recoil_jjll)[0]")


        # look at ll system
        df = df.Define("Zll_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_ll)")
        df = df.Define("Zll_p", "FCCAnalyses::ReconstructedParticle::get_p(res_ll)[0]")
        df = df.Define("Zll_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_ll)[0]")

        # look at jj system
        df = df.Define("Zjj_costheta", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(res_jj)")
        df = df.Define("Zjj_p", "FCCAnalyses::ReconstructedParticle::get_p(res_jj)[0]")
        df = df.Define("Zjj_pT", "FCCAnalyses::ReconstructedParticle::get_pt(res_jj)[0]")


        df = df.Define("dot_prod_had", "FCCAnalyses::ZHfunctions::dot_prod_had(missP, jet1, jet2)")
        df = df.Define("dot_prod_lep", "FCCAnalyses::ZHfunctions::dot_prod_lep(missP, l1, l2)")
        df = df.Define("dot_prod_ll", "FCCAnalyses::ZHfunctions::dot_prod_ll(l1, l2)")

        # Filter for orthogonality! 

        # 0) recoil_jjll against my own qqllvv analysis
        df = df.Filter("recoil_mass_jjll > 60")

        # 1) recoil_jj against llvvqq, vvllqq 
        df = df.Filter("recoil_mass < 140")
        df = df.Filter("m_recoil_ll > 140")

        # 2) llqqvv is hard to distinguish from qqllvv, so I'll leave it as a background 


        # Filter - TO THE SAME CUTS AS IN THE CUT AND COUNT APPROACH
        # like in cut and count... 

        # 1) lepton pairs (already done)

        # 2) rought cut on recoil_jj
        df = df.Filter("recoil_mass > 100 && recoil_mass < 170")

        # 3) rough cut on recoill jjll
        df = df.Filter("recoil_mass_jjll > 80 && recoil_mass_jjll < 105") 

        # 4) rough cut on missing momentum
        df = df.Filter("miss_p > 20 && miss_p < 100")

        # 5) rough cut on m_jj
        df = df.Filter("m_jj > 85 && m_jj < 105")

        # 6) rough cut on p_jj
        df = df.Filter("p_res_jj > 40 && p_res_jj < 55")

        # 7) cut on m_ll
        df = df.Filter("m_ll > 10 && m_ll < 45")

        # 8)  cut on miss p_T
        df = df.Filter("miss_pT > 10 && miss_pT < 70")

        # 9) cos theta cut between jets and neutrinos
        df = df.Filter("dot_prod_had < -0.4")

        # 10) cos theta cut between leptons and neutrinos
        df = df.Filter("dot_prod_lep < 0.95")

        # 11) on max lepton momentum
        df = df.Filter("l_p_max > 10 && l_p_max < 40")

        # 12) cos theta cut between the two leptons
        df = df.Filter("dot_prod_ll > -0.75")

        # 13) tighter cut on recoil mass
        df = df.Filter("recoil_mass > 120 && recoil_mass < 140")

        return df

    # define output branches to be saved
    def output():
        branchList = ["l1_p", "l2_p", "l1_theta", "l2_theta", "m_ll", "m_recoil_ll", "y23", "y34", "jet1_nconst_N2", "jet2_nconst_N2", "m_jj", "p_res_jj", "recoil_mass", "miss_p", "miss_e", "miss_pz", "miss_theta", "miss_pT", "recoil_mass_jjll", "Zll_costheta", "Zll_p", "Zll_pT", "Zjj_costheta", "Zjj_p", "Zjj_pT", "dot_prod_had", "dot_prod_lep", "dot_prod_ll"]
        if doInference:
            branchList.append("mva_score")
        return branchList