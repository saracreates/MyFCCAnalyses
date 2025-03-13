import os, copy

# list of processes (mandatory)
processList = {
    # cross sections given on the webpage: https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/ 
    'wzp6_ee_qqH_HZZ_ecm240':    {'fraction':1}, # 0.001409 pb -> 15200 events
    'wzp6_ee_qqH_HWW_ecm240':   {'fraction':1}, # 0.01148 pb  -> 186000 events
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
    df = df.Define(
        "muons_all",
        "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)",
    )
    # select leptons with momentum > 20 GeV
    df = df.Define(
        "muons",
        "FCCAnalyses::ReconstructedParticle::sel_p(20)(muons_all)", # 20 GeV is reasonable, I've checked on MC level
    )
    df = df.Define(
        "electrons_all",
        "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)",
    )
    df = df.Define(
        "electrons",
        "FCCAnalyses::ReconstructedParticle::sel_p(20)(electrons_all)",
    )

    # compute the muon isolation and store muons with an isolation cut of 0df = df.25 in a separate column muons_sel_iso
    df = df.Define(
        "muons_iso",
        "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(muons, ReconstructedParticles)",
    )
    df = df.Define(
        "muons_sel_iso",
        "FCCAnalyses::ZHfunctions::sel_iso(0.25)(muons, muons_iso)",
    )
    df = df.Define(
        "muons_sel_q",
        "FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)",
    )
    df = df.Define(
        "electrons_iso",
        "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(electrons, ReconstructedParticles)",
    )
    df = df.Define(
        "electrons_sel_iso",
        "FCCAnalyses::ZHfunctions::sel_iso(0.25)(electrons, electrons_iso)",
    )
    df = df.Define(
        "electrons_sel_q",
        "FCCAnalyses::ReconstructedParticle::get_charge(electrons_sel_iso)",
    )

    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))
    

    #########
    ### CUT 1: exactely two isolated leptons of opposite charge
    #########
    df = df.Filter("(muons_sel_iso.size() == 2 && Sum(muons_sel_q) == 0) || (electrons_sel_iso.size() == 2 && Sum(electrons_sel_q) == 0)")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))


    ## here cluster jets in the events but first remove muons from the list of
    ## reconstructed particles

    ## create a new collection of reconstructed particles removing muons with p>20
    df = df.Define(
        "ReconstructedParticlesNoMuons",
        "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles,muons)",
    )
    df = df.Define(
        "ReconstructedParticlesNoLeptons",
        "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticlesNoMuons,electrons)",
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
    collections_noleptons["PFParticles"] = "ReconstructedParticlesNoMuons"

    jetClusteringHelper = ExclusiveJetClusteringHelper(
        collections_noleptons["PFParticles"], 2
    )
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


    #########
    ### CUT 2: Njets = 2
    #########
    df = df.Filter("event_njet == 2")
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    ### gives 800 events, I expect 644 signal events - so need to check that qs come not from Z from higgs -> use p = 53 GeV!!

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

    """
    THE PLAN:
    Find events with
    - two jets with p around 53 GeV
    - two leptons with inv m ~ Z mass
    - missing energy
    """


    return results, weightsum