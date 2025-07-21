import sys
from array import array
from ROOT import TFile, TTree
from examples.FCCee.weaver.config import variables_pfcand, variables_jet, flavors

debug = True

if len(sys.argv) < 2:
    print(" Usage: stage2.py input_file output_file n_start n_events")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
n_start = int(sys.argv[3])
n_final = int(sys.argv[4])
n_events = n_final - n_start

# Opening the input file containing the tree (output of stage1.py)
infile = TFile.Open(input_file, "READ")

ev = infile.Get("events")
numberOfEntries = ev.GetEntries()
if debug:
    print("-> number of entries in input file: {}".format(numberOfEntries))

## basic checks
if n_final > n_start + numberOfEntries:
    print("ERROR: requesting too many events. This file only has {}".format(numberOfEntries))
    sys.exit()

branches_pfcand = ["pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", "pfcand_dptdpt", "pfcand_detadeta", "pfcand_dphidphi", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dxydz", "pfcand_dphidxy", "pfcand_dlambdadz", "pfcand_dxyc", "pfcand_dxyctgtheta", "pfcand_phic", "pfcand_phidz", "pfcand_phictgtheta", "pfcand_cdz", "pfcand_cctgtheta", "pfcand_mtof", "pfcand_dndx", "pfcand_charge", "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad", "pfcand_dxy", "pfcand_dz", "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal", "pfcand_btagSip3dSig", "pfcand_btagJetDistVal", "pfcand_btagJetDistSig", "pfcand_type", "pfcand_e", "pfcand_p"]

branches_jet = ["jet_nmu", "jet_nel", "jet_nchad", "jet_nnhad", 'jet_p_N2', 'jet_e_N2', 'jet_mass_N2', 'jet_phi_N2', 'jet_theta_N2', 'jet_nconst_N2'] # 'event_njet_N2'

if len(branches_pfcand) == 0:
    print("ERROR: branches_pfcand is empty ...")
    sys.exit()

if len(branches_jet) == 0:
    print("ERROR: branches_jet is empty ...")
    sys.exit()

# print("")
# print("-> number of events: {}".format(numberOfEntries))
# print("-> requested to run over [{},{}] event range".format(n_start, n_final))

# branches_pfcand = [branches_pfcand[0]]
# branches_jet = [branches_jet[-1]]

## define variables for output tree
maxn = 500


## output jet-wise tree
out_root = TFile(output_file, "RECREATE")
t = TTree("tree", "tree with jets")

jet_array = dict()

# number of reconstructed particles per event
n_reco_particles = array("i", [0])
t.Branch("n_reco_particles", n_reco_particles, "n_reco_particles/I")

flavors = ["U", "D", "S", "B", "C", "G", "TAU"]
for f in flavors:
    b1 = "recojet_is{}".format(f.upper())
    b2 = "score_recojet_is{}".format(f.upper())
    jet_array[b1] = array("i", [0])
    jet_array[b2] = array("f", [0.0])
    t.Branch(b1, jet_array[b1], "{}/I".format(b1))
    t.Branch(b2, jet_array[b2], "{}/F".format(b2))

    # jet clustering 
    b3 = "event_njet_N2"
    jet_array[b3] = array("i", [0])
    t.Branch(b3, jet_array[b3], "{}/I".format(b3))
for b in branches_jet:
    jet_array[b] = array("f", [0])
    t.Branch(b, jet_array[b], "{}/F".format(b))

## need this branch to define pfcand branches
jet_npfcand = array("i", [0])
t.Branch("jet_npfcand", jet_npfcand, "jet_npfcand/I")

pfcand_array = dict()
for b in branches_pfcand:
    pfcand_array[b] = array("f", maxn * [0])
    t.Branch(b, pfcand_array[b], "{}[jet_npfcand]/F".format(b))

if debug:
    for key, item in jet_array.items():
        print(key)
    for key, item in pfcand_array.items():
        print(key)

# Loop over all events
for entry in range(n_start, n_final):
    # Load selected branches with data from specified event

    # if (entry+1)%100 == 0:
    # if (entry + 1) % 1000 == 0:
    #    print(" ... processed {} events ...".format(entry + 1))

    ev.GetEntry(entry)

    njets = len(getattr(ev, branches_jet[0]))
    if debug:
        print("-> processing event {} with {} jets".format(entry, njets))

    ## loop over jets
    for j in range(njets):
        name = "event_njet_N2"
        jet_array[name][0] = getattr(ev, name)

        # number of reconstructed particles
        n_reco_particles[0] = getattr(ev, "n_reco_particles")
        if debug:
            print("-> processing event {} with {} reconstructed particles".format(entry, n_reco_particles[0]))
        

        ## fill jet-based quantities
        for f in flavors:
            name = "recojet_is{}".format(f.upper())
            jet_array[name][0] = getattr(ev, name)
            jet_array["score_" + name][0] = float(getattr(ev, "score_" + name)[j])
            if debug:
                print("   jet:", j, name, jet_array[name][0])
                print("   jet:", j, "score_" + name, jet_array["score_" + name][0])

        for name in branches_jet:
            jet_array[name][0] = getattr(ev, name)[j]
            if debug:
                print("   jet:", j, name, getattr(ev, name)[j])

        ## loop over constituents
        jet_npfcand[0] = len(getattr(ev, branches_pfcand[0])[j])
        if debug:
            print("   jet:", j, "jet_npfcand", jet_npfcand[0])
        for k in range(jet_npfcand[0]):
            for name in branches_pfcand:
                pfcand_array[name][k] = getattr(ev, name)[j][k]
                if name=="pfcand_p":
                    if debug:
                        print("       const:", k, name, getattr(ev, name)[j][k])

        ## fill tree at every jet
        t.Fill()

# write tree
t.SetDirectory(out_root)
t.Write()
