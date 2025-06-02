import os
from multiprocessing import Process
import glob

# indir = "/eos/experiment/fcc/ee/generation/DelphesStandalone/Edm4Hep/pre_winter2023_tests_v2/"
# outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/pre_winter2023_tests_v2/selvaggi_2022Nov26/"

indir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
indir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/"
outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/weavercore_v1/mgarcia_2023Jan11"
# outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_test/"
outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_test/"

outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022_MT/"


outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_20_03_2022_MCtruth/"
outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_20_03_2022_MCtruth_sf2p0Inf/"
outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022_MT_sf2p0Inf/"

outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"
outdir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_365_02_09_2024/"

#script = "examples/FCCee/weaver/analysis_inference_gen.py"
script = "examples/FCCee/weaver/analysis_inference.py"


os.system("cd build; make install -j 20; cd ..;")


# __________________________________________________________________________________
def run(cmd):
    os.system(cmd)


procs = []

# for sample in ["wzp6_ee_nunuH", "p8_ee_ZH_Znunu"]:
for sample in ["wzp6_ee_nunuH"]:
    for f in ["Hgg", "Hss", "Hcc", "Hbb", "Huu", "Hdd", "Htautau"]:
    #for f in ["Huu", "Hdd", "Htautau"]:
    #for f in ["Hgg", "Hss", "Hcc", "Hbb", "Hqq"]:
    #for f in ["Hgg", "Hss", "Hcc", "Hbb", "Hqq"]:
        #procdir = "{}_{}_ecm240".format(sample, f)
        procdir = "{}_{}_ecm365".format(sample, f)

        #infile = glob.glob("{}/{}_{}_ecm240/events_*.root".format(indir, sample, f))[0]
        #infile = glob.glob("{}/{}_{}_ecm240/events_*.root".format(indir, sample, f))
        infile = "{}/{}_{}_ecm240/events_*.root".format(indir, sample, f)
        # for f in ["Huu", "Hdd"]:
        # cmd = "fccanalysis run examples/FCCee/weaver/analysis_inference.py "
        cmd = "fccanalysis run {} ".format(script)
        # cmd = "fccanalysis run examples/FCCee/weaver/analysis_inference_zh_vvjj.py "
        cmd += "--output {}/{}/events.root ".format(outdir, procdir)
        cmd += "--files-list {} ".format(infile)
        cmd += "--ncpus 64 "
        # cmd += "--nevents 1000"

        print(cmd)

        #os.system(cmd)
        # proc = Process(target=run, args=(cmd,))
        # procs.append(proc)
        # proc.start()

# for proc in procs:
#    proc.join()
