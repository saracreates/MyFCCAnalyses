def generate_datacard(process_paths, output_file):
    with open(output_file, "w") as f:
        # Header
        f.write("imax *\njmax *\nkmax *\n---------------\n")
        
        # Shapes
        for process, path in process_paths.items():
            # f.write(f"shapes {process} * {path} recoil_mass_LHF\n")
            f.write(f"shapes {process} * {path} gamma_recoil_m_signal_cut\n")
        
        f.write("---------------\n---------------\n")
        f.write("#bin            bin1\nobservation     -1\n------------------------------\n")
        
        # Exclude data_obs from bin, process, and rate lines
        processes = [p for p in process_paths.keys() if p != "data_obs"]
        
        # Bin and process definitions
        f.write("bin          " + " ".join(["bin1"] * len(processes)) + "\n")
        f.write("process      " + " ".join(processes) + "\n")
        f.write("process      " + " ".join(map(str, range(len(processes)))) + "\n")
        f.write("rate         " + " ".join(["-1"] * len(processes)) + "\n")
        
        f.write("--------------------------------\n")
        f.write("#bkg lnU      -              1.5\n")
        f.write("#HWW_norm rateParam bin1 WW 1\n")
        f.write("#ZZ_norm rateParam bin1 ZZ 1\n")
        f.write("#Zqq_norm rateParam bin1 Zqq 1\n\n")
        
        # Systematic uncertainties (excluding data_obs and the first process, which is the signal)
        for i, process in enumerate(processes[1:]):  # Skip the first process (signal)
            # f.write(f"{process}_norm lnN " + "- " * (i + 1) + "1.05 " + "- " * (len(processes) - i - 2) + "\n")
            f.write(f"{process}_norm lnN " + "- " * (i + 1) + "1.01 " + "- " * (len(processes) - i - 2) + "\n")



# Change path here
data_path = "/afs/cern.ch/work/s/saaumill/public/FCCAnalyses/outputs/histmaker_LHF/ZZZqqllvv/"
# data_path = "/afs/cern.ch/work/s/saaumill/public/FCCAnalyses/outputs/histmaker_LHF/ZZZqqvvll/"

process_paths = {
    "HZZ": data_path + "wzp6_ee_qqH_HZZ_llvv_ecm240.root", # signal!
    "WW": data_path + "p8_ee_WW_ecm240.root",
    "ZZ": data_path + "p8_ee_ZZ_ecm240.root",
    # "nunuHZZ": data_path + "wzp6_ee_nunuH_HZZ_ecm240.root",
    "eeHZZ": data_path + "wzp6_ee_eeH_HZZ_ecm240.root",
    "mumuHZZ": data_path + "wzp6_ee_mumuH_HZZ_ecm240.root",
    # "qqHbb": data_path + "wzp6_ee_qqH_Hbb_ecm240.root",
    # "ssHbb": data_path + "wzp6_ee_ssH_Hbb_ecm240.root",
    # "ccHbb": data_path + "wzp6_ee_ccH_Hbb_ecm240.root",
    # "bbHbb": data_path + "wzp6_ee_bbH_Hbb_ecm240.root",
    "qqHtautau": data_path + "wzp6_ee_qqH_Htautau_ecm240.root",
    "ssHtautau": data_path + "wzp6_ee_ssH_Htautau_ecm240.root",
    "ccHtautau": data_path + "wzp6_ee_ccH_Htautau_ecm240.root",
    "bbHtautau": data_path + "wzp6_ee_bbH_Htautau_ecm240.root",
    "qqHWW": data_path + "wzp6_ee_qqH_HWW_ecm240.root",
    "ssHWW": data_path + "wzp6_ee_ssH_HWW_ecm240.root",
    "ccHWW": data_path + "wzp6_ee_ccH_HWW_ecm240.root",
    "bbHWW": data_path + "wzp6_ee_bbH_HWW_ecm240.root",
    "Zqq": data_path + "p8_ee_Zqq_ecm240.root",
    "data_obs": data_path + "wzp6_ee_qqH_HZZ_llvv_ecm240.root"
}

# Generate the datacard
# generate_datacard(process_paths, "datacard_qqllvv_recoil_m_jj.txt")

# Hgamma study

data_path = "/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/outputs/histmaker_fullsim/ZHgamma_btag/"
process_path = {
    "Hgamma": data_path + "p8_ee_Hgamma_ecm240.root", # signal!
    "ZH": data_path + "p8_ee_ZH_ecm240.root",
    "eegamma": data_path + "p8_ee_eegamma_ecm240.root",
    "mumugamma": data_path + "p8_ee_mumugamma_ecm240.root",
    "tautaugamma": data_path + "p8_ee_tautaugamma_ecm240.root",
    "qqgamma": data_path + "p8_ee_qqgamma_ecm240.root",
    "ccgamma": data_path + "p8_ee_ccgamma_ecm240.root",
    "bbgamma": data_path + "p8_ee_bbgamma_ecm240.root",
    "data_obs": data_path + "p8_ee_Hgamma_ecm240.root"
}
generate_datacard(process_path, "datacard_Hgamma.txt")