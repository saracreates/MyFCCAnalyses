
import uproot
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
import ROOT
import pickle


ROOT.gROOT.SetBatch(True)
# e.g. https://root.cern/doc/master/tmva101__Training_8py.html

def load_process(fIn, variables, target=0, weight_sf=1., entry_stop=1000):

    f = uproot.open(fIn)
    tree = f["events"]
    #meta = f["meta"]
    #weight = meta.values()[2]/meta.values()[1]*weight_sf

    df = tree.arrays(variables, library="pd", entry_start=0, entry_stop=entry_stop) # limit num of events for equal datasets

    num_events = len(df)

    weight = 1.0/num_events*weight_sf

    df['target'] = target # add a target column to indicate signal (1) and background (0)
    df['weight'] = weight

    print("Load {} with {} events and weight {}".format(fIn.replace(".root", ""), tree.num_entries, weight))
    return df



print("Parse inputs")

# configuration of signal, background, variables, files, ...
variables = ["l1_p", "l2_p", "l1_theta", "l2_theta", "m_ll", "m_recoil_ll", "y23", "y34", "jet1_nconst_N2", "jet2_nconst_N2", "m_jj", "p_res_jj", "recoil_mass", "miss_p", "miss_e", "miss_pz", "miss_theta", "miss_pT", "recoil_mass_jjll", "Zll_costheta", "Zll_p", "Zll_pT", "Zjj_costheta", "Zjj_p", "Zjj_pT", "dot_prod_had", "dot_prod_lep", "dot_prod_ll"]
weight_sf = 1e9
sig_df = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/wzp6_ee_qqH_HZZ_llvv_ecm240.root", variables, weight_sf=weight_sf, target=0, entry_stop=2000)
bkg_1 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_qqH_HWW_ecm240.root", variables, weight_sf=weight_sf, target=4) # not in training
bkg_2 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_ssH_HWW_ecm240.root", variables, weight_sf=weight_sf, target=4) #not
bkg_3 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_ccH_HWW_ecm240.root", variables, weight_sf=weight_sf, target=4) #not
bkg_4 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_bbH_HWW_ecm240.root", variables, weight_sf=weight_sf, target=4) #not
bkg_5 = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/p8_ee_ZZ_ecm240.root", variables, weight_sf=weight_sf, target=3) #yes
# bkg_6 = load_process("outputs/mva/ZZZ_qqllvv/preselection/p8_ee_WW_ecm240.root", variables, weight_sf=weight_sf) #yes 
# bkg_7 = load_process("outputs/mva/ZZZ_qqllvv/preselection/p8_ee_Zqq_ecm240.root", variables, weight_sf=weight_sf) #yes
# bkg_8 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_qqH_Hbb_ecm240.root", variables, weight_sf=weight_sf) #no
# bkg_9 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_ssH_Hbb_ecm240.root", variables, weight_sf=weight_sf) #no
# bkg_10 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_ccH_Hbb_ecm240.root", variables, weight_sf=weight_sf) #no
# bkg_11 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_bbH_Hbb_ecm240.root", variables, weight_sf=weight_sf) #no
bkg_12 = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/wzp6_ee_qqH_Htautau_ecm240.root", variables, weight_sf=weight_sf, target=1) # no
bkg_13 = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/wzp6_ee_ssH_Htautau_ecm240.root", variables, weight_sf=weight_sf, target=1) # no
bkg_14 = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/wzp6_ee_ccH_Htautau_ecm240.root", variables, weight_sf=weight_sf, target=1) # no
bkg_15 = load_process("outputs/mva_multi/ZZZ_qqllvv/preselection/wzp6_ee_bbH_Htautau_ecm240.root", variables, weight_sf=weight_sf, target=1) # no
# bkg_16 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_eeH_HZZ_ecm240.root", variables, weight_sf=weight_sf, target=2, entry_stop=350) #no
# bkg_17 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_mumuH_HZZ_ecm240.root", variables, weight_sf=weight_sf, target=2, entry_stop=350) #no
bkg_18 = load_process("outputs/mva/ZZZ_qqllvv/preselection/wzp6_ee_nunuH_HZZ_ecm240.root", variables, weight_sf=weight_sf, target=2, entry_stop=700) # yes

# NOTE: does is make sense to train on such unequal number of signal and background events?
# NOTE: how do I choose the best weight?




# Concatenate the dataframes into a single dataframe
data = pd.concat([sig_df, bkg_1, bkg_2, bkg_3, bkg_4, bkg_5, bkg_12, bkg_13, bkg_14, bkg_15, bkg_18], ignore_index=True)


# split data in train/test events
train_data, test_data, train_labels, test_labels, train_weights, test_weights  = train_test_split(
    data[variables], data['target'], data['weight'], test_size=0.2, random_state=42
)



# conversion to numpy needed to have default feature_names (fN), needed for conversion to TMVA
train_data = train_data.to_numpy()
test_data = test_data.to_numpy()
train_labels = train_labels.to_numpy()
test_labels = test_labels.to_numpy()
train_weights = train_weights.to_numpy()
test_weights = test_weights.to_numpy()

# set hyperparameters for the XGBoost model
params = {
    'objective': 'multi:softproba',
    'num_classes': 5,
    'eval_metric': ["merror", "mlogloss"],
    'eta': 0.1,
    'max_depth': 5,
    'subsample': 0.5,
    'colsample_bytree': 0.5,
    'seed': 42,
    'n_estimators': 350, # low number for testing purposes (default 350) ("large: 1000")
    'early_stopping_rounds': 10,
    'num_rounds': 50,
    'learning_rate': 0.20,
    'gamma': 5,
    'min_child_weight': 10,
    'max_delta_step': 0,
    'reg_lambda': 10,  # L2 regularization
    'reg_alpha': 5, 
}


# train the XGBoost model
print("Start training")
eval_set = [(train_data, train_labels), (test_data, test_labels)]
bdt = xgb.XGBClassifier(**params)
bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight=train_weights)


# export model (to ROOT and pkl)
print("Export model")
fOutName = "outputs/mva_multi/ZZZ_qqvvll/bdt_model_multi_all_bkg.root"
ROOT.TMVA.Experimental.SaveXGBoost(bdt, "bdt_model", fOutName, num_inputs=len(variables))

# append the variables
variables_ = ROOT.TList()
for var in variables:
     variables_.Add(ROOT.TObjString(var))
fOut = ROOT.TFile(fOutName, "UPDATE")
fOut.WriteObject(variables_, "variables")


save = {}
save['model'] = bdt
save['train_data'] = train_data
save['test_data'] = test_data
save['train_labels'] = train_labels
save['test_labels'] = test_labels
save['variables'] = variables
pickle.dump(save, open("outputs/mva_multi/ZZZ_qqvvll/bdt_model_multi_all_bkg.pkl", "wb"))