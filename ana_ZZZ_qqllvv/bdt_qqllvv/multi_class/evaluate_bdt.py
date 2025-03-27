import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
import pickle
import argparse
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize

parser = argparse.ArgumentParser()
# parser.add_argument("-i", "--input", type=str, default="outputs/mva_multi/ZZZ_qqllvv/bdt_model_multi.pkl", help="Input pkl file")
# parser.add_argument("-o", "--outDir", type=str, default="outputs/mva_multi/ZZZ_qqllvv/plots_training", help="Output directory")
# tick_names = ["ZZ", "Htautau", "signal"] 

# parser.add_argument("-i", "--input", type=str, default="outputs/mva_multi/ZZZ_qqllvv/bdt_model_multi_all_bkg.pkl", help="Input pkl file")
# parser.add_argument("-o", "--outDir", type=str, default="outputs/mva_multi/ZZZ_qqllvv/plots_training_all_bkg", help="Output directory")
# tick_names = ["signal", "Htautau", "llH;HZZ", "ZZ"]

parser.add_argument("-i", "--input", type=str, default="outputs/mva_multi/ZZZ_qqllvv/bdt_model_multi_all_bkg_large.pkl", help="Input pkl file")
parser.add_argument("-o", "--outDir", type=str, default="outputs/mva_multi/ZZZ_qqllvv/plots_training_all_bkg_large", help="Output directory")
tick_names = ["signal", "Htautau", "llH;HZZ", "ZZ"]

args = parser.parse_args()


def plot_roc():

    train_probs = bdt.predict_proba(train_data)
    train_preds = train_probs[:,1]
    train_fpr, train_tpr, threshold = sklearn.metrics.roc_curve(train_labels, train_preds)
    train_roc_auc = sklearn.metrics.auc(train_fpr, train_tpr)

    test_probs = bdt.predict_proba(test_data)
    test_preds = test_probs[:,1]
    test_fpr, test_tpr, threshold = sklearn.metrics.roc_curve(test_labels, test_preds)
    test_roc_auc = sklearn.metrics.auc(test_fpr, test_tpr)

    # Plot the ROC curve
    plt.figure(figsize=(8, 6))
    plt.plot(train_fpr, train_tpr, color='blue', label=f"Training ROC (AUC = {train_roc_auc:.2f})")
    plt.plot(test_fpr, test_tpr, color='red', label=f"Testing ROC (AUC = {test_roc_auc:.2f})")
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random Guess')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.grid()
    plt.savefig(f"{outDir}/roc.png")
    plt.savefig(f"{outDir}/roc.pdf")
    plt.close()

def plot_roc_multi():
    # Define class names manually
    tick_names = ["signal", "Htautau", "background"]  # Modify as needed
    n_classes = len(tick_names)

    # Binarize the labels (convert to one-hot encoding)
    y_true = label_binarize(test_labels, classes=np.arange(n_classes))
    y_pred = bdt.predict_proba(test_data)

    plt.figure(figsize=(8, 6))

    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_true[:, i], y_pred[:, i])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"{tick_names[i]} (AUC = {roc_auc:.2f})")

    # Random guess line
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random Guess')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Multiclass ROC Curve')
    plt.legend()
    plt.grid()
    plt.savefig(f"{outDir}/roc.png")
    plt.savefig(f"{outDir}/roc.pdf")
    plt.close()


def plot_score():

    train_predictions = bdt.predict_proba(train_data)[:,1]
    test_predictions = bdt.predict_proba(test_data)[:,1]

    # Separate the data into signal and background samples
    train_signal_scores = train_predictions[train_labels == 1]
    train_background_scores = train_predictions[train_labels == 0]
    test_signal_scores = test_predictions[test_labels == 1]
    test_background_scores = test_predictions[test_labels == 0]

    # Plot the BDT scores for signal and background events
    plt.figure(figsize=(8, 6))
    plt.hist(train_signal_scores, bins=50, range=(0, 1), histtype='step', label='Training Signal', color='blue', density=True)
    plt.hist(train_background_scores, bins=50, range=(0, 1), histtype='step', label='Training Background', color='red', density=True)
    plt.hist(test_signal_scores, bins=50, range=(0, 1), histtype='step', label='Testing Signal', color='blue', linestyle='dashed', density=True)
    plt.hist(test_background_scores, bins=50, range=(0, 1), histtype='step', label='Testing Background', color='red', linestyle='dashed', density=True)
    plt.xlabel('BDT Score')
    plt.ylabel('Number of Events (normalized)')
    plt.title('BDT Score Distribution')
    plt.legend()
    plt.grid()
    plt.savefig(f"{outDir}/score.png")
    plt.savefig(f"{outDir}/score.pdf")
    plt.close()

def plot_importance():

    fig, ax = plt.subplots(figsize=(12, 6))

    importance = bdt.get_booster().get_score(importance_type='weight')
    sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=False)
    sorted_indices = [int(x[0][1:]) for x in sorted_importance] # sorted indices

    # Get the sorted variable names and their corresponding importances
    sorted_vars = [variables[i] for i in sorted_indices]
    sorted_values = [x[1] for x in sorted_importance]

    # Create a DataFrame and plot the feature importances
    importance_df = pd.DataFrame({'Variable': sorted_vars, 'Importance': sorted_values})
    importance_df.plot(kind='barh', x='Variable', y='Importance', legend=None, ax=ax)
    ax.set_xlabel('BDT score')
    ax.set_title("BDT variable scores", fontsize=16)
    plt.savefig(f"{outDir}/importance.png")
    plt.savefig(f"{outDir}/importance.pdf")
    plt.close()

def plot_confusion_matrix():
    y_true = test_labels
    y_pred = np.argmax(bdt.predict_proba(test_data), axis=1)
    
    cm = confusion_matrix(y_true, y_pred)
    cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]  # Normalize per class
    
    plt.figure(figsize=(8, 6))
    plt.imshow(cm_normalized, interpolation='nearest', cmap=plt.cm.Blues)
    plt.colorbar()

    # Define class names manually
    # tick_names = ["signal", "Htautau", "background"]  # Modify as needed
    plt.xticks(np.arange(len(tick_names)), tick_names, rotation=45, ha='right')
    plt.yticks(np.arange(len(tick_names)), tick_names)
    
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, f"{cm_normalized[i, j]:.2f}", ha='center', va='center', color='white' if cm_normalized[i, j] > 0.5 else 'black')
    
    plt.xlabel("Predicted Label")
    plt.ylabel("True Label")
    plt.title("Confusion Matrix")
    plt.savefig(f"{outDir}/confusion_matrix.png")
    plt.savefig(f"{outDir}/confusion_matrix.pdf")
    plt.close()

if __name__ == "__main__":
    outDir = args.outDir

    res = pickle.load(open(args.input, "rb"))
    bdt = res['model']
    train_data = res['train_data']
    test_data = res['test_data']
    train_labels = res['train_labels']
    test_labels = res['test_labels']
    variables = res['variables']

    plot_score()
    # plot_roc_multi()
    plot_importance()
    plot_confusion_matrix()