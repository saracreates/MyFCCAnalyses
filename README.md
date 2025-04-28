# MyFCCAnalyses

Analyses I conduct at FCC-ee.

## $e^+e^- \rightarrow H(ZZ)Z \rightarrow$ 6 leptons

Analysis in **fast simultion**. This analysis contributes to the overall precision of the $ZH \rightarrow ZZZ$ cross-section. The overall precison achievable in this channel is 30~\%. To reproduce this result:

1. Clone the `pre-edm4hep1` branch from FCCAnalyses, see intructions [here](https://github.com/HEP-FCC/FCCAnalyses/blob/master/README.md#winter-2023-and-spring-2021-pre-generated-samples)
2. Set up the envirnoment with `source ./setup.sh && fccanalysis build -j 8`

Then run the selection: 

```
fccanalysis run ./ana_ZZZ_6l/histmaker.py
```

and create the plots (with the custom `do_plots.py` file in the `extras` folder): 

```
export PYTHONPATH=/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/extras:$PYTHONPATH
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:$PATH
fccanalysis plots ./ana_ZZZ_6l/plots.py 
```

Find the cutflow table in `cutFlow_cutflow.pdf` and the money-plot in `recoil_mass_cut4.pdf`. 

*Note:* Some of the data used was produced privately with `MadGraph`. For instructions on how to create the data, see my [mg_config repo](https://github.com/saracreates/mg_configs).

## $e^+e^- \rightarrow H(ZZ)Z \rightarrow qq \ell \ell \nu \nu$ and $qq \nu \nu \ell \ell$

Analysis in **fast simultion**. This analysis contributes to the overall precision of the $ZH \rightarrow ZZZ$ cross-section. The overall precison achievable in this channel is 13.5~\% for $Z_{\mathrm{prod}(qq)Z(\ell \ell)Z^*(\nu \nu)}$ and 29.7~\% for $Z_{\mathrm{prod}(qq)Z(\nu \nu)Z^*(\ell \ell)}$. To reproduce this result:

1. Clone the `pre-edm4hep1` branch from FCCAnalyses, see intructions [here](https://github.com/HEP-FCC/FCCAnalyses/blob/master/README.md#winter-2023-and-spring-2021-pre-generated-samples)
2. Set up the envirnoment with `source ./setup.sh && fccanalysis build -j 8`

Then run the selection: 

```
fccanalysis run ./ana_ZZZ_qqllvv/histmaker_LHF.py
```

The histograms `recoil_mass_LHF;1` in the output files are used for the likelihood fit with [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/). See [this section](#perform-a-likelihood-fit-with-combine) for instructions. 

Let's create the plots with: 

```
export PYTHONPATH=/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/extras:$PYTHONPATH
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:$PATH
fccanalysis plots ./ana_ZZZ_qqllvv/plots_LHF.py 
```

Again, the cutflow table can be found in `cutFlow_cutflow.pdf` and the money-plot in `recoil_mass_LHF.pdf`. 


*Note:* For the $qq \nu \nu \ell \ell$ final state do exactly the same, the files are named after the processes. 

**BDT setup**: In this code, you can also find a setup for using BDTs although they were not used for the final results presented. Here is an example on how to use the code: 

```
cd /path/to/pre-edm4hep1/FCCAnalyses
source ./setup.sh

# create training data (set doInference = False)
fccanalysis run ./ana_ZZZ_qqllvv/bdt_qqllvv/multi_class/preselection.py
fccanalysis run ./ana_ZZZ_qqvvll/bdt_qqvvll/preselection.py

# train bdt
python3 ./ana_ZZZ_qqllvv/bdt_qqllvv/multi_class/train_bdt_all_bkg.py 
python3 ./ana_ZZZ_qqvvll/bdt_qqvvll/train_bdt_all_bkg.py

# evaluate the bdt
python3 ./ana_ZZZ_qqvvll/bdt_qqvvll/evaluate_bdt.py

# then load bdt in histmaker_LHF (set doInference = True) & do everything as usual
```


## $e^+e^- \rightarrow H \gamma $

Analysis in **full simulation**. This probes the effektive $ZH\gamma$ coupling in the $H\gamma$ production mode. Lena Herrmann already did this analysis in fast simulation, see [her repo](https://github.com/herrmannlena/FCCAnalyses/blob/higgsgamma/myanalysis/plots_recoil.py). Here, we are reproducing the results in full simulation! Also to test the [new ML-based jet-flavor tagger](https://github.com/key4hep/k4MLJetTagger) in key4hep ;)

We need the **master branch** of FCCAnalyses, so simply source the nightlies `source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh` (or later, the latest stable release). If you want to set-up FCCAnalyses locally: Clone the master branch and `source ./setup.sh && fccanalysis build -j 8`. 

The background data is producted with whizard and is available in the [central production](https://fcc-physics-events.web.cern.ch/fcc-ee/full-sim/index.php). The signal sample is produced privately with MadGraph. For instructions on how to produce the data, check out my [mg_config repo](https://github.com/saracreates/mg_configs). 

Run the selection with:

```
fccanalysis run ./ana_ZHgamma/inclusive/histmaker_recoil.py
```

and create the plots:

```
export PYTHONPATH=/afs/cern.ch/work/s/saaumill/public/MyFCCAnalyses/extras:$PYTHONPATH
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:$PATH

fccanalysis plots ./ana_ZHgamma/inclusive/plots_recoil.py
```


# Perform a Likelihood fit with Combine

We use [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/) for the profile likelihood fits. The datacards used can be found in `extras/combine`

To set it up the first time do:
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmssw-cc7
export SCRAM_ARCH="slc7_amd64_gcc700"
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src/
cmsenv 
git clone -o bendavid -b tensorflowfit git@github.com:bendavid/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit 
cd HiggsAnalysis/CombinedLimit 
scram b -j 8 
```

Every time after, set up the envirnoment with

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmssw-cc7
cd CMSSW_10_6_19_patch2/src/
cmsenv 
cd HiggsAnalysis/CombinedLimit 
scram b -j 8 
```

To run the fit:

```
text2hdf5.py mydatacard.txt -o ws_datacard.hdf5
combinetf.py ws_datacard.hdf5 -o fit_output.root -t -1 --expectSignal=1 --doImpacts
```

# Citation

The analyses are published on [CDS](https://repository.cern/communities/fcc-ped-sub/records?q=Sara%20Aumiller&l=list&p=1&s=10&sort=bestmatch). If you use any of the code or results presented, please cite:

- $ZH \rightarrow 6$ leptons analysis at FCC-ee

```
@misc{aumiller_2025_aaf59-q6v03,
  author       = {Aumiller, Sara and
                  Selvaggi, Michele},
  title        = {\$ZH\rightarrow 6\$ leptons analysis at FCC-ee},
  month        = mar,
  year         = 2025,
  publisher    = {CERN},
  doi          = {10.17181/aaf59-q6v03},
  url          = {https://doi.org/10.17181/aaf59-q6v03},
}
```

- Measurements of the $Z_{\mathrm{prod}(qq)Z(\ell \ell)Z^*(\nu \nu)}$ and $Z_{\mathrm{prod}(qq)Z(\nu \nu)Z^*(\ell \ell)}$ final states at FCC-ee

```
@misc{aumiller_2025_fnwcp-qfs39,
  author       = {Aumiller, Sara and
                  Selvaggi, Michele},
  title        = {Measurements of the \$ZH\rightarrow ZZZ^*
                   \rightarrow q\bar q \ell \bar \ell \nu \bar \nu\$
                   and \$q\bar q \nu \bar \nu\ell \bar \ell \$ final
                   states at FCC-ee
                  },
  month        = mar,
  year         = 2025,
  publisher    = {CERN},
  doi          = {10.17181/fnwcp-qfs39},
  url          = {https://doi.org/10.17181/fnwcp-qfs39},
}
```

- Probing Effective HZ$\gamma$ Couplings via H$\gamma$ Production at FCC-ee

```
@misc{herrmann_2025_kd1n4-ajd66,
  author       = {Herrmann, Lena Maria and
                  Selvaggi, Michele and
                  Aumiller, Sara},
  title        = {Probing Effective HZ\$\gamma\$ Couplings via
                   H\$\gamma\$ Production at FCC-ee
                  },
  month        = mar,
  year         = 2025,
  publisher    = {CERN},
  doi          = {10.17181/kd1n4-ajd66},
  url          = {https://doi.org/10.17181/kd1n4-ajd66},
}
```
