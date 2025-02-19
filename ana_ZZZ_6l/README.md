# $e^+e^- \rightarrow Z(ll)H; H \rightarrow ZZ \rightarrow 4l$ Analysis

## Set-up

```
cd /afs/cern.ch/work/s/saaumill/public/FCCAnalyses
source ./setup.sh
mkdir build install && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 20
cd ../ana_ZZZ_6l
```

To run the analysis run

```
fccanalysis run histmaker.py
```

## Data

- Signal: `/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_llH_HZZ_llll_ecm240/`

## Resources: 
- [Tutorial](https://hep-fcc.github.io/fcc-tutorials/main/fast-sim-and-analysis/fccanalyses/doc/starterkit/FccFastSimAnalysis/Readme.html)