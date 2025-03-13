#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>
#include <tuple>
#include <algorithm>
#include <numeric>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"


// helper functions

bool haveCommonElement(const std::vector<int>& array1, const std::vector<int>& array2) {
    std::unordered_set<int> elements(array1.begin(), array1.end());

    for (int num : array2) {
        if (elements.find(num) != elements.end()) {
            return true; // Found a common element
        }
    }

    return false; // No common elements found
}

rp return_rp_from_tlv(TLorentzVector tlv) {
    rp rp_fcc;
    rp_fcc.momentum.x = tlv.Px();
    rp_fcc.momentum.y = tlv.Py();
    rp_fcc.momentum.z = tlv.Pz();
    rp_fcc.mass = tlv.M();
    return rp_fcc;
}

// FCC Analyses functions

namespace FCCAnalyses { 
    namespace MCfunctions {
        
        Vec_mc get_leptons_from_onshell_Z_decay(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers){
            // get all leptons from a onshell Z decay from Higgs H->ZZ*; on-shell Z->ll
            Vec_mc leptons_from_Z;

            // loop over all MC particles and check the daughters of the Z. If lepton, add to the list
            for(edm4hep::MCParticleData& mcp: mcparticles){ // loop over Zs, should only happen once
                if (mcp.PDG == 23) {
                    // check if parent is Higgs
                    int pab = mcp.parents_begin;
                    int pae = mcp.parents_end;
                    int size_parents = pae - pab;
                    if (size_parents != 1) continue;
                    // std::cout << "Found Z with " << size_parents << " parents" << std::endl;
                    int ind_parent = ind_parents[pab];
                    // std::cout << "Particle pdg: " << mcparticles[ind_parent].PDG << std::endl;
                    if (mcparticles[ind_parent].PDG != 25) continue;
                    // std::cout << "Found Z with parent Higgs" << std::endl;
                    // check if Z is on-shell
                    if (mcp.mass < 60 || mcp.mass > 110) continue;
                    // std::cout << "Found Z with mass " << mcp.mass << std::endl;
                    // loop over daughters
                    int pb = mcp.daughters_begin;
                    int pe = mcp.daughters_end;
                    int size_daughters = pe - pb;
                    for (int j = pb; j < pe; ++j) {
                        int ind_d = ind_daugthers[j];
                        edm4hep::MCParticleData daughter = mcparticles[ind_d];
                        if(abs(daughter.PDG) == 11 || abs(daughter.PDG) == 13 || abs(daughter.PDG) == 15){
                            // std::cout << "Found lepton from Z decay with PDG " << daughter.PDG << std::endl;
                            leptons_from_Z.push_back(daughter);
                        } 
                    }
                }
            }
            // if (leptons_from_Z.size() != 0) std::cout << "Found " << leptons_from_Z.size() << " leptons from Z decay" << std::endl;

            return leptons_from_Z;
            
        }
        
    }
}

#endif