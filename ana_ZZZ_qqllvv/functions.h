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
    namespace ZHfunctions {
        
        // STRCUT: compute the cone isolation for reco particles
        struct coneIsolation {

            // constructor
            coneIsolation(float arg_dr_min, float arg_dr_max);

            // compute deltaR
            double deltaR(double eta1, double phi1, double eta2, double phi2) { 
                return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); 
            };

            // variables
            float dr_min = 0;
            float dr_max = 0.4;

            // operator
            Vec_f operator() (Vec_rp in, Vec_rp rps) ;
        };

        // constructor - set the cone size
        coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) {};

        // operator - compute the isolation
        Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
            // in - particles for which isolation is computed
            // rps - particles used for isolation computation aka all reconstructed particles
        
            Vec_f result;
            result.reserve(in.size());

            std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
            std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
            std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

            // compute lorentz vectors for reconstructed particles
            for(size_t i = 0; i < rps.size(); ++i) {

                ROOT::Math::PxPyPzEVector tlv;
                tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
                
                if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
                else lv_charged.push_back(tlv);
            }
            // compute lorentz vectors for particles for which isolation is computed
            for(size_t i = 0; i < in.size(); ++i) {

                ROOT::Math::PxPyPzEVector tlv;
                tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
                lv_reco.push_back(tlv);
            }

            
            // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
            for (auto & lv_reco_ : lv_reco) {
            
                double sumNeutral = 0.0;
                double sumCharged = 0.0;
            
                // charged
                for (auto & lv_charged_ : lv_charged) {
                    double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
                    if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
                }
                
                // neutral
                for (auto & lv_neutral_ : lv_neutral) {
                    double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
                    if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
                }
                
                double sum = sumCharged + sumNeutral;
                double ratio= sum / lv_reco_.P();
                result.emplace_back(ratio); // isolation value for each particle: momentum sum of charged and neutral particles in the cone divided by the reco particle momentum
            }
            return result;
        }

        // STRCUT: select particles with isolation below a threshold
        struct sel_iso {
            // constructor
            sel_iso(float arg_max_iso);
            // variables
            float m_max_iso = .25;
            // operator
            Vec_rp operator() (Vec_rp in, Vec_f iso);
        };
        // constructor - set the isolation threshold
        sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};

        // operator - select particles with isolation below the threshold
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
            // in - particles for which isolation was computed
            // iso - isolation values for each particle

            Vec_rp result;
            result.reserve(in.size());
            // select particles with isolation below the threshold
            for (size_t i = 0; i < in.size(); ++i) {
                auto & p = in[i];
                if (iso[i] < m_max_iso) {
                    result.emplace_back(p);
                }
            }
            return result;
        }

        
    }
}

#endif