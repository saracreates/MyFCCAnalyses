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

        Vec_rp get_two_jets_res(const TLorentzVector &tlv1, const TLorentzVector &tlv2){
            // add the two TLorenzVector to get the resonance of the two jets

            TLorentzVector reso_lv = tlv1 + tlv2;
            rp reso = return_rp_from_tlv(reso_lv);
            Vec_rp result;
            result.emplace_back(reso);
            return result;
        }

        Vec_rp get_two_lep_res(rp l1, rp l2){

            // add the two TLorenzVector to get the resonance of the two leptons
            TLorentzVector tlv1;
            tlv1.SetPxPyPzE(l1.momentum.x, l1.momentum.y, l1.momentum.z, l1.energy);
            TLorentzVector tlv2;
            tlv2.SetPxPyPzE(l2.momentum.x, l2.momentum.y, l2.momentum.z, l2.energy);
            TLorentzVector reso_lv = tlv1 + tlv2;
            rp reso = return_rp_from_tlv(reso_lv);
            Vec_rp result;
            result.emplace_back(reso);
            return result;
        }

        // returns missing four momentum vector, based on reco particles
        Vec_rp missingParticle(float ecm, Vec_rp in, float p_cutoff = 0.0) {
            float px = 0, py = 0, pz = 0, e = 0;
            for(auto &p : in) {
                if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
                px += -p.momentum.x;
                py += -p.momentum.y;
                pz += -p.momentum.z;
                e += p.energy;
            }
            
            Vec_rp ret;
            rp res;
            res.momentum.x = px;
            res.momentum.y = py;
            res.momentum.z = pz;
            res.energy = ecm-e;
            ret.emplace_back(res);
            return ret;
        }


        Vec_rp get_recoil_jets(float ecm, const TLorentzVector &tlv1, const TLorentzVector &tlv2){
            // get the recoil of the two jets from the two leptons

            TLorentzVector tlv_jets = tlv1 + tlv2;

            // calculate recoil
            auto tlv_recoil = TLorentzVector(0, 0, 0, ecm);
            tlv_recoil -= tlv_jets;
      
            rp recoil = return_rp_from_tlv(tlv_recoil);
            Vec_rp result;
            result.emplace_back(recoil);
            return result;
        }

        Vec_rp get_recoil_from_lep_and_jets(float ecm, const TLorentzVector &jet1, const TLorentzVector &jet2, rp lep1, rp lep2){
            // get the recoil of the two jets from the two leptons

            TLorentzVector tlv_jets = jet1 + jet2;

            TLorentzVector tlv_lep1;
            tlv_lep1.SetPxPyPzE(lep1.momentum.x, lep1.momentum.y, lep1.momentum.z, lep1.energy);
            TLorentzVector tlv_lep2;
            tlv_lep2.SetPxPyPzE(lep2.momentum.x, lep2.momentum.y, lep2.momentum.z, lep2.energy);
            TLorentzVector tlv_leptons = tlv_lep1 + tlv_lep2;

            // calculate recoil
            auto tlv_recoil = TLorentzVector(0, 0, 0, ecm);
            tlv_recoil -= tlv_jets;
            tlv_recoil -= tlv_leptons;

            rp recoil = return_rp_from_tlv(tlv_recoil);
            Vec_rp result;
            result.emplace_back(recoil);
            return result;

        }

        float dot_prod(Vec_rp miss_particle, TLorentzVector jet){
            // caluclate dot product between momentum of missing particle and jet because in case of the signal, they should not be alligned

            // check size of the missing particle vector
            if (miss_particle.size() != 1){
                std::cout << "Error: missing particle vector should have size one, but got "<< miss_particle.size() << std::endl;
                exit(1);
            }
                

            rp miss_part = miss_particle.at(0);
            rp jet_rp = return_rp_from_tlv(jet);

            float dot_product = miss_part.momentum.x * jet_rp.momentum.x + miss_part.momentum.y * jet_rp.momentum.y + miss_part.momentum.z * jet_rp.momentum.z;

            // norm it via its magnitude
            float mag_miss = std::sqrt(miss_part.momentum.x * miss_part.momentum.x + miss_part.momentum.y * miss_part.momentum.y + miss_part.momentum.z * miss_part.momentum.z);
            float mag_jet = std::sqrt(jet_rp.momentum.x * jet_rp.momentum.x + jet_rp.momentum.y * jet_rp.momentum.y + jet_rp.momentum.z * jet_rp.momentum.z);

            float dot_product_norm = dot_product / (mag_miss * mag_jet);

            return dot_product_norm;
        }

        // calculate the cosine(theta) of the missing energy vector
        float get_cosTheta_miss(Vec_rp met){
            
            float costheta = 0.;
            if(met.size() > 0) {
                
                TLorentzVector lv_met;
                lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
                costheta = fabs(std::cos(lv_met.Theta()));
            }
            return costheta;
        }

        float miss_pT(Vec_rp met){
            // calculate the missing transverse energy
            float px_miss = met[0].momentum.x;
            float py_miss = met[0].momentum.y;

            float pT_miss = std::sqrt(px_miss * px_miss + py_miss * py_miss);

            return pT_miss;
        }

        // float func_neg_x(float x1){
        //     // check that x1 is smaller than 0, else throw error
        //     if (x1 > 0){
        //         std::cout << "Error: Number should be smaller than 1, but got " << x1 << std::endl;
        //         exit(1);
        //     }
        //     return -0.8 * x1 - 1;
        // }

        // float func_pos_x(float x1){
        //     // check that x1 is greater than 0, else throw error
        //     if (x1 < 0){
        //         std::cout << "Error: Number should be greater than -1, but got " << x1 << std::endl;
        //         exit(1);
        //     }
        //     return - 0.8 * x1 + 1;
        // }


        // int dot_prod_cut(float x1, float x2){
        //     // bool == 1 equals to true, 
        //     // bool == 0 equals to false

        //     // check if x1 if between -1 and 1, else through error
        //     if (x1 > 1 || x1 < -1){
        //         std::cout << "Error: dot product should be between -1 and 1, but got " << x1 << std::endl;
        //         exit(1);
        //     }

        //     bool is_in;

        //     // now define function
        //     if (x1>0){
        //         float x2_lim = func_pos_x(x1);
        //         if (x2 < x2_lim){
        //             is_in = true;
        //         } else {
        //             is_in = false;
        //         }
        //     } else if (x1 < 0){
        //         float x2_lim = func_neg_x(x1);
        //         if (x2 > x2_lim){
        //             is_in = true;
        //         } else {
        //             is_in = false;
        //         }
        //     // } else if (x1 == 0.0){
        //     //     is_in = true;
        //     } else {
        //         std::cout << "Error: dot product should be a number, but got " << x1 << std::endl;
        //         exit(1);
        //     }

        //     if (is_in){
        //         return 1;
        //     } else {
        //         std::cout << " cutting value " << std::endl;
        //         return 0;
        //     }
        // }
        
    }
}

#endif