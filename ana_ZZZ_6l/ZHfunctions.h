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

        // STRCUT: checks all lepton permuatiations and selects the best pair for the Z resonance
        struct resonanceBuilder {
            // constructor
            resonanceBuilder(float arg_resonance_mass);
            // variables
            float m_resonance_mass;
            // functions
            std::tuple<std::vector<std::vector<int>>, Vec_rp> calculate_pairs(Vec_rp legs);
            std::vector<int> find_best_combination_4(Vec_rp res, std::vector<std::vector<int> > pairs, rp other_lep_pair);
            std::vector<int> find_best_combination_6(Vec_rp res, std::vector<std::vector<int> > pairs);
            float chi_square_3Z(float m1, float m2, float m3);

            // operator
            Vec_rp operator()(Vec_rp electrons, Vec_rp muons) ;
        };

        // constructor - set the resonance mass, recoil mass, fraction of chi2 for recoil, ecm and use_MC_Kinematics flag
        resonanceBuilder::resonanceBuilder(float arg_resonance_mass) {
            m_resonance_mass = arg_resonance_mass;
        }
        // functions
        std::tuple<std::vector<std::vector<int>>, Vec_rp> resonanceBuilder::calculate_pairs(Vec_rp legs) {
            // legs - leptons/muons
            // returns all possible resonance candidates (pairs of l+l-) and their corresponding indicies
            Vec_rp result;
            result.reserve(3);

            std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
            int n = legs.size(); // number of leptons/muons
        
            if(n > 1) { // at last two leptons required
                ROOT::VecOps::RVec<bool> v(n);
                std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
                do { // calculate the resonance for all permutations (pairs of leptons)
                    std::vector<int> pair;
                    rp reso; // ReconstructedParticleData object for the resonance
                    reso.charge = 0;
                    TLorentzVector reso_lv; 
                    for(int i = 0; i < n; ++i) {
                        if(v[i]) {
                            pair.push_back(i);
                            reso.charge += legs[i].charge;

                            TLorentzVector leg_lv;
                            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);

                            reso_lv += leg_lv;
                        }
                    }

                    if(reso.charge != 0) continue; // neglect non-zero charge pairs
                    reso.momentum.x = reso_lv.Px();
                    reso.momentum.y = reso_lv.Py();
                    reso.momentum.z = reso_lv.Pz();
                    reso.mass = reso_lv.M();
                    result.emplace_back(reso);
                    pairs.push_back(pair);

                } while(std::next_permutation(v.begin(), v.end())); // loop over all permutations
            }
            else {
                std::cout << "ERROR: resonanceBuilder, at least two leptons required." << std::endl;
                exit(1);
            }

            return std::make_tuple(pairs, result);
        };

        float resonanceBuilder::chi_square_3Z(float m1, float m2, float m3){
            // mi are the masses of the potential resonances
            float Z_mass = m_resonance_mass; // 91.2 GeV
            float Z_mass_offshell = 30; // GeV

            float nat_Z_width = 2; // GeV
            float nat_Z_width_offshell = 10; // GeV

            float chi2_1 = pow((m1-Z_mass)/nat_Z_width, 2) + pow((m2-Z_mass)/nat_Z_width, 2) + pow((m3-Z_mass_offshell)/nat_Z_width_offshell, 2);
            float chi2_2 = pow((m3-Z_mass)/nat_Z_width, 2) + pow((m1-Z_mass)/nat_Z_width, 2) + pow((m2-Z_mass_offshell)/nat_Z_width_offshell, 2);
            float chi2_3 = pow((m2-Z_mass)/nat_Z_width, 2) + pow((m3-Z_mass)/nat_Z_width, 2) + pow((m1-Z_mass_offshell)/nat_Z_width_offshell, 2);

            // return the minimum chi2
            return std::min({chi2_1, chi2_2, chi2_3});
        }

        std::vector<int> resonanceBuilder::find_best_combination_6(Vec_rp res, std::vector<std::vector<int> > pairs) {
            // returns the best combination of 6 leptons
            // res - 9 resonance candidates (9 leptons paris of e^+e^- or mu^+mu^-)
            // pairs - 9 pairs of indices of leptons for each resonance candidate
            std::cout << "finding best combination out of 6" << std::endl;

            std::vector<float> masses; // will have 9 entries
            for (int i = 0; i < res.size(); ++i) {
                masses.push_back(res[i].mass);
            } 
            // which combination of 3 entries in dist_Z are the lowest? Considering that the pairs entries are unqiue (because one lepton can only be in one combination)
            std::vector<vector<int>> allowed_combinations; // {{0,1,2}, {...}, ...} - 3 pairs of l+l- (6 allowed combinations) 
            std::vector<float> chi2_values;
            for (int i = 0; i<masses.size(); ++i) {
                for (int j=i+1; j<masses.size(); ++j) {
                    for (int k=j+1; k<masses.size(); ++k) {
                        // now check if combination is allowed, check pairs
                        if (haveCommonElement(pairs[i], pairs[j]) || haveCommonElement(pairs[i], pairs[k]) || haveCommonElement(pairs[j], pairs[k])) {
                            continue;
                        } else {
                            allowed_combinations.push_back({i, j, k});
                            // float dist = masses[i] + masses[j] + masses[k]; // naive approach
                            float chi2 = chi_square_3Z(masses[i], masses[j], masses[k]);
                            std::cout << "masses: \t" << masses[i] << "\t" << masses[j] << "\t" << masses[k] << "\t chi2: \t" << chi2 << std::endl;
                            chi2_values.push_back(chi2);
                        }
                    }
                }
            }
            // std::cout << "chi2_values: " << chi2_values.size() << std::endl;  // 6
            // now find the minimum chi2 value for the allowed combinations
            int ind_best_combi = std::distance(chi2_values.begin(), 
                                   std::min_element(chi2_values.begin(), chi2_values.end()) // find the minimum chi2 value and returns an iterator pointing to the smallest element 
                                   ); // computes the index of the iterator relative to the beginning


            return allowed_combinations[ind_best_combi]; // e.g. {0, 1, 2}

        };

        std::vector<int> resonanceBuilder::find_best_combination_4(Vec_rp res, std::vector<std::vector<int> > pairs, rp other_lep_pair) {
            std::cout << "finding best combination out of 4" << std::endl;
            // returns the best combination of 4 leptons
            // res - 4 resonance candidates
            // pairs - 4 pairs of indices of leptons for each resonance candidate
            float mass_other = other_lep_pair.mass;
            std::vector<float> masses; // will have 9 entries
            for (int i = 0; i < res.size(); ++i) {
                masses.push_back(res[i].mass);
            } 
            // which combination of 2 entries in dist_Z are the lowest? Considering that the pairs entries are unqiue (because one lepton can only be in one combination)
            std::vector<vector<int>> allowed_combinations; // {{0,1}, {...}, ...}
            std::vector<float> chi2_values;
            for (int i = 0; i<masses.size(); ++i) {
                for (int j=i+1; j<masses.size(); ++j) {
                    // now check if combination is allowed, check pairs
                    if (haveCommonElement(pairs[i], pairs[j])) {
                        continue;
                    } else {
                        allowed_combinations.push_back({i, j});
                        float chi2 = chi_square_3Z(masses[i], masses[j], mass_other);
                        std::cout << "masses: \t" << masses[i] << "\t" << masses[j] << "\t" << mass_other << "\t chi2: \t" << chi2 << std::endl;
                        chi2_values.push_back(chi2);
                    }
                }
            }
            // now find the minimum chi2 value for the allowed combinations
            int ind_best_combi = std::distance(chi2_values.begin(), 
                                   std::min_element(chi2_values.begin(), chi2_values.end()) // find the minimum chi2 value and returns an iterator pointing to the smallest element 
                                   ); // computes the index of the iterator relative to the beginning

            return allowed_combinations[ind_best_combi]; // e.g. {0, 1}

        };

        // operator - checks all lepton permuatiations and selects the best pair for the Z resonance
        Vec_rp resonanceBuilder::resonanceBuilder::operator()(Vec_rp electrons, Vec_rp muons) {

            Vec_rp result_e;
            std::vector<std::vector<int>> pairs_e; 
            Vec_rp result_m;
            std::vector<std::vector<int>> pairs_m;

            // final output - 3 resonance particles
            Vec_rp resonance;

            if (electrons.size() + muons.size() != 6) {
                std::cout << "ERROR: 6 leptons required" << std::endl;
                exit(1);
            }

            // pick the best resulting pairs
            if(electrons.size() == 0){ // all leptons are muons
                auto [pairs_m, result_m] = calculate_pairs(muons);
                // now find the 3 best pairs out of 9 possible combinations
                std::vector<int> best_resonance = find_best_combination_6(result_m, pairs_m);
                for (int i = 0; i < best_resonance.size(); ++i) {
                    resonance.push_back(result_m[best_resonance[i]]);
                }
            } else if (muons.size() == 0){ // all leptons are electrons
                auto [pairs_e, result_e] = calculate_pairs(electrons);
                // now find the 3 best pairs out of 9 possible combinations
                std::vector<int> best_resonance = find_best_combination_6(result_e, pairs_e);
                for (int i = 0; i < best_resonance.size(); ++i) {
                    resonance.push_back(result_e[best_resonance[i]]);
                }
            } else if (muons.size() == 2){ // ONE muon pair and find best 2 pairs out of 4 possible electron combinations
                auto [pairs_m, result_m] = calculate_pairs(muons);
                rp muon_pair = result_m[pairs_m[0][0]]; 
                resonance.push_back(muon_pair);

                auto [pairs_e, result_e] = calculate_pairs(electrons);
                std::vector<int> best_resonance = find_best_combination_4(result_e, pairs_e, muon_pair);
                for (int i = 0; i < best_resonance.size(); ++i) {
                    resonance.push_back(result_e[best_resonance[i]]);
                }
            } else if (electrons.size() == 2){ // ONE electron pair and find best 2 pairs out of 4 possible muon combinations
                auto [pairs_e, result_e] = calculate_pairs(electrons);
                rp electron_pair = result_e[pairs_e[0][0]];
                resonance.push_back(electron_pair);

                auto [pairs_m, result_m] = calculate_pairs(muons);
                std::vector<int> best_resonance = find_best_combination_4(result_m, pairs_m, electron_pair);
                for (int i = 0; i < best_resonance.size(); ++i) {
                    resonance.push_back(result_m[best_resonance[i]]);
                }
            } else { 
                std::cout << "ERROR: 6 leptons must be electron and muon pairs" << std::endl;
                exit(1);
            }
            // check is resonance has 3 entries
            if(resonance.size() != 3) {
                std::cout << "ERROR: 3 resonance particles required" << std::endl;
                exit(1);
            }


            // sort the resonance particles by mass difference to Z mass
            std::sort(resonance.begin(), resonance.end(), [this](const auto& a, const auto& b) {
                return std::abs(a.mass - m_resonance_mass) < std::abs(b.mass - m_resonance_mass);
            });


            return resonance; // returns 3 resonance particles: on-shell, o-shell, off-shell Z
        }

        // STRUCT: compute the recoil mass
        struct recoilBuilder {
            // constructor
            recoilBuilder(float arg_resonance_mass, float arg_ecm);
            // variables
            float m_resonance_mass;
            float ecm;
            // operator
            Vec_rp operator()(Vec_rp res1, Vec_rp res2) ;
        };
        // constructor - set the resonance mass
        recoilBuilder::recoilBuilder(float arg_resonance_mass, float arg_ecm) {
            m_resonance_mass = arg_resonance_mass;
            ecm = arg_ecm;
        }
        // operator - compute the recoil mass
        Vec_rp recoilBuilder::operator()(Vec_rp res1, Vec_rp res2) {
            // resonance - 3 resonance particles: on-shell, o-shell, off-shell Z
            // if(resonance.size() != 3) {
            //     std::cout << "ERROR: 3 resonance particles required" << std::endl;
            //     exit(1);
            // }

            TLorentzVector Z1;
            Z1.SetXYZM(res1[0].momentum.x, res1[0].momentum.y, res1[0].momentum.z, res1[0].mass);
            TLorentzVector Z2;
            Z2.SetXYZM(res2[0].momentum.x, res2[0].momentum.y, res2[0].momentum.z, res2[0].mass);

            // compute the recoil mass - consider that we need to test two on-shell Zs
            

            // get recoil
            TLorentzVector recoil1;
            recoil1.SetXYZM(0, 0, 0, ecm);
            recoil1 -= Z1;
            float m_r1 = recoil1.M();

            TLorentzVector recoil2;
            recoil2.SetXYZM(0, 0, 0, ecm);
            recoil2 -= Z2;
            float m_r2 = recoil2.M();

            // check which recoil fits best with the Higgs mass 
            float chi2_1 = pow(m_r1-m_resonance_mass, 2);
            float chi2_2 = pow(m_r2-m_resonance_mass, 2);

            // select the best fit
            TLorentzVector recoil;
            if (chi2_1 < chi2_2) {
                recoil = recoil1;
            } else {
                recoil = recoil2;
            }

            // create a ReconstructedParticleData object
            rp recoil_fcc;
            recoil_fcc.momentum.x = recoil.Px();
            recoil_fcc.momentum.y = recoil.Py();
            recoil_fcc.momentum.z = recoil.Pz();
            recoil_fcc.mass = recoil.M();

            Vec_rp recoil_v;
            recoil_v.push_back(recoil_fcc);

            return recoil_v; // returns the recoil particle (MCParticleData object)

        }


        // helper functions

        // returns missing energy vector, based on reco particles
        Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
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

            std::cout << "total momentum: " << std::sqrt(px*px + py*py + pz*pz) << std::endl;

            ret.emplace_back(res);
            return ret;
        }

    }
}

#endif