#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"
#include <edm4hep/utils/ParticleIDUtils.h>
#include <edm4hep/ParticleIDCollection.h>


namespace FCCAnalyses { 

    // full sim functions
    
    Vec_rp unBoostCrossingAngle(Vec_rp in, float angle) {
        Vec_rp result;
        float ta = std::tan(angle);
        for (size_t i=0; i < in.size(); ++i) {
            auto & p = in[i];
            edm4hep::ReconstructedParticleData newp = p;
            float e = p.energy;
            float px = p.momentum.x;
            float e_prime = e * sqrt(1 + ta*ta) + px * ta;
            float px_prime = px * sqrt(1 + ta*ta) + e * ta;
            newp.momentum.x = px_prime;
            newp.energy = e_prime;
            result.push_back(newp);
        }
        return result;
    }

    Vec_mc unBoostCrossingAngle(Vec_mc in, float angle) {
        Vec_mc result;
        float ta = std::tan(angle);
        for (size_t i=0; i < in.size(); ++i) {
            auto & p = in[i];
            TLorentzVector mc_;
            mc_.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
            edm4hep::MCParticleData newp = p;
            float e = mc_.Energy();
            float px = p.momentum.x;
            float px_prime = px * sqrt(1 + ta*ta) + e * ta;
            newp.momentum.x = px_prime;
            result.push_back(newp);
        }
        return result;
    }

    Vec_rp sel_type(int type, Vec_rp in) {
        Vec_rp res;
        for(auto &p : in) {
            if(std::abs(p.PDG) == type) {
            //if(std::abs(p.PDG) == type) {
                res.push_back(p);
            }
        }
        return res;
    }
    
    Vec_mc sel_type(int type, Vec_mc in) {
        Vec_mc res;
        for(auto &p : in) {
            if(std::abs(p.PDG) == type) {
                res.push_back(p);
            }
        }
        return res;
    }

    Vec_rp shift_E_photons(Vec_rp photons){
        // shift the energy of the photons by 2.6 GeV due to wrong calibration!
        // momentum must also be modified by the appropriate factor
        Vec_rp result;
        result.reserve(photons.size());
        for (const auto & photon : photons) {
            edm4hep::ReconstructedParticleData newp = photon;
            float e = photon.energy;
            float px = photon.momentum.x;
            float py = photon.momentum.y;
            float pz = photon.momentum.z;
            float p = sqrt(px*px + py*py + pz*pz);
            float factor = (e - 2.6) / e;
            newp.energy -= 2.6;
            newp.momentum.x *= factor;
            newp.momentum.y *= factor;
            newp.momentum.z *= factor;
            result.push_back(newp);
        }

        return result;

    }

    // exclusive bb functions 
    Vec_f get_bscores(ROOT::VecOps::RVec<edm4hep::ParticleIDData> b_tags_coll) {
        // check size of the b-tag collection - must be 2
        if(b_tags_coll.size() != 2) {
            std::cout << "ERROR: Expected two b-tag collections! " << std::endl;
            exit(1);
        }

        Vec_f result;
        result.reserve(2); // two b scores 

        for (const auto b_tags : b_tags_coll) {
            // std::cout << "jet likelihood: " << b_tags.likelihood << std::endl;
            // std::cout << "jet pdg: " << b_tags.PDG << std::endl;
            result.emplace_back(b_tags.likelihood);
        }

        return result;
    }

    Vec_rp get_jet_res(Vec_rp in){
        // check if there are two jets
        if(in.size() != 2) {
            std::cout << "ERROR: Expected two jets! " << std::endl;
            exit(1);
        }
        TLorentzVector jet1;
        jet1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);
        TLorentzVector jet2;
        jet2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);
        TLorentzVector res = jet1 + jet2;

        edm4hep::ReconstructedParticleData jet_res;
        jet_res.momentum.x = res.Px();
        jet_res.momentum.y = res.Py();
        jet_res.momentum.z = res.Pz();
        jet_res.mass = res.M();

        Vec_rp result;
        result.reserve(1);
        result.emplace_back(jet_res);
        return result;
    }


    ROOT::VecOps::RVec<TLorentzVector> return_p4_jets(Vec_rp in) {
        ROOT::VecOps::RVec<TLorentzVector> result;
        result.reserve(in.size());
        for (size_t i = 0; i < in.size(); ++i) {
            auto & p = in[i];
            TLorentzVector tlv;
            tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
            result.emplace_back(tlv);
        }
        return result;
    }

    namespace fullsimtagger {

        int get_higgs_daughters_MC_pdg(Vec_mc mcparticles, Vec_i ind_daugthers){
            // get higgs
            int pdg_higgs = 25;
            edm4hep::MCParticleData higgs;
            int count_higgs = 0;
            std::vector<int> daugthers_ids;
            for(edm4hep::MCParticleData& mcp: mcparticles){
                if(mcp.PDG == pdg_higgs){
                    count_higgs++;
                    higgs = mcp;
                    int db = higgs.daughters_begin;
                    int de = higgs.daughters_end;
                    int size_daughters = de - db;
                    if (size_daughters != 2){
                        std::cout << "ERROR: Higgs should have 2 daughters" << std::endl;
                        exit(1);
                    }
                    // get PDG of daughters
                    for (int j = db; j < de; ++j) {
                        int ind_daugther = ind_daugthers[j];
                        daugthers_ids.push_back(mcparticles[ind_daugther].PDG);
                    }
                }
            }
            if (count_higgs != 1){
                std::cout << "ERROR: there should be exactly one Higgs in the event" << std::endl;
                exit(1);
            }

            // the pdgs of the daughters can be 1, 2, 3, 4, 5, 15 and -1, -2, -3, -4, -5, -15 or 21 (gluon). Check if the two values are valid and return the positive value
            int pdg1 = daugthers_ids[0];
            int pdg2 = daugthers_ids[1];
            if (std::abs(pdg1) != std::abs(pdg2)) {
                std::cout << "ERROR: the two daughters of the Higgs should have the same absolute PDG value" << std::endl;
                exit(1);
            }
            if (std::abs(pdg1) != 21) {
                if (-pdg1 != pdg2) {
                    std::cout << "ERROR: the two daughters of the Higgs should be particle-antiparticle pairs" << std::endl;
                    exit(1);
                }
            }

            return std::abs(pdg1); // return the absolute value of the PDG of the daughters
        }

        int is_of_flavor(int true_pdg, int asked_flavor){
            if (true_pdg == asked_flavor) {
                return true;
            } else {
                return false;
            }
        }

        int print_scores(Vec_f scores){
            // print the scores of the b-tagging
            std::cout << "b-tagging scores: ";
            for (size_t i = 0; i < scores.size(); ++i) {
                std::cout << scores[i] << " ";
            }
            std::cout << std::endl;
            return 0;
        }
        

    }
        
    
    
    namespace ZHfunctions {


// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass_recoil {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    Vec_rp operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) ;
};

resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

Vec_rp resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) {

    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            rp reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

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

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        Vec_rp bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }

     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}    




struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    Vec_rp operator() (Vec_rp in, Vec_f iso);
  };

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
    Vec_rp result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}

 
// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    Vec_f operator() (Vec_rp in, Vec_rp rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
  
    Vec_f result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }
    
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
        result.emplace_back(ratio);
    }
    return result;
}
 
 
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
    ret.emplace_back(res);
    return ret;
}


float  print_momentum(Vec_rp  in) {
    

    ROOT::VecOps::RVec<float> ptype = FCCAnalyses::ReconstructedParticle::get_p(in);
    auto ptype1 = FCCAnalyses::ReconstructedParticle::get_n(in);
    std::cout <<  " number: " << ptype1 << std::endl;
   

    auto max_it = std::max_element(ptype.begin(), ptype.end());
    auto index = std::distance(ptype.begin(), max_it);
    
    if (max_it != ptype.end()) {
        std::cout << "Maximum value: " << *max_it << " at index: " << index << std::endl;
    } else {
        std::cout << "The vector is empty!" << std::endl;
    }
    if(index != 0) {
        std::cout << "not sorted!" << std::endl;
        std::cout <<  " momentum: " << ptype << std::endl;
        std::cout << "Maximum value: " << *max_it << " at index: " << index << std::endl;
    }
   
    
    return 1;
}


ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sort_by_energy(Vec_rp particles) {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sorted = particles; // make a copy
    
    for (size_t i = 0; i < sorted.size(); ++i) {
        // Assume the current element is the one with maximum momentum
        size_t maxIndex = i;
        for (size_t j = i + 1; j < sorted.size(); ++j) {
            // Compare momentum values; if a particle with a higher momentum is found, record its index
            if (sorted[j].energy > sorted[maxIndex].energy) {
                maxIndex = j;
            }
        }
        // If a particle with a higher momentum was found, swap it with the current element
        if (maxIndex != i) {
            std::swap(sorted[i], sorted[maxIndex]);
        }
    }

    // Return only the highest energetic element if available
    if (!sorted.empty()) {
        return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{sorted[0]};
    }

    return sorted;
}


ROOT::VecOps::RVec<float> ee_costheta_max(ROOT::VecOps::RVec<float>  in) {

    //std::cout << "in: " << typeid(in).name() << std::endl;

    ROOT::VecOps::RVec<float> copy = in;
    
   // std::cout <<  " costheta " << in << std::endl;
    std::sort(copy.begin(), copy.end(), [](float a, float b) {
        return std::abs(a) > std::abs(b);
    });
   //std::cout << "out: " << typeid(copy).name() << std::endl;
  
    return copy ;
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




 
 

}}

#endif