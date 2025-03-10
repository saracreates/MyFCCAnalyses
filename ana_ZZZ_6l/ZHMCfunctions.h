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


// FCC Analyses functions

namespace FCCAnalyses { 
    namespace ZHMCfunctions {

        ROOT::VecOps::RVec<float> print_MC_info(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers){
            // because I need to use some info in a histogram for the function to be used, I return the pdg of MC particles but that's for the sake of the example

            ROOT::VecOps::RVec<float> pdg = MCParticle::get_pdg(mcparticles);
            auto status = MCParticle::get_genStatus(mcparticles);
            auto e = MCParticle::get_e(mcparticles);

            // print the information
            std::cout<<"--------------- MC particles: ------------ "<<std::endl;
            for (int i = 0; i < mcparticles.size(); ++i) {
                // mcparticles[i] is a MCParticleData object: https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle_data.html 
                int pb = mcparticles[i].parents_begin;
                int pe = mcparticles[i].parents_end;
                int size_parents = pe - pb;
                std::vector<int> pdg_parents;


                std::cout << "PDG: \t" << pdg[i] << "\t status: " << status[i] << " energy: \t" << e[i] << "\t parents: pdg (index) - \t";
                for (int j = pb; j < pe; ++j) {
                    int ind_par = ind_parents[j];
                    std::cout << mcparticles[ind_par].PDG << " (" << ind_par << ") ";
                } 
                std::cout << std::endl;
            }
            return pdg;
        }

        edm4hep::MCParticleData find_higgs(Vec_mc mcparticles){
            // get higgs
            int pdg_higgs = 25;
            edm4hep::MCParticleData higgs;
            int count_higgs = 0;
            for(edm4hep::MCParticleData& mcp: mcparticles){
                if(mcp.PDG == pdg_higgs){
                    count_higgs++;
                    higgs = mcp;
                }
            }
            if (count_higgs != 1){
                std::cout << "ERROR: there should be exactly one Higgs in the event" << std::endl;
                exit(1);
            }
            return higgs;
        }

        Vec_mc find_Zs_from_H(Vec_mc mcparticles, Vec_i ind_parents){
            // get Zs from H->ZZ;
            int pdg_Z = 23;
            int pdg_higgs = 25;
            Vec_mc Zs;
            for(edm4hep::MCParticleData& mcp: mcparticles){
                if(mcp.PDG == pdg_Z){
                    // check if the Z has a Higgs as a parent
                    int pb = mcp.parents_begin;
                    int pe = mcp.parents_end;
                    int size_parents = pe - pb;
                    if (size_parents != 1){
                        std::cout << "ERROR: Z should have 1 parent" << std::endl;
                        exit(1);
                    }
                    int parent_id = ind_parents[pb];
                    if (mcparticles[parent_id].PDG == pdg_higgs){
                        Zs.push_back(mcp);
                    }
                }
            }
            if (Zs.size() != 2){
                std::cout << "ERROR: there should be exactly two Zs from a Higgs in the event" << std::endl;
                exit(1);
            }

            return Zs;
        }

        Vec_mc get_leptons_from_Z(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers){
            // get the leptons from ee->ZH; Z->ll;
            // the Z is not created by whizard, but the ll come from the same parents as the H, namely e+e-
            // strategy: find higgs, get parents, find leptons with the same parents

            // get higgs
            edm4hep::MCParticleData higgs = find_higgs(mcparticles);

            // get parents of higgs
            int pb = higgs.parents_begin;
            int pe = higgs.parents_end;
            int size_parents = pe - pb;
            if (size_parents != 2){
                std::cout << "ERROR: Higgs has not 2 parents" << std::endl;
                exit(1);
            }
            // find index of Higgs parents
            std::vector<int> parent_ids;
            for (int j = pb; j < pe; ++j) {
                parent_ids.push_back(ind_parents[j]);
            }
            // loop over all mc particles and for all leptons check if they have the parents with ids "parent_ids"
            Vec_mc leptons;
            std::vector<int> pdg_leptons = {11, -11, 13, -13, 15, -15};
            for(edm4hep::MCParticleData& mcp: mcparticles){
                // check if mcp is a lepton
                if(std::find(pdg_leptons.begin(), pdg_leptons.end(), mcp.PDG) != pdg_leptons.end()){
                    // check if mcp has the same parents as the Higgs
                    int pb = mcp.parents_begin;
                    int pe = mcp.parents_end;
                    int size_parents = pe - pb;
                    if (size_parents != 2) continue;
                    // find index of lepton parents
                    std::vector<int> parent_ids_lepton;
                    for (int j = pb; j < pe; ++j) {
                        parent_ids_lepton.push_back(ind_parents[j]);
                    }
                    // check if the parents are the same
                    if (parent_ids == parent_ids_lepton){
                        leptons.push_back(mcp);
                    }
                }
            }
            if (leptons.size() != 2){
                std::cout << "ERROR: could not find the two leptons from Z decay" << std::endl;
                exit(1);
            }



            return leptons;
        }

        Vec_mc get_recoil_from_Z(Vec_mc Zleptons, int ecm){
            // get recoil from Z->ll

            // combine leptons to Z
            TLorentzVector Z;
            for(edm4hep::MCParticleData& mcp: Zleptons){
                TLorentzVector tlv;
                tlv.SetXYZM(mcp.momentum.x, mcp.momentum.y, mcp.momentum.z, mcp.mass);
                Z += tlv;
            }

            // get recoil
            TLorentzVector recoil;
            recoil.SetXYZM(0, 0, 0, ecm);
            recoil -= Z;
            // create a MCParticleData object
            auto recoil_fcc = edm4hep::MCParticleData();
            recoil_fcc.momentum.x = recoil.Px();
            recoil_fcc.momentum.y = recoil.Py();
            recoil_fcc.momentum.z = recoil.Pz();
            recoil_fcc.mass = recoil.M();

            Vec_mc recoil_v;
            recoil_v.push_back(recoil_fcc);

            return recoil_v;

        }

        Vec_mc get_Z_from_leptons(Vec_mc Zleptons){
            // get Z from two leptons
            if (Zleptons.size() != 2){
                std::cout << "ERROR: Zleptons should have 2 leptons" << std::endl;
                exit(1);
            }
            // combine leptons to Z
            TLorentzVector Z;
            for(edm4hep::MCParticleData& mcp: Zleptons){
                TLorentzVector tlv;
                tlv.SetXYZM(mcp.momentum.x, mcp.momentum.y, mcp.momentum.z, mcp.mass);
                Z += tlv;
            }
            // create a MCParticleData object
            auto Z_fcc = edm4hep::MCParticleData();
            Z_fcc.momentum.x = Z.Px();
            Z_fcc.momentum.y = Z.Py();
            Z_fcc.momentum.z = Z.Pz();
            Z_fcc.mass = Z.M();

            Vec_mc Z_v;
            Z_v.push_back(Z_fcc);

            return Z_v;
        }

        float chi_square(float m1, float m2){
            // calculate chi square: m1 on-shell Z, m2 off-shell Z
            float chi2 = pow((m1 - 91.2)/2, 2) + pow((m2 - 30)/10, 2);
            return chi2;
        }

        Vec_mc get_Z_from_H(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers, int which){
            // get the leptons from ee->ZH; H->ZZ; Z->ll;
            // so there are two Zs, so two leptons pairs. "which" can be 0: on-shell Z, 1: off-shell Z
            // strategy: find Zs from higgs
            Vec_mc Z_from_H;

            if (which != 0 && which != 1){
                std::cout << "ERROR: which should be 0 or 1 refering to the first or second Z from H->ZZ" << std::endl;
                exit(1);
            }

            // get Z from Higgs
            Vec_mc Zs = find_Zs_from_H(mcparticles, ind_parents);

            // find out which of them in on and which on-shell and which of them is off-shell
            float m1 = Zs[0].mass;
            float m2 = Zs[1].mass;
            //std::cout << "Z masses: \t" << m1 << "\t" << m2 << std::endl;
            // float chi2_1 = chi_square(m1, m2);
            // float chi2_2 = chi_square(m2, m1);
            // if (chi2_1 < chi2_2){ // m1 is on-shell Z
            //     if (which == 0){ // return on-shell Z
            //         Z_from_H.push_back(Zs[0]);
            //     }
            //     else{ // return off-shell Z
            //         Z_from_H.push_back(Zs[1]);
            //     }
            // }
            // else{ // m2 is on-shell Z
            //     if (which == 0){ // return on-shell Z
            //         Z_from_H.push_back(Zs[1]);
            //     }
            //     else{ // return off-shell Z
            //         Z_from_H.push_back(Zs[0]);
            //     }
            // }
            float Z_mass = 91.2;
            if(which == 0){ // onshell
                if(abs(m1-Z_mass) < abs(m2-Z_mass)){
                    Z_from_H.push_back(Zs[0]);
                }
                else{
                    Z_from_H.push_back(Zs[1]);
                }
            }
            else{ // offshell 
                if(abs(m1-Z_mass) < abs(m2-Z_mass)){
                    Z_from_H.push_back(Zs[1]);
                }
                else{
                    Z_from_H.push_back(Zs[0]);
                }
            }

            return Z_from_H;
        }

        Vec_mc get_leptons(Vec_mc Z, Vec_i ind_daugthers, Vec_mc mcparticles){
            // get daughters (leptons) from Z
            Vec_mc lepton_daughters;
            for(edm4hep::MCParticleData& mcp: Z){ // loop over Zs, should only happen once
                int pb = mcp.daughters_begin;
                int pe = mcp.daughters_end;
                int size_daughters = pe - pb;
                for (int j = pb; j < pe; ++j) {
                    int ind_d = ind_daugthers[j];
                    edm4hep::MCParticleData daughter = mcparticles[ind_d];
                    if(abs(daughter.PDG) == 11 || abs(daughter.PDG) == 13 || abs(daughter.PDG) == 15){
                        lepton_daughters.push_back(daughter);
                    } 
                }
            }

            // check if Z has two daughters:
            int n_daughters = 2* Z.size();
            if (lepton_daughters.size() != n_daughters){
                std::cout << "ERROR: Z(s) ("<< Z.size() <<") should have "<< n_daughters << " daughters, but found "<< lepton_daughters.size() << std::endl;
                exit(1);
            }

            return lepton_daughters;

        }

        Vec_mc create_lepton_missingE(Vec_mc ll_from_Z, Vec_mc ll_from_H1, Vec_mc ll_from_H2, Vec_i ind_daughter, Vec_mc mcparticles){
            // in: all MC particles: this doesn't work, because it takes ALL MC particles and not just the end products... 

            // fill up all 6 leptons into one vector
            Vec_mc leptons;
            for(int i=0; i<ll_from_Z.size(); i++){
                leptons.push_back(ll_from_Z[i]);
            }
            for(int i=0; i<ll_from_H1.size(); i++){
                leptons.push_back(ll_from_H1[i]);
            }
            for(int i=0; i<ll_from_H2.size(); i++){
                leptons.push_back(ll_from_H2[i]);
            }

            float px = 0, py = 0, pz = 0;
            // loop over all leptons 
            for (int i=0; i<leptons.size(); i++){
                // if lepton is e or mu add momentum directly
                if(abs(leptons[i].PDG) == 11 || abs(leptons[i].PDG) == 13){
                    px += leptons[i].momentum.x;
                    py += leptons[i].momentum.y;
                    pz += leptons[i].momentum.z;
                } else if (abs(leptons[i].PDG) == 15){ // if tau, get daughters
                    int pb = leptons[i].daughters_begin;
                    int pe = leptons[i].daughters_end;
                    int size_daughters = pe - pb;
                    for (int j = pb; j < pe; ++j) {
                        int ind_d = ind_daughter[j];
                        edm4hep::MCParticleData daughter = mcparticles[ind_d];
                        // if daughter is not neutrino, check if it has more daughters.
                        if (daughter.PDG != 12 && daughter.PDG != 14 && daughter.PDG != 16){
                            int pb_d = daughter.daughters_begin;
                            int pe_d = daughter.daughters_end;
                            int size_daughters_d = pe_d - pb_d;
                            if(size_daughters_d > 0){ // check if there are daughters
                                for (int k = pb_d; k < pe_d; ++k) {
                                    int ind_d_d = ind_daughter[k];
                                    edm4hep::MCParticleData daughter_d = mcparticles[ind_d_d];
                                    // check if daughter is neutrino
                                    if (daughter_d.PDG == 12 || daughter_d.PDG == 14 || daughter_d.PDG == 16){
                                        continue;
                                    } else { // is daughter has no daughter and is not neutrino, save momentum
                                        px += daughter_d.momentum.x;
                                        py += daughter_d.momentum.y;
                                        pz += daughter_d.momentum.z;
                                    }
                                }
                            } else { // is daughter has no daughter and is not neutrino, save momentum
                                px += daughter.momentum.x;
                                py += daughter.momentum.y;
                                pz += daughter.momentum.z;
                            }
                        } else { // is neutrino
                            continue;
                        }
                        px += daughter.momentum.x;
                        py += daughter.momentum.y;
                        pz += daughter.momentum.z;
                    }
                }   
            }

            Vec_mc ret;
            edm4hep::MCParticleData res;
            res.momentum.x = px;
            res.momentum.y = py;
            res.momentum.z = pz;

            ret.emplace_back(res);
            return ret;
        }

        Vec_mc create_lepton_missingE_emu_only(Vec_mc ll_from_Z, Vec_mc ll_from_H1, Vec_mc ll_from_H2, Vec_i ind_daughter, Vec_mc mcparticles){

            Vec_mc ret;

            // fill up all 6 leptons into one vector
            Vec_mc leptons;
            for(int i=0; i<ll_from_Z.size(); i++){
                leptons.push_back(ll_from_Z[i]);
            }
            for(int i=0; i<ll_from_H1.size(); i++){
                leptons.push_back(ll_from_H1[i]);
            }
            for(int i=0; i<ll_from_H2.size(); i++){
                leptons.push_back(ll_from_H2[i]);
            }

            float px = 0, py = 0, pz = 0;
            // loop over all leptons 
            for (int i=0; i<leptons.size(); i++){
                if(abs(leptons[i].PDG) == 15){
                    return ret; // exclude event if there is a tau in the final state
                }
            }
            for (int i=0; i<leptons.size(); i++){
                px += leptons[i].momentum.x;
                py += leptons[i].momentum.y;
                pz += leptons[i].momentum.z;
            }

            edm4hep::MCParticleData res;
            res.momentum.x = px;
            res.momentum.y = py;
            res.momentum.z = pz;

            ret.emplace_back(res);
            return ret;
        }

        Vec_mc create_higgs_from_ZZ(Vec_mc res1, Vec_mc res2){
            TLorentzVector Z1;
            Z1.SetXYZM(res1[0].momentum.x, res1[0].momentum.y, res1[0].momentum.z, res1[0].mass);
            TLorentzVector Z2;
            Z2.SetXYZM(res2[0].momentum.x, res2[0].momentum.y, res2[0].momentum.z, res2[0].mass);

            TLorentzVector Higgs = Z1 + Z2;

            Vec_mc result;
            edm4hep::MCParticleData higgs;
            higgs.momentum.x = Higgs.Px();
            higgs.momentum.y = Higgs.Py();
            higgs.momentum.z = Higgs.Pz();
            higgs.mass = Higgs.M();
            result.push_back(higgs);

            return result;

        }


        // MC functions for BACKGROUND data 
        void print_Z_info(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers){
            // find Zs
            Vec_mc Zs; 
            for (int i = 0; i < mcparticles.size(); ++i) {
                if(mcparticles[i].PDG == 23){
                    Zs.push_back(mcparticles[i]);
                }
            }
            // if no Z found, print all other mcparticles
            if (Zs.size() == 0){
                std::cout << "No Zs found, printing all mcparticles" << std::endl;
                for (int i = 0; i < mcparticles.size(); ++i) {
                    std::cout << "PDG: \t" << mcparticles[i].PDG << std::endl;
                }
            }

            // print Zs info
            std::cout<<"---- number of total Zs: "<<Zs.size()<< " ---------" <<std::endl;
            // print parents and daughtersof Zs
            for(edm4hep::MCParticleData& mcp: Zs){
                std::cout << "---- new Z ----" << std::endl;
                int pb = mcp.parents_begin;
                int pe = mcp.parents_end;
                int size_parents = pe - pb;
                for (int j = pb; j < pe; ++j) {
                    int ind_par = ind_parents[j];
                    std::cout << "Z parents: \t" << mcparticles[ind_par].PDG << " (" << ind_par << ")"<< std::endl;
                } 
                int pb_d = mcp.daughters_begin;
                int pe_d = mcp.daughters_end;
                int size_daughters = pe_d - pb_d;
                for (int j = pb_d; j < pe_d; ++j) {
                    int ind_d = ind_daugthers[j];
                    std::cout << "Z daughters: \t" << mcparticles[ind_d].PDG << " (" << ind_d << ")"<< std::endl;
                }
            }
        }

        Vec_mc get_all_Zs_bg(Vec_mc mcparticles, Vec_i ind_parents, Vec_i ind_daugthers){
            // function for BACKGROUND data ee-> ZZll; Z->ll
            // get all Zs
            int pdg_Z = 23;
            std::vector<int> pdg_leptons = {11, -11, 13, -13, 15, -15};
            Vec_mc Zs;
            for(edm4hep::MCParticleData& mcp: mcparticles){
                if(mcp.PDG == pdg_Z){
                    // check that Z has two leptons as daughters
                    int pb = mcp.daughters_begin;
                    int pe = mcp.daughters_end;
                    int size_daughters = pe - pb;
                    if (size_daughters != 2){
                        continue;
                    } else {
                        int daughter_pdg1 = mcparticles[ind_daugthers[pb]].PDG;
                        int daughter_pdg2 = mcparticles[ind_daugthers[pb+1]].PDG;
                        if(std::find(pdg_leptons.begin(), pdg_leptons.end(), daughter_pdg1) != pdg_leptons.end() && std::find(pdg_leptons.begin(), pdg_leptons.end(), daughter_pdg2) != pdg_leptons.end()){
                            std::cout << "found Z with two leptons as daughters with ID:" << daughter_pdg1 << " " << daughter_pdg2 << std::endl;
                            Zs.push_back(mcp);
                        }
                    }
                }
            }
            bool debug = true;
            if (debug){
                print_Z_info(mcparticles, ind_parents, ind_daugthers);
            }

            // check if there are 2 Zs
            if(Zs.size() != 2){
                std::cout << "ERROR: there should be exactly two Zs in the event but got "<< Zs.size() << std::endl;
                std::cout << "total number of mcparticles: " << mcparticles.size() << std::endl;
                std::cout << "Only found these Zs: " << std::endl;
                print_Z_info(mcparticles, ind_parents, ind_daugthers);
                exit(1);
            }


            return Zs;
        }


        Vec_f calculate_inv_mass_of_two_leptons(Vec_mc leptons){
            // calculate invariant mass of two leptons

            // check if leptons have size 2
            if (leptons.size() != 2){
                std::cout << "ERROR: there should be exactly two leptons as an input to this function" << std::endl;
                exit(1);
            }

            Vec_f inv_mass;

            TLorentzVector lv1;
            TLorentzVector lv2;
            lv1.SetXYZM(leptons[0].momentum.x, leptons[0].momentum.y, leptons[0].momentum.z, leptons[0].mass);
            lv2.SetXYZM(leptons[1].momentum.x, leptons[1].momentum.y, leptons[1].momentum.z, leptons[1].mass);

            TLorentzVector lv = lv1 + lv2;
            inv_mass.push_back(lv.M());

            return inv_mass;
        }

        Vec_i n_daughters(Vec_mc mcparticles){
            // caluclate number of daughters of each particle
            Vec_i n_daughters;
            for(edm4hep::MCParticleData& mcp: mcparticles){
                int pb = mcp.daughters_begin;
                int pe = mcp.daughters_end;
                int size_daughters = pe - pb;
                n_daughters.push_back(size_daughters);
            }
            return n_daughters;
        }

        Vec_i n_mcparticles(Vec_mc mcparticles){
            // return number of mc particles 
            Vec_i n_mc;
            int num_mc = mcparticles.size();
            n_mc.push_back(num_mc);

            return n_mc;
        }

        

    }

}

#endif