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


       ROOT::VecOps::RVec<float> print_MC_info(Vec_mc mcparticles){
            ROOT::VecOps::RVec<float> pdg = MCParticle::get_pdg(mcparticles);
            auto status = MCParticle::get_genStatus(mcparticles);
            auto e = MCParticle::get_e(mcparticles);

            ROOT::VecOps::RVec<int> = MCParticle::get_parentid(ROOT::VecOps::RVec<int>("partilce0", mcparticles, "particle1 continue here"));

            // now get parent ids
            // auto parent_pdg = MCParticle::get_pdg(parents);


            // print the information
            std::cout<<"--------------- MC particles: ------------ "<<std::endl;
            for (int i = 0; i < mcparticles.size(); ++i) {
                // mcparticles[i] is a MCParticleData object: https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle_data.html 
                int pb = mcparticles[i].parents_begin;
                int pe = mcparticles[i].parents_end;
                int size_parents = pe - pb;
                std::vector<int> pdg_parents;
                for (int j = pb; j < pe; ++j) {
                    pdg_parents.push_back(mcparticles[j].PDG);
                }


                std::cout << "PDG: \t" << pdg[i] << "\t status: \t" << status[i] << " energy: \t" << e[i] << "\t # parents: \t";
                for (int j = 0; j < size_parents; ++j) {
                    std::cout << pdg_parents[j] << " \t";
                } 
                std::cout << std::endl;
            }


            return pdg;
        }





    }
}

#endif