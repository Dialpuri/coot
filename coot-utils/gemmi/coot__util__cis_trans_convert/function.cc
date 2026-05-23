#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__omega_torsion/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "function.hh"

namespace coot {
namespace util {

int
cis_trans_convert_gemmi(std::pair<gemmi::Residue*, gemmi::Residue*> mol_residues,
                        std::vector<gemmi::Residue*>& trans_residues,
                        std::vector<gemmi::Residue*>& cis_residues) {
    int istatus = 0;
    
    if (!mol_residues.first || !mol_residues.second)
        return 0;
    
    // Calculate omega torsion using coot's omega_torsion_gemmi
    std::string altconf("");
    std::pair<bool, double> omega =
        coot::util::omega_torsion_gemmi(mol_residues.first, mol_residues.second, altconf);
    
    // Determine if cis or trans
    short int is_cis_flag = 0;
    std::vector<gemmi::Residue*>* cis_trans_init_match = &trans_residues;
    std::vector<gemmi::Residue*>* converted_residues = &cis_residues;
    
    if (omega.first && omega.second < 1.57 && omega.second > -1.57) {
        // CIS peptide - making it TRANS
        is_cis_flag = 1;
        cis_trans_init_match = &cis_residues;
        converted_residues = &trans_residues;
    }
    
    // Get atoms from mol_residues
    gemmi::Atom* mol_residue_CA_1 = nullptr;
    gemmi::Atom* mol_residue_C_1 = nullptr;
    gemmi::Atom* mol_residue_O_1 = nullptr;
    gemmi::Atom* mol_residue_CA_2 = nullptr;
    gemmi::Atom* mol_residue_N_2 = nullptr;
    gemmi::Atom* mol_residue_H_2 = nullptr;
    
    for (auto& atom : mol_residues.first->atoms) {
        if (atom.name == "CA") mol_residue_CA_1 = &atom;
        if (atom.name == "C") mol_residue_C_1 = &atom;
        if (atom.name == "O") mol_residue_O_1 = &atom;
    }
    
    for (auto& atom : mol_residues.second->atoms) {
        if (atom.name == "CA") mol_residue_CA_2 = &atom;
        if (atom.name == "N") mol_residue_N_2 = &atom;
        if (atom.name == "H") mol_residue_H_2 = &atom;
    }
    
    if (mol_residue_CA_1 && mol_residue_C_1 && mol_residue_O_1 &&
        mol_residue_CA_2 && mol_residue_N_2) {
        
        // Get atoms from cis_trans_init_match residues
        gemmi::Atom* match_CA_1 = nullptr;
        gemmi::Atom* match_C_1 = nullptr;
        gemmi::Atom* match_O_1 = nullptr;
        gemmi::Atom* match_CA_2 = nullptr;
        gemmi::Atom* match_N_2 = nullptr;
        
        for (auto& atom : (*cis_trans_init_match)[0]->atoms) {
            if (atom.name == "CA") match_CA_1 = &atom;
            if (atom.name == "C") match_C_1 = &atom;
            if (atom.name == "O") match_O_1 = &atom;
        }
        
        for (auto& atom : (*cis_trans_init_match)[1]->atoms) {
            if (atom.name == "CA") match_CA_2 = &atom;
            if (atom.name == "N") match_N_2 = &atom;
        }
        
        if (match_CA_1 && match_C_1 && match_O_1 && match_CA_2 && match_N_2) {
            
            // Get atoms from converted_residues
            gemmi::Atom* conv_CA_1 = nullptr;
            gemmi::Atom* conv_C_1 = nullptr;
            gemmi::Atom* conv_O_1 = nullptr;
            gemmi::Atom* conv_CA_2 = nullptr;
            gemmi::Atom* conv_N_2 = nullptr;
            
            for (auto& atom : (*converted_residues)[0]->atoms) {
                if (atom.name == "CA") conv_CA_1 = &atom;
                if (atom.name == "C") conv_C_1 = &atom;
                if (atom.name == "O") conv_O_1 = &atom;
            }
            
            for (auto& atom : (*converted_residues)[1]->atoms) {
                if (atom.name == "CA") conv_CA_2 = &atom;
                if (atom.name == "N") conv_N_2 = &atom;
            }
            
            if (conv_CA_1 && conv_C_1 && conv_O_1 && conv_CA_2 && conv_N_2) {
                // Build coordinate vectors
                std::vector<clipper::Coord_orth> current;
                std::vector<clipper::Coord_orth> cis_trans_init;
                std::vector<clipper::Coord_orth> converted;
                
                current.push_back(clipper::Coord_orth(mol_residue_CA_1->pos.x,
                                                      mol_residue_CA_1->pos.y,
                                                      mol_residue_CA_1->pos.z));
                current.push_back(clipper::Coord_orth(mol_residue_C_1->pos.x,
                                                      mol_residue_C_1->pos.y,
                                                      mol_residue_C_1->pos.z));
                current.push_back(clipper::Coord_orth(mol_residue_O_1->pos.x,
                                                      mol_residue_O_1->pos.y,
                                                      mol_residue_O_1->pos.z));
                current.push_back(clipper::Coord_orth(mol_residue_CA_2->pos.x,
                                                      mol_residue_CA_2->pos.y,
                                                      mol_residue_CA_2->pos.z));
                current.push_back(clipper::Coord_orth(mol_residue_N_2->pos.x,
                                                      mol_residue_N_2->pos.y,
                                                      mol_residue_N_2->pos.z));
                
                cis_trans_init.push_back(clipper::Coord_orth(match_CA_1->pos.x,
                                                             match_CA_1->pos.y,
                                                             match_CA_1->pos.z));
                cis_trans_init.push_back(clipper::Coord_orth(match_C_1->pos.x,
                                                             match_C_1->pos.y,
                                                             match_C_1->pos.z));
                cis_trans_init.push_back(clipper::Coord_orth(match_O_1->pos.x,
                                                             match_O_1->pos.y,
                                                             match_O_1->pos.z));
                cis_trans_init.push_back(clipper::Coord_orth(match_CA_2->pos.x,
                                                             match_CA_2->pos.y,
                                                             match_CA_2->pos.z));
                cis_trans_init.push_back(clipper::Coord_orth(match_N_2->pos.x,
                                                             match_N_2->pos.y,
                                                             match_N_2->pos.z));
                
                converted.push_back(clipper::Coord_orth(conv_CA_1->pos.x,
                                                        conv_CA_1->pos.y,
                                                        conv_CA_1->pos.z));
                converted.push_back(clipper::Coord_orth(conv_C_1->pos.x,
                                                        conv_C_1->pos.y,
                                                        conv_C_1->pos.z));
                converted.push_back(clipper::Coord_orth(conv_O_1->pos.x,
                                                        conv_O_1->pos.y,
                                                        conv_O_1->pos.z));
                converted.push_back(clipper::Coord_orth(conv_CA_2->pos.x,
                                                        conv_CA_2->pos.y,
                                                        conv_CA_2->pos.z));
                converted.push_back(clipper::Coord_orth(conv_N_2->pos.x,
                                                        conv_N_2->pos.y,
                                                        conv_N_2->pos.z));
                
                // Calculate transformation
                clipper::RTop_orth lsq_mat(cis_trans_init, current);
                
                // Apply transformation to mol_residues atoms
                clipper::Coord_orth newpos;
                
                newpos = converted[0].transform(lsq_mat);
                mol_residue_CA_1->pos.x = newpos.x();
                mol_residue_CA_1->pos.y = newpos.y();
                mol_residue_CA_1->pos.z = newpos.z();
                
                newpos = converted[1].transform(lsq_mat);
                mol_residue_C_1->pos.x = newpos.x();
                mol_residue_C_1->pos.y = newpos.y();
                mol_residue_C_1->pos.z = newpos.z();
                
                newpos = converted[2].transform(lsq_mat);
                mol_residue_O_1->pos.x = newpos.x();
                mol_residue_O_1->pos.y = newpos.y();
                mol_residue_O_1->pos.z = newpos.z();
                
                newpos = converted[3].transform(lsq_mat);
                mol_residue_CA_2->pos.x = newpos.x();
                mol_residue_CA_2->pos.y = newpos.y();
                mol_residue_CA_2->pos.z = newpos.z();
                
                newpos = converted[4].transform(lsq_mat);
                mol_residue_N_2->pos.x = newpos.x();
                mol_residue_N_2->pos.y = newpos.y();
                mol_residue_N_2->pos.z = newpos.z();
                
                if (mol_residue_H_2) {
                    // Place H on N as a riding atom
                    clipper::Coord_orth at_c_pos  = coot::co_gemmi(mol_residue_C_1);
                    clipper::Coord_orth at_n_pos  = coot::co_gemmi(mol_residue_N_2);
                    clipper::Coord_orth at_ca_pos = coot::co_gemmi(mol_residue_CA_1);
                    double bl = 0.86;
                    double angle = clipper::Util::d2rad(125.0);
                    clipper::Coord_orth H_pos(at_ca_pos, at_c_pos, at_n_pos, bl, angle, M_PI);
                    mol_residue_H_2->pos.x = H_pos.x();
                    mol_residue_H_2->pos.y = H_pos.y();
                    mol_residue_H_2->pos.z = H_pos.z();
                }
                istatus = 1;
            }
        }
    }
    return istatus;
}

}} // namespace coot::util