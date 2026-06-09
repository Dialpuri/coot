#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

#include "geometry/protein-geometry.hh"

// Helper: find chain name for a given residue pointer in a Structure
static std::string find_chain_name(gemmi::Structure& st, gemmi::Residue* res_ptr) {
    if (!st.models.empty()) {
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& ch : model.chains) {
            for (gemmi::Residue& res : ch.residues) {
                if (&res == res_ptr) {
                    return ch.name;
                }
            }
        }
    }
    return "";
}

TEST(OracleTest, missing_atoms) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Init geometry dictionary (REQUIRED)
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: without hydrogen checking (default/typical)
    {
        bool do_missing_hydrogen_atoms_flag = false;
        coot::util::missing_atom_info_gemmi mai = coot::util::missing_atoms_gemmi(&st, do_missing_hydrogen_atoms_flag, &geom);

        EXPECT_EQ(mai.residues_with_no_dictionary.size(), 0u);
        EXPECT_EQ(mai.residues_with_missing_atoms.size(), 9u);
        EXPECT_EQ(mai.atoms_in_coords_but_not_in_dict.size(), 0u);
        EXPECT_EQ(mai.residue_missing_atom_names_map.size(), 9u);

        // Verify first 5 residues with missing atoms and their missing atom names
        size_t max_print = std::min(mai.residues_with_missing_atoms.size(), size_t(5));
        std::vector<std::pair<std::string, std::string>> expected_results = {
            {"//A/20/LYS", " CG ,  CD ,  CE ,  NZ "},
            {"//A/109/PHE", " CG ,  CD1,  CD2,  CE1,  CE2,  CZ "},
            {"//A/126/ARG", " CG ,  CD ,  NE ,  CZ ,  NH1,  NH2"},
            {"//A/140/ALA", " CB "},
            {"//A/167/TRP", " CG ,  CD1,  CD2,  NE1,  CE2,  CE3,  CZ2,  CZ3,  CH2"}
        };

        for (size_t i = 0; i < max_print; i++) {
            gemmi::Residue *r = mai.residues_with_missing_atoms[i];
            std::string chain_id = find_chain_name(st, r);
            std::string cid = std::string("//") + chain_id + "/" + std::to_string(r->seqid.num.value) + "/" + std::string(r->name);
            EXPECT_EQ(cid, expected_results[i].first);

            if (mai.residue_missing_atom_names_map.count(r)) {
                std::vector<std::string> &names = mai.residue_missing_atom_names_map[r];
                std::string atom_list;
                for (size_t j = 0; j < names.size(); j++) {
                    if (j > 0) atom_list += ", ";
                    atom_list += names[j];
                }
                EXPECT_EQ(atom_list, expected_results[i].second);
            }
        }
    }

    // Case 2: with hydrogen checking — should find more missing atoms
    {
        bool do_missing_hydrogen_atoms_flag = true;
        coot::util::missing_atom_info_gemmi mai = coot::util::missing_atoms_gemmi(&st, do_missing_hydrogen_atoms_flag, &geom);

        EXPECT_EQ(mai.residues_with_missing_atoms.size(), 264u);
        EXPECT_EQ(mai.residue_missing_atom_names_map.size(), 264u);

        // Verify first 5 residues with missing atoms and their missing atom names
        size_t max_print = std::min(mai.residues_with_missing_atoms.size(), size_t(5));
        std::vector<std::pair<std::string, std::string>> expected_results = {
            {"//A/2/GLU", " H2 ,  H3 "},
            {"//A/3/ASN", " H2 ,  H3 "},
            {"//A/4/PHE", " H2 ,  H3 "},
            {"//A/5/GLN", " H2 ,  H3 "},
            {"//A/6/LYS", " H2 ,  H3 "}
        };

        for (size_t i = 0; i < max_print; i++) {
            gemmi::Residue *r = mai.residues_with_missing_atoms[i];
            std::string chain_id = find_chain_name(st, r);
            std::string cid = std::string("//") + chain_id + "/" + std::to_string(r->seqid.num.value) + "/" + std::string(r->name);
            EXPECT_EQ(cid, expected_results[i].first);

            if (mai.residue_missing_atom_names_map.count(r)) {
                std::vector<std::string> &names = mai.residue_missing_atom_names_map[r];
                std::string atom_list;
                for (size_t j = 0; j < names.size(); j++) {
                    if (j > 0) atom_list += ", ";
                    atom_list += names[j];
                }
                EXPECT_EQ(atom_list, expected_results[i].second);
            }
        }
    }

    // Case 3: null mol — should return empty results
    {
        coot::util::missing_atom_info_gemmi mai = coot::util::missing_atoms_gemmi(nullptr, false, &geom);

        EXPECT_EQ(mai.residues_with_no_dictionary.size(), 0u);
        EXPECT_EQ(mai.residues_with_missing_atoms.size(), 0u);
        EXPECT_EQ(mai.atoms_in_coords_but_not_in_dict.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
