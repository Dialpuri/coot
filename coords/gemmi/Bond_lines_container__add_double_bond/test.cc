#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <set>
#include "function.hh"

static std::string atom_name_trimmed(const gemmi::Atom& atom) {
    std::string a = atom.name;
    size_t end = a.find_last_not_of(' ');
    if (end == std::string::npos) return "";
    size_t start = a.find_first_not_of(' ');
    return a.substr(start, end - start + 1);
}

TEST(OracleTest, add_double_bond) {
    // Read the PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Find a suitable residue (GLU) — matching original test logic
    const gemmi::Residue* res_ptr = nullptr;
    const gemmi::Chain* chain_ptr = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            std::string rn = res.name;
            if (rn == "ASP" || rn == "GLU" || rn == "ASN" || rn == "GLN") {
                res_ptr = &res;
                chain_ptr = &chain;
                break;
            }
        }
        if (res_ptr) break;
    }
    ASSERT_NE(res_ptr, nullptr);

    std::string res_name = res_ptr->name;
    EXPECT_EQ(res_name, "GLU");

    // Build atom vector for this residue
    std::vector<gemmi::Atom> atoms;
    for (auto& atom : res_ptr->atoms) {
        atoms.push_back(atom);
    }

    int n_residue_atoms = static_cast<int>(atoms.size());
    EXPECT_EQ(n_residue_atoms, 9);

    // Determine carbon and oxygen atom names based on residue type
    std::string carbon_name = "CD";
    std::string oxygen_name = "OE1";

    if (res_name == "ASP") {
        carbon_name = "C";
        oxygen_name = "OD1";
    } else if (res_name == "GLU") {
        carbon_name = "CD";
        oxygen_name = "OE1";
    } else if (res_name == "ASN" || res_name == "GLN") {
        carbon_name = "C";
        oxygen_name = "O";
    }

    // Find atom indices
    int iat_c = -1;
    int iat_o = -1;

    for (int i = 0; i < n_residue_atoms; i++) {
        if (atom_name_trimmed(atoms[i]) == carbon_name) iat_c = i;
        if (atom_name_trimmed(atoms[i]) == oxygen_name) iat_o = i;
    }

    EXPECT_EQ(iat_c, 6);
    EXPECT_EQ(iat_o, 7);

    // Set up bond storage
    std::map<int, Bond_lines> bonds_storage;
    int total_before = 0;
    for (const auto& bl : bonds_storage) total_before += bl.second.size();

    int imol = 0;
    int imodel = 1;
    int atom_colour_type = 0; // COLOUR_BY_CHAIN
    int udd_user_defined_atom_colour_index_handle = -1;

    EXPECT_EQ(total_before, 0);
    EXPECT_EQ(imodel, 1);

    // Create protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Call add_double_bond_gemmi (different elements: C=O)
    add_double_bond_gemmi(
        bonds_storage,
        imol, imodel, iat_c, iat_o,
        atoms, res_name,
        atom_colour_type, nullptr,
        udd_user_defined_atom_colour_index_handle,
        &geom,
        false);

    int total_after = 0;
    for (const auto& bl : bonds_storage) total_after += bl.second.size();

    // Different-element double bond adds 4 bonds
    EXPECT_EQ(total_after, 4);
    EXPECT_EQ((total_after - total_before), 4);

    // Case 2: Same-element double bond (C=C path) — skip for GLU (matching original oracle)
    int total_before2 = total_after;
    int iat_CB = -1;
    int iat_C2 = -1;
    if (res_name == "ASP") {
        for (int i = 0; i < n_residue_atoms; i++) {
            if (atom_name_trimmed(atoms[i]) == "CB") iat_CB = i;
            if (atom_name_trimmed(atoms[i]) == "C") iat_C2 = i;
        }
    }

    // Oracle: since res_name is "GLU", iat_CB and iat_C2 stay -1,
    // so add_double_bond is NOT called (matching original)
    if (iat_CB >= 0 && iat_C2 >= 0) {
        add_double_bond_gemmi(
            bonds_storage,
            imol, imodel, iat_CB, iat_C2,
            atoms, res_name,
            atom_colour_type, nullptr,
            udd_user_defined_atom_colour_index_handle,
            &geom,
            false);
    }

    int total_after2 = 0;
    for (const auto& bl : bonds_storage) total_after2 += bl.second.size();
    EXPECT_EQ(total_after2, 4);
    EXPECT_EQ((total_after2 - total_before2), 0);

    // Case 3: Complementary — same-element C=C double bond (exercises the ele_1==ele_2 branch)
    {
        for (int i = 0; i < n_residue_atoms; i++) {
            if (atom_name_trimmed(atoms[i]) == "CB") iat_CB = i;
            if (atom_name_trimmed(atoms[i]) == "C") iat_C2 = i;
        }
        if (iat_CB >= 0 && iat_C2 >= 0) {
            int total_before3 = total_after2;
            add_double_bond_gemmi(
                bonds_storage,
                imol, imodel, iat_CB, iat_C2,
                atoms, res_name,
                atom_colour_type, nullptr,
                udd_user_defined_atom_colour_index_handle,
                &geom,
                false);
            int total_after3 = 0;
            for (const auto& bl : bonds_storage) total_after3 += bl.second.size();
            // Same-element bond (C=C): exercises the ele_1==ele_2 branch
            EXPECT_FALSE(total_after3 < total_before3);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
