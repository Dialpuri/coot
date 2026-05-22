#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <vector>
#include <string>
#include <iostream>

TEST(OracleTest, beamInLinkedResidueLsqFit) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "Failed to load PDB";

    // Get two residues for testing (chain A, residue 10 and 11)
    gemmi::Model& model = st.models[0];
    gemmi::Residue* res1 = nullptr;
    gemmi::Residue* res2 = nullptr;

    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == 10) res1 = &res;
            if (res.seqid.num.value == 11) res2 = &res;
        }
    }

    ASSERT_NE(res1, nullptr) << "Could not get residue 1 (ILE 10)";
    ASSERT_NE(res2, nullptr) << "Could not get residue 2 (GLY 11)";

    // Collect atom names from res1 (reference/template)
    std::vector<std::string> atom_names_ref;
    std::vector<std::string> atom_names_match;

    int n_atoms_ref = static_cast<int>(res1->atoms.size());
    int n_atoms_mov = static_cast<int>(res2->atoms.size());

    std::cout << "Number of atoms in res1: " << n_atoms_ref << std::endl;
    std::cout << "Number of atoms in res2: " << n_atoms_mov << std::endl;

    // Use at most 4 atoms (or fewer if residue has fewer)
    int n_use = std::min(n_atoms_ref, 4);
    for (int i = 0; i < n_use; i++) {
        const auto& atom = res1->atoms[i];
        std::string atom_name = atom.name;
        atom_names_ref.push_back(atom_name);
        atom_names_match.push_back(atom_name);
        std::cout << "Using atom: " << atom_name << std::endl;
    }

    // Get initial coordinates for res2 (mov)
    double x_before = res2->atoms[0].pos.x;
    double y_before = res2->atoms[0].pos.y;
    double z_before = res2->atoms[0].pos.z;

    std::cout << "INPUT res2 initial atom[0] coords: ("
              << x_before << ", " << y_before << ", " << z_before << ")" << std::endl;
    std::cout << "INPUT atom_names_ref count: " << atom_names_ref.size() << std::endl;
    std::cout << "INPUT atom_names_match count: " << atom_names_match.size() << std::endl;

    // Call lsq_fit_gemmi - this should return true if successful
    bool status = beam_in_linked_residue_lsq_fit_gemmi(res1, res1, res2,
                                                       atom_names_ref, atom_names_match);

    std::cout << "OUTPUT lsq_fit_status: " << (status ? "true" : "false") << std::endl;

    // Get final coordinates after transform
    double x_after = res2->atoms[0].pos.x;
    double y_after = res2->atoms[0].pos.y;
    double z_after = res2->atoms[0].pos.z;

    std::cout << "OUTPUT res2 final atom[0] coords: ("
              << x_after << ", " << y_after << ", " << z_after << ")" << std::endl;

    // Test with mismatched atom names (should return false)
    std::vector<std::string> mismatched_ref;
    mismatched_ref.push_back("CA");
    std::vector<std::string> mismatched_match;
    mismatched_match.push_back("N");
    mismatched_match.push_back("CA");
    bool status2 = beam_in_linked_residue_lsq_fit_gemmi(res1, res1, res2,
                                                        mismatched_ref, mismatched_match);
    std::cout << "OUTPUT mismatched_atom_names_status: " << (status2 ? "true" : "false") << std::endl;

    // Test with valid atoms again
    bool status3 = beam_in_linked_residue_lsq_fit_gemmi(res1, res1, res2,
                                                        atom_names_ref, atom_names_match);
    std::cout << "OUTPUT second_valid_fit_status: " << (status3 ? "true" : "false") << std::endl;

    // Verify coordinates were actually transformed by checking a different atom
    double x_transformed = res2->atoms[1].pos.x;
    double y_transformed = res2->atoms[1].pos.y;
    double z_transformed = res2->atoms[1].pos.z;
    std::cout << "OUTPUT res2 atom[1] after transform: ("
              << x_transformed << ", " << y_transformed << ", " << z_transformed << ")" << std::endl;

    // Assertions based on oracle output
    EXPECT_TRUE(status) << "lsq_fit should succeed with matching atom names";
    EXPECT_EQ(status2, false) << "lsq_fit should fail with mismatched atom names";
    EXPECT_TRUE(status3) << "second valid fit should succeed";

    // Verify coordinates after transform (from oracle output)
    EXPECT_NEAR(x_transformed, 25.487, 1e-3);
    EXPECT_NEAR(y_transformed, 13.854, 1e-3);
    EXPECT_NEAR(z_transformed, 64.707, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}