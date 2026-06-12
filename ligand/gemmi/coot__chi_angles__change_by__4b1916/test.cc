#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_contact_indices_from_restraints__a15507/gemmi/function.hh"
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include <vector>
#include <utility>
#include <cstdlib>

// Helper to find a residue by chain ID and sequence number in gemmi
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper to get insertion code string normalized for gemmi
static std::string get_ins_code(gemmi::Residue& res) {
    char ic = res.seqid.icode;
    if (ic == ' ') return "";
    return std::string(1, ic);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, chi_angles_change_by) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: VAL residue with chi1 — valid case using 7-param change_by
    {
        gemmi::Residue* residue = find_residue(model, "A", 30);
        ASSERT_NE(residue, nullptr);
        EXPECT_STREQ(residue->name.c_str(), "VAL");
        EXPECT_EQ(residue->seqid.num.value, 30);

        int nth_chi = 1;
        double diff = 180.0;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, &geom, true, false);

        EXPECT_EQ(contact_indices.size(), 16u);

        coot::atom_spec_t tree_base_atom;
        tree_base_atom.chain_id = "A";
        tree_base_atom.res_no = residue->seqid.num.value;
        tree_base_atom.ins_code = get_ins_code(*residue);
        tree_base_atom.atom_name = " N  ";

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(
            0,
            nth_chi,
            diff,
            contact_indices,
            &geom,
            tree_base_atom,
            static_cast<short int>(0),
            *residue);

        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 2.06627, 1e-4);
    }

    // Case 2: GLU (no chi angles comment in oracle but it has chi1) — edge case
    {
        gemmi::Residue* residue = find_residue(model, "A", 8);
        ASSERT_NE(residue, nullptr);
        EXPECT_STREQ(residue->name.c_str(), "GLU");
        EXPECT_EQ(residue->seqid.num.value, 8);

        int nth_chi = 1;
        double diff = 90.0;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, &geom, true, false);

        coot::atom_spec_t tree_base_atom;
        tree_base_atom.chain_id = "A";
        tree_base_atom.res_no = residue->seqid.num.value;
        tree_base_atom.ins_code = get_ins_code(*residue);
        tree_base_atom.atom_name = " N  ";

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(
            0,
            nth_chi,
            diff,
            contact_indices,
            &geom,
            tree_base_atom,
            static_cast<short int>(0),
            *residue);

        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 2.57723, 1e-4);
    }

    // Case 3: Use the 7-param change_by with atom_spec_t
    {
        gemmi::Residue* residue = find_residue(model, "A", 30);
        ASSERT_NE(residue, nullptr);
        EXPECT_STREQ(residue->name.c_str(), "VAL");
        EXPECT_EQ(residue->seqid.num.value, 30);

        int nth_chi = 1;
        double diff = 45.0;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, &geom, true, false);

        coot::atom_spec_t tree_base_atom;
        tree_base_atom.chain_id = "A";
        tree_base_atom.res_no = residue->seqid.num.value;
        tree_base_atom.ins_code = get_ins_code(*residue);
        tree_base_atom.atom_name = " N  ";

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(
            0,
            nth_chi,
            diff,
            contact_indices,
            &geom,
            tree_base_atom,
            static_cast<short int>(0),
            *residue);

        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 2.85167, 1e-4);
    }
}
