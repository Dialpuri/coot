#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, flev_attached_hydrogens_close_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model &model = st.models[0];

    // Collect all non-HOH residues from the model
    std::vector<gemmi::Residue*> env_residues;
    for (gemmi::Chain &chain : model.chains) {
        for (gemmi::Residue &res : chain.residues) {
            if (res.name != "HOH") {
                env_residues.push_back(&res);
            }
        }
    }

    // Find chain A, residue 10
    gemmi::Residue *ref_res = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue &r : chain.residues) {
                if (r.seqid.num.value == 10) {
                    ref_res = &r;
                    break;
                }
            }
            break;
        }
    }

    gemmi::Vec3 pt1(0, 0, 0);
    bool found_ref = false;
    if (ref_res) {
        if (!ref_res->atoms.empty()) {
            pt1 = ref_res->atoms[0].pos;
            found_ref = true;
        }
    }

    // Create the object (close_atoms doesn't need dictionary)
    pli::flev_attached_hydrogens_t obj;

    // --- Case 1: Point near real atoms (inside the protein) ---
    {
        EXPECT_EQ(env_residues.size(), 267);
        std::vector<gemmi::CRA> close = obj.close_atoms_gemmi(pt1, env_residues, model);
        EXPECT_EQ(close.size(), 83);
        if (!close.empty()) {
            EXPECT_EQ(close[0].atom->name, "CA");
            EXPECT_EQ(close[0].residue->name, "GLU");
            EXPECT_EQ(close[0].residue->seqid.num.value, 8);
        }
    }

    // --- Case 2: Point far away from all atoms (should return empty) ---
    {
        gemmi::Vec3 pt2(9999, 9999, 9999);
        std::vector<gemmi::CRA> close = obj.close_atoms_gemmi(pt2, env_residues, model);
        EXPECT_EQ(close.size(), 0);
    }

    // --- Case 3: Empty residue list (should return empty) ---
    {
        std::vector<gemmi::Residue*> empty_residues;
        gemmi::Vec3 pt3(0, 0, 0);
        std::vector<gemmi::CRA> close = obj.close_atoms_gemmi(pt3, empty_residues, model);
        EXPECT_EQ(close.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}