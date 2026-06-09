#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <gemmi/chemcomp.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_contact_indices_from_restraints__a15507/gemmi/function.hh"
#include "geometry/protein-geometry.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// Helper to deep copy a gemmi residue with a new chain ID
static gemmi::Residue deep_copy_residue(const gemmi::Residue& src, const std::string& new_chain) {
    gemmi::Residue res;
    res.name = src.name;
    res.seqid = src.seqid;
    for (const auto& atom : src.atoms) {
        gemmi::Atom a(atom);
        res.atoms.push_back(a);
    }
    return res;
}

TEST(OracleTest, chi_angles_change_by) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Find residue at chain A, seqnum 10
    gemmi::Model& model = st.models[0];
    gemmi::Residue* residue = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    residue = &res;
                    break;
                }
            }
            break;
        }
    }
    ASSERT_NE(residue, nullptr);

    // CASE 1: valid chi angle (ichi=1)
    {
        gemmi::Residue residue_copy = deep_copy_residue(*residue, "B");

        EXPECT_EQ(residue_copy.atoms.size(), 8u);

        int nth_chi = 1;
        float diff = 180.0f;

        std::vector<std::vector<int>> contact_indices;
        bool add_reverse_contacts = false;
        contact_indices = coot::util::get_contact_indices_from_restraints_gemmi(
            &residue_copy, &geom, true, add_reverse_contacts);

        std::pair<short int, float> result = coot::chi_angles::change_by_gemmi(
            nth_chi, diff, contact_indices, residue_copy);

        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 4.09672f, 1e-4);

        // Find CG1 atom and check coordinates
        double cg1_x = 0, cg1_y = 0, cg1_z = 0;
        bool found_cg1 = false;
        for (const auto& atom : residue_copy.atoms) {
            if (atom.name == "CG1") {
                cg1_x = atom.pos.x;
                cg1_y = atom.pos.y;
                cg1_z = atom.pos.z;
                found_cg1 = true;
                break;
            }
        }
        EXPECT_TRUE(found_cg1);
        EXPECT_NEAR(cg1_x, 29.1706, 1e-3);
        EXPECT_NEAR(cg1_y, 10.6852, 1e-3);
        EXPECT_NEAR(cg1_z, 68.1655, 1e-3);
    }

    // CASE 2: invalid chi angle (ichi=5 — out of range)
    {
        gemmi::Residue residue_copy = deep_copy_residue(*residue, "B");

        std::vector<std::vector<int>> contact_indices;
        bool add_reverse_contacts = false;
        contact_indices = coot::util::get_contact_indices_from_restraints_gemmi(
            &residue_copy, &geom, true, add_reverse_contacts);

        int ichi_bad = 5;
        float diff = 180.0f;

        std::pair<short int, float> result2 = coot::chi_angles::change_by_gemmi(
            ichi_bad, diff, contact_indices, residue_copy);

        EXPECT_EQ(result2.first, 2);
        EXPECT_NEAR(result2.second, 0.0f, 1e-4);
    }
}
