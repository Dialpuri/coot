#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_link_plane) {
    // Case 1: Valid link_type "TRANS" for adjacent peptide residues
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Get two adjacent residues from chain A
    gemmi::Residue* res1_ptr = nullptr;
    gemmi::Residue* res2_ptr = nullptr;
    gemmi::Chain* chain_ptr = nullptr;

    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& residue : chain.residues) {
                if (residue.seqid.num.value == 1) {
                    res1_ptr = &residue;
                    chain_ptr = &chain;
                }
                if (residue.seqid.num.value == 2) {
                    res2_ptr = &residue;
                }
            }
        }
    }

    ASSERT_NE(res1_ptr, nullptr);
    ASSERT_NE(res2_ptr, nullptr);

    gemmi::CRA cra1{chain_ptr, res1_ptr, nullptr};
    gemmi::CRA cra2{chain_ptr, res2_ptr, nullptr};

    std::vector<gemmi_plane_restraint_data_t> planes_out;

    // Count restraints before add_link_plane
    int n_restraints_before = static_cast<int>(planes_out.size());

    // Call add_link_plane with link_type "TRANS" for adjacent peptide residues
    int n_plane = coot::add_link_plane_gemmi("TRANS", cra1, cra2, 0, 0, geom, planes_out);

    int n_restraints_after = static_cast<int>(planes_out.size());

    EXPECT_EQ(n_plane, 1);
    EXPECT_EQ(n_restraints_before, 0);
    EXPECT_EQ(n_restraints_after, 1);

    // Case 2: non-existent link_type - should return 0
    int n_plane2 = coot::add_link_plane_gemmi("NONEXISTENT_LINK", cra1, cra2, 0, 0, geom, planes_out);

    EXPECT_EQ(n_plane2, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}