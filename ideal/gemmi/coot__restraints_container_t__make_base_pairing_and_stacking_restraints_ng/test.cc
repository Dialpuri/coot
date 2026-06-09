#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

TEST(OracleTest, make_base_pairing_and_stacking_restraints_ng) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: protein-only (example-hydrogen.pdb) — oracle case
    // Branch: no nucleic acid → no restraints added
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        EXPECT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        int nchains = static_cast<int>(model.chains.size());
        EXPECT_EQ(nchains, 2);

        gemmi::Chain* chain = nullptr;
        for (auto& ch : model.chains) {
            if (ch.name == "A") { chain = &ch; break; }
        }
        EXPECT_NE(chain, nullptr);

        int nresidues = static_cast<int>(chain->residues.size());
        EXPECT_EQ(nresidues, 267);

        std::vector<std::pair<bool, gemmi::Residue*>> residues_vec;
        for (auto& res : chain->residues) {
            if (!res.atoms.empty()) {
                residues_vec.push_back(std::make_pair(false, &res));
            }
        }
        EXPECT_EQ(residues_vec.size(), 267);

        int imol = 0;

        int restraints_before = 0;
        EXPECT_EQ(restraints_before, 0);

        int restraints_after = coot::make_base_pairing_and_stacking_restraints_ng_gemmi(
            st, model, residues_vec, imol, geom);

        EXPECT_EQ(restraints_after, 0);
        EXPECT_EQ((restraints_after - restraints_before), 0);
        EXPECT_EQ(imol, 0);
        // is_protein_only: true (no restraints were added because no nucleic acid)
        EXPECT_TRUE(restraints_after == restraints_before);
    }

    // Case 2: empty residues vector — also no restraints added
    // Branch: no atoms selected → trivially no nucleic acid → no restraints
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        EXPECT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        std::vector<std::pair<bool, gemmi::Residue*>> residues_vec;  // empty

        int imol = 0;

        int restraints_before = 0;

        EXPECT_NO_THROW(coot::make_base_pairing_and_stacking_restraints_ng_gemmi(
            st, model, residues_vec, imol, geom));

        int restraints_after = coot::make_base_pairing_and_stacking_restraints_ng_gemmi(
            st, model, residues_vec, imol, geom);
        EXPECT_EQ(restraints_after, restraints_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}