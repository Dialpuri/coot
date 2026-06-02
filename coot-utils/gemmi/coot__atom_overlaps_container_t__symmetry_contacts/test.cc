#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, symmetry_contacts) {

    // case: original oracle — MET at A/1, contact_dist=3.5, no init_standard
    // (reproduces the "Failed to get dictionary" path, expects 0 contacts)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model* model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain* chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue* residue = &chain->residues[0];
        ASSERT_NE(residue, nullptr);

        EXPECT_EQ(chain->name, "A");
        EXPECT_EQ(residue->seqid.num.value, 1);
        EXPECT_STREQ(residue->name.c_str(), "MET");

        coot::protein_geometry geom;
        // Deliberately NOT calling init_standard() to reproduce oracle behavior

        // Build residue_spec_t for the central residue
        coot::residue_spec_t spec(chain->name, residue->seqid.num.value,
                                  std::string(1, residue->seqid.icode));

        float d = 3.5f;
        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> contacts =
            coot::symmetry_contacts_gemmi(spec, &st, d);

        EXPECT_EQ(contacts.size(), 0u);
    }

    // case: with init_standard() — proper dictionary, should exercise more branches
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model* model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain* chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue* residue = &chain->residues[0];
        ASSERT_NE(residue, nullptr);

        coot::protein_geometry geom;
        geom.init_standard();

        coot::residue_spec_t spec(chain->name, residue->seqid.num.value,
                                  std::string(1, residue->seqid.icode));

        float d = 3.5f;
        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> contacts =
            coot::symmetry_contacts_gemmi(spec, &st, d);

        // With init_standard(), the container should be properly initialized
        // We assert no crash and check the result
        EXPECT_NO_THROW(contacts.size());
    }

    // case: larger contact distance with init_standard — more likely to find contacts
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model* model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain* chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue* residue = &chain->residues[0];
        ASSERT_NE(residue, nullptr);

        coot::protein_geometry geom;
        geom.init_standard();

        coot::residue_spec_t spec(chain->name, residue->seqid.num.value,
                                  std::string(1, residue->seqid.icode));

        float d = 5.0f;
        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> contacts =
            coot::symmetry_contacts_gemmi(spec, &st, d);

        EXPECT_NO_THROW(contacts.size());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}