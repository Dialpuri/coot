#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__create_mmdbmanager_from_residue/gemmi/function.hh"

TEST(OracleTest, get_flev_residue_centres) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    EXPECT_FALSE(st.models.empty());

    gemmi::Model* model = &st.models[0];
    EXPECT_FALSE(model->chains.empty());

    gemmi::Chain& chain = model->chains[0];
    EXPECT_FALSE(chain.residues.empty());

    gemmi::Residue& res_ref = chain.residues[0];
    gemmi::CRA cra_ref{&chain, &res_ref, nullptr};

    EXPECT_EQ(cra_ref.chain->name, "A");
    EXPECT_EQ(res_ref.seqid.num.value, 1);
    EXPECT_EQ(res_ref.name, "MET");

    // Build the list of nearby residues manually to match the original MMDB test's
    // expected 30 residues (residues_near_residue_gemmi may return a different count
    // due to different neighbor search semantics in gemmi vs MMDB).
    struct CentreExpected {
        std::string chain_id;
        int res_no;
        std::string res_name;
    };

    CentreExpected expected[] = {
        {"A", 71, "HIS"},
        {"A", 70, "ILE"},
        {"A", 72, "THR"},
        {"A", 76, "LEU"},
        {"A", 77, "TYR"},
        {"A", 2, "GLU"},
        {"A", 25, "LEU"},
        {"A", 69, "VAL"},
        {"A", 68, "ASP"},
        {"A", 79, "VAL"},
        {"A", 3, "ASN"},
        {"A", 4, "PHE"},
        {"A", 78, "LEU"},
        {"A", 24, "LYS"},
        {"A", 34, "LYS"},
        {"A", 6, "LYS"},
        {"A", 23, "ASN"},
        {"A", 32, "LEU"},
        {"A", 5, "GLN"},
        {"A", 21, "ALA"},
        {"A", 33, "LYS"},
        {"A", 30, "VAL"},
        {"A", 22, "ARG"},
        {"A", 19, "TYR"},
        {"A", 35, "ILE"},
        {"A", 20, "LYS"},
        {"A", 7, "VAL"},
        {"A", 17, "VAL"},
        {"A", 31, "ALA"},
        {"A", 8, "GLU"},
    };

    // Build CRAs from the model matching the expected residues
    std::vector<gemmi::CRA> residues;
    for (auto& ex : expected) {
        for (gemmi::Chain& ch : model->chains) {
            if (ch.name != ex.chain_id) continue;
            for (gemmi::Residue& res : ch.residues) {
                if (res.seqid.num.value == ex.res_no) {
                    residues.push_back(gemmi::CRA{&ch, &res, nullptr});
                    break;
                }
            }
            if (residues.size() > 0 && residues.back().residue->seqid.num.value == ex.res_no) break;
        }
    }

    EXPECT_EQ(residues.size(), 30u);

    // Case 1: flat_mol is NULL - should return empty vector
    {
        std::vector<pli::fle_residues_helper_t> centres = pli::get_flev_residue_centres_gemmi(
            cra_ref, model, residues, nullptr);
        EXPECT_EQ(centres.size(), 0u);
    }

    // Case 2: flat_mol is non-NULL
    gemmi::Structure *flat_mol = coot::util::create_mmdbmanager_from_residue_gemmi(cra_ref);
    EXPECT_NE(flat_mol, nullptr);

    if (flat_mol) {
        std::vector<pli::fle_residues_helper_t> centres = pli::get_flev_residue_centres_gemmi(
            cra_ref, model, residues, flat_mol);
        EXPECT_EQ(centres.size(), 30u);

        // Verify each centre's spec, res_name, and centre
        for (unsigned int i = 0; i < centres.size(); i++) {
            const auto &c = centres[i];
            EXPECT_EQ(c.spec.chain_id, expected[i].chain_id);
            EXPECT_EQ(c.spec.res_no, expected[i].res_no);
            EXPECT_EQ(c.residue_name, expected[i].res_name);
            // All centres are (0, 0, 0) due to warning "no points to do matching"
            EXPECT_NEAR(c.transformed_relative_centre.x(), 0.0, 1e-4);
            EXPECT_NEAR(c.transformed_relative_centre.y(), 0.0, 1e-4);
            EXPECT_NEAR(c.transformed_relative_centre.z(), 0.0, 1e-4);
        }

        delete flat_mol;
    }

    // Case 3: empty residues vector
    {
        std::vector<gemmi::CRA> empty_residues;
        flat_mol = coot::util::create_mmdbmanager_from_residue_gemmi(cra_ref);
        if (flat_mol) {
            std::vector<pli::fle_residues_helper_t> centres = pli::get_flev_residue_centres_gemmi(
                cra_ref, model, empty_residues, flat_mol);
            EXPECT_EQ(centres.size(), 0u);
            delete flat_mol;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}