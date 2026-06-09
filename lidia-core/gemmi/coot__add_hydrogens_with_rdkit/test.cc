#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <coot/api/molecules-container.hh>
#include <coot/lidia-core/rdkit-interface.hh>
#include <coot/geometry/protein-geometry.hh>
#include <cstdlib>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

// Helper to find a residue by chain ID and sequence number in a gemmi::Model
static gemmi::CRA find_residue_in_model(gemmi::Model &model, char chain_id, int seq_num) {
    for (auto &chain : model.chains) {
        if (chain.name == std::string(1, chain_id)) {
            for (auto &res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, add_hydrogens_with_rdkit) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol_pdb = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_GE(imol_pdb, 0);

    // Get restraints for ALA
    int imol_enc = mc.get_imol_enc_any();
    std::pair<bool, coot::dictionary_residue_restraints_t> restraint_pair =
        mc.get_geometry().get_monomer_restraints("ALA", imol_enc);
    ASSERT_TRUE(restraint_pair.first);

    const coot::dictionary_residue_restraints_t &restraints = restraint_pair.second;

    // Load structure with gemmi for our _gemmi function
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];

    // Case 1: PDB residue with existing hydrogens — triggers guard clause
    {
        gemmi::CRA cra = find_residue_in_model(model, 'A', 21);
        ASSERT_NE(cra.residue, nullptr);
        EXPECT_STREQ(cra.residue->name.c_str(), "ALA");
        EXPECT_EQ(cra.residue->seqid.num.value, 21);
        EXPECT_EQ(cra.residue->atoms.size(), 10u);

        std::pair<bool, std::string> result =
            coot::add_hydrogens_with_rdkit_gemmi(cra, restraints);

        EXPECT_FALSE(result.first);
        EXPECT_EQ(result.second, "Ligand contains (some) hydrogens already");
    }

    // Case 2: residue at //A/178 — also triggers guard clause if it has H
    {
        gemmi::CRA cra_hoh = find_residue_in_model(model, 'A', 178);
        if (!cra_hoh.residue) {
            // Not found in PDB — skip
            return;
        }
        EXPECT_STREQ(cra_hoh.residue->name.c_str(), "LYS");
        EXPECT_EQ(cra_hoh.residue->seqid.num.value, 178);
        EXPECT_EQ(cra_hoh.residue->atoms.size(), 22u);

        std::pair<bool, std::string> result =
            coot::add_hydrogens_with_rdkit_gemmi(cra_hoh, restraints);

        EXPECT_FALSE(result.first);
        EXPECT_EQ(result.second, "Ligand contains (some) hydrogens already");
    }

    // Case 3: nullptr residue — verifies outer guard clause
    {
        gemmi::CRA null_cra{nullptr, nullptr, nullptr};
        std::pair<bool, std::string> result =
            coot::add_hydrogens_with_rdkit_gemmi(null_cra, restraints);
        EXPECT_FALSE(result.first);
        EXPECT_EQ(result.second, "");
    }
}
