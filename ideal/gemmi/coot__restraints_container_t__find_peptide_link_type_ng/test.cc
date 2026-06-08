#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, find_peptide_link_type_ng) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1u) << "Failed to load PDB";

    // Init protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Get residues from chain A (first three residues)
    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr) << "Chain A not found";

    ASSERT_GE(chain->residues.size(), 3u) << "Not enough residues in chain A";

    const gemmi::Residue& res_ala = chain->residues[0];  // residue 1 (MET)
    const gemmi::Residue& res_val = chain->residues[1];  // residue 2 (GLU)
    const gemmi::Residue& res_pro = chain->residues[2];  // residue 3 (ASN/PRO)

    // Case 1: GLU/2 -> ASN/3 (should be TRANS)
    {
        std::string link_type = coot::restraints_container_t::find_peptide_link_type_ng_gemmi(res_val, res_pro, geom);
        EXPECT_EQ(link_type, "TRANS");
    }

    // Case 2: MET/1 -> GLU/2 (standard peptide, should be TRANS)
    {
        std::string link_type = coot::restraints_container_t::find_peptide_link_type_ng_gemmi(res_ala, res_val, geom);
        EXPECT_EQ(link_type, "TRANS");
    }

    // Case 3: MET/1 -> MET/1 (same residue, edge case)
    {
        std::string link_type = coot::restraints_container_t::find_peptide_link_type_ng_gemmi(res_ala, res_ala, geom);
        EXPECT_EQ(link_type, "CIS");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
