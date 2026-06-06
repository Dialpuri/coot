#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include "function.hh"

// Helper: find a residue by chain name and sequence number, return CRA
static gemmi::CRA find_residue(const gemmi::Model& model, const std::string& chain_name, int seqnum) {
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_name) continue;
        for (const gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == seqnum) {
                const gemmi::Atom* first_atom = nullptr;
                if (!res.atoms.empty()) {
                    first_atom = &res.atoms.front();
                }
                return gemmi::CRA{const_cast<gemmi::Chain*>(&chain),
                                  const_cast<gemmi::Residue*>(&res),
                                  const_cast<gemmi::Atom*>(first_atom)};
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

TEST(OracleTest, get_rotamer_name) {
    // Load structure using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: A typical residue with a side chain (ILE at chain A, seqnum 10)
    {
        gemmi::CRA cra = find_residue(model, "A", 10);
        ASSERT_NE(cra.chain, nullptr);
        ASSERT_NE(cra.residue, nullptr);

        std::string chain_id = cra.chain->name;
        int seqnum = cra.residue->seqid.num.value;
        std::string resname = cra.residue->name;

        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 10);
        EXPECT_EQ(resname, "ILE");

        std::string rotamer_name = coot::side_chain_densities::get_rotamer_name_gemmi(cra);

        EXPECT_EQ(rotamer_name, "pt");
    }

    // Case 2: VAL (seqnum 7)
    {
        gemmi::CRA cra = find_residue(model, "A", 7);
        ASSERT_NE(cra.chain, nullptr);
        ASSERT_NE(cra.residue, nullptr);

        std::string chain_id = cra.chain->name;
        int seqnum = cra.residue->seqid.num.value;
        std::string resname = cra.residue->name;

        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 7);
        EXPECT_EQ(resname, "VAL");

        std::string rotamer_name = coot::side_chain_densities::get_rotamer_name_gemmi(cra);

        EXPECT_EQ(rotamer_name, "t");
    }

    // Case 3: HOH water residue in chain B (non-standard — should return "none")
    {
        gemmi::CRA cra = find_residue(model, "B", 1);
        ASSERT_NE(cra.chain, nullptr);
        ASSERT_NE(cra.residue, nullptr);

        std::string chain_id = cra.chain->name;
        int seqnum = cra.residue->seqid.num.value;
        std::string resname = cra.residue->name;

        EXPECT_EQ(chain_id, "B");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(resname, "HOH");

        std::string rotamer_name = coot::side_chain_densities::get_rotamer_name_gemmi(cra);

        EXPECT_EQ(rotamer_name, "none");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}