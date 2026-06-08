#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == seqnum) {
                return &res;
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, backrub_residue_triple_t_trim_residue_atoms_generic) {
    // Setup: load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: use_keep_atom_vector=true, empty keep vector -> deletes all atoms
    {
        gemmi::Residue* res10 = find_residue(model, "A", 10);
        gemmi::Residue* res9  = find_residue(model, "A", 9);
        gemmi::Residue* res11 = find_residue(model, "A", 11);

        EXPECT_NE(res10, nullptr);
        EXPECT_NE(res9, nullptr);
        EXPECT_NE(res11, nullptr);

        coot::util::backrub_residue_triple_t* triple = new coot::util::backrub_residue_triple_t(res9, res10, res11, "");

        int n_before = static_cast<int>(res10->atoms.size());

        std::vector<std::string> empty_keep;

        triple->trim_residue_atoms_generic_gemmi(res10, empty_keep, true);

        int n_after = static_cast<int>(res10->atoms.size());

        EXPECT_EQ(n_before, 8);
        EXPECT_EQ(n_after, 0);

        delete triple;
    }

    // Case 2: use_keep_atom_vector=true, keep only CA
    {
        gemmi::Residue* res20 = find_residue(model, "A", 20);
        gemmi::Residue* res19 = find_residue(model, "A", 19);
        gemmi::Residue* res21 = find_residue(model, "A", 21);

        if (res20 && res19 && res21) {
            coot::util::backrub_residue_triple_t* triple = new coot::util::backrub_residue_triple_t(res19, res20, res21, "");

            int n_before = static_cast<int>(res20->atoms.size());

            std::vector<std::string> keep;
            keep.push_back("CA");

            triple->trim_residue_atoms_generic_gemmi(res20, keep, true);

            int n_after = static_cast<int>(res20->atoms.size());

            EXPECT_EQ(n_before, 5);
            EXPECT_EQ(n_after, 0);

            delete triple;
        }
    }

    // Case 3: nullptr residue (guarded path, no crash)
    {
        coot::util::backrub_residue_triple_t* triple = new coot::util::backrub_residue_triple_t(nullptr, nullptr, nullptr, "");
        std::vector<std::string> keep;
        keep.push_back("CA");

        EXPECT_NO_THROW(triple->trim_residue_atoms_generic_gemmi(nullptr, keep, true));

        delete triple;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
