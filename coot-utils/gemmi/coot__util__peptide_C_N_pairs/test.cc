#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, peptide_C_N_pairs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model* model = &st.models[0];
    EXPECT_NE(model, nullptr);

    // Helper: find chain by name (gemmi::Model has no find_chain_by_name)
    auto find_chain = [](gemmi::Model& m, const std::string& name) -> gemmi::Chain* {
        for (auto& ch : m.chains) {
            if (ch.name == name) return &ch;
        }
        return nullptr;
    };

    // Helper: find residue containing an atom (gemmi::Atom has no parent pointer)
    auto find_residue_of_atom = [](gemmi::Chain& chain, gemmi::Atom* atom) -> gemmi::Residue* {
        for (auto& res : chain.residues) {
            for (auto& a : res.atoms) {
                if (&a == atom) return &res;
            }
        }
        return nullptr;
    };

    // Case 1: valid chain A
    {
        gemmi::Chain* chain = find_chain(*model, "A");
        EXPECT_NE(chain, nullptr);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> pairs =
            coot::util::peptide_C_N_pairs_gemmi(chain);

        EXPECT_EQ(pairs.size(), 263u);

        // Verify some key pairs
        auto make_label = [&](gemmi::Atom* c_atom, gemmi::Atom* n_atom) {
            gemmi::Residue* c_res = find_residue_of_atom(*chain, c_atom);
            gemmi::Residue* n_res = find_residue_of_atom(*chain, n_atom);
            // Normalize insertion code: gemmi uses ' ' for no code, MMDB uses ""
            std::string c_resid = std::to_string(c_res->seqid.num.value);
            if (c_res->seqid.icode != ' ') c_resid += c_res->seqid.icode;
            std::string n_resid = std::to_string(n_res->seqid.num.value);
            if (n_res->seqid.icode != ' ') n_resid += n_res->seqid.icode;
            return std::string("C(") + c_res->name + " " + c_resid + ")-N(" + n_res->name + " " + n_resid + ")";
        };

        EXPECT_EQ(make_label(pairs[0].first, pairs[0].second), "C(MET 1)-N(GLU 2)");
        EXPECT_EQ(make_label(pairs[1].first, pairs[1].second), "C(GLU 2)-N(ASN 3)");
        EXPECT_EQ(make_label(pairs[2].first, pairs[2].second), "C(ASN 3)-N(PHE 4)");
        EXPECT_EQ(make_label(pairs[262].first, pairs[262].second), "C(ARG 297)-N(LEU 298)");

        // Check the gap around residue 36-46
        // pair_34 is C(ILE 35)-N(ARG 36), then pair_35 jumps to C(SER 46)-N(THR 47)
        EXPECT_EQ(make_label(pairs[34].first, pairs[34].second), "C(ILE 35)-N(ARG 36)");
        EXPECT_EQ(make_label(pairs[35].first, pairs[35].second), "C(SER 46)-N(THR 47)");
    }

    // Case 2: non-existent chain
    {
        gemmi::Chain* chain = find_chain(*model, "Z");
        EXPECT_EQ(chain, nullptr);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> pairs =
            coot::util::peptide_C_N_pairs_gemmi(chain);
        EXPECT_EQ(pairs.size(), 0u);
    }

    // Case 3: chain B (exists but has no peptide C-N pairs)
    {
        gemmi::Chain* chain = find_chain(*model, "B");
        EXPECT_NE(chain, nullptr);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> pairs =
            coot::util::peptide_C_N_pairs_gemmi(chain);
        EXPECT_EQ(pairs.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}