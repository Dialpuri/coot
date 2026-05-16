#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find an atom by chain, residue seqnum, and atom name.
static gemmi::Atom* find_atom(gemmi::Structure& st, const std::string& chain_id,
                              int res_seq_num, const std::string& atom_name) {
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value != res_seq_num) continue;
                for (auto& atom : res.atoms) {
                    if (atom.name == atom_name) return &atom;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, typed_distances_atom_type) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    using atom_type_t = coot::typed_distances::atom_type_t;

    // Case 1: Nitrogen atom -> element N, atom_type N
    {
        gemmi::Atom* at = find_atom(st, "A", 10, "N");
        EXPECT_NE(at, nullptr);
        EXPECT_EQ(at->element.name(), std::string("N"));
        EXPECT_EQ(coot::typed_distances::get_type_gemmi(at), atom_type_t::N);
    }

    // Case 2: Oxygen atom -> element O, atom_type O
    {
        gemmi::Atom* at = find_atom(st, "A", 10, "O");
        EXPECT_NE(at, nullptr);
        EXPECT_EQ(at->element.name(), std::string("O"));
        EXPECT_EQ(coot::typed_distances::get_type_gemmi(at), atom_type_t::O);
    }

    // Case 3: Carbon atom -> element C, atom_type C
    {
        gemmi::Atom* at = find_atom(st, "A", 10, "CA");
        EXPECT_NE(at, nullptr);
        EXPECT_EQ(at->element.name(), std::string("C"));
        EXPECT_EQ(coot::typed_distances::get_type_gemmi(at), atom_type_t::C);
    }

    // Case 4: Sulfur atom (SG in MET residue 1) -> no atom found
    {
        gemmi::Atom* at = find_atom(st, "A", 1, "SG");
        EXPECT_EQ(at, nullptr);
    }

    // Case 5: Phosphorus atom -> no atom found
    {
        gemmi::Atom* at = find_atom(st, "A", 1, "P");
        EXPECT_EQ(at, nullptr);
    }

    // Case 6: Invalid residue -> no atom found
    {
        gemmi::Atom* at = find_atom(st, "A", 9999, "XYZ");
        EXPECT_EQ(at, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}