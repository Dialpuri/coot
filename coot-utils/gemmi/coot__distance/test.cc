#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DistanceFunction) {
    // Load the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Helper to get atom from chain A, residue 10/11, atom name
    auto get_atom = [&st](const std::string& res_num, const std::string& atom_name) -> const gemmi::Atom* {
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        std::string seq_str = std::to_string(res.seqid.num.value);
                        if (seq_str == res_num) {
                            for (auto& atom : res.atoms) {
                                if (atom.name == atom_name) {
                                    return &atom;
                                }
                            }
                        }
                    }
                }
            }
        }
        return nullptr;
    };

    // case 1: two atoms in the same residue (valid atoms)
    {
        const gemmi::Atom* at_1 = get_atom("10", "CA");
        const gemmi::Atom* at_2 = get_atom("10", "CA");
        double dist = coot::distance_gemmi(at_1, at_2);
        EXPECT_EQ(dist, 0);
    }

    // case 2: two different atoms in the same residue
    {
        const gemmi::Atom* at_1 = get_atom("10", "CA");
        const gemmi::Atom* at_2 = get_atom("10", "N");
        double dist = coot::distance_gemmi(at_1, at_2);
        EXPECT_NEAR(dist, 1.45564, 1e-4);
    }

    // case 3: atoms from different residues
    {
        const gemmi::Atom* at_1 = get_atom("10", "CA");
        const gemmi::Atom* at_2 = get_atom("11", "CA");
        double dist = coot::distance_gemmi(at_1, at_2);
        EXPECT_NEAR(dist, 3.80553, 1e-4);
    }

    // case 4: null atom (should return -1)
    {
        const gemmi::Atom* at_1 = get_atom("10", "CA");
        const gemmi::Atom* at_2 = nullptr;
        double dist = coot::distance_gemmi(at_1, at_2);
        EXPECT_EQ(dist, -1);
    }

    // case 5: both atoms null
    {
        const gemmi::Atom* at_1 = nullptr;
        const gemmi::Atom* at_2 = nullptr;
        double dist = coot::distance_gemmi(at_1, at_2);
        EXPECT_EQ(dist, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}