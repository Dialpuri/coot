#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find an atom by chain name, resnum, and atom name
static const gemmi::Atom* find_atom(const gemmi::Structure& st,
                                    const std::string& chain_name,
                                    int resnum,
                                    const std::string& atom_name) {
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value != resnum) continue;
                for (const auto& atom : res.atoms) {
                    if (atom.name == atom_name) return &atom;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, is_hydrogen_p_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Case 1: Hydrogen atom (should return true/1)
    {
        std::string cid = "//A/26/H";
        const gemmi::Atom* atom = find_atom(st, "A", 26, "H");
        bool result = coot::is_hydrogen_p_gemmi(atom);
        std::cout << "INPUT  cid: " << cid << std::endl;
        std::cout << "OUTPUT is_hydrogen: " << (result ? "true" : "false") << std::endl;
        std::cout << "OUTPUT atom_element: " << (atom ? atom->element.name() : "nullptr") << std::endl;
        EXPECT_TRUE(result);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(std::string(atom->element.name()), "H");
    }

    // Case 2: Non-hydrogen atom (should return false/0)
    {
        std::string cid = "//A/26/CA";
        const gemmi::Atom* atom = find_atom(st, "A", 26, "CA");
        bool result = coot::is_hydrogen_p_gemmi(atom);
        std::cout << "INPUT  cid: " << cid << std::endl;
        std::cout << "OUTPUT is_hydrogen: " << (result ? "true" : "false") << std::endl;
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}