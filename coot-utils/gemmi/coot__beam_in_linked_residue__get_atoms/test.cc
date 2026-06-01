#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace {
// gemmi atom.name is unpadded; pad to 4 chars to match MMDB's GetAtomName format
inline std::string pad4(gemmi::Atom const& a) {
    std::string s = a.padded_name();
    s.resize(4, ' ');
    return s;
}
} // namespace

TEST(OracleTest, beam_in_linked_residue_get_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Find chain A, residue seqid 1 (MET)
    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    gemmi::Residue* res = nullptr;
    for (auto& r : chain->residues) {
        if (r.seqid.num.value == 1) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    coot::beam_in_linked_residue binder;

    // Case 1: residue with all atoms found (MET 1 in chain A has N, CA, C)
    {
        std::vector<std::string> names = {" N  ", " CA ", " C  "};
        std::vector<gemmi::Atom*> atoms = binder.get_atoms_gemmi(res, names);

        EXPECT_EQ(atoms.size(), 3u);
        EXPECT_EQ(pad4(*atoms[0]), " N  ");
        EXPECT_EQ(pad4(*atoms[1]), " CA ");
        EXPECT_EQ(pad4(*atoms[2]), " C  ");
    }

    // Case 2: residue with one non-existing atom (MET 1 has no OXT)
    {
        std::vector<std::string> names = {" N  ", "OXT ", " CA "};
        std::vector<gemmi::Atom*> atoms = binder.get_atoms_gemmi(res, names);

        EXPECT_EQ(atoms.size(), 2u);
        EXPECT_EQ(pad4(*atoms[0]), " N  ");
        EXPECT_EQ(pad4(*atoms[1]), " CA ");
    }

    // Case 3: empty atom names list
    {
        std::vector<std::string> names;
        std::vector<gemmi::Atom*> atoms = binder.get_atoms_gemmi(res, names);

        EXPECT_EQ(atoms.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}