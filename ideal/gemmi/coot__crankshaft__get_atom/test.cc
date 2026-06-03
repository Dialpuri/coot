#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: format gemmi's unpadded atom name as a 4-char PDB-style padded string
// (1 leading space for names ≤3 chars, then right-padded to 4)
static std::string pad_atom_name_pdb(const std::string &name) {
    std::string s = name;
    if (s.size() <= 3) s = ' ' + s;
    while (s.size() < 4) s += ' ';
    return s;
}

TEST(OracleTest, crankshaft_get_atom) {
    // Load PDB via gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get a residue: //A/1 (MET)
    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];
    gemmi::Residue &res = chain.residues[0];

    ASSERT_TRUE(&res != nullptr);

    std::string chain_id = chain.name;
    int seqnum = res.seqid.num.value;
    std::string resname = res.name;

    // Case 1: valid atom " CA "
    {
        std::string atom_name = " CA ";
        gemmi::Atom *atom = coot::crankshaft::get_atom_gemmi(&res, atom_name);
        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(resname, "MET");
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(pad_atom_name_pdb(atom->name), " CA ");
    }

    // Case 2: valid atom " N  "
    {
        std::string atom_name = " N  ";
        gemmi::Atom *atom = coot::crankshaft::get_atom_gemmi(&res, atom_name);
        EXPECT_NE(atom, nullptr);
        EXPECT_EQ(pad_atom_name_pdb(atom->name), " N  ");
    }

    // Case 3: non-existent atom " XXX "
    {
        std::string atom_name = " XXX ";
        gemmi::Atom *atom = coot::crankshaft::get_atom_gemmi(&res, atom_name);
        EXPECT_EQ(atom, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}