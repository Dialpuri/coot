#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include <iomanip>
#include "function.hh"

TEST(OracleTest, AtomStreamOperator) {
    // Setup: read PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain* chain = &model->chains[0];
    ASSERT_NE(chain, nullptr);

    gemmi::Residue* res = &chain->residues[0];
    ASSERT_NE(res, nullptr);

    gemmi::Atom* atom = &res->atoms[0];
    ASSERT_NE(atom, nullptr);

    // Case 1: valid atom — stream to ostringstream and capture result
    {
        EXPECT_EQ(chain->name, "A");
        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_EQ(res->name, "MET");
        // gemmi returns unpadded atom names; pad to match MMDB convention
        EXPECT_EQ(pad_atom_name(atom->name), " N  ");

        std::ostringstream os;
        os << std::setprecision(6) << std::fixed;
        gemmi::CRA cra{chain, res, atom};
        stream_atom_gemmi(os, cra, 0);
        std::string result = os.str();

        EXPECT_EQ(result, "1/A/1/MET/ N   altLoc :: pos: (12.334000,2.772000,74.951000) B-factor: 65.050003");
    }

    // Case 2: second atom in same residue
    {
        gemmi::Atom* atom2 = &res->atoms[1];
        ASSERT_NE(atom2, nullptr);

        EXPECT_EQ(chain->name, "A");
        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_EQ(pad_atom_name(atom2->name), " CA ");

        std::ostringstream os;
        os << std::setprecision(6) << std::fixed;
        gemmi::CRA cra{chain, res, atom2};
        stream_atom_gemmi(os, cra, 0);
        std::string result = os.str();

        EXPECT_EQ(result, "1/A/1/MET/ CA  altLoc :: pos: (13.401000,3.780000,74.703000) B-factor: 69.040001");
    }

    // Case 3: atom from a later residue (residue index 9)
    {
        gemmi::Residue* res10 = &chain->residues[9];
        ASSERT_NE(res10, nullptr);

        gemmi::Atom* atom3 = &res10->atoms[0];
        ASSERT_NE(atom3, nullptr);

        EXPECT_EQ(chain->name, "A");
        EXPECT_EQ(res10->seqid.num.value, 10);
        EXPECT_EQ(res10->name, "ILE");
        EXPECT_EQ(pad_atom_name(atom3->name), " N  ");

        std::ostringstream os;
        os << std::setprecision(6) << std::fixed;
        gemmi::CRA cra{chain, res10, atom3};
        stream_atom_gemmi(os, cra, 0);
        std::string result = os.str();

        EXPECT_EQ(result, "1/A/10/ILE/ N   altLoc :: pos: (26.018000,12.491000,68.556000) B-factor: 38.990002");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}