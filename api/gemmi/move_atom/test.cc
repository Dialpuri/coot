#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include "function.hh"

TEST(OracleTest, MoveAtomFunction) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case: valid_atom
    {
        // Find residue A/10
        gemmi::Residue* res = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& r : chain.residues) {
                        if (r.seqid.num.value == 10) {
                            res = &r;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr) << "Residue not found";

        // Find the CA atom before move
        gemmi::Atom* atom = nullptr;
        for (gemmi::Atom& a : res->atoms) {
            if (a.name == "CA") {
                atom = &a;
                break;
            }
        }
        ASSERT_NE(atom, nullptr) << "Atom not found";

        double orig_x = atom->pos.x;
        double orig_y = atom->pos.y;
        double orig_z = atom->pos.z;

        EXPECT_EQ(res->seqid.num.value, 10);
        EXPECT_EQ(res->name, "ILE");

        gemmi::Vec3 new_pos(100.0, 100.0, 100.0);
        bool result = coot::move_atom_gemmi("CA", res, new_pos);
        EXPECT_TRUE(result);

        // Find the CA atom after move
        gemmi::Atom* atom_after = nullptr;
        for (gemmi::Atom& a : res->atoms) {
            if (a.name == "CA") {
                atom_after = &a;
                break;
            }
        }
        ASSERT_NE(atom_after, nullptr);
        double new_x = atom_after->pos.x;
        double new_y = atom_after->pos.y;
        double new_z = atom_after->pos.z;
        EXPECT_NEAR(new_x, 100.0, 1e-4);
        EXPECT_NEAR(new_y, 100.0, 1e-4);
        EXPECT_NEAR(new_z, 100.0, 1e-4);
        EXPECT_TRUE(std::abs(new_x - orig_x) > 1.0);
    }

    // case: nonexistent_atom
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& r : chain.residues) {
                        if (r.seqid.num.value == 10) {
                            res = &r;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr) << "Residue not found";

        gemmi::Vec3 new_pos(200.0, 200.0, 200.0);
        bool result = coot::move_atom_gemmi("NONEXISTENT", res, new_pos);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}