#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomCoordRetrieval) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: valid atom - get CA atom from residue 10
    {
        gemmi::Atom* atom = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == 10)
                        for (auto& a : res.atoms)
                            if (a.name == "CA")
                                atom = &a;
        
        ASSERT_NE(atom, nullptr) << "atom not found";
        clipper::Coord_orth coord = coot::co_gemmi(atom);
        EXPECT_NEAR(coord.x(), 27.175, 1e-4);
        EXPECT_NEAR(coord.y(), 12.237, 1e-4);
        EXPECT_NEAR(coord.z(), 67.71, 1e-4);
    }

    // case 2: valid atom - get O atom from residue 20
    {
        gemmi::Atom* atom = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == 20)
                        for (auto& a : res.atoms)
                            if (a.name == "O")
                                atom = &a;
        
        ASSERT_NE(atom, nullptr) << "atom not found";
        clipper::Coord_orth coord = coot::co_gemmi(atom);
        EXPECT_NEAR(coord.x(), 24.652, 1e-4);
        EXPECT_NEAR(coord.y(), 8.415, 1e-4);
        EXPECT_NEAR(coord.z(), 73.577, 1e-4);
    }

    // case 3: another valid atom - get N atom from residue 30
    {
        gemmi::Atom* atom = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == 30)
                        for (auto& a : res.atoms)
                            if (a.name == "N")
                                atom = &a;
        
        ASSERT_NE(atom, nullptr) << "atom not found";
        clipper::Coord_orth coord = coot::co_gemmi(atom);
        EXPECT_NEAR(coord.x(), 25.871, 1e-4);
        EXPECT_NEAR(coord.y(), 1.9, 1e-4);
        EXPECT_NEAR(coord.z(), 72.423, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}