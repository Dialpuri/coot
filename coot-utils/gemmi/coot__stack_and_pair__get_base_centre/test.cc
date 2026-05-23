#define _NO_PYTHON
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, GetBaseCentre) {
    // Load the protein PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: Try with a nucleic acid residue that has enough base atoms
    {
        // Create a test residue with base atom names
        gemmi::Residue res;
        res.name = "ADE";
        res.seqid.num = 999;
        res.seqid.icode = ' ';
        
        // Add base atoms (N1, C2, N3, C4, C5, C6 - 6 atoms needed for n_centres > 3)
        // Note: coot uses padded atom names like " N1 ", " C2 " etc.
        const char* atomNames[] = {" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "};
        const char* elements[] = {"N", "C", "N", "C", "C", "C"};
        
        for (int i = 0; i < 6; i++) {
            gemmi::Atom atom;
            atom.name = atomNames[i];
            atom.element = gemmi::Element(elements[i]);
            atom.pos.x = 10.0 + i;
            atom.pos.y = 20.0 + i;
            atom.pos.z = 30.0 + i;
            atom.occ = 1.0;
            atom.b_iso = 0.0;
            res.atoms.push_back(atom);
        }
        
        std::pair<bool, clipper::Coord_orth> result = coot::stack_and_pair::get_base_centre_gemmi(&res);
        
        EXPECT_EQ(result.first, true) << "Expected success for residue with base atoms";
        if (result.first) {
            EXPECT_NEAR(result.second.x(), 12.5, 1e-4);
            EXPECT_NEAR(result.second.y(), 22.5, 1e-4);
            EXPECT_NEAR(result.second.z(), 32.5, 1e-4);
        }
    }
    
    // Case 2: Try with a residue that doesn't have base atoms (protein residue)
    // Use ILE A 10 which exists in example.pdb
    {
        const gemmi::Residue* res = nullptr;
        for (const gemmi::Model& model : st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                for (const gemmi::Residue& r : chain.residues) {
                    if (r.name == "ILE" && r.seqid.num.value == 10 && chain.name == "A") {
                        res = &r;
                        break;
                    }
                }
                if (res) break;
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr) << "Residue ILE A 10 not found";
        
        std::pair<bool, clipper::Coord_orth> result = coot::stack_and_pair::get_base_centre_gemmi(res);
        
        EXPECT_EQ(result.first, false) << "Expected failure for protein residue without base atoms";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}