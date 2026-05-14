#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include "function.hh"

std::vector<coot::molecule_t::WaterFragment> create_water_molecule() {
    std::vector<coot::molecule_t::WaterFragment> fragments;
    coot::molecule_t::WaterFragment frag;
    
    for (int i = 0; i < 3; i++) {
        coot::molecule_t::WaterResidue res;
        res.name = "HOH";
        res.seqnum = 100 + i;
        
        coot::molecule_t::WaterAtom o_atom;
        o_atom.name = "O";
        o_atom.element = "O";
        o_atom.pos = gemmi::Vec3(10.0 + i, 20.0 + i, 30.0 + i);
        o_atom.temperature_factor = 20.0;
        o_atom.altLoc = "";
        
        coot::molecule_t::WaterAtom h1_atom;
        h1_atom.name = "H1";
        h1_atom.element = "H";
        h1_atom.pos = gemmi::Vec3(11.0 + i, 21.0 + i, 31.0 + i);
        h1_atom.temperature_factor = 20.0;
        h1_atom.altLoc = "";
        
        coot::molecule_t::WaterAtom h2_atom;
        h2_atom.name = "H2";
        h2_atom.element = "H";
        h2_atom.pos = gemmi::Vec3(9.0 + i, 19.0 + i, 29.0 + i);
        h2_atom.temperature_factor = 20.0;
        h2_atom.altLoc = "";
        
        res.atoms.push_back(o_atom);
        res.atoms.push_back(h1_atom);
        res.atoms.push_back(h2_atom);
        
        frag.residues.push_back(res);
    }
    
    fragments.push_back(frag);
    return fragments;
}

TEST(OracleTest, InsertWatersIntoMolecule) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    std::vector<coot::molecule_t::WaterFragment> water_mol = create_water_molecule();
    
    int result = coot::molecule_t::insert_waters_into_molecule_gemmi(st, water_mol, "HOH");
    
    EXPECT_EQ(result, 0);
    
    int final_waters = 0;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "HOH") {
                    final_waters++;
                }
            }
        }
    }
    
    EXPECT_EQ(final_waters, 10);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}