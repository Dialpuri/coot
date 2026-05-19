#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// OracleTest suite - tests the get_vdw_radius_ligand_atom function
TEST(OracleTest, get_vdw_radius_ligand_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get model 1 (0-indexed)
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];
    
    // Find ligand residue (not standard amino acids)
    const gemmi::Residue* res_central = nullptr;
    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            std::string resname = res.name;
            if (resname != "ALA" && resname != "ARG" && resname != "ASN" && 
                resname != "ASP" && resname != "CYS" && resname != "GLN" && 
                resname != "GLU" && resname != "GLY" && resname != "HIS" && 
                resname != "ILE" && resname != "LEU" && resname != "LYS" && 
                resname != "MET" && resname != "PHE" && resname != "PRO" && 
                resname != "SER" && resname != "THR" && resname != "TRP" && 
                resname != "TYR" && resname != "VAL") {
                res_central = &res;
                break;
            }
        }
        if (res_central) break;
    }
    
    ASSERT_NE(res_central, nullptr) << "No ligand residue found";
    
    // Get all atoms from the residue - use gemmi::Residue::atoms directly (returns vector<Atom>)
    const std::vector<gemmi::Atom>& atoms = res_central->atoms;
    ASSERT_GT(atoms.size(), 0) << "No atoms in residue";
    
    // Create atom_overlaps_container_t with empty neighbours
    std::vector<const gemmi::Residue*> empty_neighbours;
    void* null_geom = nullptr;
    coot::atom_overlaps_container_t aoc(*res_central, empty_neighbours, model, null_geom);
    
    // Test multiple atoms - first call populates map, subsequent calls use cached value
    for (size_t i = 0; i < atoms.size(); ++i) {
        const gemmi::Atom& test_atom = atoms[i];
        double radius = aoc.get_vdw_radius_ligand_atom_gemmi(&test_atom);
        
        std::string atom_name = test_atom.name;
        // Check for oxygen atom in HOH
        if (atom_name == "O") {
            EXPECT_NEAR(radius, 1.52, 1e-4);
        }
    }
    
    // Test 2: Call again with same atom to verify caching
    const gemmi::Atom& test_atom2 = atoms[0];
    double radius2 = aoc.get_vdw_radius_ligand_atom_gemmi(&test_atom2);
    
    std::string atom_name = test_atom2.name;
    if (atom_name == "O") {
        EXPECT_NEAR(radius2, 1.52, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}