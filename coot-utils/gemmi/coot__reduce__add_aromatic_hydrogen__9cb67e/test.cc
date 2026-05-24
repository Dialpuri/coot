#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AddAromaticHydrogen) {
    // Create a simple PDB string for PHE residue with aromatic ring
    std::string pdb_content = 
        "ATOM      1  N   PHE A  10      10.000  10.000  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  PHE A  10      11.400  10.000  10.000  1.00 20.00           C  \n"
        "ATOM      3  C   PHE A  10      12.000  11.200  10.000  1.00 20.00           C  \n"
        "ATOM      4  CB  PHE A  10      11.800   9.500   8.600  1.00 20.00           C  \n"
        "ATOM      5  CG  PHE A  10      13.300   9.200   8.500  1.00 20.00           C  \n"
        "ATOM      6  CD1 PHE A  10      14.100   9.800   9.400  1.00 20.00           C  \n"
        "ATOM      7  CD2 PHE A  10      13.800   8.300   7.500  1.00 20.00           C  \n"
        "ATOM      8  CE1 PHE A  10      15.400   9.400   9.300  1.00 20.00           C  \n"
        "ATOM      9  CE2 PHE A  10      15.100   7.900   7.400  1.00 20.00           C  \n"
        "ATOM     10  CZ  PHE A  10      16.000   8.500   8.300  1.00 20.00           C  \n";
    
    // Write to temp file
    FILE *f = fopen("/tmp/test_phe.pdb", "w");
    if (f) {
        fwrite(pdb_content.c_str(), 1, pdb_content.size(), f);
        fclose(f);
    }
    
    // Load the PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/tmp/test_phe.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB";
    
    // Get the first model
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty()) << "No chains found";
    
    // Get chain A
    gemmi::Chain& chain = model.chains[0];
    EXPECT_EQ(chain.name, "A");
    
    // Get the PHE residue
    gemmi::Residue* res = nullptr;
    for (auto& r : chain.residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr) << "Failed to get residue with seqnum 10";
    
    // Check residue info
    EXPECT_EQ(res->name, "PHE");
    EXPECT_EQ(res->seqid.num.value, 10);
    
    // Check if aromatic ring atoms exist
    gemmi::Atom* at_cg = nullptr;
    gemmi::Atom* at_cd1 = nullptr;
    gemmi::Atom* at_ce1 = nullptr;
    gemmi::Atom* at_cz = nullptr;
    
    for (auto& atom : res->atoms) {
        if (atom.name == "CG") at_cg = &atom;
        if (atom.name == "CD1") at_cd1 = &atom;
        if (atom.name == "CE1") at_ce1 = &atom;
        if (atom.name == "CZ") at_cz = &atom;
    }
    
    EXPECT_NE(at_cg, nullptr);
    EXPECT_NE(at_cd1, nullptr);
    EXPECT_NE(at_ce1, nullptr);
    EXPECT_NE(at_cz, nullptr);
    
    // Count existing hydrogens before
    int hydrogens_before = 0;
    for (auto& atom : res->atoms) {
        if (atom.name.length() >= 1 && atom.name[0] == 'H') {
            hydrogens_before++;
        }
    }
    EXPECT_EQ(hydrogens_before, 0);
    
    // Test case 1: Add HD1 hydrogen using aromatic ring atoms
    coot::reduce::add_aromatic_hydrogen_gemmi("HD1", "CG", "CD1", "CE1", 1.08, *res);
    
    // Count hydrogens after first addition
    int hydrogens_after_1 = 0;
    for (auto& atom : res->atoms) {
        if (atom.name.length() >= 1 && atom.name[0] == 'H') {
            hydrogens_after_1++;
        }
    }
    EXPECT_EQ(hydrogens_after_1, 1);
    
    // Test case 2: Add HE1 hydrogen
    coot::reduce::add_aromatic_hydrogen_gemmi("HE1", "CD1", "CE1", "CZ", 1.08, *res);
    
    // Count hydrogens after second addition
    int hydrogens_after_2 = 0;
    for (auto& atom : res->atoms) {
        if (atom.name.length() >= 1 && atom.name[0] == 'H') {
            hydrogens_after_2++;
        }
    }
    EXPECT_EQ(hydrogens_after_2, 2);
    
    // Test case 3: Invalid neighbors (should not crash)
    EXPECT_NO_THROW(coot::reduce::add_aromatic_hydrogen_gemmi("HXXX", "NONEXIST1", "NONEXIST2", "NONEXIST3", 1.08, *res));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}