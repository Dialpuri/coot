#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, AppendToMolecule) {
    // Load the test structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Count chains before
    int chain_count_before = 0;
    for (auto& model : st.models) {
        chain_count_before += model.chains.size();
    }
    EXPECT_EQ(chain_count_before, 2) << "Initial chain count should be 2";

    // Create a new structure with a unique chain ID "W"
    gemmi::Structure water_mol;
    water_mol.name = "water";
    
    // Add a model
    gemmi::Model model1;
    
    // Add a chain "W"
    gemmi::Chain chain_w;
    chain_w.name = "W";
    
    // Add a residue HOH at position 1
    gemmi::Residue res;
    res.name = "HOH";
    res.seqid.num = 1;
    
    // Add an oxygen atom
    gemmi::Atom atom_o;
    atom_o.name = "O";
    atom_o.element = gemmi::Element("O");
    atom_o.pos = gemmi::Position(1.0, 2.0, 3.0);
    atom_o.b_iso = 20.0;
    res.atoms.push_back(atom_o);
    
    chain_w.residues.push_back(res);
    model1.chains.push_back(chain_w);
    water_mol.models.push_back(model1);

    int result = coot::append_to_molecule_gemmi(st, water_mol);
    EXPECT_EQ(result, 0) << "append_to_molecule should return 0 for new chain";

    int chain_count_after_first = 0;
    for (auto& model : st.models) {
        chain_count_after_first += model.chains.size();
    }
    EXPECT_EQ(chain_count_after_first, 3) << "Chain count should increase to 3 after adding new chain";

    // Edge case: try adding with a chain ID that already exists ("A")
    gemmi::Structure water_mol2;
    water_mol2.name = "water2";
    
    gemmi::Model model2;
    gemmi::Chain chain_a;
    chain_a.name = "A";
    
    gemmi::Residue res2;
    res2.name = "HOH";
    res2.seqid.num = 2;
    
    gemmi::Atom atom_o2;
    atom_o2.name = "O";
    atom_o2.element = gemmi::Element("O");
    atom_o2.pos = gemmi::Position(4.0, 5.0, 6.0);
    atom_o2.b_iso = 20.0;
    res2.atoms.push_back(atom_o2);
    
    chain_a.residues.push_back(res2);
    model2.chains.push_back(chain_a);
    water_mol2.models.push_back(model2);

    int result2 = coot::append_to_molecule_gemmi(st, water_mol2);
    EXPECT_EQ(result2, 1) << "append_to_molecule should return 1 for duplicate chain ID";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}