#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, ReplaceFragment) {
    // Load a molecule using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Get the original coordinates of a known atom (CA of residue 1 in chain A)
    double original_x = 0, original_y = 0, original_z = 0;
    
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    for (auto& atom : res.atoms) {
                        if (atom.name == "CA") {
                            original_x = atom.pos.x;
                            original_y = atom.pos.y;
                            original_z = atom.pos.z;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Create a new model with just one modified atom
    gemmi::Model model2;
    model2.chains.clear();
    
    gemmi::Chain chain;
    chain.name = "A";
    
    gemmi::Residue residue;
    residue.name = "ALA";
    residue.seqid.num = 1;
    residue.seqid.icode = ' ';
    
    gemmi::Atom atom;
    atom.name = "CA";
    atom.pos = gemmi::Position(original_x + 1.0, original_y + 1.0, original_z + 1.0);
    atom.occ = 1.0;
    atom.b_iso = 0.0;
    atom.element = gemmi::Element("C");
    atom.altloc = '\0';
    
    residue.atoms.push_back(atom);
    chain.residues.push_back(residue);
    model2.chains.push_back(chain);
    
    // Create atom_selection_container_t
    coot::atom_selection_container_t asc;
    asc.model = &model2;
    
    // Collect atoms from model2 with context (atom, chain_name, seq_num, ins_code)
    for (auto& chain : model2.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                asc.atoms_with_context.push_back(
                    std::make_tuple(&atom, chain.name, res.seqid.num.value, res.seqid.icode));
            }
        }
    }
    
    asc.UDDOldAtomIndexHandle = -1; // Don't use fast indexing
    
    // Call replace_fragment_gemmi
    int result = coot::replace_fragment_gemmi(model, asc);
    
    // Verify the return value
    EXPECT_EQ(result, 1);
    
    // Check if coordinates were updated
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    for (auto& atom : res.atoms) {
                        if (atom.name == "CA") {
                            EXPECT_FLOAT_EQ(atom.pos.x, original_x + 1.0);
                            EXPECT_FLOAT_EQ(atom.pos.y, original_y + 1.0);
                            EXPECT_FLOAT_EQ(atom.pos.z, original_z + 1.0);
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}