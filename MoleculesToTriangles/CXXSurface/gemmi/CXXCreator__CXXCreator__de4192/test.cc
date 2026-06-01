#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include "function.hh"

// gemmi atom.name is trimmed; MMDB GetAtomName() returns 4-char padded (PDB columns 13-16).
// PDB convention: 1-char → " X  ", 2-char → "XY  ", 3-char → "XY Z" (first letter has leading space)
static std::string pad_atom_name(const std::string& n) {
    if (n.size() == 1) return " " + n + "  ";   // " N  "
    if (n.size() == 2) return n + "  ";         // "CA  "
    if (n.size() == 3) return " " + n;           // " CA "
    return n.substr(0, 4);
}

TEST(OracleTest, CXXCreator) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Case 1: valid selection with context_selHnd=0 (default path)
    {
        // Select all atoms: "*" pattern → all atoms
        std::vector<const gemmi::Atom*> all_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    all_atoms.push_back(&atom);
                }
            }
        }

        // init() defaults
        double expected_gridSpacing = 1.2;
        double expected_probeRadius = 1.6;

        // context_selHnd == 0, so context_atoms is empty
        std::vector<const gemmi::Atom*> context_atoms;
        CXXCreator::CXXCreatorResult result = CXXCreator::CXXCreator_gemmi(model, st.cell, all_atoms, context_atoms);

        EXPECT_EQ(result.nSelAtoms, 2107);
        EXPECT_NEAR(result.gridSpacing, expected_gridSpacing, 1e-4);
        EXPECT_NEAR(result.probeRadius, expected_probeRadius, 1e-4);

        std::string first_atom_name = "nullptr";
        if (result.nSelAtoms > 0 && !result.SelAtom.empty()) {
            first_atom_name = pad_atom_name(result.SelAtom[0]->name);
        }
        EXPECT_EQ(first_atom_name, " N  ");
    }

    // Case 2: invalid selection handle (selHnd=0) → empty selection
    {
        std::vector<const gemmi::Atom*> empty_atoms;
        std::vector<const gemmi::Atom*> context_atoms;
        CXXCreator::CXXCreatorResult result = CXXCreator::CXXCreator_gemmi(model, st.cell, empty_atoms, context_atoms);

        EXPECT_EQ(result.nSelAtoms, 0);

        std::string first_atom_name = "nullptr";
        if (result.nSelAtoms > 0 && !result.SelAtom.empty()) {
            first_atom_name = pad_atom_name(result.SelAtom[0]->name);
        }
        EXPECT_EQ(first_atom_name, "nullptr");
    }

    // Case 3: valid selection with context_selHnd > 0
    {
        // Select all atoms
        std::vector<const gemmi::Atom*> all_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    all_atoms.push_back(&atom);
                }
            }
        }

        EXPECT_EQ(all_atoms.size(), 2107u);

        // context: select chain A only
        std::vector<const gemmi::Atom*> context_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    for (gemmi::Atom& atom : res.atoms) {
                        context_atoms.push_back(&atom);
                    }
                }
            }
        }

        // if (context_selHnd > 0) branch: neighbors of all_atoms, AND with context_atoms
        CXXCreator::CXXCreatorResult result = CXXCreator::CXXCreator_gemmi(model, st.cell, all_atoms, context_atoms);

        EXPECT_EQ(result.nSelAtoms, 0);

        std::string first_atom_name = "nullptr";
        if (result.nSelAtoms > 0 && !result.SelAtom.empty()) {
            first_atom_name = pad_atom_name(result.SelAtom[0]->name);
        }
        EXPECT_EQ(first_atom_name, "nullptr");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}