#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include "function.hh"

// gemmi::Atom.name is trimmed ("N") while MMDB GetAtomName() returns
// a 4-character padded string (" N  "). Format: leading space + name
// left-justified in 3 chars, padded with trailing spaces.
static std::string pad_atom_name(const std::string& s) {
    return std::string(1, ' ') + s + std::string(std::max(0, 3 - static_cast<int>(s.size())), ' ');
}

TEST(OracleTest, mol_by_symmetry) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    {
        // Case 1: Identity symmetry operation (no transformation)
        int n_models = static_cast<int>(st.models.size());
        EXPECT_EQ(n_models, 1);

        // Construct a default cell (P1, orthorhombic)
        gemmi::UnitCell cell(100.0, 100.0, 100.0, 90.0, 90.0, 90.0);

        EXPECT_NEAR(cell.a, 100.0, 1e-4);
        EXPECT_NEAR(cell.b, 100.0, 1e-4);
        EXPECT_NEAR(cell.c, 100.0, 1e-4);

        gemmi::FTransform ftrans_identity;
        ftrans_identity.mat = gemmi::Mat33(1,0,0,0,1,0,0,0,1);
        ftrans_identity.vec = gemmi::Vec3(0,0,0);
        std::vector<int> pre_shift_empty; // empty — no origin shift

        gemmi::Structure mol1 = coot::mol_by_symmetry_gemmi(st, cell, ftrans_identity, pre_shift_empty);

        // Get first atom coordinates from original
        gemmi::Atom* atom_orig = nullptr;
        if (!st.models.empty()) {
            gemmi::Model& mod = st.models[0];
            if (!mod.chains.empty()) {
                gemmi::Chain& ch = mod.chains[0];
                if (!ch.residues.empty()) {
                    gemmi::Residue& res = ch.residues[0];
                    if (!res.atoms.empty()) {
                        atom_orig = &res.atoms[0];
                    }
                }
            }
        }

        // Get first atom from transformed molecule
        gemmi::Atom* atom_new = nullptr;
        if (!mol1.models.empty()) {
            gemmi::Model& mod2 = mol1.models[0];
            if (!mod2.chains.empty()) {
                gemmi::Chain& ch2 = mod2.chains[0];
                if (!ch2.residues.empty()) {
                    gemmi::Residue& res2 = ch2.residues[0];
                    if (!res2.atoms.empty()) {
                        atom_new = &res2.atoms[0];
                    }
                }
            }
        }

        if (atom_orig && atom_new) {
            EXPECT_EQ(pad_atom_name(atom_orig->name), " N  ");
            EXPECT_NEAR(atom_orig->pos.x, 12.334, 1e-4);
            EXPECT_NEAR(atom_orig->pos.y, 2.772, 1e-4);
            EXPECT_NEAR(atom_orig->pos.z, 74.951, 1e-4);

            EXPECT_EQ(pad_atom_name(atom_new->name), " N  ");
            EXPECT_NEAR(atom_new->pos.x, 12.334, 1e-4);
            EXPECT_NEAR(atom_new->pos.y, 2.772, 1e-4);
            EXPECT_NEAR(atom_new->pos.z, 74.951, 1e-4);
        }

        EXPECT_EQ(mol1.models.size(), 1u);

        // No delete needed — gemmi::Structure is value-type
    }

    {
        // Case 2: Translation by (0.5, 0.5, 0.5) fractional
        gemmi::UnitCell cell(100.0, 100.0, 100.0, 90.0, 90.0, 90.0);

        gemmi::FTransform ftrans_trans;
        ftrans_trans.mat = gemmi::Mat33(1,0,0,0,1,0,0,0,1);
        ftrans_trans.vec = gemmi::Vec3(0.5,0.5,0.5);
        std::vector<int> pre_shift_abc;

        gemmi::Structure mol2 = coot::mol_by_symmetry_gemmi(st, cell, ftrans_trans, pre_shift_abc);

        // Get first atom from original
        gemmi::Atom* atom_orig = nullptr;
        if (!st.models.empty()) {
            gemmi::Model& mod = st.models[0];
            if (!mod.chains.empty()) {
                gemmi::Chain& ch = mod.chains[0];
                if (!ch.residues.empty()) {
                    gemmi::Residue& res = ch.residues[0];
                    if (!res.atoms.empty()) {
                        atom_orig = &res.atoms[0];
                    }
                }
            }
        }

        // Get first atom from transformed molecule
        gemmi::Atom* atom2 = nullptr;
        if (!mol2.models.empty()) {
            gemmi::Model& mod3 = mol2.models[0];
            if (!mod3.chains.empty()) {
                gemmi::Chain& ch3 = mod3.chains[0];
                if (!ch3.residues.empty()) {
                    gemmi::Residue& res3 = ch3.residues[0];
                    if (!res3.atoms.empty()) {
                        atom2 = &res3.atoms[0];
                    }
                }
            }
        }

        if (atom_orig && atom2) {
            EXPECT_EQ(pad_atom_name(atom2->name), " N  ");
            EXPECT_NEAR(atom2->pos.x, 62.334, 1e-4);
            EXPECT_NEAR(atom2->pos.y, 52.772, 1e-4);
            EXPECT_NEAR(atom2->pos.z, 124.951, 1e-4);
        }

        EXPECT_EQ(mol2.models.size(), 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}