#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(SphereSizeGemmiTest, sphere_size) {
    // Load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Create protein_geometry for vdw radii
    coot::protein_geometry geom;
    geom.init_standard();

    // Select chain A (mimics original: mol->SelectAtoms(selhnd, 0, "A", ...))
    std::vector<std::string> chain_names = {"A"};

    // Count atoms in chain A to verify selection
    int n_atoms_chain_a = 0;
    for (const auto& chain : model.chains) {
        if (chain.name == "A") {
            for (const auto& residue : chain.residues) {
                n_atoms_chain_a += static_cast<int>(residue.atoms.size());
            }
        }
    }
    EXPECT_EQ(n_atoms_chain_a, 2106);

    // Get first atom position from chain A (same as atom_selection[0] in original)
    gemmi::Vec3 first_atom_pos;
    {
        bool found = false;
        for (const auto& chain : model.chains) {
            if (chain.name == "A") {
                for (const auto& residue : chain.residues) {
                    for (const auto& atom : residue.atoms) {
                        first_atom_pos = atom.pos;
                        found = true;
                        break;
                    }
                    if (found) break;
                }
                if (found) break;
            }
        }
    }

    // case: point at atom position (sphere_size negative - inside atom)
    {
        gemmi::Vec3 pt(first_atom_pos.x, first_atom_pos.y, first_atom_pos.z);
        double ss = coot::hole::sphere_size_gemmi(model, chain_names, pt, geom);
        EXPECT_NEAR(ss, -1.55, 1e-4);
    }

    // case: point far away from all atoms (large positive sphere)
    {
        gemmi::Vec3 pt(first_atom_pos.x + 100, first_atom_pos.y + 100, first_atom_pos.z + 100);
        double ss = coot::hole::sphere_size_gemmi(model, chain_names, pt, geom);
        EXPECT_NEAR(ss, 150.528, 1e-3);
    }

    // case: empty selection — returns -1
    // Use a chain name that doesn't exist in the model
    {
        std::vector<std::string> empty_chains = {"ZZZZZ"};
        gemmi::Vec3 pt(0, 0, 0);
        double ss = coot::hole::sphere_size_gemmi(model, empty_chains, pt, geom);
        EXPECT_NEAR(ss, -1.0, 1e-4);
    }

    // Additional case: point in the middle of the structure (positive sphere)
    // This exercises the branch where largest_possible_sphere is updated multiple times
    {
        gemmi::Vec3 pt(50, 50, 50);
        double ss = coot::hole::sphere_size_gemmi(model, chain_names, pt, geom);
        // Just verify it returns a finite value (not -1, since we have atoms)
        EXPECT_GT(ss, -1.0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}