#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced) {
    int imol_no = 0;

    // Case 1: Valid bonds from real molecule
    {
        coot::instanced_mesh_t m;

        // Create a simple graphical_bonds_container with gemmi atoms from example.pdb
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        graphical_bonds_container_gemmi gbc;

        // Collect atoms from first model
        auto& model = st.models[0];
        std::vector<graphical_bonds_atom_info_t_gemmi> all_atoms;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                for (auto& atom : res.atoms) {
                    graphical_bonds_atom_info_t_gemmi info;
                    info.atom_p = &atom;
                    info.position = glm::vec3(atom.pos.x, atom.pos.y, atom.pos.z);
                    all_atoms.push_back(info);
                }
            }
        }

        // Create one consolidated_atom_centre with all atoms
        graphical_bonds_consolidated_atom_centre_gemmi centre;
        centre.num_points = all_atoms.size();
        centre.points = all_atoms;
        centre.colour = glm::vec4(1.0f, 0.5f, 0.2f, 1.0f);
        centre.radius = 1.5f;

        std::vector<graphical_bonds_consolidated_atom_centre_gemmi> centres;
        centres.push_back(centre);

        gbc.n_consolidated_atom_centres = 1;
        gbc.consolidated_atom_centres = centres;

        // Create colour table with one entry per consolidated_atom_centre
        std::vector<glm::vec4> colour_table;
        colour_table.push_back(glm::vec4(1.0f, 0.5f, 0.2f, 1.0f));

        coot::protein_geometry geom;
        geom.init_standard();

        EXPECT_EQ(m.geom.size(), 0u);

        unsigned int num_subdivisions = 3;
        coot::make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced_gemmi(
            m, gbc, num_subdivisions, colour_table, geom, imol_no);

        EXPECT_EQ(m.geom.size(), 1u);
        EXPECT_FALSE(m.geom.empty());
        // The number of instances = total atoms in the PDB
        EXPECT_EQ(m.geom[0].instancing_data_A.size(), all_atoms.size());
    }

    // Case 2: Empty graphical_bonds_container (edge case — no atoms)
    {
        coot::instanced_mesh_t m;
        graphical_bonds_container_gemmi empty_gbc;
        std::vector<glm::vec4> empty_colours;
        unsigned int num_subdivisions = 1;

        EXPECT_EQ(empty_gbc.n_consolidated_atom_centres, 0);
        EXPECT_EQ(m.geom.size(), 0u);

        coot::protein_geometry geom;
        geom.init_standard();

        coot::make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced_gemmi(
            m, empty_gbc, num_subdivisions, empty_colours, geom, 0);

        EXPECT_EQ(m.geom.size(), 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}