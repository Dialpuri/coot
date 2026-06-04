#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_instanced_graphical_bonds_hemispherical_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain &chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());
    gemmi::Residue &res = chain.residues[0];
    ASSERT_FALSE(res.atoms.empty());
    gemmi::Atom &atom1 = res.atoms[0];

    EXPECT_EQ(atom1.name, "N");

    // Build graphical_bonds_container_gemmi
    graphical_bonds_container_gemmi gbc;
    gbc.n_consolidated_atom_centres = 1;
    gbc.consolidated_atom_centres = new graphical_bonds_points_list<graphical_bonds_atom_info_t_gemmi>[1];

    gbc.consolidated_atom_centres[0].num_points = 1;
    gbc.consolidated_atom_centres[0].current_count = 1;

    graphical_bonds_atom_info_t_gemmi info;
    info.atom_p = &atom1;
    info.is_hydrogen_atom = false;
    info.is_water = false;
    info.radius_scale = 1.0f;
    info.position = coot::Cartesian(atom1.pos.x, atom1.pos.y, atom1.pos.z);
    info.atom_index = 0;
    info.model_number = 1;

    gbc.consolidated_atom_centres[0].points = new graphical_bonds_atom_info_t_gemmi[1];
    gbc.consolidated_atom_centres[0].points[0] = info;

    // Case 1: call with populated graphical_bonds_container
    {
        coot::instanced_mesh_t mesh;

        size_t n_before = mesh.geom.size();

        std::vector<glm::vec4> colour_table;
        colour_table.push_back(glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));

        float atom_radius = 1.0f;
        float bond_radius = 0.1f;
        unsigned int num_subdivisions = 1;

        EXPECT_EQ(n_before, size_t(0));
        EXPECT_EQ(gbc.n_consolidated_atom_centres, 1);
        EXPECT_EQ(colour_table.size(), size_t(1));

        make_instanced_graphical_bonds_hemispherical_atoms_gemmi(mesh, gbc,
                                                           coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS,
                                                           atom_radius, bond_radius,
                                                           num_subdivisions, colour_table);

        size_t n_after = mesh.geom.size();
        EXPECT_EQ(n_before, size_t(0));
        EXPECT_EQ(n_after, size_t(0));
        EXPECT_FALSE(n_before != n_after);
    }

    // Case 2: call with empty graphical_bonds_container
    {
        coot::instanced_mesh_t mesh;
        graphical_bonds_container_gemmi empty_gbc;
        empty_gbc.n_consolidated_atom_centres = 0;
        empty_gbc.consolidated_atom_centres = nullptr;

        size_t n_before = mesh.geom.size();

        std::vector<glm::vec4> colour_table;
        float atom_radius = 0.5f;
        float bond_radius = 0.1f;
        unsigned int num_subdivisions = 2;

        EXPECT_EQ(empty_gbc.n_consolidated_atom_centres, 0);

        make_instanced_graphical_bonds_hemispherical_atoms_gemmi(mesh, empty_gbc,
                                                           coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS,
                                                           atom_radius, bond_radius,
                                                           num_subdivisions, colour_table);

        size_t n_after = mesh.geom.size();
        EXPECT_EQ(n_before, size_t(0));
        EXPECT_EQ(n_after, size_t(0));
        EXPECT_FALSE(n_before != n_after);
    }

    // Cleanup
    delete[] gbc.consolidated_atom_centres[0].points;
    delete[] gbc.consolidated_atom_centres;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}