#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_instanced_graphical_bonds_spherical_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
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
    info.render_as_aniso = false;

    gbc.consolidated_atom_centres[0].points = new graphical_bonds_atom_info_t_gemmi[1];
    gbc.consolidated_atom_centres[0].points[0] = info;

    // Colour table
    std::vector<glm::vec4> colour_table;
    colour_table.push_back(glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));

    // Case 1: Normal case
    {
        coot::instanced_mesh_t m;
        int geom_before = m.geom.size();

        make_instanced_graphical_bonds_spherical_atoms_gemmi(
            m, gbc, coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS,
            0.12f, 0.05f,
            false, 0.5f,
            false, false,
            1, colour_table
        );

        int geom_after = m.geom.size();

        EXPECT_EQ(geom_before, 0);
        EXPECT_EQ(geom_after, 1);
    }

    // Case 2: render_aniso_atoms_as_ortep=true but atom has no aniso data
    {
        coot::instanced_mesh_t m2;
        make_instanced_graphical_bonds_spherical_atoms_gemmi(
            m2, gbc, coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS,
            0.12f, 0.05f,
            true, 0.5f,
            true, false,
            1, colour_table
        );
        EXPECT_EQ(m2.geom.size(), 1);
    }

    // Case 3: water with radius_scale=2.0 and larger base radius - tests the 0.65 cap
    {
        info.is_water = true;
        info.radius_scale = 2.0f;
        gbc.consolidated_atom_centres[0].points[0] = info;

        coot::instanced_mesh_t m3;
        make_instanced_graphical_bonds_spherical_atoms_gemmi(
            m3, gbc, coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS,
            0.5f, 0.05f,
            false, 0.5f,
            false, false,
            1, colour_table
        );
        EXPECT_EQ(m3.geom.size(), 1);
    }

    // Cleanup
    delete[] gbc.consolidated_atom_centres[0].points;
    delete[] gbc.consolidated_atom_centres;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}