#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, place_hydrogen_by_connected_2nd_neighbours) {
    // case: oracle original — MET with all hydrogens present
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model *model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain *chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue *residue = &chain->residues[0];
        ASSERT_NE(residue, nullptr);

        std::string resname = residue->name;
        EXPECT_EQ(resname, "MET");

        // Create protein geometry instance
        coot::protein_geometry geom;
        geom.init_standard();

        // Get restraints for this residue type
        std::pair<bool, coot::dictionary_residue_restraints_t> result =
            geom.get_monomer_restraints(resname, 0);

        ASSERT_TRUE(result.first);

        coot::dictionary_residue_restraints_t rest = result.second;
        EXPECT_EQ(rest.atom_info.size(), 20u);

        // Find a hydrogen atom and its neighbor with empty energy_type
        unsigned int iat = 0;
        unsigned int iat_neighb = 0;
        bool found = false;
        for (unsigned int i = 0; i < rest.atom_info.size(); i++) {
            if (rest.atom_info[i].is_hydrogen()) {
                std::vector<unsigned int> neighbs = rest.neighbours(i, false);
                if (neighbs.size() == 1) {
                    const std::string &energy_type = rest.atom_info[neighbs[0]].type_energy;
                    if (energy_type.empty()) {
                        iat = i;
                        iat_neighb = neighbs[0];
                        found = true;
                        break;
                    }
                }
            }
        }

        if (!found) {
            // Fallback: use any hydrogen with one neighbor
            for (unsigned int i = 0; i < rest.atom_info.size(); i++) {
                if (rest.atom_info[i].is_hydrogen()) {
                    std::vector<unsigned int> neighbs = rest.neighbours(i, false);
                    if (neighbs.size() == 1) {
                        iat = i;
                        iat_neighb = neighbs[0];
                        found = true;
                        break;
                    }
                }
            }
        }

        ASSERT_TRUE(found);
        EXPECT_EQ(iat, 9u);
        EXPECT_EQ(rest.atom_info[iat].atom_id_4c, " H  ");
        EXPECT_EQ(iat_neighb, 0u);
        EXPECT_EQ(rest.atom_info[iat_neighb].atom_id_4c, " N  ");
        EXPECT_EQ(rest.atom_info[iat_neighb].type_energy, "NT3");

        int atoms_before = residue->atoms.size();
        EXPECT_EQ(atoms_before, 19);

        coot::reduce::place_hydrogen_by_connected_2nd_neighbours_gemmi(iat, iat_neighb, rest, *residue);

        int atoms_after = residue->atoms.size();
        EXPECT_EQ(atoms_after, 19);
    }

    // case: complementary — strip hydrogens from residue first, then call function
    // The function is a stub so no atoms should be added, but we verify it doesn't crash
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model *model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain *chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue *residue = &chain->residues[0];
        ASSERT_NE(residue, nullptr);

        std::string resname = residue->name;
        EXPECT_EQ(resname, "MET");

        // Strip hydrogen atoms from the residue
        for (auto it = residue->atoms.begin(); it != residue->atoms.end(); ) {
            if (it->element.name() == "H") {
                it = residue->atoms.erase(it);
            } else {
                ++it;
            }
        }

        int atoms_after_strip = residue->atoms.size();
        // Should have non-H atoms left (backbone + sidechain heavy atoms)
        EXPECT_FALSE(atoms_after_strip == 0);

        // Create protein geometry instance
        coot::protein_geometry geom;
        geom.init_standard();

        // Get restraints for this residue type
        std::pair<bool, coot::dictionary_residue_restraints_t> result =
            geom.get_monomer_restraints(resname, 0);

        ASSERT_TRUE(result.first);

        coot::dictionary_residue_restraints_t rest = result.second;

        // Find a hydrogen atom and its neighbor
        unsigned int iat = 0;
        unsigned int iat_neighb = 0;
        bool found = false;
        for (unsigned int i = 0; i < rest.atom_info.size(); i++) {
            if (rest.atom_info[i].is_hydrogen()) {
                std::vector<unsigned int> neighbs = rest.neighbours(i, false);
                if (neighbs.size() == 1) {
                    const std::string &energy_type = rest.atom_info[neighbs[0]].type_energy;
                    if (energy_type.empty()) {
                        iat = i;
                        iat_neighb = neighbs[0];
                        found = true;
                        break;
                    }
                }
            }
        }

        if (!found) {
            for (unsigned int i = 0; i < rest.atom_info.size(); i++) {
                if (rest.atom_info[i].is_hydrogen()) {
                    std::vector<unsigned int> neighbs = rest.neighbours(i, false);
                    if (neighbs.size() == 1) {
                        iat = i;
                        iat_neighb = neighbs[0];
                        found = true;
                        break;
                    }
                }
            }
        }

        ASSERT_TRUE(found);

        int atoms_before = residue->atoms.size();

        EXPECT_NO_THROW(
            coot::reduce::place_hydrogen_by_connected_2nd_neighbours_gemmi(iat, iat_neighb, rest, *residue)
        );

        int atoms_after = residue->atoms.size();

        // The function is currently a stub — it does not actually place atoms.
        // Verify that atoms count is unchanged (function didn't crash, did nothing).
        EXPECT_EQ(atoms_before, atoms_after);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
