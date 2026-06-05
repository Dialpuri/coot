#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTestGemmi, place_hydrogen_by_connected_atom_energy_type) {
    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    coot::protein_geometry geom;
    geom.init_standard();

    // Get the dictionary for ALA
    std::string res_name = "ALA";
    std::pair<bool, coot::dictionary_residue_restraints_t> rest_pair =
        geom.get_monomer_restraints(res_name, coot::protein_geometry::IMOL_ENC_ANY);
    ASSERT_TRUE(rest_pair.first);
    const coot::dictionary_residue_restraints_t& rest = rest_pair.second;

    // Find a hydrogen atom with exactly one neighbor that has non-empty type_energy
    unsigned int chosen_iat = 0;
    unsigned int chosen_iat_neighb = 0;
    bool found = false;

    for (unsigned int iat = 0; iat < rest.atom_info.size(); iat++) {
        if (!rest.atom_info[iat].is_hydrogen()) continue;
        std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
        if (neighbs.size() != 1) continue;
        unsigned int iat_neighb = neighbs[0];
        const std::string& energy_type = rest.atom_info[iat_neighb].type_energy;
        if (energy_type.empty()) continue;

        chosen_iat = iat;
        chosen_iat_neighb = iat_neighb;
        found = true;
        break;
    }
    ASSERT_TRUE(found);

    // Find a real ALA residue from the PDB using gemmi
    gemmi::Residue* residue_p = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "ALA") {
                residue_p = &res;
                break;
            }
        }
        if (residue_p) break;
    }
    ASSERT_NE(residue_p, nullptr);

    // Bond length defaults from the code
    double bl_aliph = 1.09;
    double bl_arom  = 1.08;
    double bl_amino = 1.01;
    double bl_oh    = 0.96;
    double bl_sh    = 1.34;

    // Verify INPUT values match oracle
    EXPECT_EQ(chosen_iat, 6u);
    EXPECT_EQ(chosen_iat_neighb, 0u);
    EXPECT_EQ(rest.atom_info[chosen_iat].atom_id_4c, " H  ");
    EXPECT_EQ(rest.atom_info[chosen_iat_neighb].atom_id_4c, " N  ");
    EXPECT_EQ(rest.atom_info[chosen_iat_neighb].type_energy, "NT3");

    // Count atoms before
    int atoms_before = static_cast<int>(residue_p->atoms.size());
    EXPECT_EQ(atoms_before, 10);

    // Call the function
    std::vector<std::string> placed = coot::reduce::place_hydrogen_by_connected_atom_energy_type_gemmi(
        chosen_iat, chosen_iat_neighb, rest, residue_p,
        bl_aliph, bl_arom, bl_amino, bl_oh, bl_sh
    );

    int atoms_after = static_cast<int>(residue_p->atoms.size());

    // Assertions for first case
    EXPECT_EQ(atoms_after, 12);
    EXPECT_EQ(placed.size(), 3u);
    EXPECT_EQ(placed[0], " H  ");
    EXPECT_EQ(placed[1], " H2 ");
    EXPECT_EQ(placed[2], " H3 ");

    // Second case: different hydrogen on same residue type
    {
        unsigned int chosen_iat2 = 0;
        unsigned int chosen_iat_neighb2 = 0;
        bool found2 = false;

        for (unsigned int iat = 0; iat < rest.atom_info.size(); iat++) {
            if (!rest.atom_info[iat].is_hydrogen()) continue;
            if (iat == chosen_iat) continue;
            std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
            if (neighbs.size() != 1) continue;
            unsigned int iat_neighb = neighbs[0];
            const std::string& energy_type = rest.atom_info[iat_neighb].type_energy;
            if (energy_type.empty()) continue;

            chosen_iat2 = iat;
            chosen_iat_neighb2 = iat_neighb;
            found2 = true;
            break;
        }

        if (found2) {
            // Find a second ALA residue
            gemmi::Residue* residue_p2 = nullptr;
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.name == "ALA" && &res != residue_p) {
                        residue_p2 = &res;
                        break;
                    }
                }
                if (residue_p2) break;
            }

            if (residue_p2) {
                int atoms_before2 = static_cast<int>(residue_p2->atoms.size());
                std::vector<std::string> placed2 = coot::reduce::place_hydrogen_by_connected_atom_energy_type_gemmi(
                    chosen_iat2, chosen_iat_neighb2, rest, residue_p2,
                    bl_aliph, bl_arom, bl_amino, bl_oh, bl_sh
                );
                int atoms_after2 = static_cast<int>(residue_p2->atoms.size());

                EXPECT_EQ(atoms_before2, 10);
                EXPECT_EQ(atoms_after2, 12);
                EXPECT_EQ(placed2.size(), 3u);
                EXPECT_EQ(placed2[0], " H2 ");
                EXPECT_EQ(placed2[1], " H  ");
                EXPECT_EQ(placed2[2], " H3 ");
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}