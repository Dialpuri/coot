#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__residue_has_hetatms/gemmi/function.hh"

TEST(OracleTest, move_hetgroups_around_protein) {
    {
        gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        EXPECT_FALSE(st.models.empty());

        // Gather hetgroup residues before the call
        std::vector<std::string> het_residue_ids_before;
        std::vector<std::vector<double>> het_coords_before;

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& residue : chain.residues) {
                    std::string residue_name(residue.name);
                    if (residue_name == "WAT" || residue_name == "HOH") continue;
                    if (coot::util::residue_has_hetatms_gemmi(residue)) {
                        std::string res_id(chain.name);
                        res_id += "/";
                        res_id += std::to_string(residue.seqid.num.value);
                        res_id += " ";
                        res_id += residue_name;
                        het_residue_ids_before.push_back(res_id);

                        std::vector<double> coords;
                        for (gemmi::Atom& atom : residue.atoms) {
                            coords.push_back(atom.pos.x);
                            coords.push_back(atom.pos.y);
                            coords.push_back(atom.pos.z);
                        }
                        het_coords_before.push_back(coords);
                    }
                }
            }
        }

        // Oracle: INPUT hetgroup_count: 1
        EXPECT_EQ(het_residue_ids_before.size(), 1u);

        // Oracle: INPUT het_res_before(A/1299 LZA)
        EXPECT_EQ(het_residue_ids_before[0], "A/1299 LZA");

        // Verify first and last atom coordinates from oracle output
        const auto& coords = het_coords_before[0];
        EXPECT_NEAR(coords[0], 31.856, 0.001);
        EXPECT_NEAR(coords[1], 9.409, 0.001);
        EXPECT_NEAR(coords[2], 67.127, 0.001);
        EXPECT_NEAR(coords[coords.size()-3], 32.438, 0.001);
        EXPECT_NEAR(coords[coords.size()-2], 9.322, 0.001);
        EXPECT_NEAR(coords[coords.size()-1], 64.731, 0.001);

        // Call the function
        EXPECT_NO_THROW(coot::util::move_hetgroups_around_protein_gemmi(st));

        // Gather hetgroup residues after the call
        int het_count_after = 0;
        std::vector<std::string> het_residue_ids_after;
        std::vector<std::vector<double>> het_coords_after;

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& residue : chain.residues) {
                    std::string residue_name(residue.name);
                    if (residue_name == "WAT" || residue_name == "HOH") continue;
                    if (coot::util::residue_has_hetatms_gemmi(residue)) {
                        het_count_after++;
                        std::string res_id(chain.name);
                        res_id += "/";
                        res_id += std::to_string(residue.seqid.num.value);
                        res_id += " ";
                        res_id += residue_name;
                        het_residue_ids_after.push_back(res_id);

                        std::vector<double> coords;
                        for (gemmi::Atom& atom : residue.atoms) {
                            coords.push_back(atom.pos.x);
                            coords.push_back(atom.pos.y);
                            coords.push_back(atom.pos.z);
                        }
                        het_coords_after.push_back(coords);
                    }
                }
            }
        }

        // Oracle: OUTPUT hetgroup_count_after: 1
        EXPECT_EQ(het_count_after, 1);

        // Oracle: OUTPUT het_res_after(A/1299 LZA) with same coords as before
        EXPECT_EQ(het_residue_ids_after[0], "A/1299 LZA");
        EXPECT_NEAR(het_coords_after[0][0], 31.856, 0.001);
        EXPECT_NEAR(het_coords_after[0][1], 9.409, 0.001);
        EXPECT_NEAR(het_coords_after[0][2], 67.127, 0.001);
        EXPECT_NEAR(het_coords_after[0][het_coords_after[0].size()-3], 32.438, 0.001);
        EXPECT_NEAR(het_coords_after[0][het_coords_after[0].size()-2], 9.322, 0.001);
        EXPECT_NEAR(het_coords_after[0][het_coords_after[0].size()-1], 64.731, 0.001);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}