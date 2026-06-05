#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cstdio>
#include <cstdlib>
#include "function.hh"

TEST(OracleTest, side_chain_densities_write_density_box) {
    // Clean up any leftover side-chain-data from previous runs
    system("rm -rf side-chain-data 2>/dev/null");

    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1);

    gemmi::Model& model = st.models[0];
    std::string id = "test";

    // Case 1: Valid residue with filled density box
    {
        // Find residue A/3 (or fall back to A/1)
        gemmi::CRA cra{nullptr, nullptr, nullptr};
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 3) {
                        cra = gemmi::CRA{&chain, &res, nullptr};
                        break;
                    }
                }
                if (!cra.residue) {
                    // Fallback to A/1
                    for (auto& res : chain.residues) {
                        if (res.seqid.num.value == 1) {
                            cra = gemmi::CRA{&chain, &res, nullptr};
                            break;
                        }
                    }
                }
                break;
            }
        }

        EXPECT_NE(cra.residue, nullptr);

        int n_steps = 2;
        int n_per_side = 2 * n_steps + 1;
        int n_box_vol = n_per_side * n_per_side * n_per_side;

        float *density_box = new float[n_box_vol];
        for (int i = 0; i < n_box_vol; i++) {
            density_box[i] = 1.0f;
        }

        std::string res_name = (cra.residue ? cra.residue->name : "NONE");
        std::string chain_id = (cra.chain ? cra.chain->name : "?");
        int seq_num = (cra.residue ? cra.residue->seqid.num.value : -1);

        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(seq_num, 3);
        EXPECT_EQ(res_name, "ASN");
        EXPECT_EQ(n_steps, 2);
        EXPECT_EQ(id, "test");

        // Rotamer name "m120" for ASN — determined by the residue conformation
        std::string rotamer = "m120";

        EXPECT_NO_THROW(coot::side_chain_densities::write_density_box_gemmi(
            density_box, n_steps, cra, id, rotamer));

        EXPECT_EQ(rotamer, "m120");

        // Check if expected file was created (oracle says file_created: false)
        if (!rotamer.empty() && cra.residue) {
            std::string expected_file = "side-chain-data/" + res_name + "/" + rotamer + "/test-" + chain_id + "-" + std::to_string(seq_num) + ".tab";
            FILE *f = fopen(expected_file.c_str(), "r");
            EXPECT_EQ(f, nullptr); // oracle output: file_created: false
            EXPECT_EQ(expected_file, "side-chain-data/ASN/m120/test-A-3.tab");
            if (f) fclose(f);
        }

        delete[] density_box;
    }

    // Case 2: Null residue pointer — triggers the guard clause (returns immediately)
    {
        int n_steps = 2;
        int n_per_side = 2 * n_steps + 1;
        int n_box_vol = n_per_side * n_per_side * n_per_side;
        float *density_box = new float[n_box_vol];
        for (int i = 0; i < n_box_vol; i++) {
            density_box[i] = 1.0f;
        }

        gemmi::CRA cra{nullptr, nullptr, nullptr};

        EXPECT_NO_THROW(coot::side_chain_densities::write_density_box_gemmi(
            density_box, n_steps, cra, id, "m120"));

        delete[] density_box;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}