#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include <set>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "function.hh"

namespace {

bool is_flippable(const std::string& resname) {
    static const std::set<std::string> flippable = {
        "ARG", "ASP", "ASN", "CYS", "GLN", "GLU",
        "PHE", "HIS", "SER", "THR", "VAL", "TRP", "TYR"
    };
    return flippable.count(resname) > 0;
}

} // namespace

TEST(OracleTest, do_180_degree_side_chain_flip) {

    // --- Case 1: GLU at A/2 (original oracle case) ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());

        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Model& model = st.models[0];

        // Find first flippable residue (GLU at A/2)
        std::string target_resname;
        int target_res_no = -1;
        std::string target_chain_id;

        for (gemmi::Chain& chain : model.chains) {
            target_chain_id = chain.name;
            for (gemmi::Residue& res : chain.residues) {
                std::string resname = res.name;
                if (is_flippable(resname)) {
                    target_resname = resname;
                    target_res_no = res.seqid.num.value;
                    break;
                }
            }
            if (target_res_no != -1) break;
        }

        EXPECT_NE(target_res_no, -1);
        EXPECT_EQ(target_resname, "GLU");
        EXPECT_EQ(target_res_no, 2);
        EXPECT_EQ(target_chain_id, "A");

        // Find the residue pointer
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == target_chain_id) {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == target_res_no) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        EXPECT_NE(res, nullptr);

        std::string resname = res->name;
        EXPECT_EQ(resname, "GLU");

        EXPECT_EQ(res->atoms.size(), 9u);

        double com_x_before = 0, com_y_before = 0, com_z_before = 0;
        for (const auto& atom : res->atoms) {
            com_x_before += atom.pos.x;
            com_y_before += atom.pos.y;
            com_z_before += atom.pos.z;
        }
        com_x_before /= res->atoms.size();
        com_y_before /= res->atoms.size();
        com_z_before /= res->atoms.size();

        EXPECT_NEAR(com_x_before, 13.323, 1e-4);
        EXPECT_NEAR(com_y_before, 3.69133, 1e-4);
        EXPECT_NEAR(com_z_before, 79.8516, 1e-4);

        std::string alt_conf("");
        EXPECT_NO_THROW(coot::do_180_degree_side_chain_flip_gemmi(res, alt_conf, &geom));

        double com_x_after = 0, com_y_after = 0, com_z_after = 0;
        for (const auto& atom : res->atoms) {
            com_x_after += atom.pos.x;
            com_y_after += atom.pos.y;
            com_z_after += atom.pos.z;
        }
        com_x_after /= res->atoms.size();
        com_y_after /= res->atoms.size();
        com_z_after /= res->atoms.size();

        EXPECT_NEAR(com_x_after, 13.3251, 1e-4);
        EXPECT_NEAR(com_y_after, 3.69332, 1e-4);
        EXPECT_NEAR(com_z_after, 79.851, 1e-4);

        double dx = com_x_after - com_x_before;
        double dy = com_y_after - com_y_before;
        double dz = com_z_after - com_z_before;
        double total_shift = std::sqrt(dx * dx + dy * dy + dz * dz);

        EXPECT_NEAR(total_shift, 0.00294125, 1e-4);
        EXPECT_TRUE(total_shift > 0.001);
    }

    // --- Case 2: SER at A/3 (complementary — different flippable residue) ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());

        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Model& model = st.models[0];

        std::string target_chain_id = "A";
        int target_res_no = 3;

        // Find the residue pointer
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == target_chain_id) {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == target_res_no) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        EXPECT_NE(res, nullptr);

        std::string resname = res->name;

        EXPECT_GT(res->atoms.size(), 0u);

        double com_x_before = 0, com_y_before = 0, com_z_before = 0;
        for (const auto& atom : res->atoms) {
            com_x_before += atom.pos.x;
            com_y_before += atom.pos.y;
            com_z_before += atom.pos.z;
        }
        com_x_before /= res->atoms.size();
        com_y_before /= res->atoms.size();
        com_z_before /= res->atoms.size();

        std::string alt_conf("");
        EXPECT_NO_THROW(coot::do_180_degree_side_chain_flip_gemmi(res, alt_conf, &geom));

        double com_x_after = 0, com_y_after = 0, com_z_after = 0;
        for (const auto& atom : res->atoms) {
            com_x_after += atom.pos.x;
            com_y_after += atom.pos.y;
            com_z_after += atom.pos.z;
        }
        com_x_after /= res->atoms.size();
        com_y_after /= res->atoms.size();
        com_z_after /= res->atoms.size();

        double dx = com_x_after - com_x_before;
        double dy = com_y_after - com_y_before;
        double dz = com_z_after - com_z_before;
        double total_shift = std::sqrt(dx * dx + dy * dy + dz * dz);

        // SER OH flip should produce some movement
        EXPECT_GT(total_shift, 0.0);
    }

    // --- Case 3: ALA at A/1 (non-flippable — different branch) ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());

        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Model& model = st.models[0];

        std::string target_chain_id = "A";
        int target_res_no = 1;

        // Find the residue pointer
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == target_chain_id) {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == target_res_no) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        EXPECT_NE(res, nullptr);

        std::string resname = res->name;

        double com_x_before = 0, com_y_before = 0, com_z_before = 0;
        for (const auto& atom : res->atoms) {
            com_x_before += atom.pos.x;
            com_y_before += atom.pos.y;
            com_z_before += atom.pos.z;
        }
        com_x_before /= res->atoms.size();
        com_y_before /= res->atoms.size();
        com_z_before /= res->atoms.size();

        std::string alt_conf("");
        EXPECT_NO_THROW(coot::do_180_degree_side_chain_flip_gemmi(res, alt_conf, &geom));

        double com_x_after = 0, com_y_after = 0, com_z_after = 0;
        for (const auto& atom : res->atoms) {
            com_x_after += atom.pos.x;
            com_y_after += atom.pos.y;
            com_z_after += atom.pos.z;
        }
        com_x_after /= res->atoms.size();
        com_y_after /= res->atoms.size();
        com_z_after /= res->atoms.size();

        double dx = com_x_after - com_x_before;
        double dy = com_y_after - com_y_before;
        double dz = com_z_after - com_z_before;
        double total_shift = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Non-flippable residue should have negligible shift
        EXPECT_LT(total_shift, 0.001);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}