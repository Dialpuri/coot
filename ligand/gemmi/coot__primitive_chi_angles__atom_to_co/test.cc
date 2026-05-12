#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, PrimitiveChiAnglesAtomToCoGemmi) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // case 1: Test with residue A/10 (ILE)
    {
        // Find the residue in the structure
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const auto& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Failed to find chain A";

        const gemmi::Residue* res_10 = nullptr;
        for (const auto& res : chain_a->residues) {
            if (res.seqid.num.value == 10) {
                res_10 = &res;
                break;
            }
        }
        ASSERT_NE(res_10, nullptr) << "Failed to get residue A/10";

        std::cout << "INPUT residue: " << chain_a->name << " " << res_10->seqid.num.value << " " << res_10->name << std::endl;

        // Get first atom from residue
        ASSERT_GT(res_10->atoms.size(), 0) << "No atoms in residue A/10";
        const gemmi::Atom& at = res_10->atoms[0];
        std::cout << "INPUT atom_name: " << at.name << std::endl;
        std::cout << "INPUT atom_coords: (" << at.pos.x << ", " << at.pos.y << ", " << at.pos.z << ")" << std::endl;

        // Call atom_to_co_gemmi
        clipper::Coord_orth coord = coot::primitive_chi_angles::atom_to_co_gemmi(at);

        std::cout << "OUTPUT coord: (" << coord.x() << ", " << coord.y() << ", " << coord.z() << ")" << std::endl;
        
        // Check that coordinates match
        EXPECT_NEAR(at.pos.x, coord.x(), 1e-4) << "x coordinate mismatch";
        EXPECT_NEAR(at.pos.y, coord.y(), 1e-4) << "y coordinate mismatch";
        EXPECT_NEAR(at.pos.z, coord.z(), 1e-4) << "z coordinate mismatch";
        
        bool coord_match = (at.pos.x == coord.x() && at.pos.y == coord.y() && at.pos.z == coord.z());
        std::cout << "OUTPUT coord_match: " << (coord_match ? "true" : "false") << std::endl;
        EXPECT_TRUE(coord_match) << "Coordinates should match exactly";
    }

    // case 2: Test with residue A/20 (LYS)
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const auto& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr) << "Failed to find chain A";

        const gemmi::Residue* res_20 = nullptr;
        for (const auto& res : chain_a->residues) {
            if (res.seqid.num.value == 20) {
                res_20 = &res;
                break;
            }
        }
        ASSERT_NE(res_20, nullptr) << "Failed to get residue A/20";

        std::cout << "INPUT residue: " << chain_a->name << " " << res_20->seqid.num.value << " " << res_20->name << std::endl;

        ASSERT_GT(res_20->atoms.size(), 0) << "No atoms in residue A/20";
        const gemmi::Atom& at = res_20->atoms[0];
        std::cout << "INPUT atom_name: " << at.name << std::endl;
        std::cout << "INPUT atom_coords: (" << at.pos.x << ", " << at.pos.y << ", " << at.pos.z << ")" << std::endl;

        clipper::Coord_orth coord = coot::primitive_chi_angles::atom_to_co_gemmi(at);

        std::cout << "OUTPUT coord: (" << coord.x() << ", " << coord.y() << ", " << coord.z() << ")" << std::endl;
        
        // Verify coordinates match
        EXPECT_NEAR(at.pos.x, coord.x(), 1e-4);
        EXPECT_NEAR(at.pos.y, coord.y(), 1e-4);
        EXPECT_NEAR(at.pos.z, coord.z(), 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}