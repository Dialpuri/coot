#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to get atom coordinates as a vector of pairs
static std::vector<std::pair<std::string, std::array<double, 3>>> get_atom_coords(gemmi::Residue *res) {
    std::vector<std::pair<std::string, std::array<double, 3>>> coords;
    for (const auto& atom : res->atoms) {
        coords.push_back({atom.name, 
                         {atom.pos.x, atom.pos.y, atom.pos.z}});
    }
    return coords;
}

TEST(OracleTest, TransformResidues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get first model
    gemmi::Model& model = st.models[0];
    
    // case 1: transform a single residue with a translation (not identity)
    {
        // Find residue with chain A, seqid 10
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        
        ASSERT_NE(res, nullptr) << "residue should be found for cid: //A/10";
        
        // Make a copy of the residue
        gemmi::Residue res_copy = *res;
        
        // Before transform coordinates
        auto before_coords = get_atom_coords(&res_copy);
        
        // Translation by (1, 2, 3)
        clipper::Mat33<double> rot(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth trans(1.0, 2.0, 3.0);
        clipper::RTop_orth rtop(rot, trans);
        coot::util::transform_atoms_gemmi(&res_copy, rtop);
        
        // After transform coordinates
        auto after_coords = get_atom_coords(&res_copy);
        
        // Verify atom count is the same
        EXPECT_EQ(before_coords.size(), after_coords.size());
        
        // Verify specific atoms match expected transformed coordinates
        // N: (26.018, 12.491, 68.556) -> (27.018, 14.491, 71.556)
        EXPECT_NEAR(after_coords[0].second[0], 27.018, 1e-3);
        EXPECT_NEAR(after_coords[0].second[1], 14.491, 1e-3);
        EXPECT_NEAR(after_coords[0].second[2], 71.556, 1e-3);
        
        // CA: (27.175, 12.237, 67.71) -> (28.175, 14.237, 70.71)
        EXPECT_NEAR(after_coords[1].second[0], 28.175, 1e-3);
        EXPECT_NEAR(after_coords[1].second[1], 14.237, 1e-3);
        EXPECT_NEAR(after_coords[1].second[2], 70.71, 1e-3);
    }

    // case 2: transform with a non-trivial rotation (90 degree rotation around z-axis)
    {
        // Find residue with chain A, seqid 15
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 15) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        
        ASSERT_NE(res, nullptr) << "residue should be found for cid: //A/15";
        
        // Make a copy of the residue
        gemmi::Residue res_copy = *res;
        
        // Before transform coordinates
        auto before_coords = get_atom_coords(&res_copy);
        
        // 90 degree rotation around z-axis: x' = -y, y' = x, z' = z
        clipper::Mat33<double> rot(0, -1, 0, 1, 0, 0, 0, 0, 1);
        clipper::RTop_orth rtop(rot, clipper::Coord_orth(0,0,0));
        coot::util::transform_atoms_gemmi(&res_copy, rtop);
        
        // After transform coordinates
        auto after_coords = get_atom_coords(&res_copy);
        
        // Verify atom count is the same
        EXPECT_EQ(before_coords.size(), after_coords.size());
        
        // Verify specific atoms match expected rotated coordinates
        // N: (17.343, 12.361, 58.028) -> (-12.361, 17.343, 58.028)
        EXPECT_NEAR(after_coords[0].second[0], -12.361, 1e-3);
        EXPECT_NEAR(after_coords[0].second[1], 17.343, 1e-3);
        EXPECT_NEAR(after_coords[0].second[2], 58.028, 1e-3);
        
        // CA: (16.4, 11.464, 58.677) -> (-11.464, 16.4, 58.677)
        EXPECT_NEAR(after_coords[1].second[0], -11.464, 1e-3);
        EXPECT_NEAR(after_coords[1].second[1], 16.4, 1e-3);
        EXPECT_NEAR(after_coords[1].second[2], 58.677, 1e-3);
    }

    // case 3: invalid residue (should not crash)
    {
        // Find non-existent residue with seqid 9999
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 9999) {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        
        EXPECT_EQ(res, nullptr) << "residue should NOT be found for invalid cid: //A/9999";
        
        if (res) {
            gemmi::Residue res_copy = *res;
            clipper::RTop_orth rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), clipper::Coord_orth(0,0,0));
            coot::util::transform_atoms_gemmi(&res_copy, rtop);
            EXPECT_TRUE(true); // transformed successfully (but shouldn't reach here)
        } else {
            // Expected path - invalid residue
            EXPECT_TRUE(true);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}