#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, similar_normals) {
    // Load a PDB to get gemmi::Residue pointers
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GT(st.models.size(), 0) << "no models";
    
    gemmi::Model& model = st.models[0];
    EXPECT_EQ(st.models.size(), 1) << "expected 1 model";
    
    EXPECT_GT(model.chains.size(), 0) << "no chains";
    gemmi::Chain& chain = model.chains[0];
    
    EXPECT_GT(chain.residues.size(), 1) << "not enough residues";
    
    gemmi::Residue& res1 = chain.residues[0];
    gemmi::Residue& res2 = chain.residues[1];
    
    // Verify residue identities match oracle: A/1 (MET) and A/2 (GLU)
    EXPECT_EQ(chain.name, "A");
    EXPECT_EQ(res1.seqid.num.value, 1);
    
    // Build normal_map and angle_crit (23 degrees)
    std::map<gemmi::Residue *, clipper::Coord_orth> normal_map;
    double angle_crit = 23.0 * 3.14159265358979 / 180.0; // ~0.401 radians
    
    // Case 1: Identical normals (dot = 1.0, cos(23°) ≈ 0.921, 1.0 > 0.921 => true)
    {
        normal_map.clear();
        normal_map[&res1] = clipper::Coord_orth(0.0, 0.0, 1.0);
        normal_map[&res2] = clipper::Coord_orth(0.0, 0.0, 1.0);
        bool result = coot::stack_and_pair::similar_normals_gemmi(&res1, &res2, normal_map, angle_crit);
        EXPECT_TRUE(result);
    }
    
    // Case 2: Anti-parallel normals (dot = -1.0, |-1.0| > 0.921 => true)
    {
        normal_map.clear();
        normal_map[&res1] = clipper::Coord_orth(0.0, 0.0, 1.0);
        normal_map[&res2] = clipper::Coord_orth(0.0, 0.0, -1.0);
        bool result = coot::stack_and_pair::similar_normals_gemmi(&res1, &res2, normal_map, angle_crit);
        EXPECT_TRUE(result);
    }
    
    // Case 3: Orthogonal normals (dot = 0, |0| > 0.921 => false)
    {
        normal_map.clear();
        normal_map[&res1] = clipper::Coord_orth(1.0, 0.0, 0.0);
        normal_map[&res2] = clipper::Coord_orth(0.0, 1.0, 0.0);
        bool result = coot::stack_and_pair::similar_normals_gemmi(&res1, &res2, normal_map, angle_crit);
        EXPECT_FALSE(result);
    }
    
    // Case 4: res2 not in normal_map (guard clause: it2 == end, returns false)
    {
        normal_map.clear();
        normal_map[&res1] = clipper::Coord_orth(1.0, 0.0, 0.0);
        bool result = coot::stack_and_pair::similar_normals_gemmi(&res1, &res2, normal_map, angle_crit);
        EXPECT_FALSE(result);
    }
    
    // Case 5: res1 not in normal_map (guard clause: it1 == end, returns false)
    {
        normal_map.clear();
        normal_map[&res2] = clipper::Coord_orth(1.0, 0.0, 0.0);
        bool result = coot::stack_and_pair::similar_normals_gemmi(&res1, &res2, normal_map, angle_crit);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}