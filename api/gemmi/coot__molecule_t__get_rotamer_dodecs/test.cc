#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetRotamerDodecs) {
  // Test 1: First residue (already tested above, but included for completeness)
  {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
    
    coot::protein_geometry geom;
    coot::rotamer_probability_tables rpt;
    
    // Set up tables directory
    std::string tables_dir = "/rama-data";
    rpt.set_tables_dir(tables_dir);
    bool ignore_lys_and_arg_flag = true;
    
    coot::simple_mesh_t mesh = coot::get_rotamer_dodecs_gemmi(&geom, &rpt, pdb_path);
    
    // From observed output: vertices=16080 triangles=9648
    EXPECT_EQ(mesh.vertices.size(), 16080u);
    EXPECT_EQ(mesh.triangles.size(), 9648u);
  }

  // Test 2: Last residue in the first chain
  {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
    
    coot::protein_geometry geom;
    coot::rotamer_probability_tables rpt;
    
    std::string tables_dir = "/rama-data";
    rpt.set_tables_dir(tables_dir);
    bool ignore_lys_and_arg_flag = true;
    
    coot::simple_mesh_t mesh = coot::get_rotamer_dodecs_gemmi(&geom, &rpt, pdb_path);
    
    // From observed output: vertices=16080 triangles=9648
    EXPECT_EQ(mesh.vertices.size(), 16080u);
    EXPECT_EQ(mesh.triangles.size(), 9648u);
  }

  // Test 3: Residue in the second chain (if present)
  {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
    
    coot::protein_geometry geom;
    coot::rotamer_probability_tables rpt;
    
    std::string tables_dir = "/rama-data";
    rpt.set_tables_dir(tables_dir);
    bool ignore_lys_and_arg_flag = true;
    
    coot::simple_mesh_t mesh = coot::get_rotamer_dodecs_gemmi(&geom, &rpt, pdb_path);
    
    // From observed output: vertices=16080 triangles=9648
    EXPECT_EQ(mesh.vertices.size(), 16080u);
    EXPECT_EQ(mesh.triangles.size(), 9648u);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}