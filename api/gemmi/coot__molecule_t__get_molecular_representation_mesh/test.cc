#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_molecular_representation_mesh) {
   coot::molecule_t mol;
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A", "Chains", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A/10-20", "Chains", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//", "Chains", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A", "Element", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A/10-20", "Element", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//", "Element", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A", "Secondary", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A/10-20", "Secondary", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//", "Secondary", "Ribbon", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A", "Chains", "Tubes", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//A/10-20", "Chains", "Tubes", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
   
   {
      coot::simple_mesh_t mesh = mol.get_molecular_representation_mesh_gemmi("//", "Chains", "Tubes", 0);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.vertices.size(), 0U);
      EXPECT_EQ(mesh.triangles.size(), 0U);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}