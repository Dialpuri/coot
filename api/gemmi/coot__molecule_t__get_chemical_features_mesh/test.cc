#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(GetChemicalFeaturesMesh, GemmiPort) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   const gemmi::Model& model = st.models[0];
   const gemmi::Chain& chain = model.chains[0];
   const gemmi::Residue& residue = chain.residues[0];
   coot::simple_mesh_t mesh = coot::molecule_t_get_chemical_features_mesh_gemmi(residue, chain.name);
   EXPECT_EQ(mesh.vertices.size(), 0u);
   EXPECT_EQ(mesh.triangles.size(), 0u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}