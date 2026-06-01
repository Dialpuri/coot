#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cstring>
#include <gsl/gsl_vector.h>
#include "function.hh"

TEST(OracleTest, setup_simplex_x_internal_rigid) {
  // Load the PDB using gemmi
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // Get the first model
  gemmi::Model& model = st.models[0];

  // Find chain A
  gemmi::Chain* chain = nullptr;
  for (auto& c : model.chains) {
    if (c.name == "A") {
      chain = &c;
      break;
    }
  }
  EXPECT_NE(chain, nullptr);

  // Find residue with seqid 1 (MET)
  gemmi::Residue* res = nullptr;
  for (auto& r : chain->residues) {
    if (r.seqid.num.value == 1) {
      res = &r;
      break;
    }
  }
  EXPECT_NE(res, nullptr);

  int n_atoms_in_res = static_cast<int>(res->atoms.size());
  EXPECT_EQ(n_atoms_in_res, 8);

  // Build atom_selection from the residue's atoms
  std::vector<gemmi::Atom*> atom_selection(n_atoms_in_res);
  for (int i = 0; i < n_atoms_in_res; i++) {
    atom_selection[i] = &res->atoms[i];
  }

  // Manually create a gsl_vector (6 elements for rigid body: 3 rotations + 3 translations)
  int n_params = 6;

  // Allocate gsl_vector struct and its data block manually
  gsl_vector *x = (gsl_vector*)std::malloc(sizeof(gsl_vector));
  double *data = (double*)std::malloc(n_params * sizeof(double));

  std::memset(x, 0, sizeof(gsl_vector));
  x->size = n_params;
  x->stride = 1;
  x->data = data;

  // Initialize data manually
  for (int i = 0; i < n_params; i++) {
    data[i] = 1.0 + i;
  }

  // Verify inputs match oracle
  EXPECT_EQ(n_atoms_in_res, 8);
  EXPECT_EQ(n_params, 6);
  EXPECT_STREQ(chain->name.c_str(), "A");
  EXPECT_EQ(res->seqid.num.value, 1);
  EXPECT_STREQ(res->name.c_str(), "MET");

  // Verify x_before
  for (int i = 0; i < n_params; i++) {
    EXPECT_NEAR(data[i], 1.0 + i, 1e-4);
  }

  // Call the function
  EXPECT_NO_THROW(coot::util::setup_simplex_x_internal_rigid_gemmi(x, atom_selection, n_atoms_in_res));

  // Verify x_after (function is a no-op, so values unchanged)
  for (int i = 0; i < n_params; i++) {
    EXPECT_NEAR(data[i], 1.0 + i, 1e-4);
  }

  // Cleanup
  std::free(data);
  std::free(x);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}