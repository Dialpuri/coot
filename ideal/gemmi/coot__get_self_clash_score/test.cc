#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <cmath>
#include "function.hh"

// Helper: collect all atoms in a model as gemmi::CRA
static std::vector<gemmi::CRA> collect_all_atoms_cra(gemmi::Model& model) {
  std::vector<gemmi::CRA> all_cra;
  for (gemmi::Chain& chain : model.chains) {
    for (gemmi::Residue& res : chain.residues) {
      for (gemmi::Atom& atom : res.atoms) {
        all_cra.push_back(gemmi::CRA{&chain, &res, &atom});
      }
    }
  }
  return all_cra;
}

TEST(OracleTest, get_self_clash_score_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1);

    gemmi::Model& model = st.models[0];
    
    // Select all atoms as CRA
    std::vector<gemmi::CRA> atom_selection = collect_all_atoms_cra(model);

    EXPECT_EQ(atom_selection.size(), 2107);

    // Case 2: empty quads vector (all contacts counted)
    {
        std::vector<coot::torsion_atom_quad_gemmi> quads_empty;

        double score = coot::get_self_clash_score_gemmi(st, model, atom_selection, quads_empty);

        EXPECT_EQ(quads_empty.size(), 0);
        EXPECT_NEAR(score, 710.485, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}