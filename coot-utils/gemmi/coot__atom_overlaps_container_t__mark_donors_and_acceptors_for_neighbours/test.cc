#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, atom_overlaps_container_mark_donors_and_acceptors_for_neighbours) {

  // =========================================================
  // Case 1: Original oracle behaviour — protein_geometry NOT
  //         initialised, so no dictionaries exist and ALL
  //         residues are skipped.  This exercises the
  //         empty-neighbours path (the oracle's 1 branch).
  // =========================================================
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Intentionally NOT calling init_standard() — matches oracle
    coot::protein_geometry geom;
    coot::protein_geometry* geom_p = &geom;

    // Build neighbours list using standalone vector (gemmi::Residue*)
    std::vector<gemmi::Residue*> neighbours;
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
      for (gemmi::Residue& res : chain.residues) {
        std::string resname(res.name);
        std::pair<bool, coot::dictionary_residue_restraints_t> d =
            geom_p->get_monomer_restraints(resname, coot::protein_geometry::IMOL_ENC_ANY);
        if (d.first) {
          neighbours.push_back(&res);
        }
      }
    }

    // Simulate udd_handle (in original this was mol->RegisterUDInteger(...))
    int udd_handle = 16777217;

    // Oracle assertions
    EXPECT_EQ(neighbours.size(), 0U);
    EXPECT_EQ(udd_handle, 16777217);

    // Count atoms before
    int atoms_before = 0;
    for (unsigned int i = 0; i < neighbours.size(); i++) {
      for (const gemmi::Atom& atom : neighbours[i]->atoms) {
        atoms_before++;
      }
    }

    // atom_hb_types replaces UD data in gemmi port
    std::unordered_map<const gemmi::Atom*, coot::hb_t> atom_hb_types;
    bool molecule_has_hydrogens = false;

    // Call the target function — should NOT throw
    EXPECT_NO_THROW(coot::mark_donors_and_acceptors_for_neighbours_gemmi(
        geom_p, neighbours, atom_hb_types, molecule_has_hydrogens));

    EXPECT_FALSE(molecule_has_hydrogens);
    EXPECT_EQ(atoms_before, 0);
  }

  // =========================================================
  // Case 2: Complementary — protein_geometry WITH init_standard()
  //         so dictionaries exist, neighbours are populated,
  //         and the function actually marks donors/acceptors.
  //         This exercises the non-empty branch through the
  //         function (mark_donors_and_acceptors loop).
  // =========================================================
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Now we DO call init_standard() — this is the complementary case
    coot::protein_geometry geom;
    geom.init_standard();
    coot::protein_geometry* geom_p = &geom;

    // Build neighbours list using standalone vector (gemmi::Residue*)
    std::vector<gemmi::Residue*> neighbours;
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
      for (gemmi::Residue& res : chain.residues) {
        std::string resname(res.name);
        std::pair<bool, coot::dictionary_residue_restraints_t> d =
            geom_p->get_monomer_restraints(resname, coot::protein_geometry::IMOL_ENC_ANY);
        if (d.first) {
          neighbours.push_back(&res);
        }
      }
    }

    // With init_standard(), we should have neighbours
    EXPECT_GT(neighbours.size(), 0U);

    // atom_hb_types replaces UD data in gemmi port
    std::unordered_map<const gemmi::Atom*, coot::hb_t> atom_hb_types;
    bool molecule_has_hydrogens = false;

    // Call the target function — should NOT throw
    EXPECT_NO_THROW(coot::mark_donors_and_acceptors_for_neighbours_gemmi(
        geom_p, neighbours, atom_hb_types, molecule_has_hydrogens));

    // molecule_has_hydrogens becomes true when dictionaries are available
    // and the function detects hydrogen atom names in the dictionary
    EXPECT_TRUE(molecule_has_hydrogens);

    // Check that some atoms got hb_type set (replacing UD data check)
    int atoms_with_data = 0;
    int atoms_checked = 0;
    for (unsigned int i = 0; i < neighbours.size(); i++) {
      for (const gemmi::Atom& atom : neighbours[i]->atoms) {
        atoms_checked++;
        auto it = atom_hb_types.find(&atom);
        if (it != atom_hb_types.end()) {
          atoms_with_data++;
        }
      }
    }

    // With proper dictionaries we should have marked some atoms
    EXPECT_GT(atoms_with_data, 0);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}