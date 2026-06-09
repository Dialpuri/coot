#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot/geometry/protein-geometry.hh"

TEST(OracleTest, get_vdw_radius_neighb_atom) {
  // Initialize protein geometry
  coot::protein_geometry geom;
  geom.init_standard();

  // Load PDB with gemmi
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  ASSERT_TRUE(!st.models.empty());

  gemmi::Model& model = st.models[0];
  gemmi::Chain& chain = model.chains[0];
  gemmi::Residue& res = chain.residues[0];

  ASSERT_FALSE(res.atoms.empty());

  // Get first atom to test with
  const gemmi::Atom& test_atom = res.atoms[0];
  std::string res_name(res.name);

  // Get dictionary for this residue
  auto dict_result = geom.get_monomer_restraints(res_name, 0);

  // Build the dictionaries vector and cache maps (mirroring atom_overlaps_container_t)
  std::vector<coot::dictionary_residue_restraints_t> neighb_dictionaries;
  std::map<const gemmi::Atom*, double> neighbour_atoms_vdw_radius_map;
  std::map<std::string, double> type_to_vdw_radius_map;

  if (dict_result.first) {
    neighb_dictionaries.push_back(dict_result.second);
  } else {
    FAIL() << "Could not get dictionary for " << res_name;
  }

  // Case 1: first call — cache miss
  {
    double r = coot::get_vdw_radius_neighb_atom_gemmi(
        &test_atom, 0, neighb_dictionaries, &geom,
        neighbour_atoms_vdw_radius_map, type_to_vdw_radius_map);
    EXPECT_NEAR(r, 1.55, 1e-4);
  }

  // Case 2: same atom again — cache hit (should be same value)
  {
    double r2 = coot::get_vdw_radius_neighb_atom_gemmi(
        &test_atom, 0, neighb_dictionaries, &geom,
        neighbour_atoms_vdw_radius_map, type_to_vdw_radius_map);
    EXPECT_NEAR(r2, 1.55, 1e-4);
  }

  // Case 3: different atom in same residue
  if (res.atoms.size() > 1) {
    const gemmi::Atom& atom2 = res.atoms[1];
    double r3 = coot::get_vdw_radius_neighb_atom_gemmi(
        &atom2, 0, neighb_dictionaries, &geom,
        neighbour_atoms_vdw_radius_map, type_to_vdw_radius_map);
    EXPECT_NEAR(r3, 1.7, 1e-4);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
