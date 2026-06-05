#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

TEST(AddExtraTargetPositionRestraintsGemmi, EmptyRestraints) {
   restraints_container_t rc;
   rc.from_residue_vector = true;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;

   extra_restraints_t extra_restraints;
   add_extra_target_position_restraints_gemmi(rc, extra_restraints);

   EXPECT_EQ(rc.restraints_vec.size(), 0u);
}

TEST(AddExtraTargetPositionRestraintsGemmi, NonFixedAtom) {
   gemmi::Structure st = gemmi::read_pdb_file(
       "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   gemmi::Model& model = st.models[0];
   gemmi::Chain& chain = model.chains[0];
   gemmi::Residue& res = chain.residues[0];
   gemmi::Atom& atom = res.atoms[0];

   restraints_container_t rc;
   rc.from_residue_vector = true;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;

   // Register atom in atom_map with index 0
   atom_spec_t spec(chain.name, res.seqid.num.value,
                    std::string(1, res.seqid.icode),
                    atom.name, std::string(1, atom.altloc));
   rc.atom_map[spec] = 0;

   // Non-fixed atom
   extra_restraints_t extra;
   clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
   double weight = 1.0;
   extra.target_position_restraints.push_back(
       extra_restraints_t::extra_target_position_restraint_t(spec, pos, weight)
   );

   EXPECT_EQ(rc.restraints_vec.size(), 0u);
   add_extra_target_position_restraints_gemmi(rc, extra);
   EXPECT_EQ(rc.restraints_vec.size(), 1u);
}

TEST(AddExtraTargetPositionRestraintsGemmi, FixedAtom) {
   gemmi::Structure st = gemmi::read_pdb_file(
       "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   gemmi::Model& model = st.models[0];
   gemmi::Chain& chain = model.chains[0];
   gemmi::Residue& res = chain.residues[0];
   gemmi::Atom& atom = res.atoms[0];

   restraints_container_t rc;
   rc.from_residue_vector = true;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;

   // Register atom in atom_map with index 0
   atom_spec_t spec(chain.name, res.seqid.num.value,
                    std::string(1, res.seqid.icode),
                    atom.name, std::string(1, atom.altloc));
   rc.atom_map[spec] = 0;

   // Mark atom as fixed
   rc.fixed_atom_indices.insert(0);

   extra_restraints_t extra;
   clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
   double weight = 1.0;
   extra.target_position_restraints.push_back(
       extra_restraints_t::extra_target_position_restraint_t(spec, pos, weight)
   );

   EXPECT_EQ(rc.restraints_vec.size(), 0u);
   add_extra_target_position_restraints_gemmi(rc, extra);
   // Fixed atom: restraint should NOT be added
   EXPECT_EQ(rc.restraints_vec.size(), 0u);
}

TEST(AddExtraTargetPositionRestraintsGemmi, WrongAtomName) {
   gemmi::Structure st = gemmi::read_pdb_file(
       "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   gemmi::Model& model = st.models[0];
   gemmi::Chain& chain = model.chains[0];
   gemmi::Residue& res = chain.residues[0];
   gemmi::Atom& atom = res.atoms[0];

   restraints_container_t rc;
   rc.from_residue_vector = true;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;

   // Register atom in atom_map with index 0
   atom_spec_t spec(chain.name, res.seqid.num.value,
                    std::string(1, res.seqid.icode),
                    atom.name, std::string(1, atom.altloc));
   rc.atom_map[spec] = 0;

   // Use a wrong atom name "ZZZ"
   atom_spec_t wrong_spec(chain.name, res.seqid.num.value,
                          std::string(1, res.seqid.icode),
                          "ZZZ", std::string(1, atom.altloc));

   extra_restraints_t extra;
   clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
   double weight = 1.0;
   extra.target_position_restraints.push_back(
       extra_restraints_t::extra_target_position_restraint_t(wrong_spec, pos, weight)
   );

   EXPECT_EQ(rc.restraints_vec.size(), 0u);
   add_extra_target_position_restraints_gemmi(rc, extra);
   // Wrong atom name: restraint should NOT be added
   EXPECT_EQ(rc.restraints_vec.size(), 0u);
}

} // namespace coot

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}