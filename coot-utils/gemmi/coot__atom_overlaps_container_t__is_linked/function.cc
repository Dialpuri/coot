#include "function.hh"

namespace coot {

bool atom_overlaps_container_t::is_linked_gemmi(gemmi::CRA cra1,
                                                gemmi::CRA cra2,
                                                gemmi::Model* model,
                                                gemmi::Structure* structure) const {
  bool status = false;
  if (!cra1.atom) return false;
  if (!cra2.atom) return false;

  if (model && structure) {
    for (const auto& conn : structure->connections) {
      std::pair<atom_spec_t, atom_spec_t> atoms =
        coot::link_atoms_gemmi(conn, model);

      atom_spec_t spec_1(
        cra1.chain->name,
        cra1.residue->seqid.num.value,
        std::string(1, cra1.residue->seqid.icode),
        cra1.atom->name,
        std::string(1, cra1.atom->altloc)
      );
      atom_spec_t spec_2(
        cra2.chain->name,
        cra2.residue->seqid.num.value,
        std::string(1, cra2.residue->seqid.icode),
        cra2.atom->name,
        std::string(1, cra2.atom->altloc)
      );

      if (spec_1 == atoms.first) {
        if (spec_2 == atoms.second) {
          status = true;
          break;
        }
      }
      if (spec_2 == atoms.first) {
        if (spec_1 == atoms.second) {
          status = true;
          break;
        }
      }
    }
  }
  return status;
}

} // namespace coot