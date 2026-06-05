#include "function.hh"

#include <vector>
#include <string>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

// Helper: pad atom name to 4 characters matching dictionary convention
// All names get a leading space then right-padded to 4.
static std::string pad_to_4(const std::string &name) {
  std::string result = " " + name;
  if (result.size() < 4)
    result.resize(4, ' ');
  return result;
}

// Inline version of get_asc_index_gemmi — avoids including the generated
// get_asc_index port which redefines atom_spec_t and conflicts with the
// real coot::atom_spec_t from residue-and-atom-specs.hh.
static int get_cra_index(const std::vector<gemmi::CRA> &cras,
                         const coot::atom_spec_t &spec) {
  for (int i = 0; i < static_cast<int>(cras.size()); i++) {
    const gemmi::CRA &cra = cras[i];
    if (cra.atom != nullptr && cra.residue != nullptr && cra.chain != nullptr) {
      if (cra.atom->name == spec.atom_name &&
          std::string(1, cra.atom->altloc) == spec.alt_conf &&
          cra.residue->seqid.num.value == spec.res_no &&
          std::string(1, cra.residue->seqid.icode) == spec.ins_code &&
          cra.chain->name == spec.chain_id) {
        return i;
      }
    }
  }
  return -1;
}

int add_link_torsion_for_phi_psi_gemmi(
    std::vector<coot::simple_restraint> &restraints_vec,
    const std::vector<gemmi::CRA> &cras,
    std::string link_type,
    int phi_psi_restraints_type,
    gemmi::CRA first,
    gemmi::CRA second,
    short int is_fixed_first,
    short int is_fixed_second,
    const coot::protein_geometry &geom) {

  int n_torsion = 0;

  if (!first.residue || !second.residue || !first.chain || !second.chain) {
    return 0;
  }

  std::vector<bool> fixed_flag(4, false);

  for (int i = 0; i < geom.link_size(); i++) {
    if (geom.link(i).link_id == link_type) {
      for (unsigned int j = 0; j < geom.link(i).link_torsion_restraint.size(); j++) {

        const auto &torsion = geom.link(i).link_torsion_restraint[j];

        gemmi::Residue *atom_1_res = nullptr;
        gemmi::Residue *atom_2_res = nullptr;
        gemmi::Residue *atom_3_res = nullptr;
        gemmi::Residue *atom_4_res = nullptr;
        gemmi::Chain *atom_1_chain = nullptr;
        gemmi::Chain *atom_2_chain = nullptr;
        gemmi::Chain *atom_3_chain = nullptr;
        gemmi::Chain *atom_4_chain = nullptr;

        if (torsion.atom_1_comp_id == 1) {
          atom_1_res = first.residue;
          atom_1_chain = first.chain;
          fixed_flag[0] = is_fixed_first;
        } else {
          atom_1_res = second.residue;
          atom_1_chain = second.chain;
          fixed_flag[0] = is_fixed_second;
        }
        if (torsion.atom_2_comp_id == 1) {
          atom_2_res = first.residue;
          atom_2_chain = first.chain;
          fixed_flag[1] = is_fixed_first;
        } else {
          atom_2_res = second.residue;
          atom_2_chain = second.chain;
          fixed_flag[1] = is_fixed_second;
        }
        if (torsion.atom_3_comp_id == 1) {
          atom_3_res = first.residue;
          atom_3_chain = first.chain;
          fixed_flag[2] = is_fixed_first;
        } else {
          atom_3_res = second.residue;
          atom_3_chain = second.chain;
          fixed_flag[2] = is_fixed_second;
        }
        if (torsion.atom_4_comp_id == 1) {
          atom_4_res = first.residue;
          atom_4_chain = first.chain;
          fixed_flag[3] = is_fixed_first;
        } else {
          atom_4_res = second.residue;
          atom_4_chain = second.chain;
          fixed_flag[3] = is_fixed_second;
        }

        const std::string &target1 = torsion.atom_id_1_4c();
        const std::string &target2 = torsion.atom_id_2_4c();
        const std::string &target3 = torsion.atom_id_3_4c();
        const std::string &target4 = torsion.atom_id_4_4c();

        for (gemmi::Atom &a1 : atom_1_res->atoms) {
          if (pad_to_4(a1.name) == target1) {
            for (gemmi::Atom &a2 : atom_2_res->atoms) {
              if (pad_to_4(a2.name) == target2) {
                for (gemmi::Atom &a3 : atom_3_res->atoms) {
                  if (pad_to_4(a3.name) == target3) {
                    for (gemmi::Atom &a4 : atom_4_res->atoms) {
                      if (pad_to_4(a4.name) == target4) {

                        auto make_spec = [&](gemmi::Atom &atom,
                                             gemmi::Chain *ch,
                                             gemmi::Residue *res) {
                          coot::atom_spec_t spec;
                          spec.atom_name = pad_to_4(atom.name);
                          spec.alt_conf = std::string(1, atom.altloc);
                          spec.res_no = res->seqid.num.value;
                          spec.ins_code = std::string(1, res->seqid.icode);
                          spec.chain_id = ch->name;
                          return spec;
                        };

                        int index1 = get_cra_index(cras, make_spec(a1, atom_1_chain, atom_1_res));
                        int index2 = get_cra_index(cras, make_spec(a2, atom_2_chain, atom_2_res));
                        int index3 = get_cra_index(cras, make_spec(a3, atom_3_chain, atom_3_res));
                        int index4 = get_cra_index(cras, make_spec(a4, atom_4_chain, atom_4_res));

                        double target_phi = -57.82 + 360.0;
                        double target_psi = -47.0 + 360.0;
                        double esd = 5.01;

                        if (phi_psi_restraints_type ==
                            coot::restraints_container_t::LINK_TORSION_ALPHA_HELIX) {
                          target_phi = -57.82 + 360.0;
                          target_psi = -47.00 + 360.0;
                        }

                        if (phi_psi_restraints_type ==
                            coot::restraints_container_t::LINK_TORSION_BETA_STRAND) {
                          target_phi = -110.0 + 360.0;
                          target_psi = 120.0;
                          esd = 15.0;
                        }

                        if (torsion.id() == "phi") {
                          restraints_vec.push_back(coot::simple_restraint(
                              restraint_type_t::TORSION_RESTRAINT, index1, index2, index3, index4,
                              fixed_flag,
                              static_cast<float>(target_phi),
                              static_cast<float>(esd),
                              1.2f,  // weight
                              1));   // periodicity
                          n_torsion++;
                        }
                        if (torsion.id() == "psi") {
                          restraints_vec.push_back(coot::simple_restraint(
                              restraint_type_t::TORSION_RESTRAINT, index1, index2, index3, index4,
                              fixed_flag,
                              static_cast<float>(target_psi),
                              static_cast<float>(esd),
                              1.2f,  // weight
                              1));   // periodicity
                          n_torsion++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return n_torsion;
}

}  // namespace coot