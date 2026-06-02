#include "function.hh"
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

int make_header_metal_links_ng_gemmi(
    const std::vector<link_info_t> &links,
    const coot::protein_geometry &geom) {

  int n_bonds_added = 0;
  int imol = protein_geometry::IMOL_ENC_ANY;

  for (std::size_t i = 0; i < links.size(); i++) {
    const link_info_t &link = links[i];

    atom_spec_t a1(link.chainID1, link.seqNum1, link.insCode1, link.atName1, link.aloc1);
    atom_spec_t a2(link.chainID2, link.seqNum2, link.insCode2, link.atName2, link.aloc2);

    if ((a1.alt_conf == a2.alt_conf) || a1.alt_conf.empty() || a2.alt_conf.empty()) {
      std::string rn_1 = link.resName1;
      std::string rn_2 = link.resName2;

      bool is_oxygen_a1   = false;
      bool is_oxygen_a2   = false;
      bool is_nitrogen_a1 = false;
      bool is_nitrogen_a2 = false;
      bool is_sulfur_a1   = false;
      bool is_sulfur_a2   = false;

      std::pair<bool, dict_atom> da_1 = geom.get_monomer_atom_info(rn_1, a1.atom_name, imol);
      std::pair<bool, dict_atom> da_2 = geom.get_monomer_atom_info(rn_2, a2.atom_name, imol);

      if (da_1.first) {
        if (da_1.second.type_symbol == "O") is_oxygen_a1   = true;
        if (da_1.second.type_symbol == "S") is_sulfur_a1   = true;
        if (da_1.second.type_symbol == "N") is_nitrogen_a1 = true;
      }
      if (da_2.first) {
        if (da_2.second.type_symbol == "O") is_oxygen_a2   = true;
        if (da_2.second.type_symbol == "S") is_sulfur_a2   = true;
        if (da_2.second.type_symbol == "N") is_nitrogen_a2 = true;
      }

      // Check a1 -> a2 direction (O, N, S on a1, metal on a2)
      if (is_oxygen_a1) {
        auto it = geom.metal_O_map.find(rn_2);
        if (it != geom.metal_O_map.end()) {
          n_bonds_added++;
        }
      }
      if (is_nitrogen_a1) {
        auto it = geom.metal_N_map.find(rn_2);
        if (it != geom.metal_N_map.end()) {
          n_bonds_added++;
        }
      }
      if (is_sulfur_a1) {
        auto it = geom.metal_S_map.find(rn_2);
        if (it != geom.metal_S_map.end()) {
          n_bonds_added++;
        }
      }

      // Check a2 -> a1 direction (O, N, S on a2, metal on a1)
      if (is_oxygen_a2) {
        auto it = geom.metal_O_map.find(rn_1);
        if (it != geom.metal_O_map.end()) {
          n_bonds_added++;
        }
      }
      if (is_nitrogen_a2) {
        auto it = geom.metal_N_map.find(rn_1);
        if (it != geom.metal_N_map.end()) {
          n_bonds_added++;
        }
      }
      if (is_sulfur_a2) {
        auto it = geom.metal_S_map.find(rn_1);
        if (it != geom.metal_S_map.end()) {
          n_bonds_added++;
        }
      }
    }
  }

  return n_bonds_added;
}

} // namespace coot