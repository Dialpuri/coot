#include "function.hh"

#include <iostream>
#include <vector>
#include "geometry/protein-geometry.hh"

bool
coot::dict_link_info_t::check_for_order_switch_gemmi(const gemmi::Residue* residue_ref,
                                                     const gemmi::Residue* residue_new,
                                                     const std::string& link_type,
                                                     const coot::protein_geometry& geom) const {
  bool order_switch_flag = false;
  std::string comp_id_ref = residue_ref->name;
  std::string comp_id_new = residue_new->name;

  try {
    std::string group_ref = geom.get_group(comp_id_ref);
    std::string group_new = geom.get_group(comp_id_new);

    std::vector<std::pair<coot::chem_link, bool>> link_infos;
    std::vector<coot::chem_link> link_infos_f =
        geom.matching_chem_links(comp_id_ref, group_ref, comp_id_new, group_new);
    std::vector<coot::chem_link> link_infos_b =
        geom.matching_chem_links(comp_id_new, group_new, comp_id_ref, group_ref);

    for (const auto& link : link_infos_f)
      link_infos.push_back(std::make_pair(link, false));
    for (const auto& link : link_infos_b)
      link_infos.push_back(std::make_pair(link, true));

    for (unsigned int ilink = 0; ilink < link_infos.size(); ilink++) {
      if (link_infos[ilink].first.Id() == link_type) {
        order_switch_flag = link_infos[ilink].second;
        break;
      }
    }
  } catch (const std::runtime_error& rte) {
    std::cout << "WARNING:: check_for_order_switch() exception: " << rte.what()
              << std::endl;
  }
  return order_switch_flag;
}