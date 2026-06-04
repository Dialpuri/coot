#include "function.hh"

#include <algorithm>

namespace coot {

bool is_in_same_ring_gemmi(
    int imol,
    gemmi::Residue const& residue,
    std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>>& residue_ring_map_cache,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const coot::protein_geometry& geom)
{
  bool r = false;

  std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>> residue_ring_map;
  std::list<std::string> r1;
  std::list<std::string> r2;
  std::list<std::string> r3;
  std::list<std::string> r4;

  // HIS
  r1.push_back(" CG ");
  r1.push_back(" CD2");
  r1.push_back(" ND1");
  r1.push_back(" CE1");
  r1.push_back(" NE2");

  // PHE/TYR
  r2.push_back(" CG ");
  r2.push_back(" CD1");
  r2.push_back(" CD2");
  r2.push_back(" CE1");
  r2.push_back(" CE2");
  r2.push_back(" CZ ");

  // TRP
  r3.push_back(" CG ");
  r3.push_back(" CD1");
  r3.push_back(" CD2");
  r3.push_back(" CE2");
  r3.push_back(" NE1");

  r4.push_back(" CD2");
  r4.push_back(" CE2");
  r4.push_back(" CE3");
  r4.push_back(" CZ2");
  r4.push_back(" CZ3");
  r4.push_back(" CH2");

  if (residue_ring_map_cache.size() == 0) {
    r1.sort();
    r2.sort();
    r3.sort();
    r4.sort();
    residue_ring_map["HIS"].second.push_back(r1);
    residue_ring_map["PHE"].second.push_back(r2);
    residue_ring_map["TYR"].second.push_back(r2);
    residue_ring_map["TRP"].second.push_back(r3);
    residue_ring_map["TRP"].second.push_back(r4);
    residue_ring_map["HIS"].first = false;
    residue_ring_map["PHE"].first = false;
    residue_ring_map["TYR"].first = false;
    residue_ring_map["TRP"].first = false;
  }

  std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>>::const_iterator it;
  std::string res_name = residue.name; // gemmi: .name field instead of GetResName()

  it = residue_ring_map_cache.find(res_name);
  if (it != residue_ring_map_cache.end()) {

    if (it->second.first == 0) { // not looked up before and failed
      for (unsigned int i = 0; i < it->second.second.size(); i++) {
        std::list<std::string>::const_iterator it_1 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_1);
        std::list<std::string>::const_iterator it_2 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_2);
        if (it_1 != it->second.second[i].end()) {
          if (it_2 != it->second.second[i].end()) {
            r = true;
            break;
          }
        }
      }
    } else {
      // We tried to look it up before and failed
    }
  } else {

    // add it then
    std::pair<bool, dictionary_residue_restraints_t> rest = geom.get_monomer_restraints(res_name, imol);
    const auto& dict = rest.second;
    if (rest.first) {
      std::vector<std::vector<std::string>> ri = dict.get_ligand_ring_list();
      residue_ring_map_cache[res_name].first = false; // not looked up before and failed
      for (unsigned int ii = 0; ii < ri.size(); ii++) {
        std::list<std::string> l;
        for (unsigned int jj = 0; jj < ri[ii].size(); jj++)
          l.push_back(ri[ii][jj]);
        l.sort();
        residue_ring_map_cache[res_name].second.push_back(l);
      }

      std::vector<std::list<std::string>>& vl = residue_ring_map_cache[res_name].second;

      for (unsigned int ii = 0; ii < vl.size(); ii++) {
        std::list<std::string>::const_iterator it_1 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_1);
        std::list<std::string>::const_iterator it_2 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_2);
        if (it_1 != vl[ii].end()) {
          if (it_2 != vl[ii].end()) {
            r = true;
            break;
          }
        }
      }
    } else {
      // OK, the lookup failed
      std::vector<std::list<std::string>> fv;
      std::pair<bool, std::vector<std::list<std::string>>> failed_data(true, fv);
      residue_ring_map_cache[res_name] = failed_data;
    }
  }
  return r;
}

} // namespace coot