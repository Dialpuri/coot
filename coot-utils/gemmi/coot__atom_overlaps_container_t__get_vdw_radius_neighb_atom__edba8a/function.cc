#include "function.hh"

#include "coot/geometry/protein-geometry.hh"

namespace {

// Pad trimmed atom name to 4 chars using MMDB's atom_id_4c convention:
// len 1 → "X  ", len 2 → " XY", len 3 → " XYZ", len 4 → "ABCD"
std::string atom_name_4c(const std::string& name) {
  std::string trimmed = name;
  size_t start = trimmed.find_first_not_of(' ');
  size_t end = trimmed.find_last_not_of(' ');
  if (start == std::string::npos) return "    ";
  trimmed = trimmed.substr(start, end - start + 1);
  int ilen = static_cast<int>(trimmed.size());
  if (ilen >= 4) return trimmed.substr(0, 4);
  std::string leading(1, ' ');
  if (ilen == 1) return leading + trimmed + "  ";
  if (ilen == 2) return leading + trimmed + " ";
  if (ilen == 3) return leading + trimmed;
  return trimmed;
}

} // anonymous namespace

double coot::get_vdw_radius_neighb_atom_gemmi(
    const gemmi::Atom* at,
    unsigned int idx_res,
    const std::vector<dictionary_residue_restraints_t>& neighb_dictionaries,
    const protein_geometry* geom_p,
    std::map<const gemmi::Atom*, double>& neighbour_atoms_vdw_radius_map,
    std::map<std::string, double>& type_to_vdw_radius_map)
{
  double r = 1.5;

  auto it = neighbour_atoms_vdw_radius_map.find(at);
  if (it == neighbour_atoms_vdw_radius_map.end()) {
    // cache miss — look up the vdw radius
    // type_energy() compares against atom_id_4c(). gemmi atom.name preserves
    // PDB column 13-16 format (" N  ") while the dictionary uses trimmed
    // then right-padded ("N   "), so we normalize.
    std::string te = neighb_dictionaries[idx_res].type_energy(atom_name_4c(at->name));

    auto it_type = type_to_vdw_radius_map.find(te);
    if (it_type == type_to_vdw_radius_map.end()) {
      r = geom_p->get_energy_lib_atom(te).vdw_radius;
      type_to_vdw_radius_map[te] = r;
    } else {
      r = it_type->second;
    }
    neighbour_atoms_vdw_radius_map[at] = r;
  } else {
    r = it->second;
  }

  return r;
}
