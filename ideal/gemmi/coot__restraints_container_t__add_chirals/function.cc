#include "function.hh"

#include <algorithm>
#include <iostream>

using namespace coot;

// PDB atom name to 4-char dictionary format:
//   1-2 char names: right-justified in cols 14-16, col 13 is always space
//     "N"  -> " N  ",  "CA" -> " CA "
//   3-4 char names: left-justified in cols 13-16
//     "HB2" -> "HB2 "
static std::string pad_name_to_four_chars(const std::string &name) {
  std::string result = name;
  if (result.size() <= 2)
    result = " " + result;
  while (result.size() < 4)
    result.push_back(' ');
  if (result.size() > 4)
    result.resize(4);
  return result;
}

int coot::add_chirals_gemmi(
    int idr,
    std::vector<gemmi::Atom *> res_selection,
    int i_no_res_atoms,
    gemmi::CRA cra,
    const coot::protein_geometry &geom,
    std::vector<coot::simple_restraint> &restraints_vec)
{
  int n_chiral_restr = 0;
  int index1, index2, index3, indexc;

  // Build PDB-formatted 4-char atom names to match dictionary format
  std::vector<std::string> string_atom_names(i_no_res_atoms);
  for (int iat = 0; iat < i_no_res_atoms; iat++) {
    string_atom_names[iat] = pad_name_to_four_chars(res_selection[iat]->name);
  }

  for (unsigned int ic = 0; ic < geom[idr].second.chiral_restraint.size(); ic++) {
    const dict_chiral_restraint_t &dcr =
        geom[idr].second.chiral_restraint[ic];

    if (!dcr.is_a_both_restraint()) {
      for (int iat1 = 0; iat1 < i_no_res_atoms; iat1++) {
        const std::string &pdb_atom_name1 = string_atom_names[iat1];
        if (pdb_atom_name1 == dcr.atom_id_1_4c()) {

          for (int iat2 = 0; iat2 < i_no_res_atoms; iat2++) {
            const std::string &pdb_atom_name2 = string_atom_names[iat2];
            if (pdb_atom_name2 == dcr.atom_id_2_4c()) {

              for (int iat3 = 0; iat3 < i_no_res_atoms; iat3++) {
                const std::string &pdb_atom_name3 = string_atom_names[iat3];
                if (pdb_atom_name3 == dcr.atom_id_3_4c()) {

                  for (int iatc = 0; iatc < i_no_res_atoms; iatc++) {
                    const std::string &pdb_atom_namec =
                        string_atom_names[iatc];
                    if (pdb_atom_namec == dcr.atom_id_c_4c()) {

                      char alt_conf_c = res_selection[iatc]->altloc;
                      char alt_conf_1 = res_selection[iat1]->altloc;
                      char alt_conf_2 = res_selection[iat2]->altloc;
                      char alt_conf_3 = res_selection[iat3]->altloc;

                      if (((alt_conf_1 == alt_conf_c) || (alt_conf_1 == '\0')) &&
                          ((alt_conf_2 == alt_conf_c) ||
                           (alt_conf_2 == '\0')) &&
                          ((alt_conf_3 == alt_conf_c) ||
                           (alt_conf_3 == '\0'))) {

                        index1 = iat1;
                        index2 = iat2;
                        index3 = iat3;
                        indexc = iatc;

                        std::vector<bool> is_hydrogen_flags(i_no_res_atoms, false);
                        for (int i = 0; i < i_no_res_atoms; i++) {
                          if (res_selection[i]->element.is_hydrogen()) {
                            is_hydrogen_flags[i] = true;
                          }
                        }

                        int chiral_hydrogen_index =
                            get_chiral_hydrogen_index_gemmi(
                                index1, index2, index3,
                                is_hydrogen_flags);

                        if (fabs(dcr.target_volume()) < 1000.0 &&
                            fabs(dcr.target_volume()) > 0.00001) {

                          std::vector<bool> fixed_flags(4, false);

                          coot::simple_restraint sr(
                              coot::restraint_type_t::CHIRAL_VOLUME_RESTRAINT, indexc, index1,
                              index2, index3, dcr.volume_sign,
                              dcr.target_volume(), dcr.volume_sigma(),
                              fixed_flags, chiral_hydrogen_index);
                          restraints_vec.push_back(sr);
                          n_chiral_restr++;
                        } else {
                          std::cout
                              << "WARNING:: Reject chiral restraint for "
                              << res_selection[iatc]->name
                              << " with target volume "
                              << dcr.target_volume()
                              << " with volume sigma "
                              << dcr.volume_sigma()
                              << " with volume sign "
                              << dcr.volume_sign << std::endl;
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

  return n_chiral_restr;
}