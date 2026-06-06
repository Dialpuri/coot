#include <iostream>

// Gemmi port dependencies — include BEFORE function.hh
// NOTE: We skip the delete_atom port header because it defines
//   `namespace beam_in_linked_residue` which collides with our class
//   `coot::beam_in_linked_residue`.  We implement delete_atom inline below.
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__beam_in_linked_residue__lsq_fit/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__mol_from_dictionary__c7d1b3/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_first_residue/gemmi/function.hh"

#include "function.hh"

namespace coot {

// ── delete_atom inline ─────────────────────────────────────────────
static void delete_atom_inline(gemmi::Residue* res, const std::string& atom_name)
{
  if (!res) return;
  auto it = res->atoms.begin();
  for (; it != res->atoms.end(); ++it) {
    if (it->name == atom_name) {
      res->atoms.erase(it);
      return;
    }
  }
}

// ── Constructor ──────────────────────────────────────────────────────
beam_in_linked_residue::beam_in_linked_residue(
    gemmi::Residue* res, const std::string& link,
    const std::string& comp_id, coot::protein_geometry* geom)
    : residue_ref(res),
      template_res_ref(nullptr),
      template_res_mov(nullptr),
      comp_id_ref(res->name),
      comp_id_new(comp_id),
      link_type(link),
      have_template(false),
      geom_p(geom)
{
}

// ── Private helper: make_reference_atom_names ────────────────────────
static std::vector<std::string> make_reference_atom_names_impl(
    coot::protein_geometry& geom, const std::string& comp_id)
{
  std::vector<std::string> atom_names;

  int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
  gemmi::Structure* mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, 1);
  if (!mol) return atom_names;

  gemmi::CRA cra = coot::util::get_first_residue_gemmi(*mol);
  if (!cra.residue) return atom_names;

  for (auto& atom : cra.residue->atoms) {
    atom_names.push_back(atom.name);
  }

  return atom_names;
}

// ── setup_by_comp_id ────────────────────────────────────────────────
bool beam_in_linked_residue::setup_by_comp_id(
    const std::string& template_id, const std::string& reference_id)
{
  int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

  // Get reference residue from dictionary
  gemmi::Structure* ref_mol = coot::mol_from_dictionary_gemmi(*geom_p, reference_id, imol_enc, 1);
  if (!ref_mol) {
    std::cout << "ERROR:: couldn't get reference residue coords for " << reference_id << std::endl;
    return false;
  }
  gemmi::CRA cra_ref = coot::util::get_first_residue_gemmi(*ref_mol);
  if (!cra_ref.residue) {
    std::cout << "ERROR:: no residue for " << reference_id << std::endl;
    return false;
  }
  template_res_ref = cra_ref.residue;

  // Get moving residue from dictionary
  gemmi::Structure* mov_mol = coot::mol_from_dictionary_gemmi(*geom_p, template_id, imol_enc, 1);
  if (!mov_mol) {
    std::cout << "ERROR:: couldn't get moving residue coords for " << template_id << std::endl;
    return false;
  }
  gemmi::CRA cra_mov = coot::util::get_first_residue_gemmi(*mov_mol);
  if (!cra_mov.residue) {
    std::cout << "ERROR:: no residue for " << template_id << std::endl;
    return false;
  }
  template_res_mov = cra_mov.residue;

  return false; // Always return false — "no custom template loaded"
}

// ── get_residue_raw_gemmi ───────────────────────────────────────────
gemmi::Residue* beam_in_linked_residue::get_residue_raw_gemmi() const
{
  gemmi::Residue* r = nullptr;

  if (!have_template) {
    std::cout << "WARNING:: no template" << std::endl;
  } else {
    std::vector<std::string> lsq_atom_names_ref  = make_reference_atom_names_impl(*const_cast<coot::protein_geometry*>(geom_p), comp_id_ref);
    std::vector<std::string> lsq_atom_names_match = make_reference_atom_names_impl(*const_cast<coot::protein_geometry*>(geom_p), comp_id_ref);

    if (lsq_atom_names_ref.size() == 0 ||
        lsq_atom_names_ref.size() != lsq_atom_names_match.size()) {
      std::cout << "WARNING:: reference atoms for LSQing match problem" << std::endl;
    } else {
      // Fit template_res_ref -> residue_ref, move template_res_mov
      bool status = beam_in_linked_residue_lsq_fit_gemmi(
          template_res_ref, residue_ref, template_res_mov,
          lsq_atom_names_ref, lsq_atom_names_match);

      // Debug call (mirrors original)
      beam_in_linked_residue_lsq_fit_gemmi(
          template_res_ref, residue_ref, template_res_ref,
          lsq_atom_names_ref, lsq_atom_names_match);

      if (status) {
        r = template_res_mov;

        std::string r_res_name(r->name);
        if (r_res_name != comp_id_new) {
          int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
          gemmi::Structure* r_mol = coot::mol_from_dictionary_gemmi(
              *const_cast<coot::protein_geometry*>(geom_p),
              comp_id_new, imol_enc, 1);

          if (r_mol) {
            gemmi::CRA cra = coot::util::get_first_residue_gemmi(*r_mol);
            gemmi::Residue* r_new = cra.residue;

            if (!r_new) {
              std::cout << "WARNING:: couldn't get reference residue coords for "
                        << comp_id_new << " substituting " << r_res_name << std::endl;
            } else {
              bool state = beam_in_linked_residue_lsq_fit_gemmi(
                  r_new, r, r_new, lsq_atom_names_ref, lsq_atom_names_match);
              if (state)
                r = r_new;
              else
                std::cout << "WARNING:: couldn't match coords for "
                          << comp_id_new << " substituting " << r_res_name << std::endl;
            }
          }
        }
      }
    }
  }

  if (r) {
    try {
      std::pair<coot::chem_mod, coot::chem_mod> mods =
          geom_p->get_chem_mods_for_link(link_type);

      std::string res_name_ref = residue_ref->name;
      for (unsigned int i = 0; i < mods.first.atom_mods.size(); i++) {
        if (mods.first.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
          std::string atom_name = mods.first.atom_mods[i].atom_id;
          std::string at_name = geom_p->atom_id_expand(atom_name, res_name_ref, 0);
          delete_atom_inline(residue_ref, at_name);
        }
      }

      std::string res_name_new = r->name;
      for (unsigned int i = 0; i < mods.second.atom_mods.size(); i++) {
        if (mods.second.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
          std::string atom_name = mods.second.atom_mods[i].atom_id;
          std::string at_name = geom_p->atom_id_expand(atom_name, res_name_new, 0);
          delete_atom_inline(r, at_name);
        }
      }
    } catch (const std::runtime_error& rte) {
      // no chem mod for that link - that's fine
    }
  }

  return r;
}

} // namespace coot