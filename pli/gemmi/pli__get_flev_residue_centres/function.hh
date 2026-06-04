#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <iostream>

#include "coot-utils/coot-coord-utils.hh"
#include "clipper/core/coords.h"
#include "pli/flev.hh"
#include "geometry/residue-and-atom-specs.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_lsq_matrix/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_centre/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__transform_atoms/gemmi/function.hh"

namespace pli {

inline std::vector<pli::fle_residues_helper_t>
get_flev_residue_centres_gemmi(gemmi::CRA residue_ligand_3d,
                               gemmi::Model *mol_containing_residue_ligand,
                               std::vector<gemmi::CRA> residues,
                               gemmi::Structure *flat_mol) {

   std::vector<fle_residues_helper_t> centres;

   if (flat_mol) {
      // get the lsq matrix that maps the ligand in 3D onto the flat ligand
      int res_no = residue_ligand_3d.residue->seqid.num.value;
      std::string chain_id = residue_ligand_3d.chain->name;
      int every_nth = 1;
      std::vector<coot::lsq_range_match_info_t> matches;
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id, coot::lsq_t::ALL);
      matches.push_back(match);
      std::pair<short int, clipper::RTop_orth> lsq_mat =
         coot::util::get_lsq_matrix_gemmi(flat_mol->models.empty() ? nullptr : &flat_mol->models[0], mol_containing_residue_ligand, matches, every_nth, false);

      // Now make the residues
      centres.resize(residues.size());
      for (unsigned int ires=0; ires<residues.size(); ires++) {
         gemmi::Residue *res_copy = coot::util::deep_copy_this_residue_gemmi(residues[ires].residue);
         std::string res_name = residues[ires].residue->name;
         std::pair<bool, clipper::Coord_orth> absolute_centre =
            coot::util::get_residue_centre_gemmi(*res_copy);
         if (absolute_centre.first) {
            coot::util::transform_atoms_gemmi(res_copy, lsq_mat.second);
            std::pair<bool, clipper::Coord_orth> c =
               coot::util::get_residue_centre_gemmi(*res_copy);
            if (c.first) {
               // Build residue_spec_t from the CRA
               coot::residue_spec_t spec(residues[ires].chain->name,
                                         residues[ires].residue->seqid.num.value,
                                         std::string(1, residues[ires].residue->seqid.icode));
               fle_residues_helper_t fle_centre(c.second, spec, res_name);

               fle_centre.set_interaction_position(absolute_centre.second);
               centres[ires] = fle_centre;
            } else {
               coot::residue_spec_t spec(res_copy->seqid.num.value);
               std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
                         << spec << std::endl;
            }
         } else {
            coot::residue_spec_t spec(res_copy->seqid.num.value);
            std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
                      << spec << std::endl;
         }
         delete res_copy;
      }
   }
   return centres;
}

} // namespace pli