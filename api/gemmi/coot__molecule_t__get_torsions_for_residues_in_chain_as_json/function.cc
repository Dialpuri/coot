#include "function.hh"

#include <clipper/core/clipper_util.h>
#include <clipper/core/coords.h>
#include <coot/coot-utils/json.hpp>

namespace {

bool find_atom_coord(const gemmi::Residue& res, const std::string& atom_name, clipper::Coord_orth& pos) {
    for (const auto& atom : res.atoms) {
        if (atom.name == atom_name) {
            pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            return true;
        }
    }
    return false;
}

} // anonymous namespace

std::string
coot::molecule_t::get_torsions_for_residues_in_chain_as_json_gemmi(
    gemmi::Structure& structure, const std::string& chain_id) {

   auto residue_has_alt_confs = [](const gemmi::Residue& res) -> bool {
      for (const auto& atom : res.atoms) {
         if (atom.altloc != '\0')
            return true;
      }
      return false;
   };

   nlohmann::json j_residues = nlohmann::json::array();

   for (auto& model : structure.models) {
      for (auto& chain : model.chains) {
         if (chain.name != chain_id) continue;
         if (chain.residues.size() < 3) continue;

         for (size_t ires = 1; ires < (chain.residues.size() - 1); ires++) {
            gemmi::Residue& res_prev = chain.residues[ires - 1];
            gemmi::Residue& res_this = chain.residues[ires];
            gemmi::Residue& res_next = chain.residues[ires + 1];

            if (residue_has_alt_confs(res_this)) continue;

            // gemmi atom.name is unpadded ("N", "CA", "C") — NOT padded like MMDB
            clipper::Coord_orth c_prev, n_this, ca_this, c_this, n_next;
            bool ok_c_prev  = find_atom_coord(res_prev, "C", c_prev);
            bool ok_n_this  = find_atom_coord(res_this, "N", n_this);
            bool ok_ca_this = find_atom_coord(res_this, "CA", ca_this);
            bool ok_c_this  = find_atom_coord(res_this, "C", c_this);
            bool ok_n_next  = find_atom_coord(res_next, "N", n_next);

            if (!ok_c_prev || !ok_n_this || !ok_ca_this || !ok_c_this || !ok_n_next) continue;

            // Use the same torsion convention as original: ca_this.torsion(a, b, ca_this, d)
            double phi = clipper::Util::rad2d(ca_this.torsion(c_prev, n_this, ca_this, c_this));
            double psi = clipper::Util::rad2d(ca_this.torsion(n_this, ca_this, c_this, n_next));

            // Original validates peptide bonds < 2.0A
            double dist_1 = clipper::Coord_orth::length(c_prev, n_this);
            double dist_2 = clipper::Coord_orth::length(c_this, n_next);
            if (dist_1 >= 2.0 || dist_2 >= 2.0) continue;

            rama_plot::phi_psi_t pp(phi, psi);

            // Get backbone N, CA, C for tau calculation
            clipper::Coord_orth n_pos  = n_this;
            clipper::Coord_orth ca_pos = ca_this;
            clipper::Coord_orth c_pos  = c_this;
            double tau = clipper::Util::rad2d(clipper::Coord_orth::angle(n_pos, ca_pos, c_pos));

            nlohmann::json j_chis = nlohmann::json::array();
            try {
               std::vector<coot::alt_confed_chi_angles> chi_info = coot::primitive_chi_angles::get_chi_angles_gemmi(res_this);
               if (chi_info.size() == 1) {
                  for (const auto& chi_pair : chi_info[0].chi_angles) {
                     nlohmann::json j_chi;
                     j_chi["chi"] = chi_pair.first;
                     j_chi["value"] = chi_pair.second;
                     j_chis.push_back(j_chi);
                  }
               }
            } catch (const std::runtime_error&) {
               // GLY, ALA etc. have no chi angles
            }

            nlohmann::json j_res;
            j_res["residue_type"] = std::string(res_this.name);
            j_res["res_no"] = res_this.seqid.num.value;
            j_res["chain_id"] = chain_id;
            // Normalize: MMDB uses "" for "no insertion code", gemmi uses ' '
            std::string ins_code = (res_this.seqid.icode == ' ') ? std::string("") : std::string(1, res_this.seqid.icode);
            j_res["ins_code"] = ins_code;
            j_res["phi"] = pp.phi;
            j_res["psi"] = pp.psi;
            j_res["tau"] = tau;
            j_res["chi_angles"] = j_chis;

            j_residues.push_back(j_res);
         }
      }
   }

   return j_residues.dump();
}