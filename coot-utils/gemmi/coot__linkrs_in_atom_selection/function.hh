#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/seqid.hpp>
#include <vector>
#include <algorithm>
#include <string>
#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"

// gemmi-compatible bonded pair (avoids mmdb::Residue* from the original bonded_pair_t)
struct bonded_pair_gemmi {
   gemmi::Residue *res_1;
   gemmi::Residue *res_2;
   std::string link_type;
   bool is_fixed_first;
   bool is_fixed_second;

   bonded_pair_gemmi(gemmi::Residue *r1, gemmi::Residue *r2, bool isff, bool isfs, const std::string &lt)
      : res_1(r1), res_2(r2), link_type(lt), is_fixed_first(isff), is_fixed_second(isfs) {}

   bonded_pair_gemmi()
      : res_1(nullptr), res_2(nullptr), is_fixed_first(false), is_fixed_second(false) {}
};

// gemmi-compatible bonded pair container
struct bonded_pair_container_gemmi {
   std::vector<bonded_pair_gemmi> bonded_residues;
   unsigned int size() const { return bonded_residues.size(); }
};

namespace coot {

bonded_pair_container_gemmi
linkrs_in_atom_selection_gemmi(gemmi::Structure& structure,
                               std::vector<gemmi::CRA> atom_selection,
                               coot::protein_geometry* geom_p) {
   bonded_pair_container_gemmi bpc;

   // Collect unique residues from the selection
   std::vector<gemmi::CRA> unique_residues;
   for (auto& cra : atom_selection) {
      if (cra.residue) {
         bool seen = false;
         for (auto& ur : unique_residues) {
            if (ur.residue == cra.residue) { seen = true; break; }
         }
         if (!seen) unique_residues.push_back(cra);
      }
   }

   // Normalize insertion codes: MMDB uses "" for no icode, gemmi uses ' '
   auto norm_icode = [](char ic) -> std::string {
      return (ic == ' ' || ic == '\0') ? "" : std::string(1, ic);
   };

   // Iterate through connections on the Structure (gemmi equivalent of LINKRs)
   bool found = false;
   for (auto& conn : structure.connections) {
      if (conn.type != gemmi::Connection::Type::Covale) continue;

      auto& addr1 = conn.partner1;
      auto& addr2 = conn.partner2;

      // Build link specs from connection addresses
      // AtomAddress has: chain_name, res_id (which has seqid: num.value and icode)
      coot::residue_spec_t link_spec_1(addr1.chain_name, addr1.res_id.seqid.num.value, norm_icode(addr1.res_id.seqid.icode));
      coot::residue_spec_t link_spec_2(addr2.chain_name, addr2.res_id.seqid.num.value, norm_icode(addr2.res_id.seqid.icode));

      for (unsigned int i = 0; i < unique_residues.size(); i++) {
         auto& cra1 = unique_residues[i];
         std::string chain1 = (cra1.chain) ? cra1.chain->name : "";
         coot::residue_spec_t spec_1(chain1, cra1.residue->seqid.num.value, norm_icode(cra1.residue->seqid.icode));
         if (spec_1 == link_spec_1) {
            for (unsigned int j = 0; j < unique_residues.size(); j++) {
               if (i != j) {
                  auto& cra2 = unique_residues[j];
                  std::string chain2 = (cra2.chain) ? cra2.chain->name : "";
                  coot::residue_spec_t spec_2(chain2, cra2.residue->seqid.num.value, norm_icode(cra2.residue->seqid.icode));
                  if (spec_2 == link_spec_2) {
                     found = true;
                     bonded_pair_gemmi pair(cra1.residue, cra2.residue, false, false, conn.link_id);
                     bpc.bonded_residues.push_back(pair);
                     break;
                  }
               }
            }
         }
         if (found) break;
      }
      if (found) break;
   }

   return bpc;
}

} // namespace coot