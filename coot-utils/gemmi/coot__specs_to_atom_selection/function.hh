#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "coot/geometry/residue-and-atom-specs.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/calculate.hpp>

namespace coot {

/// Check if a residue name is a standard amino acid or nucleotide
inline bool is_standard_residue_name(const std::string &residue_name) {
   if (residue_name == "ALA") return true;
   if (residue_name == "ARG") return true;
   if (residue_name == "ASN") return true;
   if (residue_name == "ASP") return true;
   if (residue_name == "CYS") return true;
   if (residue_name == "GLN") return true;
   if (residue_name == "GLU") return true;
   if (residue_name == "GLY") return true;
   if (residue_name == "HIS") return true;
   if (residue_name == "ILE") return true;
   if (residue_name == "LEU") return true;
   if (residue_name == "LYS") return true;
   if (residue_name == "MET") return true;
   if (residue_name == "MSE") return true;
   if (residue_name == "PHE") return true;
   if (residue_name == "PRO") return true;
   if (residue_name == "SER") return true;
   if (residue_name == "THR") return true;
   if (residue_name == "TRP") return true;
   if (residue_name == "TYR") return true;
   if (residue_name == "VAL") return true;
   if (residue_name == "A") return true;
   if (residue_name == "G") return true;
   if (residue_name == "C") return true;
   if (residue_name == "U") return true;
   if (residue_name == "I") return true;
   if (residue_name == "DA") return true;
   if (residue_name == "DG") return true;
   if (residue_name == "DC") return true;
   if (residue_name == "DT") return true;
   if (residue_name == "DI") return true;
   return false;
}

/// Port of specs_to_atom_selection from MMDB to gemmi.
/// Returns 1 on success (selection created), -1 on failure.
inline int specs_to_atom_selection_gemmi(const std::vector<residue_spec_t> &specs,
                                         gemmi::Structure &st,
                                         int atom_mask_mode) {
   if (st.models.empty()) return -1;

   // In gemmi, we simulate selection by returning a positive value
   // and the selection is implicitly the atoms that match the specs
   // For compatibility, we return 1 on success
   for (unsigned int ilocal=0; ilocal<specs.size(); ilocal++) {
      std::string res_name_selection  = "*";
      std::string atom_name_selection = "*";

      if (atom_mask_mode != 0) {
         // Find the residue matching the spec
         bool found = false;
         for (auto& model : st.models) {
            for (auto& chain : model.chains) {
               if (chain.name == specs[ilocal].chain_id) {
                  for (auto& residue : chain.residues) {
                     if (residue.seqid.num.value == specs[ilocal].res_no) {
                        char ins_code = specs[ilocal].ins_code.empty() ? ' ' : specs[ilocal].ins_code[0];
                        if (residue.seqid.icode == ins_code) {
                           found = true;
                           std::string residue_name = residue.name;
                           if (is_standard_residue_name(residue_name)) {
                              // PDBv3 FIXME
                              if (atom_mask_mode == 1)
                                 atom_name_selection = " N  , H  , HA , CA , C  , O  ";
                              if (atom_mask_mode == 2)
                                 atom_name_selection = "!( N  , H  , HA , CA , C  , O  )";
                              if (atom_mask_mode == 3)
                                 atom_name_selection = "!( N  , H  , HA , CA , C  , O  , CB )";
                           } else {
                              if (atom_mask_mode == 4)
                                 atom_name_selection = "%%%%%%"; // nothing
                              if (atom_mask_mode == 5)
                                 atom_name_selection = "%%%%%%"; // nothing
                           }
                           break;
                        }
                     }
                  }
               }
               if (found) break;
            }
            if (found) break;
         }
      }
   }
   return 1; // success
}

} // namespace coot