#pragma once
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {

std::pair<bool, clipper::RTop_orth>
base_pair_match_matix_gemmi(gemmi::Residue *res_ref, gemmi::Residue *res_mov) {
   
   std::string res_name = res_ref->name;
   bool is_pyrimidine = false;
   bool is_purine = false;

   std::vector<std::string> purine_atom_names;
   std::vector<std::string> pyrimidine_atom_names;
   
   pyrimidine_atom_names.push_back(" N9 ");
   pyrimidine_atom_names.push_back(" C8 ");
   pyrimidine_atom_names.push_back(" N7 ");
   pyrimidine_atom_names.push_back(" C5 ");
   pyrimidine_atom_names.push_back(" C4 ");
   // 
   pyrimidine_atom_names.push_back(" N1 ");
   pyrimidine_atom_names.push_back(" C2 ");
   pyrimidine_atom_names.push_back(" N3 ");
   pyrimidine_atom_names.push_back(" C6 ");
   pyrimidine_atom_names.push_back(" N6 ");

   purine_atom_names.push_back(" N1 ");
   purine_atom_names.push_back(" C2 ");
   purine_atom_names.push_back(" N3 ");
   purine_atom_names.push_back(" C4 ");
   purine_atom_names.push_back(" C5 ");
   
   if (res_name == "G" || res_name == "Gd" || res_name == "Gr")
      is_pyrimidine = true;
   if (res_name == "A" || res_name == "Ad" || res_name == "Ar")
      is_pyrimidine = true;
   if (res_name == "DA" || res_name == "DG")
      is_pyrimidine = true;

   if (res_name == "T" || res_name == "Td" || res_name == "Tr")
      is_purine = true;
   if (res_name == "C" || res_name == "Cd" || res_name == "Cr")
      is_purine = true;
   if (res_name == "U" || res_name == "Ud" || res_name == "Ur")
      is_purine = true;
   if (res_name == "DT" || res_name == "DC" || res_name == "DT")
      is_purine = true;

   std::vector<std::string> base_atom_names;

   if (is_pyrimidine)
      base_atom_names = pyrimidine_atom_names;
   if (is_purine)
      base_atom_names = purine_atom_names;

   clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
   clipper::Coord_orth pt_dum(0,0,0);
   clipper::RTop_orth rtop(m_dum, pt_dum);
   bool rtop_is_good = false;
   if (base_atom_names.size() == 0) {
      std::cout << "  Oops neither pyrimidine nor purine" << std::endl;
   } else {

      std::vector<clipper::Coord_orth> ref_pts;
      std::vector<clipper::Coord_orth> mov_pts;
      
      for (unsigned int i=0; i<base_atom_names.size(); i++) {
         gemmi::Atom *at_ref = nullptr;
         gemmi::Atom *at_mov = nullptr;
         
         for (auto& atom : res_ref->atoms) {
            std::string atom_name = atom.name;
            if (atom_name == base_atom_names[i]) { 
               at_ref = &atom;
               break;
            } 
         } 
         
         for (auto& atom : res_mov->atoms) {
            std::string atom_name = atom.name;
            if (atom_name == base_atom_names[i]) { 
               at_mov = &atom;
               break;
            } 
         }
         
         if (at_ref && at_mov) {
            ref_pts.push_back(clipper::Coord_orth(at_ref->pos.x, at_ref->pos.y, at_ref->pos.z));
            mov_pts.push_back(clipper::Coord_orth(at_mov->pos.x, at_mov->pos.y, at_mov->pos.z));
         }
      }

      if (ref_pts.size() > 3) {
         rtop = clipper::RTop_orth(mov_pts, ref_pts);
         rtop_is_good = true;

      } else {
         std::cout << "   WARNING:: base_pair_match_matrix found only "
                   << ref_pts.size() << " matching points in the base " << std::endl;
      } 
   } 
   return std::pair<bool, clipper::RTop_orth> (rtop_is_good, rtop);
}

}