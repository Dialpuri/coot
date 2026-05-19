#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot {

std::pair<bool, clipper::Coord_orth>
centre_of_molecule_using_masses_gemmi(const gemmi::Structure& st) {
   std::map<std::string, double> pdb_element_weights = {
      {"H", 1.008},   {"HE", 4.0026}, {"LI", 6.94},   {"BE", 9.0122}, {"B", 10.81},   {"C", 12.011},
      {"N", 14.007},  {"O", 15.999},  {"F", 18.998},  {"NE", 20.180}, {"NA", 22.990}, {"MG", 24.305},
      {"AL", 26.982},  {"SI", 28.085}, {"P", 30.974},  {"S", 32.06},   {"CL", 35.45},  {"AR", 39.948},
      {"K", 39.098},   {"CA", 40.078}, {"SC", 44.956}, {"TI", 47.867}, {"V", 50.942},  {"CR", 51.996},
      {"MN", 54.938},  {"FE", 55.845}, {"CO", 58.933}, {"NI", 58.693}, {"CU", 63.546}, {"ZN", 65.38},
      {"GA", 69.723},  {"GE", 72.630}, {"AS", 74.922}, {"SE", 78.971}, {"BR", 79.904}, {"KR", 83.798},
      {"RB", 85.468},  {"SR", 87.62},  {"Y", 88.906},  {"ZR", 91.224}, {"NB", 92.906}, {"MO", 95.95},
      {"TC", 98.0},    {"RU", 101.07}, {"RH", 102.91}, {"PD", 106.42}, {"AG", 107.87}, {"CD", 112.41},
      {"IN", 114.82},  {"SN", 118.71}, {"SB", 121.76}, {"TE", 127.60}, {"I", 126.90},  {"XE", 131.29},
      {"CS", 132.91},  {"BA", 137.33}, {"LA", 138.91}, {"CE", 140.12}, {"PR", 140.91}, {"ND", 144.24},
      {"PM", 145.0},   {"SM", 150.36}, {"EU", 151.96}, {"GD", 157.25}, {"TB", 158.93}, {"DY", 162.50},
      {"HO", 164.93},  {"ER", 167.26}, {"TM", 168.93}, {"YB", 173.05}, {"LU", 174.97}, {"HF", 178.49},
      {"TA", 180.95},  {"W", 183.84},  {"RE", 186.21}, {"OS", 190.23}, {"IR", 192.22}, {"PT", 195.08},
      {"AU", 196.97},  {"HG", 200.59}, {"TL", 204.38}, {"PB", 207.2},  {"BI", 208.98}, {"PO", 209.0},
      {"AT", 210.0},   {"RN", 222.0},  {"FR", 223.0},  {"RA", 226.0},  {"AC", 227.0},  {"TH", 232.04},
      {"PA", 231.04},  {"U", 238.03},  {"NP", 237.0},  {"PU", 244.0},  {"AM", 243.0},  {"CM", 247.0},
      {"BK", 247.0},   {"CF", 251.0},  {"ES", 252.0},  {"FM", 257.0},  {"MD", 258.0},  {"NO", 259.0},
      {"LR", 262.0},   {"RF", 267.0},  {"DB", 270.0},  {"SG", 271.0},  {"BH", 270.0},  {"HS", 277.0},
      {"MT", 278.0},   {"DS", 281.0},  {"RG", 282.0},  {"CN", 285.0},  {"FL", 289.0},  {"LV", 293.0},
      {"TS", 294.0},   {"OG", 294.0}
   };

   int n_atoms = 0;
   double xs=0, ys=0, zs=0;
   double sum_weight = 0.0;
   for (const gemmi::Model& model : st.models) {
      for (const gemmi::Chain& chain : model.chains) {
         for (const gemmi::Residue& residue : chain.residues) {
            for (const gemmi::Atom& atom : residue.atoms) {
               std::string ele = atom.element.name();
               double w = 6.0;
               std::map<std::string, double>::const_iterator it;
               it = pdb_element_weights.find(ele);
               if (it != pdb_element_weights.end()) w = it->second;
               xs += w * atom.pos.x;
               ys += w * atom.pos.y;
               zs += w * atom.pos.z;
               sum_weight += w;
               n_atoms++;
            }
         }
      }
   }
   if (n_atoms > 0) {
      double x = xs / static_cast<double>(sum_weight);
      double y = ys / static_cast<double>(sum_weight);
      double z = zs / static_cast<double>(sum_weight);
      return std::make_pair(true, clipper::Coord_orth(x,y,z));
   } else {
      return std::make_pair(false, clipper::Coord_orth(0,0,0));
   }
}

} // namespace coot