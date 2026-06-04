#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <clipper/clipper.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__average_position/gemmi/function.hh"

namespace coot {
namespace secondary_structure_header_records {

class strand_relation_t {
public:
   unsigned int strand_idx;
   enum sense_t { FIRST, PARALLEL, ANTI_PARALLEL, NO_RESULT };
   sense_t sense;

   strand_relation_t(unsigned int idx, sense_t s) : strand_idx(idx), sense(s) {}

   bool operator==(const strand_relation_t &sr_in) const {
      return (sr_in.strand_idx == strand_idx);
   }

   bool operator<(const strand_relation_t &sr_in) const {
      return (sr_in.strand_idx < strand_idx);
   }

   static sense_t get_strand_sense_gemmi(const std::vector<gemmi::Residue *> &strand_1,
                                         const std::vector<gemmi::Residue *> &strand_2);
};

} // namespace secondary_structure_header_records
} // namespace coot

// ── inline definition ────────────────────────────────────────────────

namespace coot {
namespace secondary_structure_header_records {

// Helper: collect atom positions from a gemmi::Residue into
// std::vector<clipper::Coord_orth> for average_position_gemmi.
static std::vector<clipper::Coord_orth>
residue_atom_positions(const gemmi::Residue &res) {
   std::vector<clipper::Coord_orth> pts;
   pts.reserve(res.atoms.size());
   for (const gemmi::Atom &a : res.atoms) {
      pts.push_back(clipper::Coord_orth(a.pos.x, a.pos.y, a.pos.z));
   }
   return pts;
}

inline strand_relation_t::sense_t
strand_relation_t::get_strand_sense_gemmi(const std::vector<gemmi::Residue *> &strand_1,
                                          const std::vector<gemmi::Residue *> &strand_2) {
   strand_relation_t::sense_t s(NO_RESULT);

   if (strand_1.size() > 1) {
      if (strand_2.size() > 1) {
         gemmi::Residue *sr10 = strand_1[0];
         gemmi::Residue *sr1e = strand_1.back();
         gemmi::Residue *sr20 = strand_2[0];
         gemmi::Residue *sr2e = strand_2.back();

         clipper::Coord_orth sr10_pt = coot::util::average_position_gemmi(residue_atom_positions(*sr10));
         clipper::Coord_orth sr1e_pt = coot::util::average_position_gemmi(residue_atom_positions(*sr1e));
         clipper::Coord_orth sr20_pt = coot::util::average_position_gemmi(residue_atom_positions(*sr20));
         clipper::Coord_orth sr2e_pt = coot::util::average_position_gemmi(residue_atom_positions(*sr2e));

         clipper::Coord_orth v1 = sr1e_pt - sr10_pt;
         clipper::Coord_orth v2 = sr2e_pt - sr20_pt;
         clipper::Coord_orth v1u(v1.unit());
         clipper::Coord_orth v2u(v2.unit());

         double cos_theta = clipper::Coord_orth::dot(v1u, v2u);
         if (cos_theta > 0)
            s = PARALLEL;
         else
            s = ANTI_PARALLEL;
      }
   }
   return s;
}

} // namespace secondary_structure_header_records
} // namespace coot