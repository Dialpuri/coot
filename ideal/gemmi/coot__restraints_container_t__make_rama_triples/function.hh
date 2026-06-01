#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"

namespace coot {

/// gemmi-port of rama_triple_t — carries chain pointer since
/// gemmi::Residue has no parent-chain pointer (unlike mmdb::Residue).
struct rama_triple_t_gemmi {
    gemmi::Chain*  chain_1;
    gemmi::Residue* r_1;
    gemmi::Residue* r_2;
    gemmi::Residue* r_3;
    std::string     link_type;
    bool            fixed_1;
    bool            fixed_2;
    bool            fixed_3;

    rama_triple_t_gemmi(gemmi::Chain* c1, gemmi::Residue* r1, gemmi::Residue* r2,
                        gemmi::Residue* r3, const std::string& lt,
                        bool f1 = false, bool f2 = false, bool f3 = false)
        : chain_1(c1), r_1(r1), r_2(r2), r_3(r3),
          link_type(lt), fixed_1(f1), fixed_2(f2), fixed_3(f3) {}
};

/// gemmi-port of coot::restraints_container_t::make_rama_triples
///
/// Original took an MMDB selection handle; the gemmi version takes a
/// vector of CRA (chain-residue-atom) tuples that carry all the parent
/// context gemmi Residue nodes lack.
inline std::vector<rama_triple_t_gemmi> make_rama_triples_gemmi(
    const std::vector<gemmi::CRA>& residues,
    const protein_geometry& /*geom*/)
{
    std::vector<rama_triple_t_gemmi> v;
    for (size_t i = 0; i + 2 < residues.size(); i++) {
        if (residues[i].residue && residues[i + 1].residue && residues[i + 2].residue) {
            v.emplace_back(
                residues[i].chain,
                residues[i].residue,
                residues[i + 1].residue,
                residues[i + 2].residue,
                "TRANS"
            );
        }
    }
    return v;
}

} // namespace coot