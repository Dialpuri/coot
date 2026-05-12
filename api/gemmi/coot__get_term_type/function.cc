#include "function.hh"
#include <algorithm>

namespace coot {

std::string get_term_type_gemmi(const gemmi::CRA& cra) {
    std::string term_type = "not-terminal-residue"; // returned thing

    if (!cra.residue)
        return term_type;

    if (!cra.chain)
        return term_type;

    int ires_atom = cra.residue->seqid.num.value;
    int nres = static_cast<int>(cra.chain->residues.size());

    // including tests needed for single missing residue:
    short int has_up_neighb = 0;
    short int has_down_neighb = 0;
    short int has_up_up_neighb = 0;
    short int has_down_down_neighb = 0;

    // Check for neighbouring residues to the clicked atom. Don't count
    // waters as neighbours.
    for (int ires = 0; ires < nres; ++ires) {
        gemmi::Residue& res = cra.chain->residues[ires];
        
        // Skip water residues (same logic as original - don't count waters)
        if (res.name == "HOH" || res.name == "WAT")
            continue;
            
        int res_seq = res.seqid.num.value;
        if (res_seq == (ires_atom + 1))
            has_up_neighb = 1;
        if (res_seq == (ires_atom + 2))
            has_up_up_neighb = 1;
        if (res_seq == (ires_atom - 1))
            has_down_neighb = 1;
        if (res_seq == (ires_atom - 2))
            has_down_down_neighb = 1;
    }

    if ((has_up_neighb + has_down_neighb) == 1) {
        if (has_up_neighb)
            term_type = "N";
        if (has_down_neighb)
            term_type = "C";
    }

    if ((has_up_neighb == 0) && (has_down_neighb == 0))
        term_type = "singleton";

    // Now test for missing single residue, "M" (mid):
    if ((!has_up_neighb) && has_up_up_neighb)
        term_type = "MC"; // missing middle res, treat as C term

    if ((!has_down_neighb) && has_down_down_neighb)
        term_type = "MN"; // missing middle res, treat as N term

    return term_type;
}

} // namespace coot