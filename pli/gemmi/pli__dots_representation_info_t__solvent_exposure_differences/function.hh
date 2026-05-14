#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace pli {

class solvent_exposure_difference_helper_t {
public:
    coot::residue_spec_t res_spec;
    double exposure_fraction_holo;
    double exposure_fraction_apo;
    solvent_exposure_difference_helper_t(const coot::residue_spec_t& res_spec_in, double h, double a)
        : res_spec(res_spec_in), exposure_fraction_holo(h), exposure_fraction_apo(a) {}
};

namespace dots_representation_info_detail {

// Convert gemmi CRA to coot::residue_spec_t
inline coot::residue_spec_t cra_to_residue_spec(const gemmi::CRA& cra) {
    const gemmi::Chain* chain = cra.chain;
    const gemmi::Residue* res = cra.residue;
    if (!chain || !res)
        return coot::residue_spec_t();
    
    std::string chain_id = chain->name;
    int res_no = res->seqid.num.value;
    char icode = res->seqid.icode;
    // Normalize insertion code: gemmi uses ' ' for missing, coot uses ""
    std::string ins_code = (icode == ' ') ? "" : std::string(1, icode);
    return coot::residue_spec_t(chain_id, res_no, ins_code);
}

// Get residue name from gemmi residue
inline std::string get_res_name(const gemmi::Residue& res) {
    return res.name;
}

} // namespace dots_representation_info_detail

class dots_representation_info_t {
public:
    std::vector<solvent_exposure_difference_helper_t>
    solvent_exposure_differences_gemmi(const gemmi::CRA& cra_ref,
                                       const std::vector<gemmi::CRA>& near_cras) const;
};

} // namespace pli