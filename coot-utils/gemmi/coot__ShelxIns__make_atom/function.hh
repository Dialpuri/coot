#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot {

struct shelx_card_info_t {
    std::string card;
    std::vector<std::string> words;
};

class ShelxIns {
public:
    std::vector<std::string> sfac;

    // Original MMDB version ported to gemmi
    gemmi::Atom*
    make_atom_gemmi(const shelx_card_info_t& card, const std::string& altconf,
                   int udd_afix_handle_in,
                   int udd_non_riding_atom_flag_handle_in,
                   int udd_riding_atom_negative_u_value_handle_in,
                   bool have_udd_atoms, int current_afix,
                   const clipper::Cell& cell_in,
                   const std::vector<gemmi::Atom*>& atom_vector) const;

private:
    std::string make_atom_element(const std::string& name, int sfac_index) const;
    std::string make_atom_name(const std::string& name, const std::string& element) const;
    gemmi::Atom* previous_non_riding_atom(const std::vector<gemmi::Atom*>& atom_vector,
                                         int udd_non_riding_atom_flag_handle_in) const;
    static float string_to_float(const std::string& s);
};

} // namespace coot