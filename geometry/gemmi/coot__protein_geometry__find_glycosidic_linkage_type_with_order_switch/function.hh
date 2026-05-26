#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__find_glycosidic_linkage_type_by_distance/gemmi/function.hh"

namespace coot { namespace protein_geometry {

inline std::pair<std::string, bool>
find_glycosidic_linkage_type_with_order_switch_gemmi(gemmi::CRA cra1, gemmi::CRA cra2) {
    std::pair<std::string, bool> r("", false);

    std::string l = find_glycosidic_linkage_type_by_distance_gemmi(cra1, cra2);

    if (l == "") {
        l = find_glycosidic_linkage_type_by_distance_gemmi(cra2, cra1);
        if (l != "") {
            r.first = l;
            r.second = true;
        }
    } else {
        r.first = l;
        r.second = false;
    }
    return r;
}

}} // namespace coot::protein_geometry