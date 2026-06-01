#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {

struct atom_overlaps_container_t {
    bool kludge_filter_gemmi(gemmi::CRA cra_1, gemmi::CRA cra_2) const {
        bool reject = false;
        if (cra_1.chain && cra_2.chain && cra_1.residue && cra_2.residue && cra_1.atom && cra_2.atom) {
            if (cra_1.chain->name == cra_2.chain->name) {
                std::string res_name_1(cra_1.residue->name);
                if (res_name_1 == "ASN") {
                    std::string res_name_2(cra_2.residue->name);
                    if (res_name_2 == "NAG") {
                        std::string atom_name_1(cra_1.atom->name);
                        if (atom_name_1 == "NE2") {
                            std::string atom_name_2(cra_2.atom->name);
                            if (atom_name_2 == "C2") {
                                reject = true;
                            }
                        }
                    }
                }
                if (res_name_1 == "NAG") {
                    std::string res_name_2(cra_2.residue->name);
                    if (res_name_2 == "ASN") {
                        std::string atom_name_1(cra_1.atom->name);
                        if (atom_name_1 == "C2") {
                            std::string atom_name_2(cra_2.atom->name);
                            if (atom_name_2 == "NE2") {
                                reject = true;
                            }
                        }
                    }
                }
            }
        }
        return reject;
    }
};

} // namespace coot