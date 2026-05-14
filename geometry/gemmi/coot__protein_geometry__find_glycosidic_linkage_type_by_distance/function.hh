#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace protein_geometry {

struct glycosidic_distance {
    gemmi::Position at1_pos;
    gemmi::Position at2_pos;
    std::string at1_name;
    std::string at2_name;
    std::string res1_name;
    std::string res2_name;
    double distance;
    
    glycosidic_distance(const gemmi::Atom* at1, const gemmi::Atom* at2, double dist)
        : at1_pos(at1->pos),
          at2_pos(at2->pos),
          at1_name(at1->name),
          at2_name(at2->name),
          res1_name(""),
          res2_name(""),
          distance(dist) {}
    
    bool operator<(const glycosidic_distance& other) const {
        return distance < other.distance;
    }
};

std::string find_glycosidic_linkage_type_by_distance_gemmi(gemmi::CRA cra1, gemmi::CRA cra2) {
    const gemmi::Residue& res1 = *cra1.residue;
    const gemmi::Residue& res2 = *cra2.residue;
    
    double critical_dist = 2.4;
    std::vector<glycosidic_distance> close;
    
    for (const gemmi::Atom& at1 : res1.atoms) {
        gemmi::Position a1 = at1.pos;
        for (const gemmi::Atom& at2 : res2.atoms) {
            gemmi::Position a2 = at2.pos;
            gemmi::Vec3 diff = a1 - a2;
            double d = diff.length_sq();
            if (d < critical_dist * critical_dist) {
                close.push_back(glycosidic_distance(&at1, &at2, std::sqrt(d)));
            }
        }
    }
    
    std::sort(close.begin(), close.end());
    
    std::string link_type("");
    
    try {
        float smallest_link_dist = 99999.9;
        for (unsigned int i = 0; i < close.size(); i++) {
            std::string name_1 = close[i].at1_name;
            std::string name_2 = close[i].at2_name;
            
            if (name_1 == " C1 ") {
                if (name_2 == " ND2") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "pyr-ASN";
                    }
                }
            }
            
            if (name_1 == " O4 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "BETA1-4";
                    }
                }
            }
            
            if (name_1 == " O2 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "BETA1-2";
                    }
                }
            }
            
            if (name_1 == " O3 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "BETA1-3";
                    }
                }
            }
            
            if (name_1 == " C2 ") {
                if (name_2 == " O3 ") {
                    if (res1.name == "SIA") {
                        if (close[i].distance < smallest_link_dist) {
                            smallest_link_dist = close[i].distance;
                            link_type = "ALPHA2-3";
                        }
                    }
                }
            }
            
            if (name_1 == " C2 ") {
                if (name_2 == " O6 ") {
                    if (res1.name == "SIA") {
                        if (close[i].distance < smallest_link_dist) {
                            smallest_link_dist = close[i].distance;
                            link_type = "ALPHA2-6";
                        }
                    }
                }
            }
            
            if (name_1 == " O6 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "BETA1-6";
                    }
                }
            }
            
            if (name_1 == " O2 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "ALPHA1-2";
                    }
                }
            }
            
            if (name_1 == " O3 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "ALPHA1-3";
                    }
                }
            }
            
            if (name_1 == " C2 ") {
                if (name_2 == " O3 ") {
                    if (res1.name == "SIA") {
                        if (close[i].distance < smallest_link_dist) {
                            smallest_link_dist = close[i].distance;
                            link_type = "ALPHA2-3";
                        }
                    }
                }
            }
            
            if (name_1 == " O4 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "ALPHA1-4";
                    }
                }
            }
            
            if (name_1 == " O6 ") {
                if (name_2 == " C1 ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        link_type = "ALPHA1-6";
                    }
                }
            }
            
            if (name_1 == " C1 ") {
                if (name_2 == " OG ") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        if (res2.name == "SER")
                            link_type = "pyr-SER";
                    }
                }
            }
            
            if (name_1 == " C1 ") {
                if (name_2 == " OG1") {
                    if (close[i].distance < smallest_link_dist) {
                        smallest_link_dist = close[i].distance;
                        if (res2.name == "THR")
                            link_type = "pyr-THR";
                    }
                }
            }
        }
    }
    catch (const std::runtime_error &rte) {
        std::cout << "WARNING::" << rte.what() << std::endl;
    }
    
    return link_type;
}

} // namespace protein_geometry
} // namespace coot