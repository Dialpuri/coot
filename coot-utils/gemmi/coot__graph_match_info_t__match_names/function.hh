#pragma once
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

struct graph_match_info_t {
  bool success = false;
  std::vector<std::pair<std::pair<std::string, std::string>,
                        std::pair<std::string, std::string>>> matching_atom_names;

  void match_names_gemmi(gemmi::Residue& residue);
  std::string invent_new_name(const std::string& atom_name,
                              const std::string& element,
                              const std::vector<std::string>& existing_names);
};

// Implementation of invent_new_name
inline std::string
graph_match_info_t::invent_new_name(const std::string& atom_name,
                                    const std::string& element,
                                    const std::vector<std::string>& existing_names) {
    std::string base_name = element;
    if (base_name.empty())
        base_name = atom_name;

    for (int i = 1; i < 100; ++i) {
        std::string new_name = base_name + std::to_string(i);
        if (std::find(existing_names.begin(), existing_names.end(), new_name) ==
            existing_names.end()) {
            return new_name;
        }
    }
    return base_name + "_new";
}

// Implementation of match_names_gemmi
inline void
graph_match_info_t::match_names_gemmi(gemmi::Residue& residue) {
    bool debug = false;
    if (!success) {
        std::cout << "Can't do name remapping, graph match failed" << std::endl;
        return;
    }

    std::vector<std::string> residue_atom_names;
    for (const auto& atom : residue.atoms) {
        if (std::find(residue_atom_names.begin(), residue_atom_names.end(), atom.name) ==
            residue_atom_names.end())
            residue_atom_names.push_back(atom.name);
    }

    std::vector<std::string> orig_moving_atom_names_non_mapped_non_same;
    std::vector<std::string> orig_moving_atom_names_non_mapped_same;

    for (const auto& atom : residue.atoms) {
        std::string atom_name = atom.name;
        bool found_match = false;
        for (const auto& pair : matching_atom_names) {
            if (pair.first.first == atom_name) {
                found_match = true;
                break;
            }
        }

        if (!found_match) {
            bool found_match_2 = false;
            for (const auto& pair : matching_atom_names) {
                if (pair.second.first == atom_name) {
                    found_match_2 = true;
                    break;
                }
            }

            if (found_match_2)
                orig_moving_atom_names_non_mapped_non_same.push_back(atom_name);
            else
                orig_moving_atom_names_non_mapped_same.push_back(atom_name);
        }
    }

    if (debug) {
        std::cout << "Mapped atom names: " << matching_atom_names.size() << std::endl;
        for (size_t i = 0; i < matching_atom_names.size(); i++) {
            std::cout << "   " << i << " :" << matching_atom_names[i].first.first
                      << ": -> :" << matching_atom_names[i].second.first << ":"
                      << std::endl;
        }
        std::cout << "Non-mapped non-same atom names: " << orig_moving_atom_names_non_mapped_non_same.size()
                  << std::endl;
        for (size_t i = 0; i < orig_moving_atom_names_non_mapped_non_same.size(); i++) {
            std::cout << "   " << i << " :" << orig_moving_atom_names_non_mapped_non_same[i] << ":" << std::endl;
        }
        std::cout << "Non mapped same atom names: " << orig_moving_atom_names_non_mapped_same.size()
                  << std::endl;
        for (size_t i = 0; i < orig_moving_atom_names_non_mapped_same.size(); i++) {
            std::cout << "   " << i << " :" << orig_moving_atom_names_non_mapped_same[i] << ":" << std::endl;
        }
    }

    for (auto& atom : residue.atoms) {
        std::string this_atom_name = atom.name;
        bool replace_name = false;
        std::string new_atom_name = "";

        if (std::find(orig_moving_atom_names_non_mapped_non_same.begin(),
                      orig_moving_atom_names_non_mapped_non_same.end(),
                      this_atom_name) !=
            orig_moving_atom_names_non_mapped_non_same.end()) {
            std::string ele = atom.element.name();
            new_atom_name = invent_new_name(this_atom_name, ele, residue_atom_names);
            residue_atom_names.push_back(new_atom_name);
            replace_name = true;
        } else {
            if (std::find(orig_moving_atom_names_non_mapped_same.begin(),
                          orig_moving_atom_names_non_mapped_same.end(),
                          this_atom_name) !=
                orig_moving_atom_names_non_mapped_same.end()) {
                // no change to replace_name.
            } else {
                replace_name = 1;
                for (const auto& pair : matching_atom_names) {
                    if (pair.first.first == this_atom_name) {
                        if (pair.second.first == pair.first.first) {
                            replace_name = 0;
                            break;
                        } else {
                            new_atom_name = pair.second.first;
                        }
                    }
                }
            }
        }

        if (replace_name) {
            atom.name = new_atom_name;
        }
    }
}

} // namespace coot