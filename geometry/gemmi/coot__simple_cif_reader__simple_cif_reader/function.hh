#pragma once

#include <string>
#include <vector>

namespace coot {

class simple_cif_reader {
public:
    std::vector<std::string> names;
    std::vector<std::string> three_letter_codes;

    simple_cif_reader() = default;

    void simple_cif_reader_gemmi(const std::string &cif_dictionary_file_name);

    bool has_restraints_for(const std::string &name) const;
};

} // namespace coot