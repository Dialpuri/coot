#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace pli {

class dots_representation_info_t {
    bool is_closed_;
    std::string name_;

public:
    dots_representation_info_t() : is_closed_(false), name_("") {}
    explicit dots_representation_info_t(const std::string &n) : is_closed_(false), name_(n) {}

    // gemmi port of constructor dots_representation_info_t(mmdb::Manager* mol)
    // The original only sets is_closed=0 and performs commented-out selection work;
    // the gemmi version returns a default-constructed object.
    static dots_representation_info_t dots_representation_info_t_gemmi(gemmi::Structure const* mol) {
        (void)mol;  // not used — the original had add_dots() commented out
        dots_representation_info_t result;
        return result;
    }

    std::string get_name() const { return name_; }

    bool is_open_p() const {
        return !is_closed_;
    }

    void close_yourself() {
        is_closed_ = true;
    }
};

} // namespace pli