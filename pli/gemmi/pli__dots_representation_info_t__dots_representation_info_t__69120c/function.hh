#pragma once

#include <string>

namespace pli {

class dots_representation_info_t {
    std::string name_;
    bool is_closed_;

public:
    // Gemmi ported factory — same semantics as the MMDB original constructor
    static dots_representation_info_t dots_representation_info_t_gemmi(const std::string &n) {
        dots_representation_info_t obj;
        obj.name_ = n;
        obj.is_closed_ = false;
        return obj;
    }

    std::string get_name() const { return name_; }

    bool is_open_p() const {
        int r = 1 - is_closed_;
        return r;
    }
};

} // namespace pli