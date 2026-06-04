#pragma once

#include <string>

namespace pli {

class dots_representation_info_t {
    bool is_closed;
    std::string name_;

public:
    dots_representation_info_t() : is_closed(false) {}
    explicit dots_representation_info_t(const std::string &n) : name_(n), is_closed(false) {}

    // _gemmi ported constructors — static factory methods
    static dots_representation_info_t dots_representation_info_t_gemmi();
    static dots_representation_info_t dots_representation_info_t_gemmi(const std::string &n);

    std::string get_name() const { return name_; }

    bool is_open_p() const {
        int r = 1 - is_closed;
        return r;
    }

    void close_yourself() {
        is_closed = true;
    }
};

inline dots_representation_info_t dots_representation_info_t::dots_representation_info_t_gemmi() {
    return dots_representation_info_t();
}

inline dots_representation_info_t dots_representation_info_t::dots_representation_info_t_gemmi(const std::string &n) {
    return dots_representation_info_t(n);
}

} // namespace pli