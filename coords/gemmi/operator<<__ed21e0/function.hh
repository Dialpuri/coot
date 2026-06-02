#pragma once

#include <ostream>
#include <string>

struct symm_trans_t_gemmi {
    int symm_no;
    int x_shift_;
    int y_shift_;
    int z_shift_;
    std::string symm_as_string;

    symm_trans_t_gemmi(int n, int x, int y, int z)
        : symm_no(n), x_shift_(x), y_shift_(y), z_shift_(z) {}

    explicit symm_trans_t_gemmi(int idx)
        : symm_no(idx), x_shift_(0), y_shift_(0), z_shift_(0) {}

    symm_trans_t_gemmi()
        : symm_no(0), x_shift_(0), y_shift_(0), z_shift_(0) {}

    int isym() const { return symm_no; }
    int x()    const { return x_shift_; }
    int y()    const { return y_shift_; }
    int z()    const { return z_shift_; }
};

std::ostream & symm_trans_t_stream_gemmi(std::ostream &s, const symm_trans_t_gemmi &t) {
    s << "symm: " << t.symm_as_string << " (op-idx: " << t.isym() << ") trans: "
      << t.x() << " " << t.y() << " " << t.z();
    return s;
}