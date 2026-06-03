#pragma once

#include <string>

namespace coot {
namespace sequence_assignment {

enum side_chain_name_index {
    ALA, GLY, TRP, VAL
};

class side_chain_score_t {
public:
    float auto_fit_score_gemmi(const std::string &chain_id,
                               int resno,
                               const side_chain_name_index &idx);
};

} // namespace sequence_assignment
} // namespace coot