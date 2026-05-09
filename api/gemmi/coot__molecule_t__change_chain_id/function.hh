#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Change chain ID from from_chain_id to to_chain_id
// Returns: (status, message)
// status=1 on success, 0 on failure (e.g., conflict)
// message contains details
std::pair<int, std::string>
change_chain_id_gemmi(gemmi::Structure& st,
                      const std::string& from_chain_id,
                      const std::string& to_chain_id,
                      bool use_resno_range,
                      int start_resno,
                      int end_resno);

} // namespace molecule_t
} // namespace coot