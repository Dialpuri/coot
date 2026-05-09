#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot { namespace molecule_t {

bool assign_sequence_gemmi(gemmi::Model& model, const std::string& chain_id,
                          const std::string& sequence);

std::string single_letter_to_3_letter_code(char letter);

}}