#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

std::pair<int, gemmi::CRA>
find_serial_number_for_insert_gemmi(int seqnum_for_new,
                                    const std::string& ins_code_for_new,
                                    const std::string& chain_id,
                                    const gemmi::Model& model);

}