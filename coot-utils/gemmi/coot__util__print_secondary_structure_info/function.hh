#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <iostream>

namespace coot {
namespace util {

// Gemmi port of coot::util::print_secondary_structure_info
// Prints helix and sheet secondary structure info from a gemmi::Structure.
// In gemmi, helices and sheets are stored on Structure (not Model).
void print_secondary_structure_info_gemmi(const gemmi::Structure& st) {

  int nhelix = static_cast<int>(st.helices.size());
  int nsheet = static_cast<int>(st.sheets.size());

  std::cout << "INFO:: There are " << nhelix << " helices and "
            << nsheet << " sheets" << std::endl;

  if (nhelix > 0) {
    std::cout << "               Helix info: " << std::endl;
    std::cout << "------------------------------------------------\n";
    for (const auto& helix : st.helices) {
      // Original MMDB printed: serNum helixID initChainID initSeqNum endChainID endSeqNum length comment
      std::cout << helix.pdb_helix_class
                << " " << helix.start.chain_name
                << " " << helix.start.res_id.seqid.num.value
                << " " << helix.end.chain_name
                << " " << helix.end.res_id.seqid.num.value
                << " " << helix.length
                << " " << helix.start.res_id.seqid.icode << std::endl;
    }
  }

  if (nsheet > 0) {
    std::cout << "               Sheet info: " << std::endl;
    std::cout << "------------------------------------------------\n";
    for (const auto& sheet : st.sheets) {
      for (const auto& strand : sheet.strands) {
        // Original MMDB printed: sheetID strandNo initChainID initSeqNum endChainID endSeqNum
        std::cout << sheet.name
                  << " " << strand.sense
                  << " " << strand.start.chain_name
                  << " " << strand.start.res_id.seqid.num.value
                  << " " << strand.end.chain_name
                  << " " << strand.end.res_id.seqid.num.value
                  << std::endl;
      }
    }
    std::cout << "------------------------------------------------\n";
  }
}

} // namespace util
} // namespace coot