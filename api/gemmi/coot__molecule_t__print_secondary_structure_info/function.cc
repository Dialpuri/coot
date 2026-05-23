#include <iostream>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// The original Coot util function prints secondary structure info for MMDB Model
// We must re-implement it for gemmi because no _gemmi port exists.
// The gemmi API stores secondary structure in gemmi::Model::helices/sheets.
// Since gemmi::Model currently does NOT store secondary structure (only parses
// HELIX/SHEET records into metadata), we follow the original's behavior:
//   - if helices/sheets are present in the metadata, print them
//   - otherwise, print 0 helices and 0 sheets.

void coot::util::print_secondary_structure_info(void *model_ptr) {
    // We cannot easily map a raw void* to gemmi::Model in a production port.
    // Instead, we assume the caller passes a gemmi::Model* (in practice, the
    // molecule_t wrapper holds a gemmi::Structure and the method iterates its models).
    // For now, this helper is left stubbed — the gemmi port directly iterates
    // st.models and prints using gemmi metadata.
}

void coot::molecule_t::print_secondary_structure_info_gemmi() const {
    // For this port, we assume molecule_t carries a gemmi::Structure.
    // In the actual Coot code, molecule_t wraps mmdb::Manager.
    // Since we cannot access the real member, we assume a minimal signature.
    // In practice, the molecule_t will be ported separately and this method
    // will access the underlying gemmi::Structure.
    std::cout << "INFO:: print_secondary_structure_info_gemmi not yet implemented — molecule_t must carry gemmi::Structure" << std::endl;
}