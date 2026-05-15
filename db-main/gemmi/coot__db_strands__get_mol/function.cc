#include "function.hh"
#include <gemmi/pdb.hpp>

namespace coot {
namespace db_strands {

gemmi::Structure *get_mol_gemmi(const std::string &filename) {
    try {
        auto *mol = new gemmi::Structure{gemmi::read_pdb_file(filename)};
        std::cout << "Read OK: " << filename << std::endl;
        return mol;
    } catch (...) {
        std::cout << "Error reading " << filename << std::endl;
        return nullptr;
    }
}

} // namespace db_strands
} // namespace coot