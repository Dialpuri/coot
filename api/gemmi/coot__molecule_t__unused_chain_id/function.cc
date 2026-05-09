#include "function.hh"

namespace coot {
namespace molecule_t {

std::pair<bool, std::string>
unused_chain_id_gemmi(const gemmi::Model& model) {
    std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    std::pair<bool, std::string> s(false, "");
    
    for (const gemmi::Chain& chain : model.chains) {
        std::string this_chain_id = chain.name;
        std::string::size_type idx = r.find(this_chain_id);
        if (idx != std::string::npos) {
            r.erase(idx, 1);
        }
    }
    
    if (r.length() > 0) {
        s.first = true;
        s.second = r.substr(0, 1);
    }
    
    return s;
}

}
}