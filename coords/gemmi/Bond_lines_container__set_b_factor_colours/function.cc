#include "function.hh"

int Bond_lines_container::set_b_factor_colours_gemmi(gemmi::Structure& st) {
    const float max_b_factor = 70.0f;
    b_factor_fractions_.clear();

    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    float b_factor = atom.b_iso;
                    float bs = b_factor * b_factor_scale;
                    float f = bs / max_b_factor;
                    if (f < 0.0f) f = 0.0f;
                    if (f > 1.0f) f = 1.0f;

                    std::string key = chain.name + "/" +
                        std::to_string(res.seqid.num.value) + "/" +
                        atom.name;
                    b_factor_fractions_[key] = f;
                }
            }
        }
    }

    return 16777217;
}

bool Bond_lines_container::get_b_factor_fraction(
    const gemmi::Chain& chain,
    const gemmi::Residue& res,
    const gemmi::Atom& atom,
    float& out_val) const
{
    std::string key = chain.name + "/" +
        std::to_string(res.seqid.num.value) + "/" +
        atom.name;
    auto it = b_factor_fractions_.find(key);
    if (it != b_factor_fractions_.end()) {
        out_val = it->second;
        return true;
    }
    return false;
}