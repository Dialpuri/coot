#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

fragment_container_t make_overlapping_fragments_gemmi(
    const gemmi::Structure &st,
    const std::string &chain_id,
    unsigned int fragment_length) {

    fragment_container_t fc;

    if (st.models.empty()) return fc;

    const gemmi::Model &model = st.models[0];

    for (const gemmi::Chain &chain : model.chains) {
        if (chain.name == chain_id) {
            if (chain.residues.empty()) continue;

            int n_res = static_cast<int>(chain.residues.size());

            for (int i_res = 0; i_res < n_res; i_res++) {
                const gemmi::Residue &res = chain.residues[i_res];
                int resno_start = res.seqid.num.value;
                int this_run_max = static_cast<int>(fragment_length);
                if ((i_res + this_run_max) >= n_res) {
                    this_run_max = n_res - i_res - 1;
                }

                std::vector<const gemmi::Residue*> residues_running;
                residues_running.push_back(&res);

                for (int j_res = 1; j_res < this_run_max; j_res++) {
                    const gemmi::Residue &r = chain.residues[i_res + j_res];
                    int resno_this = r.seqid.num.value;
                    if ((resno_this - resno_start) == j_res) {
                        residues_running.push_back(&r);
                    }
                }

                if (residues_running.size() > 10) {
                    const gemmi::Residue *start_res = residues_running.front();
                    const gemmi::Residue *end_res = residues_running.back();

                    fragment_container_t::fragment_range_t fr(
                        chain_id,
                        residue_spec_t(chain.name, start_res->seqid.num.value, std::string(1, start_res->seqid.icode)),
                        residue_spec_t(chain.name, end_res->seqid.num.value, std::string(1, end_res->seqid.icode))
                    );
                    fc.add(fr);
                }
            }
        }
    }

    return fc;
}

} // namespace coot