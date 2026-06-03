#pragma once
#include <gemmi/model.hpp>
#include <vector>
#include <queue>
#include <map>
#include <set>

namespace coot {

class match_container_for_residues_t {
public:
    gemmi::Residue *residue_1;
    gemmi::Residue *residue_2;
    std::vector<std::pair<gemmi::Atom *, gemmi::Atom *>> atom_pairs;

    match_container_for_residues_t(gemmi::Residue *r1, gemmi::Residue *r2)
        : residue_1(r1), residue_2(r2) {}
    match_container_for_residues_t() : residue_1(nullptr), residue_2(nullptr) {}

    std::vector<gemmi::Residue*> residue_vector_from_residue_gemmi(gemmi::CRA cra) const {
        std::vector<gemmi::Residue*> v;
        if (!cra.chain || !cra.residue) return v;

        gemmi::Chain* chain_p = cra.chain;

        // Build lookup map: (seqno, icode) -> index in chain.residues
        using Key = std::pair<int, char>;
        std::map<Key, int> seqno_map;
        for (int i = 0; i < static_cast<int>(chain_p->residues.size()); i++) {
            const auto& res = chain_p->residues[i];
            seqno_map[Key{res.seqid.num.value, res.seqid.icode}] = i;
        }

        // BFS from the starting residue
        std::queue<int> q;
        std::set<Key> visited;

        Key start_key{cra.residue->seqid.num.value, cra.residue->seqid.icode};
        auto sit = seqno_map.find(start_key);
        if (sit == seqno_map.end()) return v;

        q.push(sit->second);
        visited.insert(start_key);

        while (!q.empty()) {
            int idx = q.front();
            q.pop();

            const auto& res = chain_p->residues[idx];
            int res_no = res.seqid.num.value;
            v.push_back(&chain_p->residues[idx]);

            // Previous residue (seqno - 1, no insertion code → icode = ' ')
            {
                Key prev_key{res_no - 1, ' '};
                auto pit = seqno_map.find(prev_key);
                if (pit != seqno_map.end() && visited.find(prev_key) == visited.end()) {
                    visited.insert(prev_key);
                    q.push(pit->second);
                }
            }

            // Next residue (seqno + 1, no insertion code → icode = ' ')
            {
                Key next_key{res_no + 1, ' '};
                auto nit = seqno_map.find(next_key);
                if (nit != seqno_map.end() && visited.find(next_key) == visited.end()) {
                    visited.insert(next_key);
                    q.push(nit->second);
                }
            }
        }

        return v;
    }
};

}  // namespace coot