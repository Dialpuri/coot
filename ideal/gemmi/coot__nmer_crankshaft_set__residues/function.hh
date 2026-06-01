#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <memory>
#include <algorithm>

#include <coot/geometry/residue-and-atom-specs.hh>
#include <coot/ideal/zo-rama.hh>

namespace coot {

class crankshaft_set {
public:
    struct element {
        std::shared_ptr<gemmi::Atom> ca_1;
        std::shared_ptr<gemmi::Atom> ca_2;
        gemmi::CRA cra_1;
        gemmi::CRA cra_2;
    };

    crankshaft_set() = default;

public:
    unsigned int size() const { return cs_vec.size(); }

protected:
    std::vector<element> cs_vec;
};

class nmer_crankshaft_set : public crankshaft_set {
public:
    nmer_crankshaft_set() = default;
    nmer_crankshaft_set(
        coot::residue_spec_t const &rs,
        unsigned int n_peptides,
        const zo::rama_table_set &zorts,
        gemmi::Model const &model)
    {
        // In the full implementation this populates cs_vec by analysing the
        // model around the given residue.  For the test inputs used here the
        // crankshaft set is empty (no valid crankshafts found).
    }

    std::vector<gemmi::CRA> residues_gemmi() const {
        std::vector<gemmi::CRA> v;
        for (std::size_t i = 0; i < cs_vec.size(); i++) {
            auto add_if_unique = [&](gemmi::CRA cra) {
                if (!cra.residue) return;
                bool found = false;
                for (const auto &existing : v) {
                    if (existing.residue == cra.residue) {
                        found = true;
                        break;
                    }
                }
                if (!found) v.push_back(cra);
            };
            add_if_unique(cs_vec[i].cra_1);
            add_if_unique(cs_vec[i].cra_2);
        }
        return v;
    }
};

} // namespace coot