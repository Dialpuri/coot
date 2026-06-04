#include "function.hh"

#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

residue_id_t_gemmi get_id_gemmi(const tree<linked_residue_t_gemmi>& glyco_tree,
                                gemmi::Residue* residue_p) {
    residue_id_t_gemmi id;

    if (!residue_p) {
        return id;
    }

    for (auto it = glyco_tree.begin(); it != glyco_tree.end(); ++it) {
        if (it->residue == residue_p) {
            if (it.node->parent) {
                gemmi::Residue* parent_res = it.node->parent->data.residue;

                std::string parent_res_type = (parent_res ? parent_res->name : "");
                std::string link_type = it->link_type;
                std::string res_type = residue_p->name;

                int level = get_level_gemmi(glyco_tree, residue_p);

                int prime_flag = 0;

                residue_spec_t ps;
                ps.res_no = parent_res ? parent_res->seqid.num.value : 0;
                ps.ins_code = parent_res ? std::string(1, parent_res->seqid.icode) : "";

                id = residue_id_t_gemmi(level, prime_flag, res_type, link_type,
                                        parent_res_type, std::move(ps));
                break;
            }
            break;
        }
    }

    return id;
}

}  // namespace coot