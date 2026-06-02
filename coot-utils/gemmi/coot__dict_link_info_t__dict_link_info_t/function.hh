#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include "coot/geometry/protein-geometry.hh"
#include "coot/geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__distance/gemmi/function.hh"

namespace coot {

// Lightweight mirror of the real coot::dict_link_info_t that uses gemmi types.
// We do NOT include coot-utils/dict-link-info.hh here to avoid name collision
// with the real mmdb-based class.
class dict_link_info_t {
public:
    atom_spec_t spec_ref;
    atom_spec_t spec_new;
    double dist = 0.0;

    dict_link_info_t() = default;

    static dict_link_info_t dict_link_info_t_gemmi(
        const gemmi::Chain *chain_ref,
        const gemmi::Residue *residue_ref,
        const gemmi::Chain *chain_new,
        const gemmi::Residue *residue_new,
        const std::string &link_type,
        const coot::protein_geometry &geom);
};

inline dict_link_info_t dict_link_info_t::dict_link_info_t_gemmi(
    const gemmi::Chain *chain_ref,
    const gemmi::Residue *residue_ref,
    const gemmi::Chain *chain_new,
    const gemmi::Residue *residue_new,
    const std::string &link_type,
    const coot::protein_geometry &geom)
{
    if (!residue_ref) {
        throw std::runtime_error("Null residue_ref");
    }
    if (!residue_new) {
        throw std::runtime_error("Null residue_new");
    }

    coot::dictionary_residue_link_restraints_t rr = geom.link(link_type);

    if (rr.link_id == "") {
        std::string rn1 = residue_ref->name;
        std::string rn2 = residue_new->name;
        std::string mess = "Link not found in dictionary " + link_type + " between " + rn1 + " " + rn2;
        throw std::runtime_error(mess);
    }

    // check_for_order_switch logic inline (mirrors the original private method)
    bool order_switch_flag = false;
    try {
        std::string comp_id_ref = residue_ref->name;
        std::string comp_id_new = residue_new->name;

        std::string group_ref = geom.get_group(comp_id_ref);
        std::string group_new = geom.get_group(comp_id_new);

        std::vector<coot::chem_link> link_infos_f =
            geom.matching_chem_links(comp_id_ref, group_ref, comp_id_new, group_new);
        std::vector<coot::chem_link> link_infos_b =
            geom.matching_chem_links(comp_id_new, group_new, comp_id_ref, group_ref);

        struct link_info_pair {
            coot::chem_link link;
            bool is_backwards;
        };
        std::vector<link_info_pair> link_infos;
        for (const auto &link : link_infos_f)
            link_infos.push_back({link, false});
        for (const auto &link : link_infos_b)
            link_infos.push_back({link, true});

        for (const auto &lp : link_infos) {
            if (lp.link.Id() == link_type) {
                order_switch_flag = lp.is_backwards;
                break;
            }
        }
    } catch (const std::runtime_error &) {
        // Swallow exceptions like the original
    }

    const gemmi::Residue *res_1 = residue_ref;
    const gemmi::Residue *res_2 = residue_new;
    const gemmi::Chain *ch_1 = chain_ref;
    const gemmi::Chain *ch_2 = chain_new;

    if (order_switch_flag) {
        std::swap(res_1, res_2);
        std::swap(ch_1, ch_2);
    }

    // Build residue_spec_t from gemmi residues
    // gemmi uses ' ' for no insertion code; coot uses ""
    auto make_res_spec = [](const std::string &chain_name, const gemmi::Residue &res)
        -> coot::residue_spec_t
    {
        char ic = res.seqid.icode;
        std::string ins_code;
        if (ic != ' ')
            ins_code = std::string(1, ic);
        return coot::residue_spec_t(chain_name, res.seqid.num.value, ins_code);
    };

    coot::residue_spec_t res_spec_ref = make_res_spec(ch_1->name, *res_1);
    coot::residue_spec_t res_spec_new = make_res_spec(ch_2->name, *res_2);

    bool ifound = false;
    for (unsigned int ibond = 0; ibond < rr.link_bond_restraint.size() && !ifound; ++ibond) {
        for (const gemmi::Atom &at1 : res_1->atoms) {
            if (at1.name == rr.link_bond_restraint[ibond].atom_id_1_4c()) {
                for (const gemmi::Atom &at2 : res_2->atoms) {
                    if (at2.name == rr.link_bond_restraint[ibond].atom_id_2_4c()) {
                        ifound = true;
                        break;
                    }
                }
            }
            if (ifound) break;
        }
    }

    if (!ifound) {
        throw std::runtime_error("Dictionary links atom not found in link residues");
    }

    // Fill in the matching atom specs
    dict_link_info_t result;
    for (unsigned int ibond = 0; ibond < rr.link_bond_restraint.size(); ++ibond) {
        if (!result.spec_ref.atom_name.empty()) break;
        for (const gemmi::Atom &at1 : res_1->atoms) {
            if (at1.name == rr.link_bond_restraint[ibond].atom_id_1_4c()) {
                for (const gemmi::Atom &at2 : res_2->atoms) {
                    if (at2.name == rr.link_bond_restraint[ibond].atom_id_2_4c()) {
                        result.spec_ref = coot::atom_spec_t(
                            res_spec_ref.chain_id, res_spec_ref.res_no,
                            res_spec_ref.ins_code, at1.name, "");
                        result.spec_new = coot::atom_spec_t(
                            res_spec_new.chain_id, res_spec_new.res_no,
                            res_spec_new.ins_code, at2.name, "");
                        result.dist = coot::distance_gemmi(&at1, &at2);
                        break;
                    }
                }
            }
        }
    }

    return result;
}

} // namespace coot