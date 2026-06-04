#include "function.hh"
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/coot-utils/gemmi/make_asc/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/coot-utils/gemmi/make_asc/function.cc"

namespace coot {

std::pair<std::unique_ptr<gemmi::Residue>, atom_selection_container_t>
deep_copy_this_residue_and_make_asc_gemmi(
    const gemmi::Structure& orig_st,
    gemmi::CRA cra,
    const std::string& altconf,
    short int whole_residue_flag,
    int /*atom_index_handle*/,
    int /*udd_afix_handle*/) {

    // ---- Step 1: Deep copy the residue, filtering atoms by altconf ----

    auto new_res = std::make_unique<gemmi::Residue>();
    *new_res = cra.residue->empty_copy();

    for (const auto& atom : cra.residue->atoms) {
        std::string this_atom_alt_loc;
        if (atom.altloc != '\0') {
            this_atom_alt_loc = std::string(1, atom.altloc);
        }

        if (whole_residue_flag || this_atom_alt_loc == altconf || this_atom_alt_loc == "") {
            gemmi::Atom new_atom;
            new_atom.name      = atom.name;
            new_atom.element   = atom.element;
            new_atom.pos       = atom.pos;
            new_atom.occ       = atom.occ;
            new_atom.b_iso     = atom.b_iso;
            new_atom.altloc    = atom.altloc;
            new_atom.charge    = atom.charge;
            new_atom.serial    = atom.serial;
            new_atom.calc_flag = atom.calc_flag;
            new_atom.flag      = atom.flag;
            new_res->atoms.push_back(std::move(new_atom));
        }
    }

    // ---- Step 2: Build a gemmi::Structure from the residue selection
    //    (equivalent to MMDB create_mmdbmanager_from_res_selection) ----

    gemmi::Structure result_st;
    result_st.name = orig_st.name + "_selected";
    result_st.cell = orig_st.cell;
    result_st.spacegroup_hm = orig_st.spacegroup_hm;

    gemmi::Model& new_model = result_st.models.emplace_back();
    gemmi::Chain& new_chain = new_model.chains.emplace_back();
    new_chain.name = cra.chain->name;

    std::vector<gemmi::Residue*> res_vec;
    res_vec.push_back(cra.residue);

    int atom_index_counter = 0;

    for (int ires = 0; ires < static_cast<int>(res_vec.size()); ires++) {
        bool is_flanking = (ires == 0) || (ires == static_cast<int>(res_vec.size()) - 1);
        short int whole_res_flag_internal = 0;

        // residue_from_alt_conf_split_flag = 1 (from original call),
        // so flanking residues are NOT whole (they respect altconf too)
        if (is_flanking) {
            if (!1)  // !residue_from_alt_conf_split_flag
                whole_res_flag_internal = 1;
        }

        if (altconf == "*")
            whole_res_flag_internal = 1;

        const gemmi::Residue& src_res = *res_vec[ires];

        gemmi::Residue new_res_in_mol;
        new_res_in_mol.name  = src_res.name;
        new_res_in_mol.seqid = src_res.seqid;

        for (const gemmi::Atom& src_atom : src_res.atoms) {
            bool include_atom = false;
            if (whole_res_flag_internal) {
                include_atom = true;
            } else if (!altconf.empty() && src_atom.altloc == altconf[0]) {
                include_atom = true;
            } else if (src_atom.altloc == '\0') {
                include_atom = true;
            }

            if (include_atom) {
                gemmi::Atom new_atom;
                new_atom.name      = src_atom.name;
                new_atom.element   = src_atom.element;
                new_atom.pos       = src_atom.pos;
                new_atom.occ       = src_atom.occ;
                new_atom.b_iso     = src_atom.b_iso;
                new_atom.altloc    = src_atom.altloc;
                new_atom.charge    = src_atom.charge;
                new_atom.serial    = ++atom_index_counter;
                new_atom.calc_flag = src_atom.calc_flag;
                new_atom.flag      = src_atom.flag;
                new_res_in_mol.atoms.push_back(std::move(new_atom));
            }
        }

        new_chain.residues.push_back(std::move(new_res_in_mol));
    }

    // ---- Step 3: Build atom_selection_container_t ----
    atom_selection_container_t asc = make_asc_gemmi(result_st, false);

    return std::make_pair(std::move(new_res), asc);
}

} // namespace coot