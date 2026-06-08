#include "function.hh"

namespace coot {

std::map<gemmi::Atom*, clipper::Coord_orth>
extra_restraints_t::position_point_map_gemmi(
    const gemmi::Structure* st_running,
    const gemmi::Structure* st_ref) const {

    std::map<gemmi::Atom*, clipper::Coord_orth> matching_atoms;

    if (st_running && st_ref) {
        for (const gemmi::Model& model_1 : st_running->models) {
            for (const gemmi::Model& model_2 : st_ref->models) {
                for (const gemmi::Chain& chain_1 : model_1.chains) {
                    std::string chain_1_id = chain_1.name;

                    for (const gemmi::Chain& chain_2 : model_2.chains) {
                        std::string chain_2_id = chain_2.name;
                        if (chain_1_id == chain_2_id) {
                            for (const gemmi::Residue& res_1 : chain_1.residues) {
                                int res_no_1 = res_1.seqid.num.value;
                                char ins_code_1 = res_1.seqid.icode;

                                for (const gemmi::Residue& res_2 : chain_2.residues) {
                                    int res_no_2 = res_2.seqid.num.value;
                                    char ins_code_2 = res_2.seqid.icode;

                                    if (res_no_2 == res_no_1) {
                                        if (ins_code_2 == ins_code_1) {
                                            for (const gemmi::Atom& at_1 : res_1.atoms) {
                                                std::string atom_name_1 = at_1.name;
                                                char alt_conf_1 = at_1.altloc;

                                                for (const gemmi::Atom& at_2 : res_2.atoms) {
                                                    std::string atom_name_2 = at_2.name;
                                                    char alt_conf_2 = at_2.altloc;

                                                    if (atom_name_2 == atom_name_1) {
                                                        if (alt_conf_2 == alt_conf_1) {
                                                            clipper::Coord_orth co = coot::co_gemmi(&at_2);
                                                            matching_atoms[const_cast<gemmi::Atom*>(&at_1)] = co;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                            break; // there won't be another residue that matches
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return matching_atoms;
}

} // namespace coot
