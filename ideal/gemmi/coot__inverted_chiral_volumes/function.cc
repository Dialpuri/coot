#include "function.hh"

namespace coot {

std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t>>
inverted_chiral_volumes_gemmi(int imol,
                              const gemmi::Structure& st,
                              coot::protein_geometry* geom_p,
                              int cif_dictionary_read_number) {

    std::vector<coot::atom_spec_t> v;
    std::vector<std::string> unknown_types_vec;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                if (static_cast<int>(residue.atoms.size()) <= 3) {
                    continue;
                }

                std::string residue_type(residue.name);
                if (residue_type == "UNK") {
                    residue_type = "ALA";
                }

                if (!geom_p->have_dictionary_for_residue_type(residue_type, imol, cif_dictionary_read_number)) {
                    bool irfound = false;
                    for (unsigned int ir = 0; ir < unknown_types_vec.size(); ir++) {
                        if (unknown_types_vec[ir] == residue_type) {
                            irfound = true;
                            break;
                        }
                    }
                    if (!irfound) {
                        unknown_types_vec.push_back(residue_type);
                    }
                } else {
                    std::vector<coot::dict_chiral_restraint_t> chiral_restraints =
                        geom_p->get_monomer_chiral_volumes(std::string(residue.name), imol);
                    coot::dict_chiral_restraint_t chiral_restraint;
                    for (unsigned int irestr = 0; irestr < chiral_restraints.size(); irestr++) {
                        chiral_restraint = chiral_restraints[irestr];
                        if (!chiral_restraint.is_a_both_restraint()) {
                            std::vector<std::pair<short int, coot::atom_spec_t>> c =
                                coot::is_inverted_chiral_atom_p_gemmi(
                                    residue,
                                    chiral_restraint.atom_id_c_4c(),
                                    chiral_restraint.atom_id_1_4c(),
                                    chiral_restraint.atom_id_2_4c(),
                                    chiral_restraint.atom_id_3_4c(),
                                    chiral_restraint.volume_sign);
                            for (unsigned int ibad = 0; ibad < c.size(); ibad++) {
                                if (c[ibad].first) {
                                    v.push_back(coot::atom_spec_t(
                                        chain.name,
                                        residue.seqid.num.value,
                                        std::string(1, residue.seqid.icode),
                                        chiral_restraint.atom_id_c_4c(),
                                        c[ibad].second.alt_conf));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t>>(unknown_types_vec, v);
}

} // namespace coot