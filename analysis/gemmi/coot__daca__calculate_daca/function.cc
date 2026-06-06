#include "function.hh"
#include <algorithm>
#include <cmath>
#include <map>
#include <iostream>
#include <iomanip>

// Include verified gemmi ports
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__daca__get_daca_fragments/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__daca__get_frag_to_reference_rtop/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_main_chain_p/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__daca__atom_is_close_to_a_residue_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__daca__atom_is_neighbour_mainchain/gemmi/function.hh"

// Reference data storage for REFERENCE mode
static std::map<coot::daca::daca_key_t, int> daca_data;

// Reference fragments (like in the original daca class member)
static std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> reference_fragments_gemmi;

// Helical residues (empty by default — can be filled by fill_helix_flags_gemmi)
static std::vector<const gemmi::Residue*> helical_residues_gemmi;

namespace coot {
namespace daca {

void add_to_box(
    coot::daca::mode_t mode,
    const std::string& res_name,
    bool helical_flag,
    unsigned int frag_idx,
    const box_index_t& box_index,
    const std::string& atom_type)
{
    daca_key_t key{res_name, helical_flag, frag_idx, box_index, atom_type};
    daca_data[key] += 1;
}

int get_reference_counts(
    const std::string& res_name,
    bool helical_flag,
    unsigned int frag_idx,
    const box_index_t& box_index,
    const std::string& atom_type)
{
    daca_key_t key{res_name, helical_flag, frag_idx, box_index, atom_type};
    auto it = daca_data.find(key);
    if (it != daca_data.end())
        return it->second;
    return 0;
}

int calculate_daca_gemmi(
    const gemmi::CRA& reference_cra,
    const std::vector<std::pair<gemmi::CRA, std::string>>& typed_atoms,
    coot::daca::mode_t mode)
{
    bool print_scores = true;

    double d_crit = 8.0;
    double dd_crit = d_crit * d_crit;

    int reference_counts = 0;

    if (!reference_cra.residue) return 0;

    std::string res_name(reference_cra.residue->name);
    int reference_residue_seqnum = reference_cra.residue->seqid.num.value;

    // Get fragments using the gemmi port
    std::vector<std::vector<const gemmi::Atom*>> fragments =
        coot::get_daca_fragments_gemmi(reference_cra.residue);

    if (false)
        std::cout << "debug:: fragments.size() " << fragments.size() << " "
                  << residue_spec_t_gemmi(reference_cra)
                  << " " << reference_cra.residue->name << std::endl;

    for (unsigned int ifrag = 0; ifrag < fragments.size(); ifrag++) {
        const std::vector<const gemmi::Atom*>& atom_vec = fragments[ifrag];
        std::vector<clipper::Coord_orth> reference_positions_vec;
        clipper::Coord_orth sum(0, 0, 0);

        for (const gemmi::Atom* at : atom_vec) {
            clipper::Coord_orth pos = coot::co_gemmi(at);
            reference_positions_vec.push_back(pos);
            sum += pos;
        }

        if (reference_positions_vec.size() > 2) {
            if (reference_positions_vec.size() == atom_vec.size()) {
                double m = 1.0 / static_cast<double>(reference_positions_vec.size());
                clipper::Coord_orth frag_centre(sum * m);

                // Get the RTop that transforms the fragment to a reference fragment at the origin
                std::pair<bool, clipper::RTop_orth> frag_to_reference_rtop_pair =
                    coot::daca::get_frag_to_reference_rtop_gemmi(res_name, ifrag, atom_vec, reference_fragments_gemmi);

                if (!frag_to_reference_rtop_pair.first) continue;
                const clipper::RTop_orth& frag_to_reference_rtop = frag_to_reference_rtop_pair.second;

                for (unsigned int ita = 0; ita < typed_atoms.size(); ita++) {
                    const gemmi::CRA& atom_cra = typed_atoms[ita].first;
                    const std::string& atom_type = typed_atoms[ita].second;

                    if (!atom_cra.atom) continue;

                    // don't consider atoms in this residue, of course
                    if (atom_cra.residue == reference_cra.residue)
                        continue;

                    // don't consider peptide neighbour mainchain
                    int res_no_delta = atom_cra.residue->seqid.num.value - reference_residue_seqnum;
                    if (std::abs(res_no_delta) < 2)
                        if (atom_cra.chain == reference_cra.chain)
                            if (coot::is_main_chain_p_gemmi(*atom_cra.atom, *atom_cra.residue))
                                continue;

                    double dx = atom_cra.atom->pos.x - frag_centre.x();
                    double dy = atom_cra.atom->pos.y - frag_centre.y();
                    double dz = atom_cra.atom->pos.z - frag_centre.z();
                    double dd = dx*dx + dy*dy + dz*dz;

                    if (dd < dd_crit) {
                        // Good, found something
                        if (coot::daca::atom_is_close_to_a_residue_atom_gemmi(*atom_cra.atom, *reference_cra.residue)) {
                            if (!coot::daca::atom_is_neighbour_mainchain_gemmi(atom_cra, reference_cra)) {
                                clipper::Coord_orth at_pos = coot::co_gemmi(atom_cra.atom);
                                clipper::Coord_orth transformed_pos = frag_to_reference_rtop * at_pos;
                                box_index_t box_index(transformed_pos);
                                bool helical_flag = false;
                                if (std::find(helical_residues_gemmi.begin(), helical_residues_gemmi.end(),
                                              reference_cra.residue) != helical_residues_gemmi.end())
                                    helical_flag = true;

                                if (mode == REFERENCE)
                                    add_to_box(mode, res_name, helical_flag, ifrag, box_index, atom_type);
                                if (mode == ANALYSIS) {
                                    // what reference score do we have for this box
                                    std::string box_key = res_name + "-non-helical";
                                    if (helical_flag) box_key = res_name + "-helical";
                                    int counts = get_reference_counts(res_name, helical_flag, ifrag, box_index, atom_type);
                                    if (counts > 0) {
                                        reference_counts += counts;
                                        if (print_scores)
                                            std::cout << "Score " << residue_spec_t_gemmi(reference_cra) << " "
                                                      << box_key << " " << ifrag << " " << " "
                                                      << atom_spec_t_gemmi(*atom_cra.atom) << " " << atom_type << " "
                                                      << std::setw(2) << box_index.idx_x << " "
                                                      << std::setw(2) << box_index.idx_y << " "
                                                      << std::setw(2) << box_index.idx_z << " "
                                                      << counts << "\n";
                                    } else {
                                        std::cout << "Miss " << residue_spec_t_gemmi(reference_cra) << " "
                                                  << box_key << " " << ifrag << " " << " "
                                                  << atom_spec_t_gemmi(*atom_cra.atom) << " " << atom_type << " "
                                                  << std::setw(2) << box_index.idx_x << " "
                                                  << std::setw(2) << box_index.idx_y << " "
                                                  << std::setw(2) << box_index.idx_z << " "
                                                  << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                std::cout << "OOps in atom set vs reference set size test " << std::endl;
            }
        } else {
            std::cout << "ERROR:: in calculate_daca(): This can't happen. reference positions size "
                      << reference_positions_vec.size() << " "
                      << residue_spec_t_gemmi(reference_cra) << std::endl;
        }
    }
    return reference_counts;
}

} // namespace daca
} // namespace coot