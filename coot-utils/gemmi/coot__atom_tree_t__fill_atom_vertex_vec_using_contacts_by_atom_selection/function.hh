#pragma once
#include <vector>
#include <queue>
#include <map>
#include <string>
#include <gemmi/model.hpp>
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"

namespace coot {

class atom_tree_t {
public:
    std::vector<atom_vertex> atom_vertex_vec;
    std::map<std::string, map_index_t> name_to_index;

    bool fill_atom_vertex_vec_using_contacts_by_atom_selection_gemmi(
        const std::vector<std::vector<int> > &contact_indices,
        const std::vector<gemmi::Atom>& residue_atoms,
        int base_atom_index) {
        bool debug = false;
        if (debug)
            std::cout << ":::::::::::: fill_atom_vertex_vec_using_contacts_by_atom_selection() --- start -- with n_atoms "
                      << residue_atoms.size() << std::endl;

        bool r = false;
        atom_vertex av;
        atom_vertex_vec.resize(residue_atoms.size());
        av.connection_type = atom_vertex::START;
        atom_vertex_vec[base_atom_index] = av;

        // fail to set up
        if (contact_indices.size() == 0) {
            atom_vertex_vec.clear();
            return 0;
        }

        std::queue<int> q;
        q.push(base_atom_index);
        std::vector<int> done; // things that were in the queue that are
                               // not done (on poping). So that we don't
                               // loop endlessly doing atom indices that we
                               // have already considered.

        while (q.size()) {
            int this_base_atom = q.front();
            // now what are the forward atoms of av?
            std::vector<int> av_contacts = contact_indices[this_base_atom]; // size check above

            for (size_t iav=0; iav<av_contacts.size(); iav++) {

                int i_forward = av_contacts[iav];

                // Check that this forward atom is not already in the forward
                // atoms of this_base_atom
                bool ifound_forward = 0;
                for (size_t ifo=0; ifo<atom_vertex_vec[this_base_atom].forward.size(); ifo++) {
                    if (atom_vertex_vec[this_base_atom].forward[ifo] == static_cast<size_t>(av_contacts[iav])) {
                        ifound_forward = 1;
                        break;
                    }
                }

                if (! ifound_forward) {
                    // Add the contact as a forward atom of this_base_atom but
                    // only if the forward atom does not have a forward atom
                    // which is this_base_atom (and this keeps the tree going in
                    // up and not back again).
                    bool ifound_forward_forward = 0;
                    for (size_t ifo=0; ifo<atom_vertex_vec[i_forward].forward.size(); ifo++) {
                        if (atom_vertex_vec[i_forward].forward[ifo] == static_cast<size_t>(this_base_atom)) {
                            ifound_forward_forward = 1;
                            break;
                        }
                    }
                    if (!ifound_forward_forward) {
                        atom_vertex_vec[this_base_atom].forward.push_back(av_contacts[iav]);
                    }
                }

                // add the forward atoms to the queue if they are not already
                // in the done list.
                bool in_done = 0;
                for (size_t idone=0; idone<done.size(); idone++) {
                    if (done[idone] == static_cast<size_t>(av_contacts[iav])) {
                        in_done = 1;
                        break;
                    }
                }
                if (!in_done)
                    q.push(av_contacts[iav]);
            }

            // if the forward atoms of this_atom_index do not have a back
            // atom, add one (but not, of course, if forward atom is the
            // start point of the tree).
            for (size_t iav=0; iav<av_contacts.size(); iav++) {
                if (atom_vertex_vec[av_contacts[iav]].backward.size() == 0) {
                    if (atom_vertex_vec[av_contacts[iav]].connection_type != atom_vertex::START)
                        atom_vertex_vec[av_contacts[iav]].backward.push_back(this_base_atom);
                }
            }
            q.pop();
            done.push_back(this_base_atom);
            r = 1; // return success
        }

        return r;
    }
};

} // namespace coot