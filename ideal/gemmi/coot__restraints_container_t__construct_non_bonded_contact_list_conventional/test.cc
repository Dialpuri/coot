#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include "function.hh"

#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <algorithm>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

    // Build bonded_atom_indices from covalent connections in the structure
static std::vector<std::vector<int>> build_bonds_from_structure(
    const gemmi::Structure& st,
    const gemmi::Model& model,
    const std::vector<gemmi::Atom*>& selected_atoms
) {
    int n = static_cast<int>(selected_atoms.size());
    std::vector<std::set<int>> bond_set(n);

    // Build a lookup: (chain_name, residue_seqnum, atom_name) -> index
    // Atoms don't have parent pointers, so iterate model to get chain info
    std::unordered_map<std::string, int> atom_key_to_index;
    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            for (const gemmi::Atom& atom : res.atoms) {
                const gemmi::Atom* ptr = &atom;
                auto it = std::find(selected_atoms.begin(), selected_atoms.end(), ptr);
                if (it != selected_atoms.end()) {
                    int idx = static_cast<int>(it - selected_atoms.begin());
                    std::string key = chain.name + "|" +
                        std::to_string(res.seqid.num.value) + "|" +
                        atom.name;
                    atom_key_to_index[key] = idx;
                }
            }
        }
    }

    // Process covalent connections
    for (const auto& conn : st.connections) {
        if (conn.type != gemmi::Connection::Type::Covale) continue;

        auto make_key = [](const gemmi::AtomAddress& addr) -> std::string {
            return addr.chain_name + "|" +
                std::to_string(addr.res_id.seqid.num.value) + "|" +
                addr.atom_name;
        };

        std::string key1 = make_key(conn.partner1);
        std::string key2 = make_key(conn.partner2);

        auto it1 = atom_key_to_index.find(key1);
        auto it2 = atom_key_to_index.find(key2);

        if (it1 != atom_key_to_index.end() && it2 != atom_key_to_index.end()) {
            int idx1 = it1->second;
            int idx2 = it2->second;
            bond_set[idx1].insert(idx2);
            bond_set[idx2].insert(idx1);
        }
    }

    // Convert to vector<vector<int>>
    std::vector<std::vector<int>> bonded_indices(n);
    for (int i = 0; i < n; i++) {
        for (int idx : bond_set[i]) {
            bonded_indices[i].push_back(idx);
        }
    }
    return bonded_indices;
}

TEST(OracleTest, construct_non_bonded_contact_list_conventional) {
    {
        // Load the PDB with gemmi
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model& model = st.models[0];

        // Collect atoms from residues 10-14 in chain A (the selected range)
        std::vector<gemmi::Atom*> gemmi_atoms;
        std::vector<gemmi::Residue*> selected_residues;
        std::vector<gemmi::Residue*> active_residues; // just GLY(10)

        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (const gemmi::Residue& res : chain.residues) {
                int seq = res.seqid.num.value;
                if (seq >= 10 && seq <= 14) {
                    selected_residues.push_back(const_cast<gemmi::Residue*>(&res));
                    for (const gemmi::Atom& atom : res.atoms) {
                        gemmi_atoms.push_back(const_cast<gemmi::Atom*>(&atom));
                    }
                    // Only GLY (residue 10) is active
                    if (seq == 10) {
                        active_residues.push_back(const_cast<gemmi::Residue*>(&res));
                    }
                }
            }
        }

        ASSERT_FALSE(selected_residues.empty());

        int n_atoms = static_cast<int>(gemmi_atoms.size());

        // Build bonded_atom_indices from covalent connections
        std::vector<std::vector<int>> bonded_indices =
            build_bonds_from_structure(st, model, gemmi_atoms);

        // Create restraints container
        coot::restraints_container_t restraints;
        restraints.n_atoms = n_atoms;
        restraints.gemmi_atoms = gemmi_atoms;
        restraints.bonded_atom_indices = bonded_indices;

        // Active residues organized by chain - only 1 active (GLY 10)
        std::vector<std::vector<gemmi::Residue*>> active_residues_vec;
        active_residues_vec.resize(1);
        active_residues_vec[0] = active_residues;
        restraints.active_residues_vec = active_residues_vec;
        restraints.have_oxt_flag = false;

        // OUTPUT BEFORE
        std::size_t before_size = restraints.filtered_non_bonded_atom_indices.size();
        std::size_t before_total = 0;
        for (const auto &vec : restraints.filtered_non_bonded_atom_indices) {
            before_total += vec.size();
        }
        EXPECT_EQ(before_size, 0u);
        EXPECT_EQ(before_total, 0u);

        // Call the target function
        restraints.construct_non_bonded_contact_list_conventional_gemmi();

        // OUTPUT AFTER - check size first
        std::size_t after_size = restraints.filtered_non_bonded_atom_indices.size();
        std::size_t after_total = 0;
        for (const auto &vec : restraints.filtered_non_bonded_atom_indices) {
            after_total += vec.size();
        }

        // Debug: print what we got
        std::cout << "DEBUG: n_atoms=" << n_atoms
                  << " active_residues=" << active_residues.size()
                  << " after_size=" << after_size
                  << " after_total=" << after_total << std::endl;

        EXPECT_EQ(after_size, static_cast<std::size_t>(n_atoms));
        // Total non-bonded contacts within 8A
        EXPECT_GT(after_total, 0u);
    }
}
