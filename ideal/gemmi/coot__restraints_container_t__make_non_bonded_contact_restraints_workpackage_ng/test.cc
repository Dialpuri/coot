#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <set>
#include <map>
#include <atomic>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "function.hh"

TEST(NbcTest, make_non_bonded_contact_restraints_workpackage_ng_full_range) {
    // Case 1: Full range
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // Flatten atoms from first chain
    std::vector<gemmi::Atom*> atoms;
    std::vector<gemmi::Residue*> atom_residues;
    std::vector<gemmi::Chain*> atom_chains;
    std::vector<int> atom_residue_indices;
    std::vector<bool> atom_is_hydrogen;
    std::vector<bool> atom_is_metal;

    size_t atom_idx = 0;
    for (size_t ci = 0; ci < model.chains.size(); ci++) {
        gemmi::Chain& chain = model.chains[ci];
        for (size_t ri = 0; ri < chain.residues.size(); ri++) {
            gemmi::Residue& res = chain.residues[ri];
            for (gemmi::Atom& atom : res.atoms) {
                atoms.push_back(&atom);
                atom_residues.push_back(&res);
                atom_chains.push_back(&chain);
                atom_residue_indices.push_back(static_cast<int>(ci * chain.residues.size() + ri));
                atom_is_hydrogen.push_back(atom.element.is_hydrogen());
                atom_is_metal.push_back(atom.element.is_metal());
                ++atom_idx;
            }
        }
    }

    int n_atoms = static_cast<int>(atoms.size());
    ASSERT_GT(n_atoms, 0);

    // Build vcontacts using gemmi NeighborSearch
    double max_radius = 8.0;
    gemmi::NeighborSearch ns(model, st.cell, max_radius);
    ns.populate(true);  // include hydrogen atoms

    std::vector<std::set<unsigned int> > vcontacts(n_atoms);
    for (int i = 0; i < n_atoms; i++) {
        if (!atoms[i]) continue;
        if (atoms[i]->name.empty()) continue;

        std::vector<gemmi::NeighborSearch::Mark*> hits =
            ns.find_atoms(atoms[i]->pos, '\0', 0.0, max_radius);
        for (auto* m : hits) {
            gemmi::CRA cra = m->to_cra(model);
            if (cra.atom) {
                for (int k = 0; k < n_atoms; k++) {
                    if (atoms[k] == cra.atom) {
                        vcontacts[i].insert(static_cast<unsigned int>(k));
                        break;
                    }
                }
            }
        }
    }

    // Empty bonded_atom_indices (no bonds to skip)
    std::vector<std::set<int> > bonded_atom_indices(n_atoms);

    // angle_atom_indices: -1 for all atoms (no angle restraints means is_1_4 always true)
    std::vector<int> angle_atom_indices(n_atoms, -1);

    // Empty fixed_atom_indices
    std::set<int> fixed_atom_indices;

    // Energy types - use placeholder
    std::vector<std::string> energy_type_for_atom(n_atoms, "C  ");

    // H_atom_parent_atom_is_donor_vec and atom_is_acceptor_vec
    std::vector<bool> H_atom_parent_atom_is_donor_vec(n_atoms, false);
    std::vector<bool> atom_is_acceptor_vec(n_atoms, false);

    // Output
    std::vector<std::set<int> > non_bonded_contacts_atom_indices(n_atoms);
    std::vector<coot::simple_restraint> nbc_restraints_fragment;
    std::atomic<unsigned int> done_count(0);

    std::pair<unsigned int, unsigned int> range(0, static_cast<unsigned int>(n_atoms));

    EXPECT_EQ(nbc_restraints_fragment.size(), 0u);

    // Call the workpackage function
    coot::make_non_bonded_contact_restraints_workpackage_ng_gemmi(
        0,
        0,  // imol
        coot::protein_geometry(),  // default geometry
        bonded_atom_indices,
        angle_atom_indices,
        vcontacts,
        range,
        fixed_atom_indices,
        energy_type_for_atom,
        false,  // extended_atom_mode
        atoms,
        atom_residues,
        atom_chains,
        atom_residue_indices,
        atom_is_metal,
        atom_is_hydrogen,
        H_atom_parent_atom_is_donor_vec,
        atom_is_acceptor_vec,
        &non_bonded_contacts_atom_indices,
        &nbc_restraints_fragment,
        done_count
    );

    EXPECT_GT(nbc_restraints_fragment.size(), 0u);
    EXPECT_EQ(done_count.load(), 1u);
}

TEST(NbcTest, make_non_bonded_contact_restraints_workpackage_ng_partial_range) {
    // Case 2: Restricted range
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    std::vector<gemmi::Atom*> atoms;
    std::vector<gemmi::Residue*> atom_residues;
    std::vector<gemmi::Chain*> atom_chains;
    std::vector<int> atom_residue_indices;
    std::vector<bool> atom_is_hydrogen;
    std::vector<bool> atom_is_metal;

    size_t atom_idx = 0;
    for (size_t ci = 0; ci < model.chains.size(); ci++) {
        gemmi::Chain& chain = model.chains[ci];
        for (size_t ri = 0; ri < chain.residues.size(); ri++) {
            gemmi::Residue& res = chain.residues[ri];
            for (gemmi::Atom& atom : res.atoms) {
                atoms.push_back(&atom);
                atom_residues.push_back(&res);
                atom_chains.push_back(&chain);
                atom_residue_indices.push_back(static_cast<int>(ci * chain.residues.size() + ri));
                atom_is_hydrogen.push_back(atom.element.is_hydrogen());
                atom_is_metal.push_back(atom.element.is_metal());
                ++atom_idx;
            }
        }
    }

    int n_atoms = static_cast<int>(atoms.size());
    ASSERT_GT(n_atoms, 0);

    // Build vcontacts
    double max_radius = 8.0;
    gemmi::NeighborSearch ns(model, st.cell, max_radius);
    ns.populate(true);

    std::vector<std::set<unsigned int> > vcontacts(n_atoms);
    for (int i = 0; i < n_atoms; i++) {
        if (!atoms[i]) continue;
        if (atoms[i]->name.empty()) continue;

        std::vector<gemmi::NeighborSearch::Mark*> hits =
            ns.find_atoms(atoms[i]->pos, '\0', 0.0, max_radius);
        for (auto* m : hits) {
            gemmi::CRA cra = m->to_cra(model);
            if (cra.atom) {
                for (int k = 0; k < n_atoms; k++) {
                    if (atoms[k] == cra.atom) {
                        vcontacts[i].insert(static_cast<unsigned int>(k));
                        break;
                    }
                }
            }
        }
    }

    std::vector<std::set<int> > bonded_atom_indices(n_atoms);
    std::vector<int> angle_atom_indices(n_atoms, -1);
    std::set<int> fixed_atom_indices;
    std::vector<std::string> energy_type_for_atom(n_atoms, "C  ");
    std::vector<bool> H_atom_parent_atom_is_donor_vec(n_atoms, false);
    std::vector<bool> atom_is_acceptor_vec(n_atoms, false);

    std::vector<std::set<int> > non_bonded_contacts_atom_indices(n_atoms);
    std::vector<coot::simple_restraint> nbc_restraints_fragment;
    std::atomic<unsigned int> done_count(0);

    // Restricted range: only first half
    unsigned int half = static_cast<unsigned int>(n_atoms / 2);
    std::pair<unsigned int, unsigned int> range(0, half);

    EXPECT_EQ(nbc_restraints_fragment.size(), 0u);

    EXPECT_NO_THROW(
        coot::make_non_bonded_contact_restraints_workpackage_ng_gemmi(
            0,
            0,
            coot::protein_geometry(),
            bonded_atom_indices,
            angle_atom_indices,
            vcontacts,
            range,
            fixed_atom_indices,
            energy_type_for_atom,
            false,
            atoms,
            atom_residues,
            atom_chains,
            atom_residue_indices,
            atom_is_metal,
            atom_is_hydrogen,
            H_atom_parent_atom_is_donor_vec,
            atom_is_acceptor_vec,
            &non_bonded_contacts_atom_indices,
            &nbc_restraints_fragment,
            done_count
        )
    );

    EXPECT_GT(nbc_restraints_fragment.size(), 0u);
    EXPECT_EQ(done_count.load(), 1u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}