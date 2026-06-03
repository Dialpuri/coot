#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <cmath>
#include "function.hh"

// Pad a gemmi atom name back to PDB 4-char format
static std::string pad4(const std::string& name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    if (name.size() >= 3) return name.substr(0, 3) + " ";
    return "    ";
}

TEST(OracleTest, mark_neighbours) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1u);

    // Collect all atoms into flat vectors (file order)
    std::vector<gemmi::Atom*> atom_ptrs;
    std::vector<std::string> atom_names;
    std::vector<std::string> res_names;

    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                atom_ptrs.push_back(&atom);
                atom_names.push_back(pad4(atom.name));
                res_names.push_back(res.name);
            }
        }
    }

    int n_atoms = static_cast<int>(atom_names.size());
    EXPECT_EQ(n_atoms, 2107);
    ASSERT_GT(n_atoms, 0);

    // Build directed contacts list (each pair stored as both (i,j) and (j,i))
    // This matches MMDB SeekContacts which returns directed pairs
    struct Contact { int id1; int id2; };
    std::vector<Contact> contacts;
    double max_dist_sq = 4.0 * 4.0;
    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double dx = atom_ptrs[i]->pos.x - atom_ptrs[j]->pos.x;
            double dy = atom_ptrs[i]->pos.y - atom_ptrs[j]->pos.y;
            double dz = atom_ptrs[i]->pos.z - atom_ptrs[j]->pos.z;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 <= max_dist_sq) {
                contacts.push_back({i, j});
                contacts.push_back({j, i});
            }
        }
    }

    int ncontacts = static_cast<int>(contacts.size());
    EXPECT_EQ(ncontacts, 24846);
    ASSERT_GT(ncontacts, 0);

    // Build neighbours adjacency list from contacts (bidirectional push from each directed contact)
    // This matches: neighbours[contact[ic].id1].push_back(contact[ic].id2);
    //              neighbours[contact[ic].id2].push_back(contact[ic].id1);
    std::vector<std::vector<int>> neighbours(n_atoms);
    for (int ic = 0; ic < ncontacts; ic++) {
        neighbours[contacts[ic].id1].push_back(contacts[ic].id2);
        neighbours[contacts[ic].id2].push_back(contacts[ic].id1);
    }

    // Initialize group array (replaces MMDB UDData)
    std::vector<int> group(n_atoms, -1);

    // Pick first atom
    int iatom_start = 0;
    std::string atom_name = atom_names[iatom_start];

    EXPECT_EQ(iatom_start, 0);
    EXPECT_EQ(atom_name, " N  ");
    EXPECT_EQ(neighbours[iatom_start].size(), 12u);

    // Count marked before
    int marked_before = 0;
    for (int i = 0; i < n_atoms; i++) {
        if (group[i] != -1) marked_before++;
    }
    EXPECT_EQ(marked_before, 0);

    // Call mark_neighbours_gemmi
    coot::high_res::mark_neighbours_gemmi(iatom_start, 0, atom_name, neighbours, atom_names, group);

    // Count marked after
    int marked_after = 0;
    std::vector<int> marked_indices;
    for (int i = 0; i < n_atoms; i++) {
        if (group[i] != -1) {
            marked_after++;
            if (marked_indices.size() < 5) marked_indices.push_back(i);
        }
    }
    EXPECT_EQ(marked_after, 36);

    // Verify first 5 marked atoms
    EXPECT_EQ(marked_indices.size(), 5u);

    {
        int idx = marked_indices[0];
        EXPECT_EQ(idx, 0);
        EXPECT_EQ(atom_names[idx], " N  ");
        EXPECT_EQ(res_names[idx], "MET");
    }
    {
        int idx = marked_indices[1];
        EXPECT_EQ(idx, 8);
        EXPECT_EQ(atom_names[idx], " N  ");
        EXPECT_EQ(res_names[idx], "GLU");
    }
    {
        int idx = marked_indices[2];
        EXPECT_EQ(idx, 17);
        EXPECT_EQ(atom_names[idx], " N  ");
        EXPECT_EQ(res_names[idx], "ASN");
    }
    {
        int idx = marked_indices[3];
        EXPECT_EQ(idx, 25);
        EXPECT_EQ(atom_names[idx], " N  ");
        EXPECT_EQ(res_names[idx], "PHE");
    }
    {
        int idx = marked_indices[4];
        EXPECT_EQ(idx, 36);
        EXPECT_EQ(atom_names[idx], " N  ");
        EXPECT_EQ(res_names[idx], "GLN");
    }

    // Case 2: re-call with same starting atom — should be no-op
    {
        int igroup2 = 99;
        int marked_before_2 = marked_after;
        coot::high_res::mark_neighbours_gemmi(iatom_start, igroup2, atom_name, neighbours, atom_names, group);

        int marked_after_2 = 0;
        for (int i = 0; i < n_atoms; i++) {
            if (group[i] != -1) marked_after_2++;
        }
        EXPECT_EQ(igroup2, 99);
        EXPECT_EQ(marked_before_2, 36);
        EXPECT_EQ(marked_after_2, 36);
    }

    // Case 3: start from an unmarked atom with a different atom name (CA)
    {
        int ca_start = -1;
        for (int i = 0; i < n_atoms; i++) {
            if (atom_names[i] == " CA ") { ca_start = i; break; }
        }
        if (ca_start >= 0) {
            int marked_before_3 = 0;
            for (int i = 0; i < n_atoms; i++) {
                if (group[i] != -1) marked_before_3++;
            }

            std::string ca_name = " CA ";

            EXPECT_EQ(marked_before_3, 36);

            coot::high_res::mark_neighbours_gemmi(ca_start, 5, ca_name, neighbours, atom_names, group);

            int marked_after_3 = 0;
            for (int i = 0; i < n_atoms; i++) {
                if (group[i] != -1) marked_after_3++;
            }
            EXPECT_GT(marked_after_3, 36);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}