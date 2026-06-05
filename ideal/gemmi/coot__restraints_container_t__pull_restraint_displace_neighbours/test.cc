#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <cmath>
#include <vector>
#include <string>
#include <set>

// Helper: find atom by CID string like "//A/10/CA"
static gemmi::Atom* find_atom_by_cid(gemmi::Model& model, const std::string& cid) {
    size_t pos = 0;
    while (pos < cid.size() && cid[pos] == '/') pos++;
    size_t start = pos;

    size_t slash1 = cid.find('/', start);
    if (slash1 == std::string::npos) return nullptr;
    std::string chain_id = cid.substr(start, slash1 - start);

    start = slash1 + 1;
    size_t slash2 = cid.find('/', start);
    if (slash2 == std::string::npos) return nullptr;
    std::string seqnum_str = cid.substr(start, slash2 - start);
    int seqnum = std::stoi(seqnum_str);

    start = slash2 + 1;
    std::string atom_name = cid.substr(start);

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == atom_name || atom.name.find(atom_name) == 0) {
                            return &atom;
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, pull_restraint_displace_neighbours) {
    {
        // Case 1: Displace neighbours with a large delta and radius
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];

        // Collect all atoms in order (like MMDB atom array)
        std::vector<gemmi::Atom*> atom_ptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    atom_ptr.push_back(&atom);
                }
            }
        }

        // No fixed atoms
        std::set<unsigned int> fixed_atom_indices;

        // Pick a pull atom: CA of residue 10 in chain A
        std::string cid = "//A/10/CA";
        gemmi::Atom* pull_atom = find_atom_by_cid(model, cid);
        EXPECT_NE(pull_atom, nullptr);

        float radius = 5.0f;
        // New target position: move 3 Angstroms in X
        clipper::Coord_orth new_target_pos(pull_atom->pos.x + 3.0, pull_atom->pos.y, pull_atom->pos.z);

        // Record positions of atoms near the pull atom BEFORE
        std::vector<gemmi::Atom*> nearby_atoms;
        float r_squared = radius * radius;
        for (unsigned int i = 0; i < atom_ptr.size(); i++) {
            gemmi::Atom* at = atom_ptr[i];
            float d2 = (at->pos.x - pull_atom->pos.x) * (at->pos.x - pull_atom->pos.x) +
                       (at->pos.y - pull_atom->pos.y) * (at->pos.y - pull_atom->pos.y) +
                       (at->pos.z - pull_atom->pos.z) * (at->pos.z - pull_atom->pos.z);
            if (d2 < r_squared) {
                nearby_atoms.push_back(at);
            }
        }

        EXPECT_EQ(nearby_atoms.size(), 22u);

        // Store original positions
        struct AtomPos { double x, y, z; };
        std::vector<AtomPos> orig_positions;
        for (gemmi::Atom* at : nearby_atoms) {
            orig_positions.push_back(AtomPos{at->pos.x, at->pos.y, at->pos.z});
        }

        // Call the function
        coot::pull_restraint_displace_neighbours_gemmi(
            pull_atom, atom_ptr, fixed_atom_indices, new_target_pos, radius);

        // Check how many atoms actually moved
        int n_moved = 0;
        double total_displacement = 0.0;
        for (size_t i = 0; i < nearby_atoms.size(); i++) {
            gemmi::Atom* at = nearby_atoms[i];
            double dx = at->pos.x - orig_positions[i].x;
            double dy = at->pos.y - orig_positions[i].y;
            double dz = at->pos.z - orig_positions[i].z;
            double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d > 1e-6) {
                n_moved++;
                total_displacement += d;
            }
        }

        EXPECT_EQ(n_moved, 18);
        EXPECT_NEAR(total_displacement, 14.5503, 1e-4);

        // Check first 3 moved atoms' names and displacements
        int printed = 0;
        std::vector<std::pair<std::string, double>> moved_atom_info;
        for (size_t i = 0; i < nearby_atoms.size() && printed < 3; i++) {
            gemmi::Atom* at = nearby_atoms[i];
            double dx = at->pos.x - orig_positions[i].x;
            double dy = at->pos.y - orig_positions[i].y;
            double dz = at->pos.z - orig_positions[i].z;
            double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d > 1e-6) {
                // Pad atom name to PDB-style 4-char field (MMDB GetAtomName format)
                // PDB convention: single-char element names get a leading space: " C  "
                std::string padded_name = at->name;
                if (padded_name.size() == 1) {
                    padded_name = " " + padded_name + "  ";
                } else {
                    while (padded_name.size() < 4) padded_name += " ";
                }
                moved_atom_info.push_back({padded_name, d});
                printed++;
            }
        }

        EXPECT_EQ(moved_atom_info.size(), 3u);
        EXPECT_EQ(moved_atom_info[0].first, " C  ");
        EXPECT_NEAR(moved_atom_info[0].second, 0.197413, 1e-4);
        EXPECT_EQ(moved_atom_info[1].first, " O  ");
        EXPECT_NEAR(moved_atom_info[1].second, 0.422013, 1e-4);
        EXPECT_EQ(moved_atom_info[2].first, " N  ");
        EXPECT_NEAR(moved_atom_info[2].second, 0.393878, 1e-4);
    }

    {
        // Case 2: Zero radius — no atoms should move
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model& model = st.models[0];

        std::vector<gemmi::Atom*> atom_ptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    atom_ptr.push_back(&atom);
                }
            }
        }

        std::set<unsigned int> fixed_atom_indices;

        std::string cid2 = "//A/50/CA";
        gemmi::Atom* pull_atom2 = find_atom_by_cid(model, cid2);
        EXPECT_NE(pull_atom2, nullptr);

        clipper::Coord_orth target2(pull_atom2->pos.x + 3.0, pull_atom2->pos.y, pull_atom2->pos.z);
        float radius2 = 0.0f;

        EXPECT_NO_THROW(coot::pull_restraint_displace_neighbours_gemmi(
            pull_atom2, atom_ptr, fixed_atom_indices, target2, radius2));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}