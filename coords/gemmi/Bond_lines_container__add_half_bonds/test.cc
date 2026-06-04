#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::CRA find_atom_by_cid(const gemmi::Structure &st, const std::string &cid) {
    // CID format: //chain_id/seq_num/atom_name
    // e.g. //A/1/N
    std::string s = cid;
    if (s.size() >= 2 && s.substr(0, 2) == "//") {
        s = s.substr(2);
    }
    
    // Find first /
    auto p1 = s.find('/');
    if (p1 == std::string::npos) return gemmi::CRA{nullptr, nullptr, nullptr};
    std::string chain_id = s.substr(0, p1);
    
    // Find second /
    auto p2 = s.find('/', p1 + 1);
    if (p2 == std::string::npos) return gemmi::CRA{nullptr, nullptr, nullptr};
    int seq_num = std::stoi(s.substr(p1 + 1, p2 - p1 - 1));
    std::string atom_name = s.substr(p2 + 1);

    for (const gemmi::Model &model : st.models) {
        for (const gemmi::Chain &chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (const gemmi::Residue &res : chain.residues) {
                if (res.seqid.num.value != seq_num) continue;
                for (const gemmi::Atom &atom : res.atoms) {
                    if (atom.name == atom_name) {
                        return gemmi::CRA{const_cast<gemmi::Chain*>(&chain),
                                          const_cast<gemmi::Residue*>(&res),
                                          const_cast<gemmi::Atom*>(&atom)};
                    }
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

TEST(OracleTest, add_half_bonds) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get two atoms from MET 1 in chain A: N and CA
    gemmi::CRA cra_n = find_atom_by_cid(st, "//A/1/N");
    gemmi::CRA cra_ca = find_atom_by_cid(st, "//A/1/CA");
    ASSERT_NE(cra_n.atom, nullptr);
    ASSERT_NE(cra_ca.atom, nullptr);

    // Build Bond_lines_container_gemmi (inherits from Bond_lines_container, exposes protected members)
    Bond_lines_container_gemmi blc;

    // Positions from the atoms
    coot::Cartesian pos_n(cra_n.atom->pos.x, cra_n.atom->pos.y, cra_n.atom->pos.z);
    coot::Cartesian pos_ca(cra_ca.atom->pos.x, cra_ca.atom->pos.y, cra_ca.atom->pos.z);

    // Record bond count before
    int total_bonds_before = blc.total_bond_count();

    // --- Case 1: add_half_bonds for N--CA bond (different elements: N vs C) ---
    {
        int model_number = 1;
        int atom_index_1 = 0;
        int atom_index_2 = 1;
        int atom_colour_type = 0;  // default colour type
        int udd_handle = -1;       // no UDD handle
        my_atom_colour_map_t *atom_colour_map_p = nullptr;
        bool add_begin_end_cap = false;
        bool add_end_end_cap = false;
        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

        add_half_bonds_gemmi(blc, pos_n, pos_ca, cra_n, cra_ca, cc,
                             model_number, atom_index_1, atom_index_2,
                             atom_colour_type, udd_handle, atom_colour_map_p,
                             add_begin_end_cap, add_end_end_cap);

        int total_bonds_after = blc.total_bond_count();

        EXPECT_EQ(total_bonds_before, 0);
        EXPECT_EQ(total_bonds_after, 2);
        EXPECT_EQ((total_bonds_after - total_bonds_before), 2);
    }

    // --- Case 2: add_half_bonds with end caps for C--O bond ---
    {
        gemmi::CRA cra_c = find_atom_by_cid(st, "//A/1/C");
        gemmi::CRA cra_o = find_atom_by_cid(st, "//A/1/O");

        ASSERT_NE(cra_c.atom, nullptr);
        ASSERT_NE(cra_o.atom, nullptr);

        coot::Cartesian pos_c(cra_c.atom->pos.x, cra_c.atom->pos.y, cra_c.atom->pos.z);
        coot::Cartesian pos_o(cra_o.atom->pos.x, cra_o.atom->pos.y, cra_o.atom->pos.z);

        int total_bonds_before2 = blc.total_bond_count();

        int model_number = 1;
        int atom_index_1 = 2;
        int atom_index_2 = 3;
        int atom_colour_type = 0;
        int udd_handle = -1;
        my_atom_colour_map_t *atom_colour_map_p = nullptr;
        bool add_begin_end_cap = true;
        bool add_end_end_cap = true;
        graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;

        add_half_bonds_gemmi(blc, pos_c, pos_o, cra_c, cra_o, cc,
                             model_number, atom_index_1, atom_index_2,
                             atom_colour_type, udd_handle, atom_colour_map_p,
                             add_begin_end_cap, add_end_end_cap);

        int total_bonds_after2 = blc.total_bond_count();

        EXPECT_EQ(total_bonds_before2, 2);
        EXPECT_EQ(total_bonds_after2, 4);
        EXPECT_EQ((total_bonds_after2 - total_bonds_before2), 2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}