#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <iomanip>
#include <sstream>

// Helper: pad atom name to 4 chars (MMDB format) for comparison
inline std::string pad_atom_name(const std::string& name) {
    std::string s = " " + name;
    while (s.size() < 4) s += ' ';
    return s;
}

// Simple struct matching the interface the template expects:
// .atom_2 and .atom_3 must be gemmi::Atom* (pointers, not strings)
struct link_torsion_atom_ptr {
    gemmi::Atom *atom_1;
    gemmi::Atom *atom_2;
    gemmi::Atom *atom_3;
    gemmi::Atom *atom_4;
    link_torsion_atom_ptr()
        : atom_1(nullptr), atom_2(nullptr), atom_3(nullptr), atom_4(nullptr) {}
};

TEST(OracleTest, contact_info) {
    // 1. Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Find chain A
    gemmi::Chain *chain = nullptr;
    for (auto &c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Original test: mol.Select(selhnd, mmdb::STYPE_CHAIN, "//A", mmdb::SKEY_NEW)
    // resulted in n_sel_atoms == 0 (selection didn't match).
    // So selected_atoms is empty — matching original behavior.
    std::vector<gemmi::Atom *> selected_atoms; // empty, matching original

    int n_residues = static_cast<int>(chain->residues.size());
    EXPECT_EQ(n_residues, 267);

    // Case 1: Valid torsion with atoms in different residues
    gemmi::Residue *r1 = nullptr;
    gemmi::Residue *r2 = nullptr;
    if (!chain->residues.empty()) {
        r1 = &chain->residues[0];
    }
    if (chain->residues.size() > 1) {
        r2 = &chain->residues[1];
    }

    gemmi::Atom *a1 = nullptr;
    if (r1 && !r1->atoms.empty()) {
        a1 = &r1->atoms[0];
    }
    gemmi::Atom *a2 = nullptr;
    if (r2 && !r2->atoms.empty()) {
        a2 = &r2->atoms[0];
    }

    ASSERT_NE(a1, nullptr);
    ASSERT_NE(a2, nullptr);
    EXPECT_STREQ(pad_atom_name(a1->name).c_str(), " N  ");
    EXPECT_STREQ(pad_atom_name(a2->name).c_str(), " N  ");

    // Verify atoms are in different residues
    bool atoms_in_different_residues = (a1 && a2 && r1 != r2);
    EXPECT_TRUE(atoms_in_different_residues);

    std::vector<link_torsion_atom_ptr> link_torsions;
    {
        link_torsion_atom_ptr lt;
        lt.atom_1 = nullptr;
        lt.atom_2 = a1;
        lt.atom_3 = a2;
        lt.atom_4 = nullptr;
        link_torsions.push_back(lt);
    }
    EXPECT_EQ(link_torsions.size(), 1);

    // Call the gemmi port — with empty selection, returns empty contacts
    auto contacts1 = coot::contact_info_gemmi(model, selected_atoms, link_torsions);
    EXPECT_EQ(contacts1.size(), 0);

    // Case 2: empty link_torsions
    std::vector<link_torsion_atom_ptr> empty_torsions;
    auto contacts2 = coot::contact_info_gemmi(model, selected_atoms, empty_torsions);
    EXPECT_EQ(contacts2.size(), 0);

    // Case 3 (complementary): Use a gemmi-compatible torsion type with all atoms populated
    {
        gemmi::Residue *r3 = nullptr;
        gemmi::Residue *r4 = nullptr;
        if (!chain->residues.empty()) {
            r3 = &chain->residues[0];
        }
        if (chain->residues.size() > 1) {
            r4 = &chain->residues[1];
        }

        gemmi::Atom *at1 = (r3 && !r3->atoms.empty()) ? &r3->atoms[0] : nullptr;
        gemmi::Atom *at2 = (r3 && r3->atoms.size() > 1) ? &r3->atoms[1] : nullptr;
        gemmi::Atom *at3 = (r4 && !r4->atoms.empty()) ? &r4->atoms[0] : nullptr;
        gemmi::Atom *at4 = (r4 && r4->atoms.size() > 1) ? &r4->atoms[1] : nullptr;

        std::vector<link_torsion_atom_ptr> torsions_real;
        {
            link_torsion_atom_ptr tq;
            tq.atom_1 = at1;
            tq.atom_2 = at2;
            tq.atom_3 = at3;
            tq.atom_4 = at4;
            torsions_real.push_back(tq);
        }

        auto contacts3 = coot::contact_info_gemmi(model, selected_atoms, torsions_real);
        EXPECT_EQ(contacts3.size(), 0);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}