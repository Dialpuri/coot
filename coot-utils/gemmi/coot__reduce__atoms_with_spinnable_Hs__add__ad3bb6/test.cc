#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <algorithm>
#include "function.hh"

static std::string strip_spaces(std::string s) {
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return s;
}

static gemmi::Atom* find_atom_by_name(gemmi::Residue& res, const char* name) {
    std::string target = strip_spaces(name);
    for (gemmi::Atom& a : res.atoms) {
        if (strip_spaces(a.name) == target) return &a;
    }
    return nullptr;
}

TEST(OracleTest, atoms_with_spinnable_Hs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    auto count_entries = [](coot::reduce::atoms_with_spinnable_Hs& s) -> size_t {
        size_t total = 0;
        for (auto const& kv : s.typed_atoms) {
            total += kv.second.size();
        }
        return total;
    };

    // Case 1: HYDROXYL type — SER 46 OG with HG hydrogen
    {
        coot::reduce::atoms_with_spinnable_Hs spinables;

        gemmi::Atom* og_atom = nullptr;
        gemmi::Atom* hg_atom = nullptr;

        for (gemmi::Residue& r : chain.residues) {
            if (r.name == "SER" && r.seqid.num.value == 46) {
                og_atom = find_atom_by_name(r, "OG");
                hg_atom = find_atom_by_name(r, "HG");
                break;
            }
        }

        EXPECT_NE(og_atom, nullptr);
        EXPECT_NE(hg_atom, nullptr);

        size_t size_before = count_entries(spinables);
        EXPECT_EQ(size_before, 0);

        if (og_atom && hg_atom) {
            spinables.add_gemmi(og_atom, coot::reduce::atom_with_attached_Hs::hydrogen_t::HYDROXYL, hg_atom);
        }

        size_t size_after = count_entries(spinables);
        EXPECT_EQ(size_after, 1);
    }

    // Case 2: METHYL type — find residue with CB and HB atoms
    {
        coot::reduce::atoms_with_spinnable_Hs spinables;

        gemmi::Atom* cb_atom = nullptr;
        std::vector<gemmi::Atom*> h_atoms;

        for (gemmi::Residue& r : chain.residues) {
            cb_atom = find_atom_by_name(r, "CB");
            if (cb_atom) {
                gemmi::Atom* hb1 = find_atom_by_name(r, "HB1");
                gemmi::Atom* hb2 = find_atom_by_name(r, "HB2");
                gemmi::Atom* hb3 = find_atom_by_name(r, "HB3");
                if (hb1) h_atoms.push_back(hb1);
                if (hb2) h_atoms.push_back(hb2);
                if (hb3) h_atoms.push_back(hb3);
                break;
            }
        }

        EXPECT_NE(cb_atom, nullptr);
        EXPECT_EQ(h_atoms.size(), 2);

        size_t size_before = count_entries(spinables);
        EXPECT_EQ(size_before, 0);

        if (cb_atom && !h_atoms.empty()) {
            spinables.add_gemmi(cb_atom, coot::reduce::atom_with_attached_Hs::hydrogen_t::METHYL, h_atoms);
        }

        size_t size_after = count_entries(spinables);
        EXPECT_EQ(size_after, 1);
    }

    // Case 3: SULFHYDRYL — CYS SG with HG
    {
        coot::reduce::atoms_with_spinnable_Hs spinables;

        gemmi::Atom* sg_atom = nullptr;
        gemmi::Atom* hg_atom = nullptr;

        for (gemmi::Residue& r : chain.residues) {
            if (r.name == "CYS") {
                sg_atom = find_atom_by_name(r, "SG");
                hg_atom = find_atom_by_name(r, "HG");
                if (sg_atom) break;
            }
        }

        EXPECT_NE(sg_atom, nullptr);
        EXPECT_NE(hg_atom, nullptr);

        size_t size_before = count_entries(spinables);
        EXPECT_EQ(size_before, 0);

        if (sg_atom && hg_atom) {
            spinables.add_gemmi(sg_atom, coot::reduce::atom_with_attached_Hs::hydrogen_t::SULFHYDRYL, hg_atom);
        }

        size_t size_after = count_entries(spinables);
        EXPECT_EQ(size_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}