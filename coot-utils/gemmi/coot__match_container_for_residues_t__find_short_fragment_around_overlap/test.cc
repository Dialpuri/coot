#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a specific atom by chain, residue seqnum, and atom name in a gemmi model
static gemmi::Atom* find_atom(gemmi::Model& model, const std::string& chain_id, int seqnum, const std::string& atom_name) {
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name != chain_id) continue;
        for (gemmi::Residue& res : ch.residues) {
            if (res.seqid.num.value != seqnum) continue;
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.name == atom_name) return &atom;
            }
        }
    }
    return nullptr;
}

// Helper: find a residue by chain and seqnum
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name != chain_id) continue;
        for (gemmi::Residue& res : ch.residues) {
            if (res.seqid.num.value == seqnum) return &res;
        }
    }
    return nullptr;
}

// Helper: build a selection of atoms (N, CA, C, O) for residues in [first, last] inclusive
static std::vector<gemmi::Atom*> build_selection(gemmi::Model& model, const std::string& chain_id, int first, int last) {
    std::vector<gemmi::Atom*> sel;
    const char* names[] = {"N", "CA", "C", "O"};
    for (int seqnum = first; seqnum <= last; seqnum++) {
        for (const char* aname : names) {
            gemmi::Atom* a = find_atom(model, chain_id, seqnum, aname);
            if (a) sel.push_back(a);
        }
    }
    return sel;
}

TEST(OracleTest, find_short_fragment_around_overlap) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Build selection 1: atoms from residues 10-20 in chain A
    std::vector<gemmi::Atom*> sel1 = build_selection(model, "A", 10, 20);
    ASSERT_FALSE(sel1.empty());

    // Build selection 2: atoms from residues 15-25 in chain A (overlaps 15-20)
    std::vector<gemmi::Atom*> sel2 = build_selection(model, "A", 15, 25);
    ASSERT_FALSE(sel2.empty());

    // --- Case 1: Upstream short fragment in first selection ---
    {
        coot::match_container_for_residues_t mcr;
        mcr.residue_1 = find_residue(model, "A", 14);
        mcr.residue_2 = find_residue(model, "A", 21);
        gemmi::Atom* ca14 = find_atom(model, "A", 14, "CA");
        gemmi::Atom* ca21 = find_atom(model, "A", 21, "CA");
        if (ca14 && ca21) {
            mcr.atom_pairs.push_back(std::make_pair(ca14, ca21));
        }

        coot::delete_a_chain_t dac = mcr.find_short_fragment_around_overlap_gemmi(sel1, sel2);

        EXPECT_TRUE(dac.chains_were_mergeable);
        EXPECT_TRUE(dac.short_fragment_is_in_first_selection);
        EXPECT_TRUE(dac.short_fragment_is_upstream_fragment);
        EXPECT_EQ(static_cast<int>(dac.delete_type), 0);
    }

    // --- Case 2: Downstream short fragment in first selection ---
    {
        coot::match_container_for_residues_t mcr;
        mcr.residue_1 = find_residue(model, "A", 16);
        mcr.residue_2 = find_residue(model, "A", 19);
        gemmi::Atom* ca16 = find_atom(model, "A", 16, "CA");
        gemmi::Atom* ca19 = find_atom(model, "A", 19, "CA");
        if (ca16 && ca19) {
            mcr.atom_pairs.push_back(std::make_pair(ca16, ca19));
        }

        coot::delete_a_chain_t dac = mcr.find_short_fragment_around_overlap_gemmi(sel1, sel2);

        EXPECT_TRUE(dac.chains_were_mergeable);
        EXPECT_FALSE(dac.short_fragment_is_in_first_selection);
        EXPECT_TRUE(dac.short_fragment_is_upstream_fragment);
        EXPECT_EQ(static_cast<int>(dac.delete_type), 0);
    }

    // --- Case 3: Non-mergeable (both sides equally long) ---
    {
        coot::match_container_for_residues_t mcr;
        mcr.residue_1 = find_residue(model, "A", 15);
        mcr.residue_2 = find_residue(model, "A", 20);
        gemmi::Atom* ca15 = find_atom(model, "A", 15, "CA");
        gemmi::Atom* ca20 = find_atom(model, "A", 20, "CA");
        if (ca15 && ca20) {
            mcr.atom_pairs.push_back(std::make_pair(ca15, ca20));
        }

        coot::delete_a_chain_t dac = mcr.find_short_fragment_around_overlap_gemmi(sel1, sel2);

        EXPECT_FALSE(dac.chains_were_mergeable);
        EXPECT_FALSE(dac.short_fragment_is_in_first_selection);
        EXPECT_TRUE(dac.short_fragment_is_upstream_fragment);
        EXPECT_EQ(static_cast<int>(dac.delete_type), 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}