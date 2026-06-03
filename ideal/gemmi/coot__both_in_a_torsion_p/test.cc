#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find an atom in the structure by CID-style string "//chain/res_seq/atom_name"
static gemmi::Atom *find_atom_by_cid(const std::string &cid, gemmi::Structure &st) {
    // Parse "//A/1/N" -> chain_id='A', res_seq=1, atom_name="N"
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/')
        return nullptr;
    size_t pos = 2;
    size_t slash1 = cid.find('/', pos);
    if (slash1 == std::string::npos) return nullptr;
    std::string chain_id = cid.substr(pos, slash1 - pos);

    pos = slash1 + 1;
    size_t slash2 = cid.find('/', pos);
    if (slash2 == std::string::npos) return nullptr;
    std::string res_seq_str = cid.substr(pos, slash2 - pos);
    int res_seq = std::stoi(res_seq_str);

    std::string atom_name = cid.substr(slash2 + 1);

    gemmi::Model &model = st.models[0];
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue &res : chain.residues) {
            if (res.seqid.num.value != res_seq) continue;
            for (gemmi::Atom &atom : res.atoms) {
                if (atom.name == atom_name) {
                    return &atom;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, both_in_a_torsion_p) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // --- Case 1: Two atoms that ARE in the same torsion quad ---
    {
        gemmi::Atom *n1  = find_atom_by_cid("//A/1/N", st);
        gemmi::Atom *ca1 = find_atom_by_cid("//A/1/CA", st);
        gemmi::Atom *c1  = find_atom_by_cid("//A/1/C", st);
        gemmi::Atom *n2  = find_atom_by_cid("//A/2/N", st);

        ASSERT_NE(n1, nullptr);
        ASSERT_NE(ca1, nullptr);
        ASSERT_NE(c1, nullptr);
        ASSERT_NE(n2, nullptr);

        coot::torsion_atom_quad_gemmi quad(n1, ca1, c1, n2, 0.0, 0.0, 1);
        std::vector<coot::torsion_atom_quad_gemmi> quads;
        quads.push_back(quad);

        gemmi::Atom *at_1 = n1;
        gemmi::Atom *at_2 = ca1;

        bool result = coot::both_in_a_torsion_p_gemmi(at_1, at_2, quads);

        EXPECT_EQ(quads.size(), 1u);
        EXPECT_TRUE(result);
    }

    // --- Case 2: Two atoms that are NOT together in any torsion quad ---
    {
        gemmi::Atom *n1  = find_atom_by_cid("//A/1/N", st);
        gemmi::Atom *ca1 = find_atom_by_cid("//A/1/CA", st);
        gemmi::Atom *c1  = find_atom_by_cid("//A/1/C", st);
        gemmi::Atom *n2  = find_atom_by_cid("//A/2/N", st);

        coot::torsion_atom_quad_gemmi quad(n1, ca1, c1, n2, 0.0, 0.0, 1);
        std::vector<coot::torsion_atom_quad_gemmi> quads;
        quads.push_back(quad);

        gemmi::Atom *cb1 = find_atom_by_cid("//A/1/CB", st);
        gemmi::Atom *ca3 = find_atom_by_cid("//A/3/CA", st);

        ASSERT_NE(cb1, nullptr);
        ASSERT_NE(ca3, nullptr);

        bool result = coot::both_in_a_torsion_p_gemmi(cb1, ca3, quads);

        EXPECT_EQ(quads.size(), 1u);
        EXPECT_FALSE(result);
    }

    // --- Case 3: Empty quads vector ---
    {
        gemmi::Atom *at_1 = find_atom_by_cid("//A/1/N", st);
        gemmi::Atom *at_2 = find_atom_by_cid("//A/1/CA", st);

        std::vector<coot::torsion_atom_quad_gemmi> quads; // empty

        bool result = coot::both_in_a_torsion_p_gemmi(at_1, at_2, quads);

        EXPECT_EQ(quads.size(), 0u);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}