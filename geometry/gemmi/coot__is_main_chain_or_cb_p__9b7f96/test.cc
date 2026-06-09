#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Minimal CID-to-Atom lookup for testing: parses "//chain/resid/atom_name"
static gemmi::Atom* cid_to_atom(const std::string& cid, gemmi::Structure& st) {
    if (cid.size() < 2 || cid[0] != '/' || cid[1] != '/') return nullptr;
    std::string rest = cid.substr(2);
    size_t p1 = rest.find('/');
    if (p1 == std::string::npos) return nullptr;
    std::string chain_id = rest.substr(0, p1);
    std::string rest2 = rest.substr(p1 + 1);
    size_t p2 = rest2.find('/');
    std::string res_str, atom_name;
    if (p2 == std::string::npos) { res_str = rest2; atom_name = ""; }
    else { res_str = rest2.substr(0, p2); atom_name = rest2.substr(p2 + 1); }
    int res_num = std::stoi(res_str);
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value != res_num) continue;
                if (atom_name.empty()) {
                    if (!res.atoms.empty()) return &res.atoms[0];
                } else {
                    for (auto& atom : res.atoms) {
                        std::string trimmed = atom.name;
                        auto s = trimmed.find_first_not_of(' ');
                        if (s != std::string::npos) trimmed = trimmed.substr(s);
                        while (!trimmed.empty() && trimmed.back() == ' ') trimmed.pop_back();
                        if (trimmed == atom_name) return &atom;
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, is_main_chain_or_cb_p) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: CA atom (main chain) — should be true
    {
        std::string cid = "//A/1/CA";
        gemmi::Atom *atom = cid_to_atom(cid, st);
        ASSERT_NE(atom, nullptr);
        bool result = coot::is_main_chain_or_cb_p_gemmi(*atom);
        EXPECT_TRUE(result);
    }

    // Case 2: CB atom (allowed as special case) — should be true
    {
        std::string cid = "//A/1/CB";
        gemmi::Atom *atom = cid_to_atom(cid, st);
        ASSERT_NE(atom, nullptr);
        bool result = coot::is_main_chain_or_cb_p_gemmi(*atom);
        EXPECT_TRUE(result);
    }

    // Case 3: CE atom (side chain, MET) — should be false
    {
        std::string cid = "//A/1/CE";
        gemmi::Atom *atom = cid_to_atom(cid, st);
        ASSERT_NE(atom, nullptr);
        bool result = coot::is_main_chain_or_cb_p_gemmi(*atom);
        EXPECT_FALSE(result);
    }

    // Case 4: test the string-based overload directly
    {
        std::string atom_name = "CD1";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
