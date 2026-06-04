#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

// Helper: find atom pointer and its CRA from a CID like "//A/1/CA"
static std::pair<gemmi::Atom*, gemmi::CRA> cid_to_cra(const std::string &cid, gemmi::Structure &st) {
    std::string rest = cid;
    if (rest.size() >= 2 && rest[0] == '/' && rest[1] == '/') {
        rest = rest.substr(2);
    }
    size_t pos1 = rest.find('/');
    if (pos1 == std::string::npos) return {nullptr, gemmi::CRA{}};
    std::string chain_id = rest.substr(0, pos1);
    std::string rest2 = rest.substr(pos1 + 1);
    size_t pos2 = rest2.find('/');
    std::string res_seq_num_str, atom_name;
    if (pos2 == std::string::npos) {
        res_seq_num_str = rest2;
        atom_name = "";
    } else {
        res_seq_num_str = rest2.substr(0, pos2);
        atom_name = rest2.substr(pos2 + 1);
    }
    int res_seq_num = 0;
    try { res_seq_num = std::stoi(res_seq_num_str); } catch (...) { return {nullptr, gemmi::CRA{}}; }

    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (auto &residue : chain.residues) {
                if (residue.seqid.num.value != res_seq_num) continue;
                if (atom_name.empty()) {
                    if (!residue.atoms.empty()) {
                        return {&residue.atoms[0], gemmi::CRA{&chain, &residue, &residue.atoms[0]}};
                    }
                } else {
                    for (auto &atom : residue.atoms) {
                        std::string trimmed = atom.name;
                        auto start = trimmed.find_first_not_of(' ');
                        if (start != std::string::npos) trimmed = trimmed.substr(start);
                        while (!trimmed.empty() && trimmed.back() == ' ') trimmed.pop_back();
                        if (trimmed == atom_name) {
                            return {&atom, gemmi::CRA{&chain, &residue, &atom}};
                        }
                    }
                }
            }
        }
    }
    return {nullptr, gemmi::CRA{}};
}

TEST(OracleTest, PAtom_stream_operator) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: valid atom — pick a known atom from chain A
    {
        std::string cid = "//A/1/CA";
        auto [atom, cra] = cid_to_cra(cid, st);

        EXPECT_NE(atom, nullptr);

        std::ostringstream oss;
        stream_cra_gemmi(oss, cra, st.models[0].num);
        EXPECT_EQ(oss.str(),
            "1/A/1 {MET}/ CA  altLoc :: segid :: pos: (13.401,3.78,74.703) B-factor: 69.04");
    }

    // case 2: invalid CID (nonexistent residue)
    {
        std::string cid = "//A/9999/CA";
        auto [atom, cra] = cid_to_cra(cid, st);

        EXPECT_EQ(atom, nullptr);

        std::ostringstream oss;
        stream_cra_gemmi(oss, cra, st.models[0].num);
        EXPECT_EQ(oss.str(), "NULL");
    }

    // case 3: atom from chain B (water molecule)
    {
        std::string cid = "//B/1/O";
        auto [atom, cra] = cid_to_cra(cid, st);

        EXPECT_NE(atom, nullptr);

        std::ostringstream oss;
        stream_cra_gemmi(oss, cra, st.models[0].num);
        EXPECT_EQ(oss.str(),
            "1/B/1 {HOH}/ O   altLoc :: segid :: pos: (30.744,8.475,36.748) B-factor: 33.29");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}