#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: parse CID "//A/1/CA" -> (chain_name, res_seqnum, inscode, atom_name)
static bool parse_cid(const std::string& cid, std::string& chain_name, int& res_seq,
                      char& icode, std::string& atom_name) {
    if (cid.size() < 5 || cid[0] != '/' || cid[1] != '/') return false;
    size_t p = cid.find('/', 2);
    if (p == std::string::npos) return false;
    chain_name = cid.substr(2, p - 2);

    size_t p2 = cid.find('/', p + 1);
    if (p2 == std::string::npos) return false;
    std::string res_part = cid.substr(p + 1, p2 - p - 1);

    if (!res_part.empty() && (res_part.back() >= 'A' && res_part.back() <= 'Z' ||
                              res_part.back() >= 'a' && res_part.back() <= 'z')) {
        icode = res_part.back();
        res_part.pop_back();
    } else {
        icode = ' ';
    }

    try {
        res_seq = std::stoi(res_part);
    } catch (...) {
        return false;
    }

    atom_name = cid.substr(p2 + 1);
    return true;
}

// Helper: find an atom by CID in a gemmi Structure
static gemmi::Atom* find_atom_by_cid(const std::string& cid, gemmi::Structure& st) {
    std::string chain_name, atom_name;
    int res_seq;
    char icode;
    if (!parse_cid(cid, chain_name, res_seq, icode, atom_name)) return nullptr;

    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value != res_seq) continue;
                if (res.seqid.icode != icode && icode != ' ') continue;
                for (auto& atom : res.atoms) {
                    if (atom.name == atom_name) return &atom;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, atom_to_co) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "load failed";

    // Case 1: valid atom
    {
        std::string cid = "//A/1/CA";
        gemmi::Atom* atom = find_atom_by_cid(cid, st);
        ASSERT_NE(atom, nullptr) << "atom not found for cid: " << cid;
        clipper::Coord_orth co = coot::atom_to_co_gemmi(*atom);
        EXPECT_NEAR(co.x(), 13.401, 1e-3);
        EXPECT_NEAR(co.y(), 3.78, 1e-3);
        EXPECT_NEAR(co.z(), 74.703, 1e-3);
    }

    // Case 2: another valid atom (different chain)
    {
        std::string cid = "//B/1/O";
        gemmi::Atom* atom = find_atom_by_cid(cid, st);
        if (!atom) {
            cid = "//B/1/OW";
            atom = find_atom_by_cid(cid, st);
        }
        ASSERT_NE(atom, nullptr) << "atom not found for cid: " << cid;
        clipper::Coord_orth co = coot::atom_to_co_gemmi(*atom);
        EXPECT_NEAR(co.x(), 30.744, 1e-3);
        EXPECT_NEAR(co.y(), 8.475, 1e-3);
        EXPECT_NEAR(co.z(), 36.748, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}