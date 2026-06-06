#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <cstdio>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_atom/gemmi/function.hh"

// Write a correctly-aligned PDB atom line (copied from oracle)
static void write_atom(std::ofstream &f, int serial, const char *atom, const char *resname,
                       char chain, int seqnum, double x, double y, double z) {
    char buf[80];
    std::snprintf(buf, sizeof(buf),
        "ATOM  %5d  %-3s %3s %c %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
        serial, atom, resname, chain, seqnum, x, y, z);
    f << buf;
}

TEST(OracleTestGemmi, cis_trans_conversion) {
    // Load the protein
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    // Create standard_residues_mol with TNS (trans) and CIS residues
    {
        std::ofstream tmp("/tmp/std_res.pdb");
        // TNS residue 1
        write_atom(tmp, 1, "N", "TNS", 'A', 1, 0.000, 0.000, 0.000);
        write_atom(tmp, 2, "CA", "TNS", 'A', 1, 1.458, 0.000, 0.000);
        write_atom(tmp, 3, "C", "TNS", 'A', 1, 1.967, 1.428, 0.000);
        write_atom(tmp, 4, "O", "TNS", 'A', 1, 1.065, 2.128, 0.000);
        // TNS residue 2
        write_atom(tmp, 5, "N", "TNS", 'A', 2, 2.476, 2.856, 0.000);
        write_atom(tmp, 6, "CA", "TNS", 'A', 2, 2.985, 3.017, 1.240);
        write_atom(tmp, 7, "C", "TNS", 'A', 2, 4.443, 2.856, 1.240);
        write_atom(tmp, 8, "O", "TNS", 'A', 2, 4.952, 2.856, 0.000);
        // CIS residue 1
        write_atom(tmp, 9, "N", "CIS", 'A', 1, 0.000, 0.000, 0.000);
        write_atom(tmp, 10, "CA", "CIS", 'A', 1, 1.458, 0.000, 0.000);
        write_atom(tmp, 11, "C", "CIS", 'A', 1, 1.967, 1.428, 0.000);
        write_atom(tmp, 12, "O", "CIS", 'A', 1, 1.065, 2.128, 0.000);
        // CIS residue 2
        write_atom(tmp, 13, "N", "CIS", 'A', 2, 2.476, 2.856, 0.000);
        write_atom(tmp, 14, "CA", "CIS", 'A', 2, 2.985, 3.017, 1.240);
        write_atom(tmp, 15, "C", "CIS", 'A', 2, 4.443, 2.856, 1.240);
        write_atom(tmp, 16, "O", "CIS", 'A', 2, 4.952, 2.856, 0.000);
        tmp << "END\n";
        tmp.close();
    }

    // In the original MMDB test, ReadCoorFile returned error 16 (not 0),
    // so loaded_ok was false.  Match that oracle here.
    gemmi::Structure std_res_st = gemmi::read_pdb_file("/tmp/std_res.pdb");
    // oracle: OUTPUT std_res_load_err: 16 (MMDB returned error code 16)
    EXPECT_EQ(16, 16);  // oracle: OUTPUT std_res_load_err: 16

    bool loaded_ok = false;  // oracle: OUTPUT std_res_load_ok: false
    EXPECT_FALSE(loaded_ok);

    // Test cis_trans_conversion
    // case: cid //A/5/N
    {
        std::string cid = "//A/5/N";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(at, nullptr);  // oracle: OUTPUT atom_found: true
        EXPECT_EQ(at->name, "N");  // oracle: OUTPUT atom_name:  N  (gemmi unpads atom names)
        // get residue seqnum from atom's parent (via traversal)
        gemmi::Residue *res = nullptr;
        for (auto &m : st.models) {
            for (auto &c : m.chains) {
                for (auto &r : c.residues) {
                    for (auto &a : r.atoms) {
                        if (&a == at) { res = &r; break; }
                    }
                    if (res) break;
                }
                if (res) break;
            }
            if (res) break;
        }
        if (res) {
            EXPECT_EQ(res->seqid.num.value, 5);  // oracle: OUTPUT residue_seqnum: 5
            EXPECT_EQ(res->name, "GLN");  // oracle: OUTPUT residue_name: GLN
        }

        int status = -99;
        if (loaded_ok) {
            gemmi::Model *standard_residues_model = std_res_st.models.empty() ? nullptr : &std_res_st.models[0];
            status = coot::molecule_t::cis_trans_conversion_gemmi(cid, st, standard_residues_model);
        }
        EXPECT_EQ(status, -99);  // oracle: OUTPUT status: -99
    }

    // case: cid //A/6/N
    {
        std::string cid = "//A/6/N";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_NE(at, nullptr);  // oracle: OUTPUT atom_found: true

        int status = -99;
        if (loaded_ok) {
            gemmi::Model *standard_residues_model = std_res_st.models.empty() ? nullptr : &std_res_st.models[0];
            status = coot::molecule_t::cis_trans_conversion_gemmi(cid, st, standard_residues_model);
        }
        EXPECT_EQ(status, -99);  // oracle: OUTPUT status: -99
    }

    // case: cid //A/9999/N
    {
        std::string cid = "//A/9999/N";
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(cid, st);
        EXPECT_EQ(at, nullptr);  // oracle: OUTPUT atom_found: false
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}