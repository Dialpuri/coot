#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <string>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Helper: find an atom by CID string like "//A/60/CD2"
static gemmi::CRA find_atom_by_cid(gemmi::Model& model, const std::string& cid) {
    // Parse CID: //chain/resseq/atomname
    std::string c = cid;
    size_t pos = 0;
    while (pos < c.size() && c[pos] == '/') pos++;

    // Extract chain ID
    size_t start = pos;
    while (pos < c.size() && c[pos] != '/') pos++;
    std::string chain_id = c.substr(start, pos - start);

    // Skip slash
    if (pos < c.size()) pos++;

    // Extract residue sequence number
    start = pos;
    while (pos < c.size() && c[pos] != '/') pos++;
    std::string res_seq_str = c.substr(start, pos - start);
    int res_seq = std::stoi(res_seq_str);

    // Skip slash
    if (pos < c.size()) pos++;

    // Extract atom name
    std::string atom_name = c.substr(pos);
    atom_name.erase(atom_name.find_last_not_of(" \t\r\n") + 1);

    // Search the model
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == res_seq) {
                    for (gemmi::Atom& atom : res.atoms) {
                        std::string aname = atom.name;
                        aname.erase(aname.find_last_not_of(" \t\r\n") + 1);
                        if (aname == atom_name) {
                            return gemmi::CRA{&chain, &res, &atom};
                        }
                    }
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

TEST(OracleTest, atom_overlaps_container_t_in_same_ring) {
    // Initialize protein geometry (required before dictionary lookups)
    coot::protein_geometry geom;
    geom.init_standard();

    // Read PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GT(st.models.size(), 0u);

    gemmi::Model& model = st.models[0];

    std::map<std::string, std::vector<std::vector<std::string>>> ring_list_map;

    // === Case 1: Two atoms in HIS ring (CD2 and CE1 are both in the imidazole ring) ===
    {
        std::string cid1 = "//A/60/CD2";
        std::string cid2 = "//A/60/CE1";

        gemmi::CRA cra1 = find_atom_by_cid(model, cid1);
        gemmi::CRA cra2 = find_atom_by_cid(model, cid2);

        EXPECT_NE(cra1.atom, nullptr);
        EXPECT_NE(cra2.atom, nullptr);

        bool result1 = false;
        if (cra1.atom && cra2.atom) {
            result1 = coot::in_same_ring_gemmi(&geom, cra1, cra2, ring_list_map);
        }

        EXPECT_TRUE(result1);
    }

    // === Case 2: Two atoms from different residues -> guard clause returns false ===
    {
        std::string cid1 = "//A/60/CD2";
        std::string cid3 = "//A/1/N";

        gemmi::CRA cra1 = find_atom_by_cid(model, cid1);
        gemmi::CRA cra3 = find_atom_by_cid(model, cid3);

        EXPECT_NE(cra1.atom, nullptr);
        EXPECT_NE(cra3.atom, nullptr);

        bool result2 = false;
        if (cra1.atom && cra3.atom) {
            result2 = coot::in_same_ring_gemmi(&geom, cra1, cra3, ring_list_map);
        }

        EXPECT_FALSE(result2);
    }

    // === Case 3: Two ring atoms from PHE (CD1 and CZ are in the benzene ring) ===
    {
        std::string cid4 = "//A/4/CD1";
        std::string cid5 = "//A/4/CZ";

        gemmi::CRA cra4 = find_atom_by_cid(model, cid4);
        gemmi::CRA cra5 = find_atom_by_cid(model, cid5);

        EXPECT_NE(cra4.atom, nullptr);
        EXPECT_NE(cra5.atom, nullptr);

        bool result3 = false;
        if (cra4.atom && cra5.atom) {
            result3 = coot::in_same_ring_gemmi(&geom, cra4, cra5, ring_list_map);
        }

        EXPECT_TRUE(result3);
    }

    // === Case 4: Two atoms in same HIS residue but NOT in same ring (CD2 and CB) ===
    {
        std::string cid1 = "//A/60/CD2";
        std::string cid6 = "//A/60/CB";

        gemmi::CRA cra1 = find_atom_by_cid(model, cid1);
        gemmi::CRA cra6 = find_atom_by_cid(model, cid6);

        EXPECT_NE(cra1.atom, nullptr);
        EXPECT_NE(cra6.atom, nullptr);

        bool result4 = false;
        if (cra1.atom && cra6.atom) {
            result4 = coot::in_same_ring_gemmi(&geom, cra1, cra6, ring_list_map);
        }

        EXPECT_FALSE(result4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}