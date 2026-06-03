#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "function.hh"
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" ");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" ");
    return s.substr(start, end - start + 1);
}

TEST(OracleTest, atom_tree) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    gemmi::Residue* res = nullptr;
    for (auto& r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    int n_atoms = static_cast<int>(res->atoms.size());

    std::vector<std::string> raw_names;
    for (auto& atom : res->atoms) {
        raw_names.push_back(atom.name);
    }

    auto find_atom = [&](const std::string& trimmed) {
        for (auto& rn : raw_names) {
            if (trim(rn) == trimmed) return rn;
        }
        return std::string("");
    };

    coot::dictionary_residue_restraints_t rest;
    rest.residue_info.comp_id = res->name;

    std::string N   = find_atom("N");
    std::string CA  = find_atom("CA");
    std::string C   = find_atom("C");
    std::string O   = find_atom("O");
    std::string CB  = find_atom("CB");

    coot::dict_chem_comp_tree_t entry;

    if (!N.empty() && !CA.empty()) {
        entry.atom_id = N;
        entry.atom_back = "";
        entry.atom_forward = CA;
        entry.connect_type = "START";
        rest.tree.push_back(entry);
    }

    if (!CA.empty() && !C.empty()) {
        entry.atom_id = CA;
        entry.atom_back = N;
        entry.atom_forward = C;
        entry.connect_type = "";
        rest.tree.push_back(entry);
    }

    if (!C.empty() && !O.empty()) {
        entry.atom_id = C;
        entry.atom_back = CA;
        entry.atom_forward = O;
        entry.connect_type = "";
        rest.tree.push_back(entry);
    }

    if (!O.empty()) {
        entry.atom_id = O;
        entry.atom_back = C;
        entry.atom_forward = "";
        entry.connect_type = "END";
        rest.tree.push_back(entry);
    }

    if (!CB.empty()) {
        entry.atom_id = CB;
        entry.atom_back = CA;
        entry.atom_forward = "";
        entry.connect_type = "END";
        rest.tree.push_back(entry);
    }

    EXPECT_EQ(rest.comp_id(), "ILE");
    EXPECT_EQ(rest.tree.size(), 5u);

    // Case 1: Valid tree with matching atoms
    {
        std::string altconf = "";

        std::map<std::string, int> name_to_index_int;
        std::vector<coot::atom_vertex> atom_vertex_vec;
        atom_tree_t_gemmi({}, 0, res, altconf, name_to_index_int, atom_vertex_vec);

        std::map<std::string, coot::map_index_t> name_to_index;
        for (auto& p : name_to_index_int) {
            name_to_index[p.first] = coot::map_index_t(p.second);
        }

        fill_atom_vertex_vec_gemmi(rest, res, altconf, name_to_index, atom_vertex_vec);

        EXPECT_EQ(atom_vertex_vec.size(), 8u);
        EXPECT_EQ(name_to_index.size(), 9u);

        // vertex[0]: forward=[1] backward=[] type=START
        {
            const coot::atom_vertex &v = atom_vertex_vec[0];
            EXPECT_EQ(v.forward.size(), 1u);
            EXPECT_EQ(v.forward[0], 1);
            EXPECT_EQ(v.backward.size(), 0u);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::START);
        }
        // vertex[1]: forward=[2,4] backward=[0] type=STANDARD
        {
            const coot::atom_vertex &v = atom_vertex_vec[1];
            EXPECT_EQ(v.forward.size(), 2u);
            EXPECT_EQ(v.forward[0], 2);
            EXPECT_EQ(v.forward[1], 4);
            EXPECT_EQ(v.backward.size(), 1u);
            EXPECT_EQ(v.backward[0], 0);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::STANDARD);
        }
        // vertex[2]: forward=[3] backward=[1] type=STANDARD
        {
            const coot::atom_vertex &v = atom_vertex_vec[2];
            EXPECT_EQ(v.forward.size(), 1u);
            EXPECT_EQ(v.forward[0], 3);
            EXPECT_EQ(v.backward.size(), 1u);
            EXPECT_EQ(v.backward[0], 1);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::STANDARD);
        }
        // vertex[3]: forward=[] backward=[2] type=END
        {
            const coot::atom_vertex &v = atom_vertex_vec[3];
            EXPECT_EQ(v.forward.size(), 0u);
            EXPECT_EQ(v.backward.size(), 1u);
            EXPECT_EQ(v.backward[0], 2);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::END);
        }
        // vertex[4]: forward=[] backward=[1] type=END
        {
            const coot::atom_vertex &v = atom_vertex_vec[4];
            EXPECT_EQ(v.forward.size(), 0u);
            EXPECT_EQ(v.backward.size(), 1u);
            EXPECT_EQ(v.backward[0], 1);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::END);
        }
        // vertex[5]: forward=[] backward=[] type=NONE
        {
            const coot::atom_vertex &v = atom_vertex_vec[5];
            EXPECT_EQ(v.forward.size(), 0u);
            EXPECT_EQ(v.backward.size(), 0u);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::NONE);
        }
        // vertex[6]: forward=[] backward=[] type=NONE
        {
            const coot::atom_vertex &v = atom_vertex_vec[6];
            EXPECT_EQ(v.forward.size(), 0u);
            EXPECT_EQ(v.backward.size(), 0u);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::NONE);
        }
        // vertex[7]: forward=[] backward=[] type=NONE
        {
            const coot::atom_vertex &v = atom_vertex_vec[7];
            EXPECT_EQ(v.forward.size(), 0u);
            EXPECT_EQ(v.backward.size(), 0u);
            EXPECT_EQ(v.connection_type, coot::atom_vertex::NONE);
        }
    }

    // Case 2: Empty tree — no START
    {
        coot::dictionary_residue_restraints_t rest_empty;
        rest_empty.residue_info.comp_id = "FAKE";
        EXPECT_EQ(rest_empty.tree.size(), 0u);

        std::string altconf = "";
        std::map<std::string, int> name_to_index_int;
        std::vector<coot::atom_vertex> atom_vertex_vec;
        atom_tree_t_gemmi({}, 0, res, altconf, name_to_index_int, atom_vertex_vec);

        std::map<std::string, coot::map_index_t> name_to_index;
        for (auto& p : name_to_index_int) {
            name_to_index[p.first] = coot::map_index_t(p.second);
        }

        fill_atom_vertex_vec_gemmi(rest, res, altconf, name_to_index, atom_vertex_vec);
        size_t old_size = atom_vertex_vec.size();

        bool result = fill_atom_vertex_vec_gemmi(rest_empty, res, altconf, name_to_index, atom_vertex_vec);

        EXPECT_FALSE(result);
        EXPECT_EQ(atom_vertex_vec.size(), 8u);
        EXPECT_EQ(old_size, atom_vertex_vec.size());
    }

    // Case 3: Tree with non-matching atom names
    {
        coot::dictionary_residue_restraints_t rest_fake;
        rest_fake.residue_info.comp_id = "FAKE";

        coot::dict_chem_comp_tree_t fake_entry;
        fake_entry.atom_id = " ZZZ";
        fake_entry.atom_back = "";
        fake_entry.atom_forward = " YYY";
        fake_entry.connect_type = "START";
        rest_fake.tree.push_back(fake_entry);

        EXPECT_EQ(rest_fake.tree.size(), 1u);

        std::string altconf = "";
        std::map<std::string, int> name_to_index_int;
        std::vector<coot::atom_vertex> atom_vertex_vec;
        atom_tree_t_gemmi({}, 0, res, altconf, name_to_index_int, atom_vertex_vec);

        std::map<std::string, coot::map_index_t> name_to_index;
        for (auto& p : name_to_index_int) {
            name_to_index[p.first] = coot::map_index_t(p.second);
        }

        fill_atom_vertex_vec_gemmi(rest, res, altconf, name_to_index, atom_vertex_vec);
        size_t old_size_4 = atom_vertex_vec.size();

        bool result_fake = fill_atom_vertex_vec_gemmi(rest_fake, res, altconf, name_to_index, atom_vertex_vec);

        EXPECT_FALSE(result_fake);
        EXPECT_EQ(atom_vertex_vec.size(), 8u);
        EXPECT_EQ(old_size_4, atom_vertex_vec.size());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}