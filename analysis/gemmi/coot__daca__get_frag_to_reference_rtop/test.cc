#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::vector<const gemmi::Atom*> get_atoms_from_residue(const gemmi::Residue& res) {
    std::vector<const gemmi::Atom*> atoms;
    for (const gemmi::Atom& atom : res.atoms)
        atoms.push_back(&atom);
    return atoms;
}

// Build a simple reference fragments map for testing
static std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>>
make_test_reference_fragments() {
    std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> ref;
    // Add a dummy entry for ALA with one fragment of 4 atoms
    std::vector<clipper::Coord_orth> frag1;
    frag1.push_back(clipper::Coord_orth(1.0, 2.0, 3.0));
    frag1.push_back(clipper::Coord_orth(4.0, 5.0, 6.0));
    frag1.push_back(clipper::Coord_orth(7.0, 8.0, 9.0));
    frag1.push_back(clipper::Coord_orth(10.0, 11.0, 12.0));
    ref["ALA"].push_back(frag1);
    return ref;
}

TEST(OracleTest, GetFragToReferenceRtop) {
    // Case 1: Valid fragment index - should succeed
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0];
        const gemmi::Residue& res = chain.residues[9]; // residue 10

        std::vector<const gemmi::Atom*> frag_atoms = get_atoms_from_residue(res);
        auto ref_frags = make_test_reference_fragments();

        // Scale the fragment atoms to match the reference (for a valid transformation)
        // We will just verify that the function runs without crashing and returns success when sizes match
        std::vector<const gemmi::Atom*> simple_frag;
        simple_frag.push_back(&res.atoms[0]);
        simple_frag.push_back(&res.atoms[1]);
        simple_frag.push_back(&res.atoms[2]);
        simple_frag.push_back(&res.atoms[3]);

        std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> simple_ref;
        std::vector<clipper::Coord_orth> rf;
        rf.push_back(clipper::Coord_orth(1.0, 2.0, 3.0));
        rf.push_back(clipper::Coord_orth(4.0, 5.0, 6.0));
        rf.push_back(clipper::Coord_orth(7.0, 8.0, 9.0));
        rf.push_back(clipper::Coord_orth(10.0, 11.0, 12.0));
        simple_ref["ALA"].push_back(rf);

        std::pair<bool, clipper::RTop_orth> result =
            coot::daca::get_frag_to_reference_rtop_gemmi("ALA", 0, simple_frag, simple_ref);

        EXPECT_TRUE(result.first) << "Valid fragment index should return success";

        if (result.first) {
            clipper::RTop_orth rtop = result.second;
            // rtop is valid
        }
    }

    // Case 2: Invalid fragment index - should fail
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0];
        const gemmi::Residue& res = chain.residues[9];

        std::vector<const gemmi::Atom*> simple_frag;
        simple_frag.push_back(&res.atoms[0]);
        simple_frag.push_back(&res.atoms[1]);

        std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> simple_ref;
        std::vector<clipper::Coord_orth> rf;
        rf.push_back(clipper::Coord_orth(1.0, 2.0, 3.0));
        rf.push_back(clipper::Coord_orth(4.0, 5.0, 6.0));
        simple_ref["ALA"].push_back(rf);

        unsigned int frag_idx = 999; // out of bounds
        std::pair<bool, clipper::RTop_orth> result =
            coot::daca::get_frag_to_reference_rtop_gemmi("ALA", frag_idx, simple_frag, simple_ref);

        EXPECT_FALSE(result.first) << "Out-of-bounds fragment index should return failure";
    }

    // Case 3: Invalid residue name - should fail
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0];
        const gemmi::Residue& res = chain.residues[9];

        std::vector<const gemmi::Atom*> simple_frag;
        simple_frag.push_back(&res.atoms[0]);
        simple_frag.push_back(&res.atoms[1]);

        std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> simple_ref;
        std::vector<clipper::Coord_orth> rf;
        rf.push_back(clipper::Coord_orth(1.0, 2.0, 3.0));
        rf.push_back(clipper::Coord_orth(4.0, 5.0, 6.0));
        simple_ref["ALA"].push_back(rf);

        std::pair<bool, clipper::RTop_orth> result =
            coot::daca::get_frag_to_reference_rtop_gemmi("XXX", 0, simple_frag, simple_ref);

        EXPECT_FALSE(result.first) << "Invalid residue name should return failure";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}