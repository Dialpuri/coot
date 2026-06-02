#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, symmetry_move_atoms) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u) << "Failed to load PDB";

    // Construct cell and spacegroup directly (PDB has no crystallography)
    clipper::Cell cell(clipper::Cell_descr(78.0, 78.0, 78.0, 90.0, 90.0, 90.0));
    clipper::Spacegroup spacegroup(clipper::Spgr_descr("P 21 21 21"));

    // Collect protein coords and water atoms from the molecule
    std::vector<clipper::Coord_orth> protein_coords;
    std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth> > water_atoms_in;

    {
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::string residue_name = residue.name;
                for (gemmi::Atom& at : residue.atoms) {
                    if (residue_name == "WAT" || residue_name == "HOH") {
                        if (at.altloc == '\0') {
                            clipper::Coord_orth c(at.pos.x, at.pos.y, at.pos.z);
                            water_atoms_in.push_back(std::make_pair(&at, c));
                        }
                    } else {
                        if (at.altloc == '\0') {
                            std::string ele = at.element.name();
                            if (ele != "C") {
                                clipper::Coord_orth pt(at.pos.x, at.pos.y, at.pos.z);
                                protein_coords.push_back(pt);
                            }
                        }
                    }
                }
            }
        }
    }

    // Case 1: Original oracle case — water far from protein, should NOT be moved
    {
        EXPECT_EQ(protein_coords.size(), 742u);
        EXPECT_EQ(water_atoms_in.size(), 1u);
        EXPECT_NEAR(cell.a(), 78.0, 1e-4);
        EXPECT_NEAR(cell.b(), 78.0, 1e-4);
        EXPECT_NEAR(cell.c(), 78.0, 1e-4);
        EXPECT_EQ(spacegroup.symbol_laue(), "mmm");
        EXPECT_NEAR(water_atoms_in[0].second.x(), 30.744, 1e-4);
        EXPECT_NEAR(water_atoms_in[0].second.y(), 8.475, 1e-4);
        EXPECT_NEAR(water_atoms_in[0].second.z(), 36.748, 1e-4);

        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth> > water_atoms_out =
            coot::util::symmetry_move_atoms_gemmi(protein_coords, water_atoms_in, cell, spacegroup);

        EXPECT_EQ(water_atoms_out.size(), 1u);
        EXPECT_EQ(water_atoms_out[0].first, nullptr);  // not moved
        EXPECT_NEAR(water_atoms_out[0].second.x(), 30.744, 1e-4);
        EXPECT_NEAR(water_atoms_out[0].second.y(), 8.475, 1e-4);
        EXPECT_NEAR(water_atoms_out[0].second.z(), 36.748, 1e-4);

        // Count moved waters
        int n_moved = 0;
        for (unsigned int i = 0; i < water_atoms_out.size(); i++) {
            if (water_atoms_out[i].first != nullptr) {
                n_moved++;
            }
        }
        EXPECT_EQ(n_moved, 0);
    }

    // Case 2: Water placed on top of a protein atom
    // Because closest_distance=0 => closest_distance*0.8=0 => dist<0 is never true,
    // the identity position always passes, so water stays in place
    {
        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth> > water_atoms_clash;
        clipper::Coord_orth water_at_protein(protein_coords[0].x(), protein_coords[0].y(), protein_coords[0].z());

        gemmi::Atom* real_atom_ptr = nullptr;
        {
            gemmi::Model& model = st.models[0];
            if (!model.chains.empty()) {
                gemmi::Chain& chain = model.chains[0];
                if (!chain.residues.empty()) {
                    gemmi::Residue& res = chain.residues[0];
                    if (!res.atoms.empty()) {
                        real_atom_ptr = &res.atoms[0];
                    }
                }
            }
        }
        EXPECT_NE(real_atom_ptr, nullptr);
        water_atoms_clash.push_back(std::make_pair(real_atom_ptr, water_at_protein));

        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth> > water_atoms_out =
            coot::util::symmetry_move_atoms_gemmi(protein_coords, water_atoms_clash, cell, spacegroup);

        EXPECT_EQ(water_atoms_out.size(), 1u);
        // closest_distance=0, so closest_distance*0.8=0, dist<0 always false => stays
        EXPECT_EQ(water_atoms_out[0].first, nullptr);
        EXPECT_NEAR(water_atoms_out[0].second.x(), water_at_protein.x(), 1e-4);
        EXPECT_NEAR(water_atoms_out[0].second.y(), water_at_protein.y(), 1e-4);
        EXPECT_NEAR(water_atoms_out[0].second.z(), water_at_protein.z(), 1e-4);
    }

    // Case 3: Empty protein_coords throws std::runtime_error (no atoms to compute median)
    {
        std::vector<clipper::Coord_orth> empty_protein;
        std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth> > water_atoms_empty;
        water_atoms_empty.push_back(std::make_pair(static_cast<gemmi::Atom*>(nullptr), clipper::Coord_orth(10.0, 10.0, 10.0)));

        EXPECT_THROW(
            coot::util::symmetry_move_atoms_gemmi(empty_protein, water_atoms_empty, cell, spacegroup),
            std::runtime_error
        );
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}