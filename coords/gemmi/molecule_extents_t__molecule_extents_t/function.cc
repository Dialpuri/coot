#include "function.hh"

molecule_extents_t molecule_extents_t::molecule_extents_t_gemmi(
    const gemmi::Model& model, float expansion_size)
{
    molecule_extents_t result;
    result.expansion_size_ = expansion_size;

    float atom_x, atom_y, atom_z;
    float max_x, max_y, max_z, min_x, min_y, min_z;

    max_x = -99999999.9f;
    max_y = -99999999.9f;
    max_z = -99999999.9f;

    min_x = 99999999.9f;
    min_y = 99999999.9f;
    min_z = 99999999.9f;

    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            for (const gemmi::Atom& atom : res.atoms) {
                atom_x = atom.pos.x;
                atom_y = atom.pos.y;
                atom_z = atom.pos.z;

                if (atom_x > max_x) {
                    max_x = atom_x;
                    result.right = coot::Cartesian(atom_x, atom_y, atom_z);
                }

                if (atom_x < min_x) {
                    min_x = atom_x;
                    result.left = coot::Cartesian(atom_x, atom_y, atom_z);
                }

                if (atom_y > max_y) {
                    max_y = atom_y;
                    result.top = coot::Cartesian(atom_x, atom_y, atom_z);
                }

                if (atom_y < min_y) {
                    min_y = atom_y;
                    result.bottom = coot::Cartesian(atom_x, atom_y, atom_z);
                }

                if (atom_z > max_z) {
                    max_z = atom_z;
                    result.back = coot::Cartesian(atom_x, atom_y, atom_z);
                }

                if (atom_z < min_z) {
                    min_z = atom_z;
                    result.front = coot::Cartesian(atom_x, atom_y, atom_z);
                }
            }
        }
    }

    float mid_x, mid_y, mid_z;
    mid_x = (result.left.get_x() + result.right.get_x()) * 0.5f;
    mid_y = (result.bottom.get_y() + result.top.get_y()) * 0.5f;
    mid_z = (result.front.get_z() + result.back.get_z()) * 0.5f;

    // Adjust the extents so that they are on the midpoints of the other axes.
    result.left   = coot::Cartesian(result.left.get_x() - expansion_size, mid_y, mid_z);
    result.right  = coot::Cartesian(result.right.get_x() + expansion_size, mid_y, mid_z);
    result.front  = coot::Cartesian(mid_x, mid_y, result.front.get_z() - expansion_size);
    result.back   = coot::Cartesian(mid_x, mid_y, result.back.get_z() + expansion_size);
    result.bottom = coot::Cartesian(mid_x, result.bottom.get_y() - expansion_size, mid_z);
    result.top    = coot::Cartesian(mid_x, result.top.get_y() + expansion_size, mid_z);

    // now make the centre for the above coordinates
    result.centre = result.front + result.back + result.left + result.right + result.top + result.bottom;
    result.centre *= 0.16666666f;

    return result;
}