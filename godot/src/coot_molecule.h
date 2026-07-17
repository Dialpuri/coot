/*
 * godot/src/coot_molecule.h
 *
 * GDExtension bindings for libcootapi: a CootMolecule node that turns
 * Coot's simple_mesh_t into a Godot ArrayMesh at runtime.
 */
#ifndef COOT_GODOT_COOT_MOLECULE_H
#define COOT_GODOT_COOT_MOLECULE_H

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/vector3.hpp>

// libcootapi
#include "api/molecules-container.hh"
#include "coot-utils/simple-mesh.hh"

class CootMolecule : public godot::RefCounted {
    GDCLASS(CootMolecule, godot::RefCounted)

    // The whole of libcootapi lives behind this one object.
    molecules_container_t mc{false}; // false = not verbose

    // Recentre offset (in Angstroms), computed from the first mesh generated so
    // that every representation of this molecule shares one origin and stays
    // registered. World position = (pos - centre) * scale.
    godot::Vector3 centre;
    bool centre_valid = false;

    // Angstrom -> Godot world-unit scale. ~0.1 is comfortable for room-scale VR.
    float scale = 0.1f;

    // Last human-readable message (surfaced from Coot on failure).
    godot::String message;

protected:
    static void _bind_methods();

public:
    CootMolecule();
    ~CootMolecule();

    //! Load a PDB or mmCIF coordinate file. Returns the molecule index (imol),
    //! or -1 on failure (see get_last_message()).
    int load_pdb(const godot::String &path);

    void set_scale(float s);
    float get_scale() const;

    //! Ball-and-stick / liquorice / VDW bonds mesh.
    //! mode: "COLOUR-BY-CHAIN-AND-DICTIONARY", "CA+LIGANDS" or "VDW-BALLS".
    //! bond_width ~0.12, ratio 1.0 (liquorice) .. 3.0 (ball & stick), smoothness 1..3.
    godot::Ref<godot::ArrayMesh> build_bonds_mesh(int imol, const godot::String &mode,
                                                  float bond_width, float ratio,
                                                  int smoothness);

    //! Gaussian molecular surface (one surface per chain). Defaults: sigma 4.4,
    //! contour 4.0.
    godot::Ref<godot::ArrayMesh> build_gaussian_surface(int imol, float sigma,
                                                        float contour_level);

    //! Cartoon / ribbon. colour_scheme: "colorRampChainsScheme",
    //! "colorBySecondaryScheme" or "Chain". style: "Ribbon" or "MolecularSurface".
    godot::Ref<godot::ArrayMesh> build_cartoon(int imol, const godot::String &cid,
                                               const godot::String &colour_scheme,
                                               const godot::String &style);

    godot::String get_last_message() const;

private:
    // The core converter: coot::simple_mesh_t -> godot::ArrayMesh.
    godot::Ref<godot::ArrayMesh> to_array_mesh(const coot::simple_mesh_t &m);
    // Compute the recentre offset from the first non-empty mesh we see.
    void ensure_centre(const coot::simple_mesh_t &m);
};

#endif // COOT_GODOT_COOT_MOLECULE_H
