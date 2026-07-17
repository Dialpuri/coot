/*
 * godot/src/coot_molecule.cpp
 */
#include "coot_molecule.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/packed_color_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/classes/mesh.hpp>

using namespace godot;

CootMolecule::CootMolecule() {}
CootMolecule::~CootMolecule() {}

void CootMolecule::_bind_methods() {
    ClassDB::bind_method(D_METHOD("load_pdb", "path"), &CootMolecule::load_pdb);

    ClassDB::bind_method(D_METHOD("set_scale", "s"), &CootMolecule::set_scale);
    ClassDB::bind_method(D_METHOD("get_scale"), &CootMolecule::get_scale);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "scale"), "set_scale", "get_scale");

    ClassDB::bind_method(D_METHOD("build_bonds_mesh", "imol", "mode",
                                  "bond_width", "ratio", "smoothness"),
                         &CootMolecule::build_bonds_mesh);
    ClassDB::bind_method(D_METHOD("build_gaussian_surface", "imol", "sigma",
                                  "contour_level"),
                         &CootMolecule::build_gaussian_surface);
    ClassDB::bind_method(D_METHOD("build_cartoon", "imol", "cid",
                                  "colour_scheme", "style"),
                         &CootMolecule::build_cartoon);

    ClassDB::bind_method(D_METHOD("get_last_message"), &CootMolecule::get_last_message);
}

int CootMolecule::load_pdb(const String &path) {
    std::string p = path.utf8().get_data();
    int imol = mc.read_pdb(p);
    if (imol < 0)
        message = String("Failed to read coordinates: ") + path;
    else
        message = String();
    return imol;
}

void CootMolecule::set_scale(float s) { scale = s; }
float CootMolecule::get_scale() const { return scale; }

String CootMolecule::get_last_message() const { return message; }

void CootMolecule::ensure_centre(const coot::simple_mesh_t &m) {
    if (centre_valid) return;
    if (m.vertices.empty()) return;
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (const auto &v : m.vertices) {
        sx += v.pos.x;
        sy += v.pos.y;
        sz += v.pos.z;
    }
    const double n = static_cast<double>(m.vertices.size());
    centre = Vector3(sx / n, sy / n, sz / n);
    centre_valid = true;
}

Ref<ArrayMesh> CootMolecule::to_array_mesh(const coot::simple_mesh_t &m) {
    Ref<ArrayMesh> am;
    am.instantiate();

    if (m.status == 0) {
        message = String("Coot reported a bad mesh (status 0)");
        return am;
    }
    if (m.vertices.empty() || m.triangles.empty()) {
        message = String("Empty mesh");
        return am;
    }

    ensure_centre(m);

    PackedVector3Array verts;
    PackedVector3Array normals;
    PackedColorArray   colors;
    PackedInt32Array   indices;

    const int nv = static_cast<int>(m.vertices.size());
    verts.resize(nv);
    normals.resize(nv);
    colors.resize(nv);
    {
        Vector3 *vp = verts.ptrw();
        Vector3 *np = normals.ptrw();
        Color   *cp = colors.ptrw();
        for (int i = 0; i < nv; ++i) {
            const auto &v = m.vertices[i];
            vp[i] = (Vector3(v.pos.x, v.pos.y, v.pos.z) - centre) * scale;
            np[i] = Vector3(v.normal.x, v.normal.y, v.normal.z);
            cp[i] = Color(v.color.r, v.color.g, v.color.b, v.color.a);
        }
    }

    const int nt = static_cast<int>(m.triangles.size());
    indices.resize(nt * 3);
    {
        int32_t *ip = indices.ptrw();
        for (int t = 0; t < nt; ++t) {
            ip[t * 3 + 0] = static_cast<int32_t>(m.triangles[t][0]);
            ip[t * 3 + 1] = static_cast<int32_t>(m.triangles[t][1]);
            ip[t * 3 + 2] = static_cast<int32_t>(m.triangles[t][2]);
        }
    }

    Array arrays;
    arrays.resize(Mesh::ARRAY_MAX);
    arrays[Mesh::ARRAY_VERTEX] = verts;
    arrays[Mesh::ARRAY_NORMAL] = normals;
    arrays[Mesh::ARRAY_COLOR]  = colors;
    arrays[Mesh::ARRAY_INDEX]  = indices;

    am->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);
    return am;
}

Ref<ArrayMesh> CootMolecule::build_bonds_mesh(int imol, const String &mode,
                                              float bond_width, float ratio,
                                              int smoothness) {
    std::string m = mode.utf8().get_data();
    coot::simple_mesh_t mesh =
        mc.get_bonds_mesh(imol, m, /*against_a_dark_background*/ true,
                          bond_width, ratio, smoothness);
    return to_array_mesh(mesh);
}

Ref<ArrayMesh> CootMolecule::build_gaussian_surface(int imol, float sigma,
                                                    float contour_level) {
    // box_radius 5.0, grid_scale 0.7, b_factor 100.0 are Coot's defaults.
    coot::simple_mesh_t mesh =
        mc.get_gaussian_surface(imol, sigma, contour_level, 5.0f, 0.7f, 100.0f);
    return to_array_mesh(mesh);
}

Ref<ArrayMesh> CootMolecule::build_cartoon(int imol, const String &cid,
                                           const String &colour_scheme,
                                           const String &style) {
    std::string c = cid.utf8().get_data();
    std::string s = colour_scheme.utf8().get_data();
    std::string st = style.utf8().get_data();
    // secondary_structure_usage_flag: 2 = CALC_SECONDARY_STRUCTURE.
    coot::simple_mesh_t mesh =
        mc.get_molecular_representation_mesh(imol, c, s, st, 2);
    return to_array_mesh(mesh);
}
