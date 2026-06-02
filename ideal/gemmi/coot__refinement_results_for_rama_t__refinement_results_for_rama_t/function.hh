#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "geometry/residue-and-atom-specs.hh"
#include <sstream>

namespace coot {

static std::string pad_atom_name_4(const std::string &name) {
    std::string padded = " " + name;
    while (padded.size() < 4) padded += ' ';
    return padded;
}

struct refinement_results_for_rama_t_gemmi {
    float distortion;
    coot::atom_spec_t atom_spec_CA;
    float ball_pos_x;
    float ball_pos_y;
    float ball_pos_z;

    refinement_results_for_rama_t_gemmi()
        : distortion(0),
          atom_spec_CA(),
          ball_pos_x(0),
          ball_pos_y(0),
          ball_pos_z(0) {}

    refinement_results_for_rama_t_gemmi(
        gemmi::Model &model,
        const gemmi::Atom *at_1,
        const gemmi::Atom *at_2,
        const gemmi::Atom *at_3,
        const gemmi::Atom *at_4,
        const gemmi::Atom *at_5,
        float distortion_in)
        : distortion(distortion_in),
          atom_spec_CA(make_spec(model, at_3)),
          ball_pos_x(0),
          ball_pos_y(0),
          ball_pos_z(0) {

        if (at_3) {
            ball_pos_x = at_3->pos.x + 0.5;
            ball_pos_y = at_3->pos.y;
            ball_pos_z = at_3->pos.z;
        }

        if (at_1 && at_2 && at_3 && at_4 && at_5) {
            clipper::Coord_orth p2 = coot::co_gemmi(at_2);
            clipper::Coord_orth p3 = coot::co_gemmi(at_3);
            clipper::Coord_orth p4 = coot::co_gemmi(at_4);
            clipper::Coord_orth v1(p3 - p2);
            clipper::Coord_orth v2(p3 - p4);
            clipper::Coord_orth v3(p4 - p2);
            clipper::Coord_orth v1_uv(v1.unit());
            clipper::Coord_orth v2_uv(v2.unit());
            clipper::Coord_orth v3_uv(v3.unit());
            clipper::Coord_orth v4(clipper::Coord_orth::cross(v2_uv, v1_uv));
            clipper::Coord_orth p2p24_mid_point(0.5 * (p4 + p2));
            clipper::Coord_orth mid_point_to_CA(p3 - p2p24_mid_point);
            clipper::Coord_orth delta = 0.2 * mid_point_to_CA + 0.4 * v4;
            ball_pos_x = delta.x() + at_3->pos.x;
            ball_pos_y = delta.y() + at_3->pos.y;
            ball_pos_z = delta.z() + at_3->pos.z;
        }
    }

  private:
    static coot::atom_spec_t make_spec(gemmi::Model &model, const gemmi::Atom *at) {
        coot::atom_spec_t spec;
        if (!at) return spec;
        for (auto &chain : model.chains) {
            for (auto &res : chain.residues) {
                for (auto &atom : res.atoms) {
                    if (&atom == at) {
                        spec.model_number = model.num;
                        spec.chain_id = chain.name;
                        spec.res_no = res.seqid.num.value;
                        spec.ins_code = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
                        spec.atom_name = pad_atom_name_4(atom.name);
                        spec.alt_conf = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
                        spec.int_user_data = -1;
                        spec.float_user_data = -1;
                        return spec;
                    }
                }
            }
        }
        return spec;
    }
};

} // namespace coot