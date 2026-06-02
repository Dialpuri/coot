#pragma once
#include <ostream>

struct matrix4x4 {
    double data[4][4];

    // Default constructor — identity matrix
    matrix4x4() {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                data[i][j] = (i == j) ? 1.0 : 0.0;
    }

    // Constructor from 4×4 array
    matrix4x4(const double (&m)[4][4]) {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                data[i][j] = m[i][j];
    }

    // From quaternion (x, y, z, w)
    static matrix4x4 matrixWithQuaternion(const double quat[4]) {
        double x = quat[0], y = quat[1], z = quat[2], w = quat[3];
        matrix4x4 m;
        m.data[0][0] = 1.0 - 2.0 * (y * y + z * z);
        m.data[0][1] = 2.0 * (x * y - w * z);
        m.data[0][2] = 2.0 * (x * z + w * y);
        m.data[0][3] = 0.0;
        m.data[1][0] = 2.0 * (x * y + w * z);
        m.data[1][1] = 1.0 - 2.0 * (x * x + z * z);
        m.data[1][2] = 2.0 * (y * z - w * x);
        m.data[1][3] = 0.0;
        m.data[2][0] = 2.0 * (x * z - w * y);
        m.data[2][1] = 2.0 * (y * z + w * x);
        m.data[2][2] = 1.0 - 2.0 * (x * x + y * y);
        m.data[2][3] = 0.0;
        m.data[3][0] = 0.0;
        m.data[3][1] = 0.0;
        m.data[3][2] = 0.0;
        m.data[3][3] = 1.0;
        return m;
    }
};

inline std::ostream &print_matrix_gemmi(std::ostream &out, const matrix4x4 &m) {
    out << "(";
    for (int i = 0; i < 4; i++) {
        out << " [ " << m.data[i][0] << " " << m.data[i][1] << " "
            << m.data[i][2] << " " << m.data[i][3] << " ]";
        out << (i == 3 ? ")\n" : "\n");
    }
    return out;
}