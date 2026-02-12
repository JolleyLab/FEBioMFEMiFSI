#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

#include "MFEMFSI/ImmersedFSI/NodalAverageOutputs.h"

#include <FECore/mat3d.h>
#include <FECore/vec3d.h>

static bool nearly_equal(double a, double b, double tol)
{
    double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= tol * scale;
}

static void sort_desc(double &x, double &y, double &z)
{
    if (x < y) std::swap(x, y);
    if (y < z) std::swap(y, z);
    if (x < y) std::swap(x, y);
}

int main()
{
    std::mt19937_64 rng(12345ULL);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    const int trials = 500; // keep quick
    int fails = 0;

    for (int t = 0; t < trials; ++t)
    {
        // Generate random symmetric components
        double a11 = dist(rng);
        double a22 = dist(rng);
        double a33 = dist(rng);
        double a12 = dist(rng);
        double a13 = dist(rng);
        double a23 = dist(rng);

        // Our analytic eigenvalues
        double l1, l2, l3;
        nodal_avg::SymmetricEigenvalues(a11, a22, a33, a12, a13, a23, l1, l2, l3);

        // FEBio mat3ds expects components (xx, yy, zz, xy, yz, xz)
        mat3ds m(a11, a22, a33, a12, a23, a13);
        double d[3];
        vec3d v[3];
        m.eigen(d, v);

        double f1 = d[0], f2 = d[1], f3 = d[2];
        sort_desc(f1, f2, f3);

        bool ok = nearly_equal(l1, f1, 1e-12) &&
                  nearly_equal(l2, f2, 1e-12) &&
                  nearly_equal(l3, f3, 1e-12);

        if (!ok)
        {
            ++fails;
            std::cout << "[Trial " << t << "]\n"
                      << "  A: [" << a11 << ", " << a12 << ", " << a13 << "; "
                      << a12 << ", " << a22 << ", " << a23 << "; "
                      << a13 << ", " << a23 << ", " << a33 << "]\n"
                      << "  ours:    " << l1 << ", " << l2 << ", " << l3 << "\n"
                      << "  FEBio:   " << f1 << ", " << f2 << ", " << f3 << "\n"
                      << std::endl;
        }
    }

    if (fails == 0)
    {
        std::cout << "Random symmetric compare: PASSED (" << trials << " trials)" << std::endl;
        return 0;
    }
    else
    {
        std::cerr << fails << " trials FAILED out of " << trials << std::endl;
        return 1;
    }
}
