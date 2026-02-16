#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "MFEMFSI/ImmersedFSI/NodalAverageOutputs.h"

static bool nearly_equal(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

static bool check_case(const char* name,
                       double a11, double a22, double a33,
                       double a12, double a13, double a23,
                       double e1, double e2, double e3,
                       double tol = 1e-12)
{
    double l1, l2, l3;
    nodal_avg::SymmetricEigenvalues(a11, a22, a33, a12, a13, a23, l1, l2, l3);

    // Sort expected to descending order
    std::vector<double> expected{e1, e2, e3};
    std::sort(expected.begin(), expected.end(), std::greater<double>());

    bool ok = nearly_equal(l1, expected[0], tol) &&
              nearly_equal(l2, expected[1], tol) &&
              nearly_equal(l3, expected[2], tol);

    std::cout << "[" << name << "]\n"
              << "  computed: " << l1 << ", " << l2 << ", " << l3 << "\n"
              << "  expected: " << expected[0] << ", " << expected[1] << ", " << expected[2] << "\n"
              << (ok ? "  RESULT: PASS\n" : "  RESULT: FAIL\n")
              << std::endl;
    return ok;
}

int main()
{
    int fails = 0;

    // Case 1: Diagonal matrix — straightforward eigenvalues
    fails += !check_case("Diagonal", 1.0, 2.0, 3.0, 0.0, 0.0, 0.0,
                         1.0, 2.0, 3.0);

    // Case 2: Symmetric with 2x2 block [2 1; 1 2] and a separate 1
    // Block eigenvalues are 3 and 1; plus 1 => {3, 1, 1}
    fails += !check_case("Block2x2", 2.0, 2.0, 1.0, 1.0, 0.0, 0.0,
                         3.0, 1.0, 1.0);

    // Case 3: Another symmetric with 2x2 block [[4,2],[2,3]]
    // Eigenvalues of block: (7 ± sqrt(17))/2 ≈ 5.5615528128 and 1.4384471872, plus 1
    fails += !check_case("Block2x2+1", 4.0, 3.0, 1.0, 2.0, 0.0, 0.0,
                         (7.0 + std::sqrt(17.0))/2.0, (7.0 - std::sqrt(17.0))/2.0, 1.0);

    // Case 4: Identity matrix — triple eigenvalue
    fails += !check_case("Identity", 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                         1.0, 1.0, 1.0);

    // Case 5: Repeated off-diagonal symmetric
    // For matrix with diagonal a and all off-diagonals b, eigenvalues are {a+2b, a-b, a-b}
    fails += !check_case("OffDiag", 3.0, 3.0, 3.0, 0.5, 0.5, 0.5,
                         4.0, 2.5, 2.5, 1e-10);

    // Case 6: Negative off-diagonals (same structure)
    // a = 3, b = -0.5 -> {a+2b, a-b, a-b} = {2.0, 3.5, 3.5}
    fails += !check_case("OffDiagNegative", 3.0, 3.0, 3.0, -0.5, -0.5, -0.5,
                         3.5, 3.5, 2.0, 1e-10);

    // Case 7: Complete graph with zero diagonal and ones off-diagonal
    // Matrix [[0,1,1],[1,0,1],[1,1,0]] has eigenvalues {2,-1,-1}
    fails += !check_case("ZeroDiag_OnesOff", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                         2.0, -1.0, -1.0, 1e-12);

    // Case 8: Rank-1 matrix u u^T with u=(1,2,3) -> eigenvalues {||u||^2, 0, 0} = {14,0,0}
    // Entries: a11=1, a22=4, a33=9, a12=2, a13=3, a23=6
    fails += !check_case("Rank1_uuT", 1.0, 4.0, 9.0, 2.0, 3.0, 6.0,
                         14.0, 0.0, 0.0, 1e-12);

    // Case 9: 2x2 block with negative off-diagonal and separate diagonal
    // Block [[2,-1],[-1,2]] -> eigenvalues {3,1}; with a33=4 -> {4,3,1}
    fails += !check_case("Block2x2Neg+4", 2.0, 2.0, 4.0, -1.0, 0.0, 0.0,
                         4.0, 3.0, 1.0, 1e-12);

    // Case 10: Scaling test — scale OffDiag case by 1e8 (eigenvalues scale accordingly)
    fails += !check_case("Scaled_OffDiag", 3.0e8, 3.0e8, 3.0e8, 0.5e8, 0.5e8, 0.5e8,
                         4.0e8, 2.5e8, 2.5e8, 1e-6);

    // Case 11: Zero matrix — triple zero eigenvalue
    fails += !check_case("ZeroMatrix", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 1e-12);

    if (fails == 0)
    {
        std::cout << "All tests PASSED" << std::endl;
        return 0;
    }
    else
    {
        std::cerr << fails << " tests FAILED" << std::endl;
        return 1;
    }
}
