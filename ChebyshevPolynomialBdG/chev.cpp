#include "chev.h"
#include <cmath>    // for std::acos
#include <vector>   // for std::vector
#include <iostream>
#include <Eigen/EigenValues>

namespace chebyshevpolynomialbdg {
    Chev::Chev()
        :   vec_ai_(Chev::NC),
            vec_delta_(Chev::LN_2, Chev::LN_2)
    {
        init_delta();
        calc_A();
        //for (auto ii = 0; ii < Chev::LN; ii++) {
        //    for (auto jj = 0; jj < Chev::LN; jj++) {
        //        std::cout << A_.coeff(ii, jj) << " ";
        //    }
        //    std::cout << std::endl;
        //}
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A_);
        //auto w = es.eigenvalues();
        //auto v = es.eigenvectors();
        //std::cout << w << std::endl;
        //for (auto ii = 0; ii < v.rows(); ii++) {
        //    for (auto jj = 0; jj < v.cols(); jj++) {
        //        std::cout << v.coeff(ii, jj) << " ";
        //    }
        //    std::cout << std::endl;
        //}

        //calc_polynomials(1, 1 + Chev::LN_2);
    }

    void Chev::iteration(bool full)
    {

    }

    void Chev::calc_A()
    {
        A_.resize(Chev::LN, Chev::LN);

        for (auto ix = 0; ix < Chev::NX; ix++) {
            for (auto iy = 0; iy < Chev::NY; iy++) {
                // A_.setdiag(-mu)
                auto const ii = xy2i(ix, iy);
                auto jx = ix;
                auto jy = iy;
                auto jj = xy2i(jx, jy);
                A_.coeffRef(ii, jj) = -Chev::MYU;

                // +1 in x direction
                jx = ix + 1;
                
                if (jx == Chev::NX) {
                    jx = 0;
                }

                jy = iy;
                jj = xy2i(jx, jy);
                A_.coeffRef(ii, jj) = -1.0;

                // -1 in x direction
                jx = ix - 1;

                if (jx == -1) {
                    jx = Chev::NX - 1;
                }

                jy = iy;
                jj = xy2i(jx, jy);
                
                A_.coeffRef(ii, jj) = -1.0;

                // + 1 in y direction
                jx = ix;
                jy = iy + 1;
                
                if (jy == Chev::NY) {
                    jy = 0;
                }
                
                jj = xy2i(jx, jy);
                A_.coeffRef(ii, jj) = -1.0;

                // -1 in y direction
                jx = ix;
                jy = iy - 1;
                if (jy == -1) {
                    jy = Chev::NY - 1;
                }
                jj = xy2i(jx, jy);
                A_.coeffRef(ii, jj) = -1.0;

                for (auto i = 0; i < Chev::LN_2; i++) {
                    for (auto j = 0; j < Chev::LN_2; j++) {
                        A_.coeffRef(i + Chev::LN_2, j + Chev::LN_2) = -A_.coeff(i, j);
                        A_.coeffRef(i, j + Chev::LN_2) = vec_delta_.coeff(i, j);
                        A_.coeffRef(i + Chev::LN_2, j) = vec_delta_.coeff(j, i);
                    }
                }
            }
        }

        A_ /= Chev::AA;
    }

    template <>
    void Chev::calc_meanfields<true>()
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(A_ * Chev::AA);
        auto const w = saes.eigenvalues();
        auto const v = saes.eigenvectors();

        for (auto ix = 0; ix < Chev::NX; ix++) {
            for (auto iy = 0; iy < Chev::NY; iy++) {
                auto const ii = xy2i(ix, iy);
                auto const jj = ii + Chev::LN_2;
                
                auto delta = 0.0;
                for (auto i = 0; i < Chev::LN; i++) {
                    if (w[i] <= 0.0) {
                        if (std::abs(w[i]) <= Chev::OMEGAC) {
                            delta += v.coeff(ii, i) * v.coeff(jj, i);
                            vec_delta_.coeffRef(ii, ii) = delta;
                        }
                    }
                }
            }
        }
    }

    template <>
    void Chev::calc_meanfields<false>()
    {
        for (auto ix = 0; ix < Chev::NX; ix++) {
            for (auto iy = 0; iy < Chev::NY; iy++) {
                auto ii = xy2i(ix, iy);
                auto jj = ii + Chev::LN_2;
                auto const right_j = jj;
                auto const left_i = ii;

                calc_polynomials(left_i, right_j);

                vec_delta_.coeffRef(ii, ii) = calc_meanfield();
            }
        }
    }
    
    double Chev::calc_meanfield() const
    {
        auto const ba = std::acos(-Chev::BB / Chev::AA);
        auto const omeb = std::acos(-(Chev::OMEGAC + Chev::BB) / Chev::AA);

        auto density = 0.0;
        for (auto j = 0; j < Chev::NC - 1; j++) {
            auto i = j + 1;
            density += vec_ai_[i] * (std::sin(static_cast<double>(i) * omeb) - std::sin(static_cast<double>(i) * ba)) / static_cast<double>(i);
            density += vec_ai_[0] * (omeb - ba) / 2.0;
        }

        return density * 2.0 / Chev::PI;
    }

    void Chev::calc_polynomials(std::int32_t left_i, std::int32_t right_j)
    {
        Eigen::VectorXd vec_jn(Chev::LN), vec_jnm(Chev::LN), vec_jnmm(Chev::LN);
        vec_jn.fill(0.0);
        vec_jnm.fill(0.0);
        vec_jnmm.fill(0.0);

        vec_jn(right_j) = 1.0;
        vec_ai_.fill(0.0);

        for (auto n = 0; n < Chev::NC; n++) {
            auto v = 0.0;
            switch (n) {
            case 0:
                vec_jnm.resize(1);
                vec_jnm.fill(0.0);
                vec_jnmm.resize(1);
                vec_jnmm.fill(0.0);
                vec_jn[right_j] = 1.0;
                break;

            case 1:
                vec_jn = A_ * vec_jn;
                break;

            default:
                vec_jn = 2.0 * A_ * vec_jnm - vec_jnmm;
                break;
            }

            vec_ai_[n] = vec_jn[left_i];
            vec_jnmm = vec_jnm;
            vec_jnm = vec_jn;
        }
    }

    void Chev::init_delta()
    {
        std::vector<Eigen::Triplet<double>> A(Chev::LN_2);
        for (auto i = 0; i < Chev::LN_2; i++) {
            A[i] = Eigen::Triplet<double>(i, i, Chev::DELTA);
        }

        vec_delta_.setFromTriplets(A.begin(), A.end());
    }

    std::int32_t Chev::xy2i(std::int32_t ix, std::int32_t iy) const
    {
        return iy * Chev::NX + ix;
    }

    template <> void Chev::calc_meanfields<true>();
    template <> void Chev::calc_meanfields<false>();
}
