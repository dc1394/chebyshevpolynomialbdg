/*! \file foelement.cpp
    \brief Bogoliubov-de Gennes方程式を解くクラスの実装
    Copyright ©  2016 @dc1394 All Rights Reserved.
    (but this is originally adapted by cometscome for Chev.py from https://github.com/cometscome/ChebyshevPolynomialBdG )
    This software is released under the BSD 2-Clause License.
*/

#include "chev.h"
#include <cmath>                // for std::acos, for std::sin
#include <iomanip>              // for std::setprecision
#include <iostream>             // for std::cout
#include <vector>               // for std::vector
#include <Eigen/EigenValues>    // for Eigen::SelfAdjointEigenSolver

namespace chebyshevpolynomialbdg {
    // #region コンストラクタ

    Chev::Chev()
        :   vec_ai_(Chev::NC),
            vec_delta_(Chev::LN_2, Chev::LN_2)
    {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    void Chev::iteration(bool full)
    {
        std::cout << std::setprecision(15);
        
        init_delta();
        calc_A();
        Eigen::SparseMatrix<double> vec_delta_old(vec_delta_);

        for (auto ite = 0; ite < ITERMAX; ite++) {
            if (full) {
                calc_meanfields<true>();
            }
            else {
                calc_meanfields<false>();
            }

            vec_delta_ = vec_delta_ * Chev::U;

            calc_A2();

            auto eps = 0.0;
            auto nor = 0.0;
            for (auto i = 0; i < Chev::LN_2; i++) {
                eps += sqr(vec_delta_.coeff(i, i) - vec_delta_old.coeff(i, i));
                nor += sqr(vec_delta_old.coeff(i, i));
            }

            eps /= nor;

            std::cout << "ite = " << ite << ", eps = " << eps << '\n';
            if (eps <= Chev::EPSTHRESHOLD) {
                std::cout << "End " << vec_delta_.coeff(Chev::NX / 2, Chev::NY / 2) << std::endl;
                break;
            }
            
            vec_delta_old = vec_delta_;
        }
    }

    // #region publicメンバ関数

    // #region privateメンバ関数

    void Chev::calc_A()
    {
        A_ = Eigen::SparseMatrix<double>(Chev::LN, Chev::LN);

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

    void Chev::calc_A2()
    {
        A_ *= Chev::AA;

        for (auto i = 0; i < Chev::LN_2; i++) {
            for (auto j = 0; j < Chev::LN_2; j++) {
                A_.coeffRef(i, j + Chev::LN_2) = vec_delta_.coeff(i, j);
                A_.coeffRef(i + Chev::LN_2, j) = vec_delta_.coeff(j, i);
            }
        }
        
        A_ /= Chev::AA;
    }
    
    double Chev::calc_meanfield() const
    {
        auto const ba = std::acos(-Chev::BB / Chev::AA);
        auto const omeb = std::acos(-(Chev::OMEGAC + Chev::BB) / Chev::AA);

        auto density = 0.0;
        for (auto j = 0; j < Chev::NC - 1; j++) {
            auto const i = j + 1;
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

        vec_jn.coeffRef(right_j) = 1.0;
        vec_ai_.fill(0.0);

        auto A = A_.toDense();
        for (auto n = 0; n < Chev::NC; n++) {
            switch (n) {
            case 0:
                vec_jnm.resize(1);
                vec_jnm.fill(0.0);
                vec_jnmm.resize(1);
                vec_jnmm.fill(0.0);
                vec_jn.coeffRef(right_j) = 1.0;
                break;

            case 1:
                vec_jn = A * vec_jn;
                break;

            default:
                vec_jn = 2.0 * A * vec_jnm - vec_jnmm;
                break;
            }

            vec_ai_[n] = vec_jn.coeff(left_i);
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

    // #endregion privateメンバ関数
}
