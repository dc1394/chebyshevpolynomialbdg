#include "chev.h"
#include <cmath>                // for std::acos, for std::sin
#include <vector>               // for std::vector
#include <iomanip>              // for std::setprecision
#include <iostream>             // for std::cout
#include <Eigen/EigenValues>    // for Eigen::SelfAdjointEigenSolver

namespace chebyshevpolynomialbdg {
    // #region コンストラクタ

    Chev::Chev()
        :   vec_ai_(Chev::NC),
            vec_delta_(Chev::LN_2, Chev::LN_2)
    {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        
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
            
            //for (auto ii = 0; ii < Chev::LN_2; ii++) {
            //    for (auto jj = 0; jj < Chev::LN_2; jj++) {
            //        std::cout << vec_delta_.coeff(ii, jj) << " ";
            //    }
            //    std::cout << std::endl;
            //}

            vec_delta_ = vec_delta_ * Chev::U;

            calc_A2();

            //for (auto ii = 0; ii < Chev::LN; ii++) {
            //    for (auto jj = 0; jj < Chev::LN; jj++) {
            //        std::cout << A_.coeff(ii, jj) << " ";
            //    }
            //    std::cout << std::endl;
            //}

            //auto i = 0;

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

    // #endregion privateメンバ関数
}
