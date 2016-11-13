/*! \file chev.h
    \brief Bogoliubov-de Gennes方程式を解くクラスの宣言
    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _CHEV_H_
#define _CHEV_H_

#pragma once

#include <cstdint>          // for std::int32_t
#include <Eigen/Dense>      // for Eigen::SparseMatrix
#include <Eigen/Sparse>     // for Eigen::VectorXd

#undef EIGEN_NO_DEBUG
#undef NDEBUG

namespace chebyshevpolynomialbdg {
    // #region フリー関数の宣言

    template <typename T>
    constexpr T sqr(T x);

    // #endregion フリー関数の宣言 

    //! A class.
    /*!
        Bogoliubov-de Gennes方程式を解くクラス
    */
    class Chev final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
        */
        Chev();

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Chev() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
        */
        void iteration(bool full);

        // #endregion publicメンバ関数
        
        // #region privateメンバ関数

    private:

        //! A private member function.
        /*!
        */
        void calc_A();

        //! A private member function.
        /*!
        */
        void calc_A2();

        //! A private member function.
        /*!
        */
        double calc_meanfield() const;

        //! A private member function.
        /*!
        */
        template <bool F>
        void calc_meanfields();

        //! A private member function.
        /*!
        */
        void calc_polynomials(std::int32_t left_i, std::int32_t right_j);

        //! A private member function.
        /*!
        */
        void init_delta();

        //! A private member function.
        /*!
        */
        std::int32_t xy2i(std::int32_t ix, std::int32_t iy) const;

        // #endregion privateメンバ関数

        // #region メンバ変数

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr AA = 10.0;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr BB = 0.0;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr DELTA = 0.1;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr EPSTHRESHOLD = 1.0E-6;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr ITERMAX = 2;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr MYU = 1.0E-12;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NC = 1000;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NX = 20;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NY = 20;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr LN = Chev::NX * Chev::NY * 2;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr LN_2 = Chev::NX * Chev::NY;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr OMEGAC = 10.0;
        
        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr PI = 3.14159265359;
        
        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr U = -2.0;

        //! A private member variable.
        /*!
        */
        Eigen::SparseMatrix<double> A_;

        //! A private member variable.
        /*!
        */
        Eigen::VectorXd vec_ai_;

        //! A private member variable.
        /*!
        */
        Eigen::SparseMatrix<double> vec_delta_;

        // #region 禁止されたコンストラクタ・メンバ関数

        private:

       //! A private copy constructor (deleted).
       /*!
           コピーコンストラクタ（禁止）
           \param コピー元のオブジェクト（未使用）
       */
       Chev(Chev const &) = delete;

       //! A private member function (deleted).
       /*!
           operator=()の宣言（禁止）
           \param コピー元のオブジェクト（未使用）
           \return コピー元のオブジェクト
       */
       Chev & operator=(Chev const &) = delete;

       // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    // #region templateメンバ関数の実装

    template <>
    inline void Chev::calc_meanfields<true>()
    {
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double, Eigen::RowMajor>> saes(A_ * Chev::AA);
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
    inline void Chev::calc_meanfields<false>()
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

    // #endregion templateメンバ関数の実装

    // #region templateフリー関数の実装

    template <typename T>
    constexpr T sqr(T x)
    {
        return x * x;
    }

    // #endregion templateフリー関数の実装
}

#endif  // _CHEV_H_
