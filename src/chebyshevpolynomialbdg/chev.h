#include <cstdint>          // for std::int32_t
#include <Eigen/Dense>      // for Eigen::SparseMatrix
#include <Eigen/Sparse>     // for Eigen::VectorXd

namespace chebyshevpolynomialbdg {
    // #region �t���[�֐��̐錾

    template <typename T>
    constexpr T sqr(T x);

    // #endregion �t���[�֐��̐錾 

    //! A class.
    /*!
        Bogoliubov-de Gennes�������������N���X
    */
    class Chev final {
        // #region �R���X�g���N�^�E�f�X�g���N�^

    public:
        //! A constructor.
        /*!
        */
        Chev();

        //! A destructor.
        /*!
            �f�t�H���g�f�X�g���N�^
        */
        ~Chev() = default;

        // #endregion �R���X�g���N�^�E�f�X�g���N�^

        // #region public�����o�֐�

        //! A public member function.
        /*!
        */
        void iteration(bool full);

        // #endregion public�����o�֐�
        
        // #region private�����o�֐�

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

        // #endregion private�����o�֐�

        // #region �����o�ϐ�

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
        static auto constexpr ITERMAX = 100;

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

        // #region �֎~���ꂽ�R���X�g���N�^�E�����o�֐�

        private:

       //! A private copy constructor (deleted).
       /*!
           �R�s�[�R���X�g���N�^�i�֎~�j
           \param �R�s�[���̃I�u�W�F�N�g�i���g�p�j
       */
       Chev(Chev const &) = delete;

       //! A private member function (deleted).
       /*!
           operator=()�̐錾�i�֎~�j
           \param �R�s�[���̃I�u�W�F�N�g�i���g�p�j
           \return �R�s�[���̃I�u�W�F�N�g
       */
       Chev & operator=(Chev const &) = delete;

       // #endregion �֎~���ꂽ�R���X�g���N�^�E�����o�֐�
    };

    // #region template�����o�֐��̎���

    template <>
    inline void Chev::calc_meanfields<true>()
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

    // #endregion template�����o�֐��̎���

    // #region template�t���[�֐��̎���

    template <typename T>
    constexpr T sqr(T x)
    {
        return x * x;
    }

    // #endregion template�t���[�֐��̎���
}
