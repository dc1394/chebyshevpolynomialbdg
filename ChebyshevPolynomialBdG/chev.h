#include <cstdint>          // for std::int32_t
#include <Eigen/Dense>      // for Eigen::SparseMatrix
#include <Eigen/Sparse>     // for Eigen::VectorXd

namespace chebyshevpolynomialbdg {
    class Chev final {
        // #region �R���X�g���N�^�E�f�X�g���N�^

    public:
        //! A constructor.
        /*!
            �B��̃R���X�g���N�^
            \param func_ ��ϕ��֐�
            \param n simpson�̌����̕�����
            \param x1 �ϕ��̉���
            \param x2 �ϕ��̏��
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
        static auto constexpr MYU = 1.0E-12;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NC = 10;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NX = 2;

        //! A private static member variable (constant expression).
        /*!
        */
        static auto constexpr NY = 2;

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
    };
}