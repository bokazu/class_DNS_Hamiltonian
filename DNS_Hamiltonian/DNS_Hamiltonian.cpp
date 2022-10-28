#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <random>

#include "DNS_Hamiltonian.hpp"
#include "cblas.h"

using namespace std;

//コピーコンストラクタ
DNS_Hamiltonian::DNS_Hamiltonian(const DNS_Hamiltonian& h)
    : jset_filename(h.jset_filename),
      tot_site_num(h.tot_site_num),
      mat_dim(h.mat_dim),
      J(h.J),
      Eig(h.Eig)
{
    H = new double[mat_dim * mat_dim];
    cblas_dcopy(mat_dim * mat_dim, h.H, 1, H, 1);
    std::cout << "copy has done.\n";
}

//代入演算子
DNS_Hamiltonian& DNS_Hamiltonian::operator=(const DNS_Hamiltonian& h)
{
    cblas_dcopy(mat_dim * mat_dim, h.H, 1, H, 1);
    return *this;
}

//行列の初期化
void DNS_Hamiltonian::init()
{
    for (int i = 0; i < mat_dim * mat_dim; i++) H[i] = 0.0;
}

//ハミルトニアン行列を求める
void DNS_Hamiltonian::hamiltonian()
{
    std::cout << "hamiltonian()\n";
    for (int m = 0; m < mat_dim; m++)
    {
        for (int site_i = 0; site_i < tot_site_num; site_i++)
        {
            spin(m, site_i);
        }
    }
}

//スピン演算子と状態ベクトルの計算を行う
void DNS_Hamiltonian::spin(int m, int site_i)
{
    boost::dynamic_bitset<> ket_m(tot_site_num, m);  // LH->1N_4
    bool bit_check0, bit_check1;
    int j_line = 0;
    for (int l = j_line; l < J.get_line(); l++)
    {
        if (J.index(0, l) == site_i)
        {
            int site_j = J.index(1, l);
            bit_check0 = ket_m.test(site_i);  // LH->1N->2Y_1
            bit_check1 = ket_m.test(site_j);  // LH->1N->2Y_1

            // LH->1N->3 (注目している2サイトについての場合分け)
            if (bit_check0 == bit_check1)
            {
                // S^z_{i}S^z_{i+1}|1_{i+1} 1_{i}> or S^z_{i}S^z_{i+1}|0_{i+1}
                // 0_{i}>
                H[m + m * mat_dim] += 0.25 * J.val(l);
            }
            else
            {
                boost::dynamic_bitset<> ket_m1(tot_site_num, m);  // LH->1N_4
                // LH->1N->3Y->1
                ket_m1.flip(site_j);        // LH->1N->3Y->1Y_1
                ket_m1.flip(site_i);        // LH->1N->3Y->1Y_1
                int n = ket_m1.to_ulong();  // LH->1N->3N_1

                // S^-_{i}S^+_{i+1} or S^+_{i}S^-_{i+1}
                H[m + n * mat_dim] += 0.5 * J.val(l);
                // S^z_{i}S^z_{i+1}|0_{i+1} 1_{i}>
                H[m + m * mat_dim] -= 0.25 * J.val(l);
            }
            j_line++;
        }
    }
}

//ベクトル規格化用の関数(非メンバ関数)
void sdz(int dim, double* vec){
    double a = 1./cblas_dnrm2(dim, vec, 1);
    cblas_dscal(dim, a, vec, 1);
}


//Hamiltonianの基底状態の固有値と固有ベクトルを計算する
void DNS_Hamiltonian::dns_lanczos(int tri_mat_dim , char c, char info_ls){
    int ls_count = 0; //lanczos法で要した反復回数を記録する
    double eps = 1.0; //step間の誤差を代入するようの変数
    double err = 1.0e-16; //要求精度
    bool err_checker = true; //誤差が要求精度の範囲内に収まっているかを確認するためのflag

    if(c == 'V') Eig.evec_elem(mat_dim);

    /*==============初期ベクトルの用意==============*/
    double **u = new double *[2];
    for(int k = 0;k<2;k++){
        u[k] = new double[mat_dim];
    }

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 1);
    for(int k = 0; k < mat_dim; k++){
        u[0][k] = rand1(mt);
        u[1][k] = 0.0;
    }
    sdz(mat_dim,u[0]);

    if(c == 'V') cblas_dcopy(mat_dim, u[0], 1, Eig.data(),1);
    /*============================================*/
    //三重対角行列の主対角成分
    double *alpha = new double[tri_mat_dim];
    vec_init(tri_mat_dim, alpha);
    
    //三重対角行列の次対角成分
    double *beta = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1, beta);

    //ls = 偶数stepでの近似固有値
    double *eval_even = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_even);

    //ls = 奇数stepでの近似固有値
    double *eval_odd = new double[tri_mat_dim];
    vec_init(tri_mat_dim, eval_odd);

    //LAPACKに三重対角行列の主対角成分を渡す用の配列
    double *diag = new double[tri_mat_dim];
    vec_init(tri_mat_dim, diag);

    //LAPACKに三重対角行列の主対角成分を渡す用の配列
    double *sub_diag = new double[tri_mat_dim - 1];
    vec_init(tri_mat_dim - 1, sub_diag);

    //LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
    double *tri_diag_evec;
    if(c == 'V') //固有ベクトルを計算する場合は配列を確保する
    {
        tri_diag_evec = new double[tri_mat_dim];
        vec_init(tri_mat_dim, tri_diag_evec);
    }

    /*==============================lanczos Algorithm=======================================*/
    for(int ls = 0; ls < tri_mat_dim; ls++){
        if(err_checker)
        {
            ls_count = ls;
            /*省メモリのためのlanczosベクトルの更新*/
            if(ls > 0)
            {
                if(ls % 2 == 0) cblas_dscal(mat_dim, -beta[ls - 1], u[1], 1);
                else cblas_dscal(mat_dim , -beta[ls - 1] , u[0], 1);
            }

            if(ls % 2 == 0)
            {
                //ls = 偶数step
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim, mat_dim, 1.0, H, mat_dim, u[0], 1 , 1.0, u[1], 1);
                alpha[ls] = cblas_ddot(mat_dim, u[1], 1, u[0], 1);
                cblas_daxpy(mat_dim, -alpha[ls], u[0], 1, u[1], 1);
                if(ls != tri_mat_dim - 1)
                {
                    beta[ls] = cblas_dnrm2(mat_dim, u[1],1);
                    cblas_dscal(mat_dim, 1./beta[ls], u[1],1);
                }
            }
            else{
                //ls = 奇数step
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim, mat_dim, 1.0, H, mat_dim, u[1], 1 , 1.0 , u[0], 1);
                alpha[ls] = cblas_ddot(mat_dim, u[1], 1, u[0], 1);
                cblas_daxpy(mat_dim,-alpha[ls],u[1],1,u[0],1);
                if(ls != tri_mat_dim - 1)
                {
                    beta[ls] = cblas_dnrm2(mat_dim, u[0], 1);
                    cblas_dscal(mat_dim, 1./beta[ls], u[0], 1);
                }
            }
            /*===================================================================================*/

            /*===========================三重対角行列の数値対角(LAPACK)=============================*/
            vec_init(tri_mat_dim, diag);
            vec_init(tri_mat_dim - 1 , sub_diag);
            int info = 0;
             if(ls % 2 == 0){
                 //偶数step
                 cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
                 cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

                 if (ls < tri_mat_dim - 1)
                 {
                     sub_diag[ls] = 0.;
                     if(c == 'N')
                     {
                         //固有値のみを計算する場合
                         info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                           sub_diag, tri_diag_evec, ls + 1);
                     }
                     else
                     {//固有ベクトルも計算する場合
                         info =
                     LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                                           sub_diag, tri_diag_evec, ls + 1);
                     }
                     
                     if (info != 0)
                     {
                         std::cout << "@ls = " << ls
                                       << " , LAPACKE_detev errored." << std::endl;
                     }
                 }
                 else
                 {
                     if(c == 'N'){
                         //固有値のみを計算する場合
                         info =LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                                           sub_diag, tri_diag_evec, ls + 1);
                     }
                     else{
                         //固有ベクトルを計算する場合
                         info =LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                                           sub_diag, tri_diag_evec, ls + 1);
                     }
                     
                     if (info != 0)
                     {
                        std::cout << "@ls = " << ls
                                  << " , LAPACKE_detev errored." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_even[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    std::cout << "@ls = " << ls << std::endl;
                }
            }
            else{
                cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
                cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

                if (ls < tri_mat_dim - 1)
                {
                    sub_diag[ls] = 0.;
                    if(c == 'N'){
                        //固有値のみを計算する場合
                        info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag, sub_diag, tri_diag_evec, ls + 1);
                    }
                    else{
                        //固有ベクトルのみを計算する場合
                        info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag, sub_diag, tri_diag_evec, ls + 1);
                    }

                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                else
                {
                    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag, sub_diag, tri_diag_evec, ls + 1);
                    if (info != 0)
                    {
                        std::cout << "@ls = " << ls << " , LAPACKE_detev's error." << std::endl;
                    }
                }
                cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
                if (info_ls == 'y')
                {
                    std::cout << "@ls = " << ls << " : eigen value = " << eval_odd[0] << std::endl;
                }
                else if (info_ls == 's')
                {
                    cout << "@ls = " << ls << endl;
                }
            }
            /*======================================================================*/

            /*============================収束状況の確認==============================*/
            if (ls > 0)
            {
                eps = abs(eval_even[0] - eval_odd[0]);
                if (info_ls == 'y')
                {
                    cout << "eps = " << std::setprecision(16) << eps << endl;
                }

                if (eps > err) err_checker = true;
                else err_checker = false;     
            }
            /*=====================================================================*/
        }
        else
        {
            cout << "break @" << ls_count << endl;
        }
    }
    /*========================基底状態の固有値===========================*/
    if(ls_count % 2 == 0) Eig.set_eval(eval_even[0]);
    else Eig.set_eval(eval_odd[0]);
    /*========================配列リソースのリリース part1===================*/
    delete[] eval_even;
    delete[] eval_odd;

    //固有ベクトルを計算する
    if(c == 'V'){
        vec_init(mat_dim, u[0]);
        vec_init(mat_dim, u[1]);
        cblas_dcopy(mat_dim, Eig.data(), 1 , u[0], 1);

        for(int ls = 0; ls < ls_count + 2; ls++)
        {
            if (ls % 2 == 0)
            {
                if (ls == 0) cblas_dscal(mat_dim, tri_diag_evec[ls], Eig.data(), 1);
                else cblas_daxpy(mat_dim, tri_diag_evec[ls], u[0], 1, Eig.data(), 1);
            }
            else
            {
                cblas_daxpy(mat_dim, tri_diag_evec[ls], u[1], 1, Eig.data(), 1);
            }

            if (ls > 0)
            {
                if (ls % 2 == 0) cblas_dscal(mat_dim, -beta[ls - 1], u[1], 1);
                else cblas_dscal(mat_dim, -beta[ls - 1], u[0], 1); 
            }

            if (ls % 2 == 0)
            {
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim, mat_dim, 1.0, H,
                        mat_dim, u[0], 1, 1.0, u[1], 1);
                cblas_daxpy(mat_dim, -alpha[ls], u[0], 1, u[1], 1);
                if (ls != tri_mat_dim - 1) cblas_dscal(mat_dim, 1. / beta[ls], u[1], 1);
            }
            else
            {
                cblas_dgemv(CblasColMajor, CblasNoTrans, mat_dim, mat_dim, 1.0, H,
                        mat_dim, u[1], 1, 1.0, u[0], 1);
                cblas_daxpy(mat_dim, -alpha[ls], u[1], 1, u[0], 1);
                if (ls != tri_mat_dim - 1) cblas_dscal(mat_dim, 1. / beta[ls], u[0], 1);
            }
        }

        sdz(mat_dim, Eig.data());
        
    }



    /*==========================配列リソースのリリース part2 ====================*/
    for (int i = 0; i < 2; i++)
    {
        delete[] u[i];
    }
    delete[] u;
    delete[] alpha;
    delete[] beta;

    delete[] diag;
    delete[] sub_diag;
    if(c == 'V') delete[] tri_diag_evec;
}   



// Hamiltonianの第[i][j]要素を標準出力する
void DNS_Hamiltonian::print(int i, int j) const
{
    cout << "H[" << i << "]"
         << "[" << j << "] = " << scientific << setprecision(print_precision)
         << H[i + j * mat_dim] << "\n";
}

// Hamiltonian行列を標準出力する
void DNS_Hamiltonian::print() const
{
    double mtmp;
    for (int row_num = 0; row_num < mat_dim; row_num++)
    {
        cout << "[";
        for (int col_num = 0; col_num < mat_dim; col_num++)
        {
            mtmp = H[col_num + mat_dim * row_num];
            // printf("%5.8e", mtmp);
            cout << scientific << setprecision(print_precision) << mtmp;
            if (col_num < mat_dim - 1)
            {
                cout << "  ";
            }
        }
        if (row_num < mat_dim - 1)
        {
            cout << "];" << endl;
        }
        else
        {
            cout << "]";
        }
    }
    cout << "]" << endl;
}
