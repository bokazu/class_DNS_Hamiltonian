#ifndef ___Class_Hamiltonian
#define ___Class_Hamiltonian

#include <iomanip>
#include <string>

#include "../Jset/Jset.hpp"
#include "../EIGEN/EIGEN.hpp"

class DNS_Hamiltonian
{
        private:
    std::string jset_filename;
    int tot_site_num;
    int mat_dim;
    double* H;
    Jset J;
    void set_J() { J.set(); }
    void spin(int m, int site_i);

        public:
        EIGEN Eig;

    //明示的コンストラクタ
    explicit DNS_Hamiltonian(std::string filename, int site)
        : jset_filename(filename),
          tot_site_num(site),
          mat_dim(1 << tot_site_num),
          J(filename),
          Eig(1)
    {
        //相互作用の情報をsetする
        set_J();
        //配列を動的確保し、0で初期化する
        H = new double[mat_dim * mat_dim];
        for (int i = 0; i < mat_dim * mat_dim; i++)
        {
            H[i] = 0.0;
        }
        std::cout << "DNS_Hamiltonian::constructed.\n";
    }

    //コピーコンストラクタ
    DNS_Hamiltonian(const DNS_Hamiltonian& h);
    //デストラクタ
    ~DNS_Hamiltonian()
    {
        delete[] H;
        std::cout << "DNS_Hamiltonian::destructed.\n";
    }

    // Hamiltonian行列の先頭要素へのポインタを返す
    double* data() { return H; }

    //系のサイト数を返す
    int site() const { return tot_site_num; }

    //行列の次元を返す
    int dim() const { return mat_dim; }

    // Hamiltonianの第(i,j)成分を取得する
    double at(int i, int j) { return H[j + mat_dim * i]; }

    // jsetの情報を書き込んだファイルの名前を返す
    std::string jsetfile() const { return jset_filename; }

    //代入演算子=
    DNS_Hamiltonian& operator=(const DNS_Hamiltonian& h);

    //添字演算子[]
    double& operator[](int i) { return H[i]; }

    // const版添字演算子
    const double& operator[](int i) const { return H[i]; }

    // Hamiltonianを0で初期化
    void init();

    // Hamiltonianの行列要素を計算する
    void hamiltonian();

    //Hamiltonian行列の基底状態の固有値を計算する
    void dns_lanczos(int tri_mat_dim , char c = 'N', char info_ls = 'n');

    //=================Hamiltonianの標準出力関係======================//
    int print_precision = 5;  //標準出力する際の行列要素の桁数を指定する
    void print(int i,
               int j) const;  // Hamiltonian行列の第(i,j)成分を標準出力する
    void print() const;  // Hamiltonian行列を行列の形で標準出力する
    //================================================================//

    // jset_filenameと系の相互作用の情報を出力する
    void jset_print() const { J.print(); }
};



void sdz(int mat_dim, double* vec);

template <typename T>
void vec_init(int dim, T *vec){
    for(int i = 0;i < dim;i++){
        vec[i] = 0;  
    }
}


#endif
