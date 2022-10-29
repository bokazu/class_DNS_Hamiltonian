#include <iostream>

#include "DNS_Hamiltonian/DNS_Hamiltonian.hpp"

using namespace std;

int main()
{
    DNS_Hamiltonian H1("jset1.txt", 2); //オブジェクトの構築(Hamiltonian行列は0で初期化される)

    H1.hamiltonian(); //ハミルトニアンの行列要素の計算と代入
   
    H1.dns_lanczos(4,'N'); //Lanczos法により固有値または固有ベクトルを計算する
   
    cout << H1 << endl; //オブジェクトの情報を出力する
   
}
