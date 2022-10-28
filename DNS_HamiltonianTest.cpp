#include <iostream>

#include "DNS_Hamiltonian/DNS_Hamiltonian.hpp"

using namespace std;

int main()
{
    DNS_Hamiltonian H1("jset1.txt", 2);
    std::cout << "H1's site = " << H1.site() << std::endl;
    std::cout << "H1's dim = " << H1.dim() << std::endl;
    std::cout << "H1's jset_filename = " << H1.jsetfile() << std::endl;
    std::cout << "\nH1's jset information\n";
    H1.jset_print();

    H1.hamiltonian();
    H1.print();
    std::cout << std::endl;

    H1.dns_lanczos(4,'V');
    cout << H1.Eig << endl;
    // for (int i = 0; i < H1.dim(); i++)
    // {
    //     for (int j = 0; j < H1.dim(); j++)
    //     {
    //         H1.print(i, j);
    //     }
    // }

    // std::cout << std::endl;
    // H1.print_precision = 3;  //デフォルトは5ケタ
    // H1.print();
    // std::cout << std::endl;
    // for (int i = 0; i < H1.dim(); i++)
    // {
    //     for (int j = 0; j < H1.dim(); j++)
    //     {
    //         H1.print(i, j);
    //     }
    // }

    // std::cout << std::endl;
    // H1.init();
    // H1.print();

    // H1.Eig.evec_elem(H1.dim());
}
