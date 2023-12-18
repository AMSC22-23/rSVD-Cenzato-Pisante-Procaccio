#include "QR_Decomposition.hpp"


int main(){
size_t m=60;
size_t n=60;



Matrix A(m,n);
for(size_t i=0;i<m;i++){
    A(i,i)=4;
    for (size_t j=0;j<i;j++){
        A(i,j)=A(i,i)-i+j;
    }
    for (size_t j=i+1;j<n;j++){
        A(i,j)=A(i,i)+i-j;
    }
}

    QR_Decomposition QR_A;

    /**
     * Serial execution with Givens
    */
    
    auto start_givens = std::chrono::high_resolution_clock::now();
    auto [Qg,Rg]=QR_A.Givens_solve(A);
    auto end_givens = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_g = end_givens  - start_givens;
    /*
    std::cout<<"Serial Givens:"<<std::endl;
    std::cout<<"R="<<std::endl;
    #ifdef EIGEN
        std::cout<<Rg<<std::endl;
    #else
        Rg.print(std::cout);
    #endif
    std::cout<<"Q="<<std::endl;
    #ifdef EIGEN
        std::cout<<Qg<<std::endl;
    #else
        Qg.print(std::cout);
    #endif
    */
    /**
     * Serial execution with HouseHolder
    */

   
    auto start_serial = std::chrono::high_resolution_clock::now();
    auto [Q,R]=QR_A.HouseHolder_solve(A);
    auto end_serial = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_s = end_serial - start_serial;

    /*
    std::cout<<"Serial HouseHolder:"<<std::endl;
    std::cout<<"R="<<std::endl;
    #ifdef EIGEN
        std::cout<<R<<std::endl;
    #else
        R.print(std::cout);
    #endif
    std::cout<<"Q="<<std::endl;
    #ifdef EIGEN
        std::cout<<Q<<std::endl;
    #else
        Q.print(std::cout);
    #endif
    */
    

   

    /**
     * Serial execution with Givens
    */

    auto start_givensp = std::chrono::high_resolution_clock::now();
    auto [Qgp,Rgp]=QR_A.Givens_solve_parallel(A);
    auto end_givensp = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_gp = end_givensp  - start_givensp;

    /*
    std::cout<<std::endl;
    std::cout<<"Parallel Givens:"<<std::endl;
    std::cout<<"R="<<std::endl;
    #ifdef EIGEN
        std::cout<<Rgp<<std::endl;
    #else
        Rgp.print(std::cout);
    #endif
    std::cout<<"Q="<<std::endl;
    #ifdef EIGEN
        std::cout<<Qgp<<std::endl;
    #else
        Qgp.print(std::cout);
    #endif
    */
    
    /**
     * Parallel execution of HouseHolder on OpenMP
    */

    std::cout<<std::endl;
    
   
    auto start_parallel = std::chrono::high_resolution_clock::now();
    auto [Qp,Rp]=QR_A.HouseHolder_solve_parallel(A);
    auto end_parallel = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_p = end_parallel - start_parallel;
       
    /*
    std::cout<<"Parallel HouseHolder:"<<std::endl;
    std::cout<<"R="<<std::endl;
    #ifdef EIGEN
        std::cout<<Rp<<std::endl;
    #else
        Rp.print(std::cout);
    #endif
    std::cout<<"Q="<<std::endl;
    #ifdef EIGEN
        std::cout<<Qp<<std::endl;
    #else
        Qp.print(std::cout);
    #endif
    */

    std::cout << "Time of execution serial givens: " << duration_g.count() << " s" << std::endl;
    std::cout << "Time of execution serial HouseHolder: " << duration_s.count() << " s" << std::endl;
    std::cout << "Time of execution parallel givens: " << duration_gp.count() << " s" << std::endl;
    std::cout << "Time of execution parallel: " << duration_p.count() << " s" << std::endl;

    double SpeedUpG=duration_g.count()/duration_gp.count();
    double SpeedUpHH=duration_s.count()/duration_p.count();
    std::cout << "Speed Up Givens: "<<SpeedUpG << std::endl;
    std::cout << "Speed Up HouseHolder: "<<SpeedUpHH << std::endl;

	return 0;
}
