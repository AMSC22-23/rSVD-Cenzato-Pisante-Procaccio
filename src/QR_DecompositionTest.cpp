#include "QR_Decomposition.hpp"


int main(){
size_t m=4;
size_t n=4;


    // Assegnazione di valori agli elementi
 // Test matrix
/* 
 Matrix A(4,4);
    A<< 4,3,2,1,
3,4,3,2,
2,3,4,3,
1,2,3,4;
*/  
 /*
 Matrix A(4,4);
 A<<121.0991,   15.0583,    1.3624,   -0.3454,
   15.0583,   13.4599,   1.7818,  -0.6114,
    1.3624,    1.7818,   0.8663,   -0.3217,
   -0.3454,   -0.6114,   -0.3217,    0.5747;  
   */
/* 
Matrix A(5,3);
A<<0.8147, 0.0975, 0.1576,
0.9058, 0.2785, 0.9706,
0.1270, 0.5469, 0.9572,
0.9134, 0.9575, 0.4854,
0.6324, 0.9649, 0.8003;
*/ 

/*
Matrix A(3,3);
    A<<12, -51, 4,
6, 167, -68,
-4, 24, -41;
*/
/*
    Matrix A(5,5);
    A << 2,-1,0,0,-1.e-12,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        1.e-12,-1.e-12,-1,2,-1,
        1.e-12,1.e-12,1.e-12,-1,2;
 */   
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
std::cout<<A<<std::endl;

    QR_Decomposition QR_A;
    /**
     * Serial execution with Givens
    */
    std::cout<<"Serial:"<<std::endl;
    auto start_givens = std::chrono::high_resolution_clock::now();
    auto [Qg,Rg]=QR_A.Givens_solve(A);
    auto end_givens = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_g = end_givens  - start_givens;
    std::cout<<"R Givens="<<std::endl;
    #ifdef EIGEN
        std::cout<<Rg<<std::endl;
    #else
        Rg.print(std::cout);
    #endif
    std::cout<<"Q Givens="<<std::endl;
    #ifdef EIGEN
        std::cout<<Qg<<std::endl;
    #else
        Qg.print(std::cout);
    #endif

    /**
     * Serial execution with HouseHolder
    */
    std::cout<<"Serial:"<<std::endl;
    auto start_serial = std::chrono::high_resolution_clock::now();
    auto [Q,R]=QR_A.HouseHolder_solve(A);
    auto end_serial = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_s = end_serial - start_serial;
    std::cout<<"R_Serial="<<std::endl;
    #ifdef EIGEN
        std::cout<<R<<std::endl;
    #else
        R.print(std::cout);
    #endif
    std::cout<<"Q_Serial="<<std::endl;
    #ifdef EIGEN
        std::cout<<Q<<std::endl;
    #else
        Q.print(std::cout);
    #endif

    

    std::cout<<std::endl;
    std::cout<<"Parallel:"<<std::endl;
    
    /**
     * Parallel execution on OpenMP
    */
   
    auto start_parallel = std::chrono::high_resolution_clock::now();
    auto [Qp,Rp]=QR_A.HouseHolder_solve_parallel(A);
    auto end_parallel = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_p = end_parallel - start_parallel;
       
    
    std::cout<<"R_parallel="<<std::endl;
    #ifdef EIGEN
        std::cout<<Rp<<std::endl;
    #else
        Rp.print(std::cout);
    #endif
    std::cout<<"Q_parallel="<<std::endl;
    #ifdef EIGEN
        std::cout<<Qp<<std::endl;
    #else
        Qp.print(std::cout);
    #endif
    
    std::cout << "Time of execution serial givens: " << duration_g.count() << " secondi" << std::endl;
    std::cout << "Time of execution serial: " << duration_s.count() << " secondi" << std::endl;
    std::cout << "Time of execution parallel: " << duration_p.count() << " s" << std::endl;
    double SpeedUp=duration_s.count()/duration_p.count();
    std::cout << "Speed Up: "<<SpeedUp << std::endl;
/* 
    auto [Qh,Rh]=QR_A.HouseHolder_solve(A);
    std::cout<<"R_Serial="<<std::endl;
    std::cout<<Rh<<std::endl;
    std::cout<<"Q_Serial="<<std::endl;
    std::cout<<Qh<<std::endl;
*/    


	return 0;
}
