#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h> 
#include <random>
#include <cstdio>	
#include <chrono>  // for high_resolution_clock
#include "/home/amir-prof/Dropbox/Research/Nima/Codes/Eigen-header/Eigen"
#include "/home/amir-prof/Dropbox/Research/Nima/Codes/Eigen-header/Dense"

using namespace std; // this is a using directive telling the compiler to check the std namespace when resolving identifiers with no prefix
using namespace Eigen;

const int nodes              = 500;                              //Nodes Number of Network
const int temp_iteration     = 1;                               //Number of Temprature Iteration
const int ensemble_iteration = 1;                               //Number of Enemble Iteration
const int iteration          = pow (nodes , 3);                  //Number of Iterations


int main () {
	FILE *ENERGY, *TRIANGLE, *NETWORK;
	std::random_device rd;                                     // Only used once to initialise (seed) engine
	std::mt19937 rng(rd());                                    // Random-number engine used (Mersenne-Twister in this case)
	std::uniform_int_distribution<int> uni1(1,1000);           // Guaranteed unbiased
	std::uniform_int_distribution<int> uni2(0, nodes*nodes-1);
	std::uniform_real_distribution<> rand(0, 1);
	MatrixXd adj(nodes,nodes), adj3(nodes,nodes);
	ENERGY   = fopen ("E-T.txt"                  ,"w");
	TRIANGLE = fopen ("r-T.txt"                  ,"w");
	NETWORK  = fopen ("Adjacency-Matrix.txt"     ,"w");
//Varibles
	int    initial_edge, final_edge, delta_E, energy_normal;
	int    m, n, triangle_num, edge_num, edge_energy, s, w;
	double alpha, theta, T, triangles, energy, r;
//Some prints
	cout << "\033[1;31m Number of Nodes = \033[0m"            << nodes          << "\n";
	cout << "\033[1;32m Number of iteration = \033[0m"        << iteration      << "\n";
	cout << "\033[1;33m Number of thermal iteration = \033[0m"<< temp_iteration << "\n";
//------------------------- Hamiltonian Parameters -------------------//   
    alpha = 6./nodes;
    theta = 2.3;	
//------------------------- Number of Configuration -------------------// 
    edge_num      = nodes*(nodes-1)/2;  
	triangle_num  = nodes*(nodes-1)*(nodes-2)/6;
	energy_normal = theta*edge_num - alpha*triangle_num ;
//--------------------------------------------------------------------//
//------------------------- Thermal Loop -----------------------------//

	for ( int i = 0; i < temp_iteration ;i++ ) {
// Record start time
	    auto start = std::chrono::high_resolution_clock::now();
		T = 0;
		cout << "\033[1;34m Temprature is =  \033[0m"<< T << "\n";
		cout << "\033[1;33m ****************************************************** \033[0m" << "\n";  
//-----------------------------Initial Configuration------------------//   
		for ( int i = 0; i < nodes ;i++ ) {
			adj(i,i) = 0;
			for ( int j = 0; j < i ;j++ ) {
				adj(i,j) = int(rand(rng)*2);
				adj(j,i) = adj(i,j);
				//cout << adj(i,j) << "\n";  
			} 
		}
//--------------------------------------------------------------------//

//----------- Initial Value for Energy, Edge Mean and Triangles -----//
		adj3      = adj*adj*adj;
		triangles = adj3.trace()/6;
		energy    = theta*adj.sum()/2 - alpha*triangles;
		cout << "Initial Energy = " <<  energy/energy_normal  << " \n";
//-------------------------------------------------------------------//

//---------------------------Monte Carlo loop------------------------//
		for ( int i = 1; i <= iteration ;i++ ) {
			s = uni2(rng);
			w = uni2(rng);
			m = floor (s / nodes) ;
			n = w % nodes;
			if (m!=n) {
				initial_edge = adj(n,m);
				if (initial_edge == 0) {
					final_edge = 1;
				}	
				else {
					final_edge = 0;
				}	
			//cout << "aa= " <<  final_edge  << " \n";
//---------------------------Energy of One Edge----------------------//
				edge_energy = 0;
				for ( int k = 0; k < nodes ;k++ ) { 
					edge_energy += adj(m,k)*adj(n,k);
				}             
				delta_E =  (final_edge - initial_edge)*(theta - alpha*edge_energy);
//-------------------------------------------------------------------//

//--------------------------Accept or Reject--------------------------//
				if (delta_E < 0) {
					adj(m,n) = final_edge;
					adj(n,m) = adj(m,n);
				}
				//else {
					//r = rand(rng);
					//if ( r < exp(-delta_E/T) ) {
						//adj(m,n) = final_edge;
						//adj(n,m) = adj(m,n);
					//}
				//}
			}
//-------------------------------------------------------------------//
		}	
//---------------------End of Iteration Loop------------------------//
	    adj3        = adj*adj*adj;
		triangles   = adj3.trace()/6;
		energy      = theta*adj.sum()/2 - alpha*triangles;
	    // Record end time
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		cout << "Elapsed time           : " << elapsed.count()        << " s\n";
		cout << "Final Triangles        = " << triangles/triangle_num << " \n";
		cout << "Final Energy Iteration = " << energy/energy_normal   << " \n";
		//write in files
		fprintf (ENERGY,     "%f \t %f \n", T, energy     /energy_normal);    
		fprintf (TRIANGLE,   "%f \t %f \n", T, triangles  /triangle_num ); 
    }
	for (int i=0; i<nodes; i++){
		for (int j=0; j<nodes; j++){
			if (j==nodes-1){	
				fprintf (NETWORK, "%d", int(adj(i,j)) );
			}	
			else {
				fprintf (NETWORK, "%d,", int(adj(i,j)) );
			}		
		}			
		fprintf (NETWORK,"\n");
	}
//---------------------End of Thermal Loop---------------------------//
return 0;
}
