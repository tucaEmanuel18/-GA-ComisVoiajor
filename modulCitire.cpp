#include <iostream>
#include <fstream>
#include <string>
using namespace std;


int main(int argc, char* argv[])
{

	if(argc != 2)
	{
		printf("Please get the name of test file!\n");
		exit(0);
	}
	else
	{
		printf("Name of testFile: %s\n", argv[1]);
		ifstream infile;
		infile.open(argv[1]);
		int nr;
		infile >> nr;
		const int dim = nr;
		printf("dim = %d\n", dim);
		int A[dim + 1][dim + 1];


		for(int i = 0; i < dim; i++)
		{
			for(int j = 0; j < dim; j++)
			{
				infile >> nr;
				A[i][j] = nr;
			}		
		}
		for(int i = 0; i < dim; i++)
		{
			for(int j = 0; j < dim; j++)
			{
				printf("%d ", A[i][j]);
			}
			printf("\n");
		}
	}
	
	return 0;
}