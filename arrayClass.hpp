#include <vector>

class Array{

    public:

        Array(int,int,int,int,int,double); // Constructor for 5d arrays
        Array(int,int,int,int,double); // Constructor for 4d arrays
        Array(int,int,int,double); // Constructor for 3d arrays

        int A, B, C, D, E; // Sizes of each dimension
        std::vector<double> array; //Empty vector container for the array

        double& operator()(int a, int b, int c, int d, int e); // Overloads operator to allow easy indexing of 5d arrays

        double& operator()(int a, int b, int c, int d); // Same as above for 4d arrays

        double& operator()(int a, int b, int c); // Same as above for 3d arrays

};