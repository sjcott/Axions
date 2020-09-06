#include <vector>

class Array{

    public:

        Array(int,int,int,int,int,double); // Constructor for 5d arrays
        Array(int,int,int,int,double); // Constructor for 4d arrays
        Array(int,int,int,double); // Constructor for 3d arrays
        Array(int,int,double); //Constructor for 2d arrays
        Array(int,double); //Constructor for 1d arrays

        int A, B, C, D, E, dim; // Sizes of each dimension and #dimensions of the array (for indexing error messages)
        std::vector<double> array; //Empty vector container for the array

        double& operator()(int a, int b, int c, int d, int e); // Overloads operator to allow easy indexing of 5d arrays

        double& operator()(int a, int b, int c, int d); // Same as above for 4d arrays

        double& operator()(int a, int b, int c); // Same as above for 3d arrays

        double& operator()(int a, int b);

        double& operator()(int a);

};