#include "array.hpp"

Array::Array(int a, int b, int c, int d, int e, double v){

    A = a;
    B = b;
    C = c;
    D = d;
    E = e;

    dim = 5;

    array.assign(A*B*C*D*E,v);

}

Array::Array(int a, int b, int c, int d, double v){

    A = a;
    B = b;
    C = c;
    D = d;

    dim = 4;

    array.assign(A*B*C*D, v);

}

Array::Array(int a, int b, int c, double v){

    A = a;
    B = b;
    C = c;

    dim = 3;

    array.assign(A*B*C,v);

}

Array::Array(int a, int b, double v){

    A = a;
    B = b;

    dim = 2;

    array.assign(A*B,v);

}

Array::Array(int a, double v){

    A = a;

    dim = 1;

    array.assign(A,v);

}

Array::Array(){}


double& Array::operator()(int a, int b, int c, int d, int e){

    return array[E*(D*(C*(B*a + b) + c) + d) + e];

}

double& Array::operator()(int a, int b, int c, int d){

    return array[D*(C*(B*a + b) + c) + d];

}

double& Array::operator()(int a, int b, int c){

    return array[C*(B*a + b) + c];

}

double& Array::operator()(int a, int b){

    return array[B*a + b];

}

double& Array::operator()(int a){

    return array[a];

}

void Array::clear(){

    clearArray.swap(array);

}

//         double& operator()(int a, int b, int c, int d, int e){ // Overloads operator to allow easy indexing of 5d arrays

//             return array[E*(D*(C*(B*a + b) + c) + d) + e]; 

//         }

//         double& operator()(int a, int b, int c, int d){ // Same as above for 4d arrays

//             return array[D*(C*(B*a + b) + c) + d];

//         }

//         double& operator()(int a, int b, int c){ // Same as above for 3d arrays

//             return array[C*(B*a + b) + c];

//         }