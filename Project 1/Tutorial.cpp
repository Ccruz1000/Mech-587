#include <iostream>

int main() {
double a = 0.3;
double * aptr;
aptr = &a;

std::cout << "a is " << a << std::endl;
std::cout << "The address of a is " << aptr << std::endl;

double* b;


return 0;
}
