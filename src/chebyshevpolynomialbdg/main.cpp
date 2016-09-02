#include "chev.h"
#include <chrono>   // for std::chrono
#include <iomanip>  // for std::setprecision
#include <iostream> // for std::cout

int main()
{
    using namespace std::chrono;

    chebyshevpolynomialbdg::Chev chev;

    //auto const begtrue = high_resolution_clock::now();
    //chev.iteration(true);
    //auto const endtrue = high_resolution_clock::now();
    //
    //std::cout << std::setprecision(3);
    //std::cout << "elapsed_time = " << (duration_cast<duration<double>>(endtrue - begtrue)).count() << " (sec)\n";

    auto const begfalse = high_resolution_clock::now();
    chev.iteration(false);
    auto const endfalse = high_resolution_clock::now();
    
    std::cout << std::setprecision(3);
    std::cout << "elapsed_time = " << (duration_cast<duration<double>>(endfalse - begfalse)).count() << " (sec)" << std::endl;
    
    return 0;
}