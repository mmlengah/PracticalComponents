#include "../PracticalComponents/PCM.h"
#include <chrono>
#include <iostream>

using namespace PC;

int main() {
    Matrix4x4int a, b; // Assume these are initialized properly

    int amount = 10000000;
    // Time first function
    auto start1 = std::chrono::high_resolution_clock::now();
    Matrix4x4int res1 = a;
    for (int i = 0; i < amount; i++) {
        res1 = res1 * b;
    }    
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff1 = end1 - start1;
    std::cout << "Time to multiply matrices using operator* : " << diff1.count() << " s\n";

    // Time second function
    auto start2 = std::chrono::high_resolution_clock::now();
    Matrix4x4int res2 = a;
    for (int i = 0; i < amount; i++) {
        res2 = res2 * b;
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to multiply matrices using Multiplication : " << diff2.count() << " s\n";
}