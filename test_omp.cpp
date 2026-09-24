//
// Created by Peter Zimmerman on 30.11.25.
//
#include <omp.h>
#include <iostream>

int main() {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp critical
        std::cout << "Hello from thread " << tid << "\n";
    }
    return 0;
}