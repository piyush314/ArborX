#pragma once
#include <string>
#include <chrono>
#include <cstdio>
// #include "cuda_runtime.h"

// #define TIMING_BARRIER 1


class mytimer_t
{
protected:
    double wtime;

    std::chrono::high_resolution_clock::time_point tStart;
    int isStarted;

public:
    std::string name;
    mytimer_t()
    {
        wtime = 0;
        isStarted = 0;
    };
    mytimer_t(std::string str)
    {
        wtime = 0;
        isStarted = 0;
        name = str;
    };
    void start()
    {
        isStarted = 1;
        tStart = std::chrono::high_resolution_clock::now();
    }; // start the timer
    void end()
    {
        auto elapsed = std::chrono::high_resolution_clock::now() - tStart;
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        wtime += 1e-6 * (double)microseconds;
        isStarted = 0;

    }; // end the timer
    double getWtime() { return wtime; };
    void reset()
    {
        wtime = 0;
        isStarted = 0;
    }; // resets the timer

    void print()
    {
        std::cout << name;
        printf("\t\t:  \t%.3g \n", wtime);
    }
};