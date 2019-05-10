#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>

class Timer
{
    public:
	Timer();

	std::chrono::steady_clock::duration duration();

	std::string durationToString();

    private:
	std::chrono::steady_clock::time_point lastTime;
};


#endif
