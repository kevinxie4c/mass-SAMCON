#include <iostream>
#include "timer.h"

Timer::Timer(): lastTime(std::chrono::steady_clock::now()) {}

std::chrono::steady_clock::duration Timer::duration()
{
    auto result = std::chrono::steady_clock::now() - lastTime;
    lastTime = std::chrono::steady_clock::now();
    return result;
}

std::string Timer::durationToString()
{
    std::chrono::steady_clock::duration d = duration();
    return std::to_string(std::chrono::duration_cast<std::chrono::seconds>(d).count()) + "s = " +
	std::to_string(std::chrono::duration_cast<std::chrono::minutes>(d).count()) + "m = " +
	std::to_string(std::chrono::duration_cast<std::chrono::hours>(d).count()) + "h";
}
