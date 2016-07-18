#pragma once
#ifndef _TIMER
	#define _TIMER
#endif

#include <chrono>
#include <string>
#include <iostream>

namespace mytimer
{
class timer
{
private:
	std::chrono::high_resolution_clock::time_point starttime;
	std::string startstring;

public:
	// Constructors

	timer () {};

	timer( const std::string& s )
	{
		start(s);
	};

	// Destructor
	~timer() {}

	void start(const std::string& s)
	{
		starttime = std::chrono::high_resolution_clock::now();
		startstring = s;
		std::cout << startstring << std::flush;
	}

	void stop() const
	{
		std::chrono::high_resolution_clock::time_point endtime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_elapsed = endtime - starttime;
		std::cout << "DONE!" << " (time taken: " << time_elapsed.count() << " seconds)" << std::endl;
	}
};
} //end namespace timer
