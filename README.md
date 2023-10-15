## Problem
I found that calculating the concave hull of a set of data-points was a primary
source of latency in a system I was working on. Given that this calculation was in the 
hot-path, I set about looking for ways to optimise the calculation. 

In previous experiences, I had found that in some cases, simply running the algorithm
in C++ and using Pybind to call it from Python resulted in a significant reduction
in calculation time for a given function call. As such, I was curious if porting the 
python based concave hull calculation over to C++ would result in any improvements.
At the same time, I was planning to use this project as a way to learn some low-latency
techniques by simply seeing how fast I could get it to run. 

## Objectives
- Implement the C++ version of the python concave hull calculation
- Provide an easy to use C++ based concave hull calculation library
- Use Modern C++ (C++20+)
- Attempt to minimise latency using low-latency techniques and tools

## Progress
The first goal was to port as much as of the Python implementation over to C++. 
This proved to be a somewhat challenging task, simply due to the lack of certain
conveniences that Python and certain libraries such as Numpy and Shapely provides. 
Using a test set of points, I was able to successfully implement code that replicated 
the output of the Python based version. A unit test has been written that 
tests the output of the hull calculation against the known result. The next step
will be to cover the code with more unit tests before moving to performane measuring.
