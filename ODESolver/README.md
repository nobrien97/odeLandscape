# Build instructions

You will need a C++ compiler (I used GCC 9.4.0), the [Ascent](https://github.com/AnyarInc/Ascent/releases/tag/v0.7.1) library, and Cmake (I used 3.16.3).
This software also requires header libraries [stats](https://github.com/kthohr/stats/tree/f8dcb15ae51cce7142b239805745a0de56aa509f) and [gcem](https://github.com/kthohr/gcem/tree/012ae73c6d0a2cb09ffe86475f5c6fba3926e200).

1. Copy the Ascent library into include/ascent and stats and gcem into include/
2. Run `mkdir build && cd build`
3. Run `cmake ../`
4. Run `sudo make install ODELandscaper`
5. Run `ODELandscaper -h` for usage instructions