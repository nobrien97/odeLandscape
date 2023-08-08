# Build instructions

You will need a C++ compiler (I used GCC 9.4.0), the [Ascent](https://github.com/AnyarInc/Ascent/releases/tag/v0.7.1) library, and Cmake (I used 3.16.3)

1. Copy the Ascent library into include/ascent
2. Run `mkdir build && cd build`
3. Run `cmake ../`
4. Run `sudo make install ODELandscaper`
5. Run `ODELandscaper -h` for usage instructions