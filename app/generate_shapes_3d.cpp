#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include "ra/helpers.hpp"

int main(int argc, char* argv[]) {
    for (int i = 0; i < 1000; i++) {
        std::string filename = "shape" + std::to_string(i);
        ra::helpers::create_random_shape_3d(filename);
    }
}