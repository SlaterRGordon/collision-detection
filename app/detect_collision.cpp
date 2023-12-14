#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include "ra/shape.hpp"
#include "ra/helpers.hpp"

int main(int argc, char* argv[]) {
    // Check if the number of command-line arguments is correct
    if (argc != 7) {
        std::cout << "Usage: detect_collision --shape1 $path_to_shape1 --shape2 $path_to_shape2 --out $path_to_output\n";
        return 1;
    }

    // Read the command-line arguments
    std::string path_to_shape1;
    std::string path_to_shape2;
    std::string path_to_output;

    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (arg == "--shape1") {
            path_to_shape1 = argv[i + 1];
        } else if (arg == "--shape2") {
            path_to_shape2 = argv[i + 1];
        } else if (arg == "--out") {
            path_to_output = argv[i + 1];
        } else {
            std::cout << "Invalid argument: " << arg << "\n";
            return 1;
        }
    }

    ra::shapes::ShapeBase* shape1;
    ra::shapes::ShapeBase* shape2;

    if(ra::helpers::is_2d(path_to_shape1) && ra::helpers::is_2d(path_to_shape2)) {
        shape1 = new ra::shapes::Shape2D(path_to_shape1);
        shape2 = new ra::shapes::Shape2D(path_to_shape2);
    } else if (!ra::helpers::is_2d(path_to_shape1) && !ra::helpers::is_2d(path_to_shape2)) {
        shape1 = new ra::shapes::Shape3D(path_to_shape1);
        shape2 = new ra::shapes::Shape3D(path_to_shape2);
    } else {
        std::cout << "Error: Cannot compare 2D and 3D shapes.\n";
        return 1;
    }

    // Check if the shapes intersect
    bool intersects = shape1->minimum_box_intersects(*shape2);
    if (intersects) {
        intersects = shape1->gjk_intersects(*shape2);
    }

    // Deallocate the memory for the shapes
    delete shape1;
    delete shape2;

    // Write the result to the output file
    std::ofstream output_file(path_to_output);
    if (!output_file) {
        std::cout << "Error: Failed to open the output file.\n";
        return 1;
    }

    // Write the result (true or false) to the output file
    output_file << path_to_shape1 << " " << path_to_shape2 << " " << std::boolalpha << intersects;
    output_file.close();

    return 0;
}