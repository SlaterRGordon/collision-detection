#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <thread>
#include "ra/shape.hpp"
#include "ra/helpers.hpp"

int main(int argc, char* argv[])
{
     if (argc != 5) {
        std::cout << "Usage: test_shapes_3d --root_dir $path_to_root_directory --out $path_to_output\n";
        return 1;
    }
    std::string path_to_root_dir;
    std::string path_to_output;

    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (arg == "--root_dir") {
            path_to_root_dir = argv[i + 1];
        } else if (arg == "--out") {
            path_to_output = argv[i + 1];
        } else {
            std::cout << "Invalid argument: " << arg << "\n";
            return 1;
        }
    }

    // Write the result to the output file
    std::ofstream output_file(path_to_output);
    if (!output_file) {
        std::cerr << "Error: Failed to open the output file.\n";
        return 1;
    }

    // generate two random 3d shapes
    for (int i = 0; i < 100; i++) {
        std::string off_file_name_1 = "shape1_" + std::to_string(i);
        ra::helpers::create_random_shape_3d(off_file_name_1, path_to_root_dir);
        std::string off_file_name_2 = "shape2_" + std::to_string(i);
        ra::helpers::create_random_shape_3d(off_file_name_2, path_to_root_dir);

        std::string off_file_3d_1 = path_to_root_dir + "input_shapes_3d/" + off_file_name_1 + ".off";
        std::string off_file_3d_2 = path_to_root_dir + "input_shapes_3d/" + off_file_name_2 + ".off";

        ra::shapes::ShapeBase* shape1 = new ra::shapes::Shape3D(off_file_3d_1);
        ra::shapes::ShapeBase* shape2 = new ra::shapes::Shape3D(off_file_3d_2);

        // Check if the shapes intersect
        bool intersects = shape1->minimum_box_intersects(*shape2);
        if (intersects) {
            intersects = shape1->gjk_intersects(*shape2);
        }

        // Deallocate the memory for the shapes
        delete shape1;
        delete shape2;

        // Write the result (true or false) to the output file
        output_file << off_file_3d_1 << " " << off_file_3d_2 << ": " << std::boolalpha << intersects << "\n";
    }

    output_file.close();

    return 0;
}