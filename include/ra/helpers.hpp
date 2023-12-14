#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <random>
#include <eigen3/Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Object.h>
#include <CGAL/convexity_check_3.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Simple_cartesian.h>


// The type used to represent the Kernel with exact predicates
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>::HalfedgeDS HalfedgeDS;

// The type used to represent points in two dimensions.
typedef Kernel::Point_2 Point_2;

typedef CGAL::Random_points_in_square_2<Point_2, CGAL::Creator_uniform_2< double, Point_2>>Point_generator;
typedef Kernel::Point_3 Point_3;


typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Kernel::Segment_3 Segment_3;

// The type used to represent vectors in two dimensions.
typedef CGAL::Polygon_2<Kernel> Polygon_2;

typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;

namespace ra {
    namespace helpers {

        bool is_2d(std::string& filename) {
            // read in shape from off file and set vertices and polygon
            std::ifstream input(filename);
            assert(input);

            std::string signature;
            assert(input >> signature && signature == "OFF");

            int num_vertices, num_faces, num_edges;
            assert(input >> num_vertices >> num_faces >> num_edges);

            double x, y, z;
            for (int i = 0; i < num_vertices; ++i) {
                input >> x >> y >> z;
                if (z != 0) {
                    return false;
                }
            }

            return true;
        }

        // Function to read the OFF file and store its contents in a string
        void create_random_shape_2d(std::string new_filename) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(3, 50);
            int random_n = dis(gen);

            Polygon_2 random_polygon;

            CGAL::random_convex_set_2(random_n, std::back_inserter(random_polygon), Point_generator(10.0));

            std::ofstream off_file(new_filename);
            if (!off_file) {
                std::cerr << "Error: Failed to open the output file.\n";
                return;
            }

            std::vector<Eigen::Vector3d> vertices;
            vertices.reserve(random_n);
            for (Polygon_2::Vertex_iterator it = random_polygon.vertices_begin(); it != random_polygon.vertices_end(); ++it) {
                Point_2 point = *it;
                vertices.push_back(Eigen::Vector3d(point.x(), point.y(), 0));
            }

            std::uniform_real_distribution<> translation_dist(-18.0, 18.0); // Adjust the range as needed
            // Generate random translation and rotation values
            double tx = translation_dist(gen);
            double ty = translation_dist(gen);
            double tz = translation_dist(gen);

             // Create the transformation matrix
            Eigen::Affine3d transform = Eigen::Affine3d::Identity();
            transform.translation() << tx, ty, tz;

            std::vector<Eigen::Vector3d> new_vertices;
            new_vertices.reserve(random_n);
            for (const Eigen::Vector3d& vertex : vertices) {
                new_vertices.push_back(transform * vertex);
            }

            off_file << "OFF\n";
            off_file << random_polygon.size() << " 1 " << random_polygon.size() << "\n";
            for (const Eigen::Vector3d& vertex : new_vertices) {
                off_file << vertex.x() << " " << vertex.y() << " 0.0\n";
            }
            off_file << random_polygon.size() << " "; // Number of vertices in the single face (the entire polygon)
            for (std::size_t i = 0; i < random_polygon.size(); ++i) {
                off_file << i << " "; // Write the indices of the vertices forming the face
            }
            off_file << "\n";
            off_file.close();
        }

        // Function to read the OFF file and store its contents in a string
        void create_random_shape_3d(std::string new_filename, std::string root_dir) {
            Polyhedron_3 polyhedron;
            
            std::random_device rd;
            std::mt19937 gen(rd());

            while(true) {
                std::uniform_int_distribution<> dis(5, 50);
                int random_n = dis(gen);

                CGAL::Random_points_in_sphere_3<Point_3, PointCreator> shape_gen(10.0);

                std::vector<Point_3> points;
                CGAL::copy_n( shape_gen, random_n, std::back_inserter(points) );
                
                // define object to hold convex hull 
                CGAL::Object ch_object;

                // compute convex hull 
                CGAL::convex_hull_3(points.begin(), points.end(), ch_object);

                if ( CGAL::assign (polyhedron, ch_object) ) {
                    break;
                }
            }

            std::ofstream untransformed_shape_file(root_dir + "input_shapes_3d/untransformed_shape.off");
            untransformed_shape_file << polyhedron;

            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::vector<size_t>> faces;
            std::ifstream input(root_dir + "input_shapes_3d/untransformed_shape.off");
            if (!input) {
                std::cerr << "Error: Unable to open file\n";
                return;
            }

            std::string signature;
            if (!(input >> signature) || signature != "OFF") {
                std::cerr << "not OFF format\n";
                return;
            }

            int num_vertices, num_faces, num_edges;
            if (!(input >> num_vertices >> num_faces >> num_edges)) {
                std::cerr << "cannot get number of vertices/faces/edges\n";
                return;
            }

            vertices.reserve(num_vertices);
            double x, y, z;
            for (int i = 0; i < num_vertices; ++i) {
                input >> x >> y >> z;
                vertices.push_back(Eigen::Vector3d(x, y, z));
            }
            
            int num_vertices_in_face;
            for (int i = 0; i < num_faces; ++i) {
                input >> num_vertices_in_face;
                std::vector<size_t> face_vertices;
                face_vertices.reserve(num_vertices_in_face);
                int vertex_index;
                for (int i = 0; i < num_vertices_in_face; ++i) {
                    input >> vertex_index;
                    face_vertices.push_back(vertex_index);
                }
                faces.push_back(face_vertices);
            }

            std::uniform_real_distribution<> translation_dist(-18.0, 18.0);
            std::uniform_real_distribution<> rotation_dist(-M_PI, M_PI);

            // Generate random translation and rotation values
            double tx = translation_dist(gen);
            double ty = translation_dist(gen);
            double tz = translation_dist(gen);
            double rx = rotation_dist(gen);
            double ry = rotation_dist(gen);
            double rz = rotation_dist(gen);

            // Create the transformation matrix
            Eigen::Affine3d transform = Eigen::Affine3d::Identity();
            transform.translation() << tx, ty, tz;
            transform.rotate(Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()));
            transform.rotate(Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()));
            transform.rotate(Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ()));

            std::vector<Eigen::Vector3d> new_vertices;
            new_vertices.reserve(num_vertices);
            for (const Eigen::Vector3d& vertex : vertices) {
                new_vertices.push_back(transform * vertex);
            }

            // Write random transformed shape to a new file
            std::string output_path = root_dir + "input_shapes_3d/" + new_filename + ".off";
            std::ofstream file_output(output_path);
            if (!file_output) {
                std::cerr << "Error opening file\n";
                return;
            }

            file_output << signature << "\n";
            file_output << num_vertices << " " << num_faces << " " << num_edges << "\n";
            for (const Eigen::Vector3d& vertex : new_vertices) {
                file_output << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
            }

            for (const std::vector<size_t>& face : faces) {
                file_output << face.size() << " ";
                for (const size_t& vertex_index : face) {
                    file_output << vertex_index << " ";
                }
                file_output << "\n";
            }

            file_output.close();
        }
    }
}