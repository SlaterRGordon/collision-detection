#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/OFF_to_nef_3.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/convexity_check_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <limits>
#include <cassert>

// The type used to represent the Kernel with exact predicates
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>::HalfedgeDS HalfedgeDS;

// The type used to represent points in two dimensions.
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_2 Vector_2;

// The type used to represent points in three dimensions.
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_3 Vector_3;

// The type used to represent vectors in two dimensions.
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

// The type used to represent vectors in three dimensions.
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron_3;

struct BoundingBox2D {
    Point_2 lower_left;  // Bottom-left point of the bounding box
    Point_2 upper_right; // Top-right point of the bounding box
};

struct BoundingBox3D {
    Point_3 lower_left;  // Bottom-left point of the bounding box
    Point_3 upper_right; // Top-right point of the bounding box
};

namespace ra
{
    namespace shapes
    {   
        struct Simplex {
        private:
            std::array<Point_3, 4> m_points;
            unsigned m_size;

        public:
            Simplex()
                : m_points({ Point_3(), Point_3(), Point_3(), Point_3() })
                , m_size(0)
            {}

            Simplex& operator=(std::vector<Point_3> list) {
                for (auto v = list.begin(); v != list.end(); v++) {
                    m_points[std::distance(list.begin(), v)] = *v;
                }
                m_size = list.size();

                return *this;
            }

            void push_front(Point_3 point) {
                m_points = { point, m_points[0], m_points[1], m_points[2] };
                m_size = std::min(m_size + 1, 4u);
            }

            Point_3& operator[](unsigned i) { return m_points[i]; }
            unsigned size() const { return m_size; }

            auto begin() const { return m_points.begin(); }
            auto end()   const { return m_points.end() - (4 - m_size); }
        };

        // Base class for Shape
        class ShapeBase {
        public:
            ShapeBase() {}
            virtual ~ShapeBase() {}

            // Common methods for both 2D and 3D shapes go here
            virtual bool minimum_box_intersects(const ShapeBase& other) = 0;
            virtual bool gjk_intersects(const ShapeBase& other) = 0;
        };

        class Shape2D : public ShapeBase {
        public:
            Shape2D(const std::string& off_file) {
                read_from_off(off_file);
                compute_bounding_box();
                set_centroid();
            }

            bool minimum_box_intersects(const ShapeBase& other) override {
                const Shape2D* other_shape = dynamic_cast<const Shape2D*>(&other);
                
                // Check if the bounding boxes intersect first
                if (bounding_box.lower_left.x() > other_shape->bounding_box.upper_right.x() || bounding_box.upper_right.x() < other_shape->bounding_box.lower_left.x() ||
                    bounding_box.lower_left.y() > other_shape->bounding_box.upper_right.y() || bounding_box.upper_right.y() < other_shape->bounding_box.lower_left.y()) {
                    return false;
                }

                return true;
            }

            bool gjk_intersects(const ShapeBase& other) override {
                const Shape2D* other_shape = dynamic_cast<const Shape2D*>(&other);

                // Pick starting direction
                Point_2 center_point(other_shape->centroid.x() - centroid.x(), other_shape->centroid.y() - centroid.y());
                direction = normalize(center_point);

                simplex.push_back(support(other_shape));

                // change direction towards origin
                direction = Vector_2(-direction.x(), -direction.y());

                while(true) {
                    Point_2 support_point = support(other_shape);

                    // if support point is not past origin in direction d
                    if (dot_product(support_point, direction) <= 0) {
                        return false;
                    }

                    simplex.push_back(support_point);

                    // check if simplex contains origin
                    if (handle_simplex()) {
                        return true;
                    }
                }
            }
        private:
            std::vector<Point_2> vertices;
            Polygon_2 polygon;
            BoundingBox2D bounding_box;
            Point_2 centroid;

            Polygon_2 simplex;
            Vector_2 direction;

            // read in shape from off file and set vertices and polygon
            void read_from_off(const std::string& filename) {
                std::ifstream input(filename);
                if (!input) {
                    std::cerr << "Error: Unable to open file " << filename << std::endl;
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
                    vertices.emplace_back(x, y);
                }

                int num_vertices_in_face;
                input >> num_vertices_in_face;
                std::vector<Point_2> face_vertices;
                face_vertices.reserve(num_vertices_in_face);
                int vertex_index;
                for (int i = 0; i < num_vertices_in_face; ++i) {
                    input >> vertex_index;
                    face_vertices.push_back(vertices[vertex_index]);
                }

                polygon = Polygon_2(face_vertices.begin(), face_vertices.end());

                if (!polygon.is_convex()) {
                    std::cerr << "Shapes must be convex to detect collision\n";
                    return;
                }
            }

            // compute the minimum bounding box around the shape
            void compute_bounding_box() {
                double min_lower_left_x = vertices[0].x();
                double min_lower_left_y = vertices[0].y();
                double max_upper_right_x = vertices[0].x();
                double max_upper_right_y = vertices[0].y();

                // find the minimum and maximum extents along the x and y axes
                for (const Point_2& p : vertices) {
                    min_lower_left_x = std::min(min_lower_left_x, p.x());
                    min_lower_left_y = std::min(min_lower_left_y, p.y());
                    max_upper_right_x = std::max(max_upper_right_x, p.x());
                    max_upper_right_y = std::max(max_upper_right_y, p.y());
                }

                // set the bounding box
                bounding_box.lower_left = Point_2(min_lower_left_x, min_lower_left_y);
                bounding_box.upper_right = Point_2(max_upper_right_x, max_upper_right_y);
            }

            // compute the centroid of the shape
            void set_centroid() {
                int num_vertices = static_cast<int>(vertices.size());
                assert(num_vertices > 0);

                double sumX = 0.0;
                double sumY = 0.0;

                for (const Point_2& vertex : vertices) {
                    sumX += vertex.x();
                    sumY += vertex.y();
                }

                centroid = Point_2(sumX / num_vertices, sumY / num_vertices);
            }

            bool handle_simplex() {
                assert(simplex.size() > 0);

                if (simplex.size() == 2) {
                    return handle_line();
                }
                else if (simplex.size() == 3) {
                    return handle_triangle();
                }

                return false;
            }

            bool handle_line() {
                Point_2 a = simplex[1];
                Point_2 b = simplex[0];

                Vector_2 ab = Vector_2(b.x() - a.x(), b.y() - a.y());
                Vector_2 ao = Vector_2(-a.x(), -a.y()); // origin - a

                // get direction perpendicular to ab in direction of origin
                direction = triple_cross_product(ab, ao, ab);

                return false;
            }

            bool handle_triangle() {
                Point_2 c = simplex[2];
                Point_2 b = simplex[1];
                Point_2 a = simplex[0];

                Vector_2 ab = Vector_2(b.x() - a.x(), b.y() - a.y());
                Vector_2 ac = Vector_2(c.x() - a.x(), c.y() - a.y());
                Vector_2 ao = Vector_2(-a.x(), -a.y()); // origin - a
                
                Vector_2 ab_perp = triple_cross_product(ac, ab, ab);
                Vector_2 ac_perp = triple_cross_product(ab, ac, ac);

                // check if origin is in region ab
                if (dot_product(ab_perp, ao) > 0) {
                    simplex = Polygon_2(); // remove point c
                    simplex.push_back(a);
                    simplex.push_back(b);
                    direction = ab_perp;

                    return false;
                } else if (dot_product(ac_perp, ao) > 0) {
                    simplex = Polygon_2(); // remove point b
                    simplex.push_back(a);
                    simplex.push_back(c);
                    direction = ac_perp;
                    
                    return false;
                } else {
                    // contains origin
                    return true;
                }
            }
            
            // compute cross product of v1 and v2 and then cross product of result and v1
            Vector_2 triple_cross_product(const Vector_2& v1, const Vector_2& v2, const Vector_2& v3) const {
                double new_new_x = -(v1.x() * v2.y() - v1.y() * v2.x()) * v3.y();
                double new_new_y = (v1.x() * v2.y() - v1.y() * v2.x()) * v3.x();

                return Vector_2(new_new_x, new_new_y);
            }

            // calculate the furthest points in direction d and -d and return the difference
            Point_2 support(const Shape2D* other) const {
                Point_2 furthest_point_v1 = furthest_point(direction);
                Point_2 furthest_point_v2 = other->furthest_point(Vector_2(-direction.x(), -direction.y()));

                return Point_2(furthest_point_v1.x() - furthest_point_v2.x(), furthest_point_v1.y() - furthest_point_v2.y());
            }

            // calculate furthest point in direction d
            Point_2 furthest_point(const Vector_2& d) const {
                Point_2 furthest_point = vertices[0];
                double max_dot_product = dot_product(vertices[0], d);
                for (const Point_2& vertex : vertices) {
                    double new_dot_product = dot_product(vertex, d);
                    if (new_dot_product > max_dot_product) {
                        max_dot_product = new_dot_product;
                        furthest_point = vertex;
                    }
                }

                return furthest_point;
            }
            
            // normalize vector v
            Point_2 normalize(Point_2& v) const {
                double length = std::sqrt(v.x() * v.x() + v.y() * v.y());
                assert(length > 0.0);

                return Point_2(v.x() / length, v.y() / length);
            }

            // Mag A * Mag B * cos(theta)
            double dot_product(const Vector_2& v1, const Point_2& v2) const {
                double maginutude_a = std::sqrt(v1.x() * v1.x() + v1.y() * v1.y());
                double maginutude_b = std::sqrt(v2.x() * v2.x() + v2.y() * v2.y());

                if(maginutude_a == 0 || maginutude_b == 0) {
                    return 0;
                }

                double theta = std::acos((v1.x() * v2.x() + v1.y() * v2.y()) / (maginutude_a * maginutude_b));
                double dot_product = maginutude_a * maginutude_b * std::cos(theta);

                return dot_product;
            }

            bool is_convex() const {
                return polygon.is_convex();
            }
        };

        class Shape3D : public ShapeBase {
        public:
            Shape3D(const std::string& off_file) {
                read_from_off(off_file);
                compute_bounding_box();
                set_centroid();
            }

            bool minimum_box_intersects(const ShapeBase& other) override {
                const Shape3D* other_shape = dynamic_cast<const Shape3D*>(&other);

                if (is_same(*other_shape)) {
                    return true;
                }

                // Check if the bounding boxes intersect first
                if (bounding_box.lower_left.x() > other_shape->bounding_box.upper_right.x() ||
                    bounding_box.upper_right.x() < other_shape->bounding_box.lower_left.x() ||
                    bounding_box.lower_left.y() > other_shape->bounding_box.upper_right.y() ||
                    bounding_box.upper_right.y() < other_shape->bounding_box.lower_left.y() ||
                    bounding_box.lower_left.z() > other_shape->bounding_box.upper_right.z() ||
                    bounding_box.upper_right.z() < other_shape->bounding_box.lower_left.z()) {
                    return false;
                }

                return true;
            }

            bool gjk_intersects(const ShapeBase& other) override {
                const Shape3D* other_shape = dynamic_cast<const Shape3D*>(&other);

                if (is_same(*other_shape)) {
                    return true;
                }

                // Pick starting direction
                Point_3 center_point(other_shape->centroid.x() - centroid.x(), other_shape->centroid.y() - centroid.y(), other_shape->centroid.z() - centroid.z());
                direction = normalize(center_point);
                
                // Add initial support point to simplex
                simplex = Simplex();
                simplex.push_front(support(other_shape));

                // change direction towards origin
                direction = Vector_3(-direction.x(), -direction.y(), -direction.z());

                while(true) {
                    Point_3 support_point = support(other_shape);

                    // if support point is not past origin in direction d
                    if (dot_product(support_point, direction) <= 0) {
                        return false;
                    }
                    
                    simplex.push_front(support_point);

                    // check if simplex contains origin
                    if (handle_simplex()) {
                        return true;
                    }
                }
            }

            std::vector<Point_3> get_vertices() const {
                return vertices;
            }

        private:
            std::vector<Point_3> vertices;
            std::vector<std::vector<size_t>> faces;
            BoundingBox3D bounding_box;
            Point_3 centroid;
            Simplex simplex;
            Vector_3 direction;
            Polyhedron_3 polyhedron;
            Nef_polyhedron_3 nef_polyhedron;

            bool is_same(const Shape3D& other) const {
                if (vertices == other.get_vertices()) {
                    return true;
                } else {
                    return false;
                }
            }

            // read in shape from off file and set vertices and polygon
            void read_from_off(const std::string& filename) {
                std::ifstream input(filename);
                if (!input) {
                    std::cerr << "Error: Unable to open file " << filename << std::endl;
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
                    vertices.push_back(Point_3(x, y, z));
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
                input.close();

                std::ifstream polyhedronInput(filename);
                if (!polyhedronInput || !(polyhedronInput >> polyhedron) || polyhedron.empty()) {
                    std::cerr << "Error: Unable to read the OFF file or the polyhedron is empty.\n";
                    return;
                }
            }

            // compute the minimum bounding box around the shape
            void compute_bounding_box() {
                double min_lower_left_x = vertices[0].x();
                double min_lower_left_y = vertices[0].y();
                double min_lower_left_z = vertices[0].z();
                double max_upper_right_x = vertices[0].x();
                double max_upper_right_y = vertices[0].y();
                double max_upper_right_z = vertices[0].z();

                // find the minimum and maximum extents along the x and y axes
                for (const Point_3& p : vertices) {
                    min_lower_left_x = std::min(min_lower_left_x, p.x());
                    min_lower_left_y = std::min(min_lower_left_y, p.y());
                    min_lower_left_z = std::min(min_lower_left_z, p.z());
                    max_upper_right_x = std::max(max_upper_right_x, p.x());
                    max_upper_right_y = std::max(max_upper_right_y, p.y());
                    max_upper_right_z = std::max(max_upper_right_z, p.z());
                }

                // set the bounding box
                bounding_box.lower_left = Point_3(min_lower_left_x, min_lower_left_y, min_lower_left_z);
                bounding_box.upper_right = Point_3(max_upper_right_x, max_upper_right_y, max_upper_right_z);
            }

            // compute the centroid of the shape
            void set_centroid() {
                int num_vertices = static_cast<int>(vertices.size());
                assert(num_vertices > 0);

                double sumX = 0.0;
                double sumY = 0.0;
                double sumZ = 0.0;

                for (const Point_3& vertex : vertices) {
                    sumX += vertex.x();
                    sumY += vertex.y();
                    sumZ += vertex.z();
                }

                centroid = Point_3(sumX / num_vertices, sumY / num_vertices, sumZ / num_vertices);
            }

            bool handle_simplex() {
                assert(simplex.size() > 0);
                if (simplex.size() == 2) {
                    return handle_line();
                }
                else if (simplex.size() == 3) {
                    return handle_triangle();
                } else if (simplex.size() == 4) {
                    return handle_tetrahedron();
                }

                return false;
            }

            bool handle_line() {
                Point_3 a = simplex[0];
                Point_3 b = simplex[1];

                Vector_3 ab = Vector_3(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
                Vector_3 ao = Vector_3(-a.x(), -a.y(), -a.z()); // origin - a

                // get direction perpendicular to ab in direction of origin
                Vector_3 ab_perp = triple_cross_product(ab, ao, ab);

                if(dot_product(ab_perp, ao) > 0.0) {
                    direction = ab_perp;
                } else { // Not in ab region pick new point
                    simplex = { a };
                    direction = ao;
                }

                return false;
            }

            bool handle_triangle() {
                Point_3 a = simplex[0];
                Point_3 b = simplex[1];
                Point_3 c = simplex[2];

                Vector_3 ab = Vector_3(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
                Vector_3 ac = Vector_3(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
                Vector_3 ao = Vector_3(-a.x(), -a.y(), -a.z()); // origin - a

                Vector_3 abc_perp = cross_product(ab, ac);
                
                // get direction perpendicular to ab in direction of origin
                Vector_3 ab_perp = triple_cross_product(ac, ab, ab);
                Vector_3 ac_perp = triple_cross_product(ab, ac, ac);
            
                if (dot_product(ac_perp, ao) > 0.0) {
                    if (dot_product(ac, ao) > 0.0) {
                        simplex = { a, c };
                        direction = ac_perp;
                    } else {
                        simplex = { a, b };
                        return handle_line();
                    }
                } else if (dot_product(ab_perp, ao) > 0.0) { // check if origin is in region ab
                    simplex = { a, b };
                    return handle_line();
                } else if (dot_product(abc_perp, ao) > 0.0) { // is origin above triangle
                    direction = abc_perp;
                } else { // origin is below triangle
                    simplex = { a, c, b };
                    direction = Vector_3(-abc_perp.x(), -abc_perp.y(), -abc_perp.z());
                }

                return false;
            }

            bool handle_tetrahedron() {
                Point_3 a = simplex[0];
                Point_3 b = simplex[1];
                Point_3 c = simplex[2];
                Point_3 d = simplex[3];

                Vector_3 ab = Vector_3(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
                Vector_3 ac = Vector_3(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
                Vector_3 ad = Vector_3(d.x() - a.x(), d.y() - a.y(), d.z() - a.z());
                Vector_3 ao = Vector_3(-a.x(), -a.y(), -a.z()); // origin - a

                // vectors perpendicular to faces of tetrahedron
                Vector_3 abc = cross_product(ab, ac);
                Vector_3 acd = cross_product(ac, ad);
                Vector_3 adb = cross_product(ad, ab);

                if (dot_product(abc, ao) > 0.0) {
                    simplex = { a, b, c };
                    return handle_triangle();
                } else if (dot_product(acd, ao) > 0.0) {
                    simplex = { a, c, d };
                    return handle_triangle();
                } else if (dot_product(adb, ao) > 0.0) {
                    simplex = { a, d, b };
                    return handle_triangle();
                } else {
                    return true;
                }
            }

            Vector_3 cross_product(const Vector_3& v1, const Vector_3& v2) {
                double new_x = (v1.y() * v2.z() - v1.z() * v2.y());
                double new_y = -(v1.x() * v2.z() - v1.z() * v2.x());
                double new_z = (v1.x() * v2.y() - v1.y() * v2.x());

                return Vector_3(new_x, new_y, new_z);
            }
            
            // compute cross product of v1 and v2 and then cross product of result and v1
            Vector_3 triple_cross_product(const Vector_3& v1, const Vector_3& v2, const Vector_3& v3) const {
                double new_x = (v1.y() * v2.z() - v1.z() * v2.y());
                double new_y = -(v1.x() * v2.z() - v1.z() * v2.x());
                double new_z = (v1.x() * v2.y() - v1.y() * v2.x());

                double new_new_x = (new_y * v3.z() - new_z * v3.y());
                double new_new_y = -(new_x * v3.z() - new_z * v3.x());
                double new_new_z = (new_x * v3.y() - new_y * v3.x());

                return Vector_3(new_new_x, new_new_y, new_new_z);
            }

            // calculate the furthest points in direction d and -d and return the difference
            Point_3 support(const Shape3D* other) const {
                Point_3 furthest_point_v1 = furthest_point(direction);
                Point_3 furthest_point_v2 = other->furthest_point(Vector_3(-direction.x(), -direction.y(), -direction.z()));

                return Point_3(furthest_point_v1.x() - furthest_point_v2.x(), furthest_point_v1.y() - furthest_point_v2.y(), furthest_point_v1.z() - furthest_point_v2.z());
            }

            // calculate furthest point in direction d
            Point_3 furthest_point(const Vector_3& d) const {
                Point_3 furthest_point = vertices[0];
                double max_dot_product = dot_product(vertices[0], d);
                for (const Point_3& vertex : vertices) {
                    double new_dot_product = dot_product(vertex, d);
                    if (new_dot_product > max_dot_product) {
                        max_dot_product = new_dot_product;
                        furthest_point = vertex;
                    }
                }

                return furthest_point;
            }
            
            // normalize vector v
            Point_3 normalize(Point_3& v) const {
                double length = std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
                assert(length > 0.0);

                return Point_3(v.x() / length, v.y() / length, v.z() / length);
            }

            // Mag A * Mag B * cos(theta)
            double dot_product(const Vector_3& v1, const Vector_3& v2) const {
                double maginutude_a = std::sqrt(v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z());
                double maginutude_b = std::sqrt(v2.x() * v2.x() + v2.y() * v2.y() + v2.z() * v2.z());

                if(maginutude_a == 0 || maginutude_b == 0) {
                    return 0;
                }

                double theta = std::acos((v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z()) / (maginutude_a * maginutude_b));
                double dot_product = maginutude_a * maginutude_b * std::cos(theta);

                return dot_product;
            }
        };
    }
}