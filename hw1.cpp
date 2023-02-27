#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
typedef CGAL::Simple_cartesian<double>  K;
typedef K::Plane_3                  Plane;
struct FaceInfo {
  bool interior;
  FaceInfo() {
    interior = false;
  }
};
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;

const std::string input_file = "../station.hw1";
const std::string output_file = "../station.obj";

struct Vertex {
  int id;
  double x, y, z;
};

struct Face {
  int id;
  std::list<int> outer_ring;
  std::list<std::list<int>> inner_rings;
  Kernel::Plane_3 best_plane;
  Triangulation triangulation;
};

int main(int argc, const char * argv[]) {
  
  std::map<int, Vertex> vertices;
  std::map<int, Face> faces;
  
  // Read file
  std::ifstream input_stream;
  input_stream.open(input_file);
  if (input_stream.is_open()) {
    std::string line;
    
    // Read vertex header
    getline(input_stream, line);
    std::cout << "Vertex header: " << line << std::endl;
    std::istringstream vertex_header_stream(line);
    int number_of_vertices;
    vertex_header_stream >> number_of_vertices;
    std::cout << "Parsing " << number_of_vertices << " vertices..." << std::endl;
    
    // Read vertices
    for (int i = 0; i < number_of_vertices; ++i) {
      getline(input_stream, line);
      std::istringstream line_stream(line);
      int id;
      double x, y, z;
      line_stream >> id >> x >> y >> z;
      std::cout << "Vertex " << id << ": (" << x << ", " << y << ", " << z << ")" << std::endl;
      vertices[id].x = x;
      vertices[id].y = y;
      vertices[id].z = z;
    }

    // Read faces
    getline(input_stream, line);
    std::cout << "Face header: " << line << std::endl;
    std::istringstream face_header_stream(line);
    int number_of_faces;
    face_header_stream >> number_of_faces;
    std::cout << "Parsing " << number_of_faces << " faces..." << std::endl;

    for (int i = 0; i < number_of_faces; ++i) {
        faces[i].id = i;
        getline(input_stream, line);
        std::istringstream face_stream(line);
        int outer_num, inner_num;
        face_stream >> outer_num >> inner_num;
        std::cout << "Face " << i << " has inner ring number:" << inner_num << std::endl;

        // read outer ring
        getline(input_stream, line);
        std::istringstream outer_stream(line);
        int o_vertices;
        outer_stream >> o_vertices;
        std::cout << "Outer ring "  << " has " << o_vertices << " vertices:" << std::endl;
        while (outer_stream >> o_vertices) {
            std::cout << o_vertices << " ";
            faces[i].outer_ring.push_back(o_vertices);
        }
        std::cout << std::endl;

        // read inner ring
        for (int k = outer_num; k < inner_num + 1; ++k) {
            getline(input_stream, line);
            std::istringstream inner_stream(line);
            int i_vertices;
            inner_stream >> i_vertices;
            int innering_id = k - outer_num;
            std::cout << "Inner_ring " << innering_id << " has " << i_vertices << " vertices" << std::endl;

            std::list<int> inner_ring;

            while (inner_stream >> i_vertices) {
                std::cout << i_vertices << " ";
                inner_ring.push_back(i_vertices);
            }
            faces[i].inner_rings.push_back(inner_ring);
            std::cout << std::endl;
        }

        std::cout << std::endl;


//        std::list<int> outer_ring;
//        std::list<std::list<int>> inner_rings;
    }

    //



    // TO DO
  }

  // fitting plane
    for (auto const &face: faces) {
        std::cout << "Face " << face.first << ": " << std::endl;
        std::cout << "\touter:";
        for (auto const &vertex: face.second.outer_ring) std::cout << " " << vertex;
        std::cout << std::endl;
        for (auto const &ring: face.second.inner_rings) {
            std::cout << "\tinner:";
            for (auto const &vertex: ring) std::cout << " " << vertex;
            std::cout << std::endl;
        }
    }

//  for (auto const &face: faces)
//  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());

  std::cout << "debug" << std::endl;
  return 0;
}
