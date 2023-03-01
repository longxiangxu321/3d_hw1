#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Plane_3 Plane;

struct FaceInfo {
  bool interior;
  bool processed;
  FaceInfo() {
    interior = false;
  }
};

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;
typedef Triangulation::Face_handle Face_handle;

const std::string input_file = "../station.hw1";
const std::string output_file = "../station.obj";

//struct Vertex {
//  int id;
//  double x, y, z;
//};

struct Face {
  int id;
  std::list<int> outer_ring;
  std::list<std::list<int>> inner_rings;
  Plane best_plane;
  Triangulation triangulation;
};




void label_triangles(Triangulation &triangulation) {
    std::list<Triangulation::Face_handle> to_check;
    triangulation.infinite_face()->info().processed = true;
    to_check.push_back(triangulation.infinite_face());
    while (!to_check.empty()) {
//        CGAL_assertion(to_check.front()->info().processed == true);
        for (int neighbour = 0; neighbour < 3; ++neighbour) {
            if (to_check.front()->neighbor(neighbour)->info().processed) {
            } else {
                to_check.front()->neighbor(neighbour)->info().processed = true;
//                CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
                    to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                    to_check.push_back(to_check.front()->neighbor(neighbour));
                } else {
                    to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                    to_check.push_back(to_check.front()->neighbor(neighbour));
                }
            }
        } to_check.pop_front();
    }
}



int main(int argc, const char * argv[]) {
  
  std::map<int, Point> vertices;
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
      vertices[id] = Point(x, y, z);
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


    }

    //


    // TO DO
  }

  // fitting plane
  for (auto &face: faces) {
    std::vector<Point> pts;
    for (auto const &vertex: face.second.outer_ring) pts.push_back(vertices[vertex]);
    Plane plane;
    linear_least_squares_fitting_3(pts.begin(),pts.end(), face.second.best_plane,CGAL::Dimension_tag<0>());
    std::cout << face.second.best_plane << std::endl;
  }

  //constrained delaunay triangulation

  // creating constraints
  for (auto &face: faces) {
    std::vector<std::pair<int, int>> segments;
    // obtaining vertices
    std::list<int> vertice_ids;

    auto it = face.second.outer_ring.begin();
    for (int i=0; i<face.second.outer_ring.size(); i++) {
        vertice_ids.push_back(*it);
        int first = *it;
        it ++;
        if (it == face.second.outer_ring.end()) it = face.second.outer_ring.begin();
        int second = *it;
        segments.push_back(std::make_pair(first, second));
    }

    for (auto const &ring: face.second.inner_rings) {
        auto it = ring.begin();
        for (int i=0; i<ring.size(); i++) {
            vertice_ids.push_back(*it);
            int first = *it;
            it ++;
            if (it == ring.end()) it = ring.begin();
            int second = *it;
            segments.push_back(std::make_pair(first, second));

        }
    }

    Triangulation cdt;
    std::vector< std::pair<Point_2, int> > subset;
    for (auto id : vertice_ids) {
        Point_2 vertex_to_insert = face.second.best_plane.to_2d(vertices[id]);
        auto to_insert = std::make_pair(vertex_to_insert, id);
        subset.push_back(to_insert);
    }
    cdt.insert(subset.begin(), subset.end());


    for (const auto& s : segments) {
      Point_2 s_1 = face.second.best_plane.to_2d(vertices[s.first]);
      Point_2 s_2 = face.second.best_plane.to_2d(vertices[s.second]);
      cdt.insert_constraint(s_1, s_2);
    }

    face.second.triangulation = cdt;

    label_triangles(face.second.triangulation);

    std::ifstream output_stream;
    output_stream.open(output_file);
//    for (auto vit = face.second.triangulation.finite_vertices_begin();
//    vit != face.second.triangulation.finite_vertices_end(); ++vit) {
//        auto temp_pt = face.second.best_plane.to_3d(vit->point());
//        auto temp_pt_pair = std::make_pair(temp_pt, vit->info());
//    }
  }


  return 0;
}
