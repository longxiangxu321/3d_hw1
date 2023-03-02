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
    int nesting_level;
    bool in_domain(){
        return nesting_level%2 == 1;
    }
};

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;
typedef Triangulation::Face_handle Face_handle;

const std::string input_file = "../stilted_house.hw1";
const std::string output_file = "../stilted_house.obj";

//struct Vertex {
//  int id;
//  double x, y, z;
//};
void
mark_domains(Triangulation & ct,
             Face_handle start,
             int index,
             std::list<Triangulation::Edge>& border )
{
    if(start->info().nesting_level != -1){
        return;
    }
    std::list<Face_handle> queue;
    queue.push_back(start);
    while(! queue.empty()){
        Face_handle fh = queue.front();
        queue.pop_front();
        if(fh->info().nesting_level == -1){
            fh->info().nesting_level = index;
            for(int i = 0; i < 3; i++){
                Triangulation::Edge e(fh,i);
                Face_handle n = fh->neighbor(i);
                if(n->info().nesting_level == -1){
                    if(ct.is_constrained(e)) border.push_back(e);
                    else queue.push_back(n);
                }
            }
        }
    }
}


void
mark_domains(Triangulation& cdt)
{
    for(Triangulation::Face_handle f : cdt.all_face_handles()){
        f->info().nesting_level = -1;
    }
    std::list<Triangulation::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while(! border.empty()){
        Triangulation::Edge e = border.front();
        border.pop_front();
        Face_handle n = e.first->neighbor(e.second);
        if(n->info().nesting_level == -1){
            mark_domains(cdt, n, e.first->info().nesting_level+1, border);
        }
    }
}

struct Face {
    int id;
    std::list<int> outer_ring;
    std::list<std::list<int>> inner_rings;
    Plane best_plane;
    Triangulation triangulation;
    std::vector< std::pair<Point_2,int>> points_2d;
    std::vector<Point_2> outer_points;
    std::vector<std::vector<Point_2>> inner_points;
    int count;
    Face(int id, std::list<int> outer_ring, std::list<std::list<int>> inner_rings) :
    id(id), outer_ring(outer_ring), inner_rings(inner_rings), best_plane(), triangulation(), points_2d() {}


    void fitting_plane(std::map<int, Point> vertices) {
        std::vector<Point> pts;
        for (auto const &vertex: outer_ring) pts.push_back(vertices[vertex]);
        for (const auto& inner: inner_rings) {
            for (auto const &vi :inner) pts.push_back(vertices[vi]);
        }
        linear_least_squares_fitting_3(pts.begin(), pts.end(), best_plane, CGAL::Dimension_tag<0>());
//        std::cout << best_plane << std::endl;
    }

    void cdt(std::map<int, Point> vertices) {
        std::vector<std::pair<int, int>> segments;
        // obtaining vertices
        int counter = 0;

        for (int vertexIndex_o : outer_ring) {
            Point_2 vertex1 = best_plane.to_2d(vertices[vertexIndex_o]);
            outer_points.push_back(vertex1);
//            auto vt = std::make_pair(vertex1, counter);
            Triangulation::Vertex_handle v1 = triangulation.insert(vertex1);
            v1->info() = counter;
            counter++;
        }

        for (const auto& ring:inner_rings) {
            std::vector<Point_2> in_pt;
            for (int vertexIndex_i : ring) {
                Point_2 vert1 = best_plane.to_2d(vertices[vertexIndex_i]);
                in_pt.push_back(vert1);
//                points_2d.push_back(std::make_pair(vert1, counter));
                Triangulation::Vertex_handle v1 = triangulation.insert(vert1);
                v1->info() = counter;
                counter++;
            }
            inner_points.push_back(in_pt);
        }

        for (std::size_t i = 0; i < outer_points.size(); ++i) {
            triangulation.insert_constraint(outer_points[i], outer_points[(i + 1) % outer_points.size()]);
        }

        // Insert inner polygons as constraints
        for (const auto& inner_poly : inner_points) {
            for (std::size_t i = 0; i < inner_poly.size(); ++i) {
                triangulation.insert_constraint(inner_poly[i], inner_poly[(i + 1) % inner_poly.size()]);
            }
        }

//        for (auto it = outer_ring.begin(); it != outer_ring.end(); ++it) {
//            auto nextIt = std::next(it);
//            if (nextIt == outer_ring.end()) {
//                nextIt = outer_ring.begin();
//            }
//            auto vertexIndex1 = best_plane.to_2d(vertices[*it]);
//            auto vertexIndex2 = best_plane.to_2d(vertices[*nextIt]);
//            triangulation.insert_constraint(vertexIndex1, vertexIndex2);
//        }
//
//        for (auto inner: inner_rings) {
//            for (auto it = inner.begin(); it != inner.end(); ++it) {
//                auto nextIt = std::next(it);
//                if (nextIt == inner.end()) {
//                    nextIt = inner.begin();
//                }
//                auto vertexIndex1 = best_plane.to_2d(vertices[*it]);
//                auto vertexIndex2 = best_plane.to_2d(vertices[*nextIt]);
//                triangulation.insert_constraint(vertexIndex1, vertexIndex2);
//            }
//        }


        count = counter;

//        triangulation.insert(points_2d.begin(), points_2d.end());

    };
};


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
            getline(input_stream, line);
            std::istringstream face_stream(line);
            int outer_num, inner_num;
            face_stream >> outer_num >> inner_num;
            std::cout << "Face " << i << " has inner ring number:" << inner_num << std::endl;
            std::list<int> outer_ring;
            std::list<std::list<int>> inner_rings;

            // read outer ring
            getline(input_stream, line);
            std::istringstream outer_stream(line);
            int o_vertices;
            outer_stream >> o_vertices;
            std::cout << "Outer ring " << " has " << o_vertices << " vertices:" << std::endl;
            while (outer_stream >> o_vertices) {
                std::cout << o_vertices << " ";
                outer_ring.push_back(o_vertices);
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
                inner_rings.push_back(inner_ring);
                std::cout << std::endl;
            }

            Face face(i, outer_ring, inner_rings);
            face.fitting_plane(vertices);
            face.cdt(vertices);
            mark_domains(face.triangulation);
            faces.insert(std::make_pair(i, face));
        }

    }

    std::ofstream output_stream;
    output_stream.open(output_file);


    int record = 1;
    for (auto &face: faces) {
        int count = 0;
        std::vector<Kernel::Point_2> point_vector;
        for (auto v = face.second.triangulation.finite_vertices_begin(); v != face.second.triangulation.finite_vertices_end(); ++v) {
            Point v_3d= face.second.best_plane.to_3d(v->point());
            point_vector.push_back(v->point());
            output_stream << "v " << v_3d.x()<< " " <<
                          v_3d.y() << " " << v_3d.z() << std::endl;
            count++;
        }

        for (auto t = face.second.triangulation.all_faces_begin();
             t != face.second.triangulation.all_faces_end(); ++t) {
            if (t->info().in_domain()) {
                output_stream <<"f ";
                for (int i = 0; i < 3; ++i) {
                int index;
                Point_2 p = t->vertex(i)->point();
                auto it = std::find(point_vector.begin(), point_vector.end(), p);
                index = int(std::distance(point_vector.begin(), it));
                output_stream <<index  + record << " ";}
                output_stream << "\n";
            }
        }
    record = record + count;
    }

    output_stream << std::endl;
    output_stream.close();


  return 0;
}