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

const std::string input_file = "../station.hw1";
const std::string output_file = "../t1.obj";

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
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
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

    void fitting_plane(std::map<int, Point> vertices) {
        std::vector<Point> pts;
        for (auto const &vertex: outer_ring) pts.push_back(vertices[vertex]);
        Plane plane;
        linear_least_squares_fitting_3(pts.begin(), pts.end(), best_plane, CGAL::Dimension_tag<0>());
        std::cout << best_plane << std::endl;
    }

    void cdt(std::map<int, Point> vertices) {
        std::vector<std::pair<int, int>> segments;
        // obtaining vertices
        std::list<int> vertice_ids;
        std::list<Triangulation::Vertex_handle> v_handles;

        Triangulation cdt;

//        std::list<Triangulation::Vertex_handle> outer_segments;
//        std::list<Triangulation::Vertex_handle> inner_segments;

        for (auto it = outer_ring.begin(); it != outer_ring.end(); ++it) {
            auto nextIt = std::next(it);
            if (nextIt == outer_ring.end()) {
                nextIt = outer_ring.begin();
            }
            auto vertexIndex1 = best_plane.to_2d(vertices[*it]);
            auto vertexIndex2 = best_plane.to_2d(vertices[*nextIt]);
            cdt.insert_constraint(vertexIndex1, vertexIndex2);
        }

        for (auto inner: inner_rings) {
            for (auto it = inner.begin(); it != inner.end(); ++it) {
                auto nextIt = std::next(it);
                if (nextIt == inner.end()) {
                    nextIt = inner.begin();
                }
                auto vertexIndex1 = best_plane.to_2d(vertices[*it]);
                auto vertexIndex2 = best_plane.to_2d(vertices[*nextIt]);
                cdt.insert_constraint(vertexIndex1, vertexIndex2);
            }
        }

        for (int vertexIndex_o : outer_ring) {
            Point_2 vertex1 = best_plane.to_2d(vertices[vertexIndex_o]);
            Triangulation::Vertex_handle v1 = cdt.insert(vertex1);
            v1->info() = vertexIndex_o;
//            outer_segments.push_back(v1);
        }

        for (auto ring:inner_rings) {
            for (int vertexIndex_i : ring) {
                Point_2 vert1 = best_plane.to_2d(vertices[vertexIndex_i]);
                Triangulation::Vertex_handle vi1 = cdt.insert(vert1);
                vi1->info() = vertexIndex_i;
//                inner_segments.push_back(vi1);
            }
        }

        triangulation = cdt;

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
            std::cout << "Outer ring " << " has " << o_vertices << " vertices:" << std::endl;
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

    }

    std::ofstream output_stream;
    output_stream.open(output_file);


    std::map<int, Point> points;

    for (auto face: faces) {
        face.second.fitting_plane(vertices);
        face.second.cdt(vertices);
        mark_domains(face.second.triangulation);
        for (auto v = face.second.triangulation.finite_vertices_begin(); v != face.second.triangulation.finite_vertices_end(); ++v) {
            Point v_3d= face.second.best_plane.to_3d(v->point());
            points[v->info()] = v_3d;
        }
        for (auto t = face.second.triangulation.finite_faces_begin();
             t != face.second.triangulation.finite_faces_end(); ++t) {
            if (t->info().in_domain()) {
                output_stream << "f " << t->vertex(0)->info() << " " <<
                              t->vertex(1)->info() << " " << t->vertex(2)->info() << std::endl;
            }
        }
    }

    for (int i = 0; i < points.size(); i++) {
        output_stream << "v " << points[i].x()<< " " <<
                      points[i].y() << " " << points[i].z() << std::endl;
    }

    output_stream << std::endl;

//    for (auto face: faces) {
//
//
//    }


  return 0;
}