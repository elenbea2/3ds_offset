// libraries
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h> // represents polygonal solids
#include <CGAL/Nef_polyhedron_3.h> // structure for boolean operations 
#include <CGAL/Surface_mesh.h> // structure for 3d meshes
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/boost/graph/copy_face_graph.h> // copies one mesh type to another
#include <boost/property_map/property_map.hpp>
#include <CGAL/convex_decomposition_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/squared_distance_3.h>
#include <cstdlib>
#include <CGAL/Polygon_mesh_processing/self_intersections.h> // detects self intersections
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // Para fechar buracos
#include <CGAL/Polygon_mesh_processing/repair.h> // Outras funções de reparo
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Surface_mesh<Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;



void repair_mesh(Mesh& mesh) {
    // Costura bordas para unir componentes desconectados
    PMP::stitch_borders(mesh);

    // Identifica e preenche buracos na malha
    std::vector<Mesh::halfedge_index> border_cycles;
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

    for (const auto& cycle : border_cycles) {
        PMP::triangulate_hole(mesh, cycle);
    }

    // Remove faces quase degeneradas
    PMP::remove_almost_degenerate_faces(mesh);

    // Remove componentes conectados de tamanho negligenciável
    PMP::remove_connected_components_of_negligible_size(mesh);
}


// function to remove thin regions
// min thickness: determines the minimum distance between the faces
void remove_small_faces(Mesh& mesh, double min_area) {
    std::vector<Mesh::Face_index> faces_to_remove;
    
    for (auto f : mesh.faces()) {
        std::vector<Point_3> face_points;

        for (auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            face_points.push_back(mesh.point(v));
        }

        if (face_points.size() < 3) continue; // Garante que seja um polígono válido

        // Calcula a área da face usando a triangulação
        double area = 0.0;
        for (size_t i = 1; i + 1 < face_points.size(); ++i) {
            area += std::sqrt(CGAL::to_double(CGAL::squared_area(face_points[0], face_points[i], face_points[i + 1])));
        }

        if (area < min_area) {
            faces_to_remove.push_back(f);
        }
    }

    for (auto f : faces_to_remove) {
        CGAL::Euler::remove_face(mesh.halfedge(f), mesh);
    }

    std::cout << "Small regions removed: " << faces_to_remove.size() << " faces." << std::endl;
}


// apply the offset to the mesh
void offset_mesh(Mesh& mesh, double step, double offset_value, double min_thickness) {
    std::map<Mesh::Vertex_index, Vector_3> vertex_normals;
    
    
    std::cout << "Calculating vertex normal..." << std::endl;
    
    // calculates the vertex normal
    for (auto v : mesh.vertices()) {
        vertex_normals[v] = PMP::compute_vertex_normal(v, mesh);
    }
    
    // initializes the offset number
    double applied_offset = 0.0;
    
    // limit to avoid unnecessary applications
    const double epsilon = 1e-10;
    
    // apply the offset iteratively
    while (std::abs(applied_offset - offset_value) > epsilon) {

        double current_step = std::min(step, std::abs(offset_value - applied_offset));
        
        // case for the negative offset
        if (offset_value < 0) { current_step = -current_step; }
        
        // avoid unnecessary applications
        if (std::abs(current_step) < epsilon) break; 
        
        std::cout << "Applying offset: " << current_step << std::endl;
        
        // apply the offset for each vertex
        for (auto v : mesh.vertices()) {
            Vector_3 normal = vertex_normals[v];
  
            double length = std::sqrt(CGAL::to_double(normal.squared_length()));
            
            // normalization
            if (length > 0) {
                normal = normal / length;
            } else {
                std::cerr << "Zero-length normal, ignoring normalization.\n";
                continue;
            }
            
            // calculates the new position of the point
            Vector_3 offset_direction = (current_step > 0) ? normal : -normal;
            mesh.point(v) = mesh.point(v) + offset_direction * std::abs(current_step);
        }
        
        // checks for thin regions 
        remove_small_faces(mesh, min_thickness);
        
        // update the iterator
        applied_offset += current_step;
      
        repair_mesh(mesh);
    }
    //repair_mesh(mesh);
}

// verifies if the polyhedron is valid
bool is_valid_polyhedron(const Polyhedron& poly) {
    return poly.is_valid() && !poly.empty();
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <offset> <model_file>" << std::endl;
        return 1;
    }

    double offset_value;
    try {
        offset_value = std::stod(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "Error: invalid offset value!" << std::endl;
        return 1;
    }
    
    std::string input_filename = argv[2];
    std::ifstream input_model(input_filename);
    Polyhedron P;
    if (!input_model || !(input_model >> P) || P.empty()) {
        std::cerr << "Error loading the model: " << input_filename << std::endl;
        return 1;
    }
    input_model.close();
    std::cout << "Model successfully loaded!" << std::endl;

    Nef_polyhedron nef_P(P);
    Polyhedron poly_offset;
    Mesh mesh_part;
    repair_mesh(mesh_part);
    CGAL::copy_face_graph(P, mesh_part);
    

    double step = 0.001; 
    double min_thickness = 0.01; 
    offset_mesh(mesh_part, step, offset_value, min_thickness);


    CGAL::copy_face_graph(mesh_part, poly_offset);

    if (!is_valid_polyhedron(poly_offset)) {
        std::cerr << "Error: Invalid polyhedron after the offset!" << std::endl;
        return 1;
    }

    std::ofstream debug_file("offset.off");
    debug_file << poly_offset;
    debug_file.close();
    std::cout << "'offset.off' file successfully saved!" << std::endl;

    return 0;
}
