#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }
        
        double sinus_fct ( vertex v )
        {
        	double pi = 3.14159265;
        	return 2*(pi*pi)*sin(pi*v.x)*sin(pi*v.y);
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            
            Mesh mesh;
            mesh.load(mesh_filename);
            
            std::vector<double> F_globale(mesh.nb_vertices(), 0.);
            SparseMatrix K_globale(mesh.nb_vertices());

            for ( int a = 0; a < mesh.nb_triangles(); ++a ) {
            	DenseMatrix Ke;
            	ElementMapping mapping(mesh, false, a);
            	Quadrature quadrature = Quadrature::get_quadrature(2);
            	ShapeFunctions shapefunctions(2, 1);
            	assemble_elementary_matrix(mapping, shapefunctions, quadrature, unit_fct, Ke);
            	local_to_global_matrix(mesh, a, Ke, K_globale);
            }
            
            std::vector< bool > attribute_dirichlet(2, false);
            attribute_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector < double > imposed_v(mesh.nb_vertices());
            
            for ( int b = 0; b<mesh.nb_vertices(); ++b) {
            	imposed_v[b] = xy_fct(mesh.get_vertex(b));
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribute_dirichlet, imposed_v, K_globale, F_globale);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_globale, F_globale, u);
            std::string export_name = "square_pure_dirichlet";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name +".bb");
            
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
        }
        
        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            
            Mesh mesh;
            mesh.load(mesh_filename);
            
            std::vector<double> F_globale(mesh.nb_vertices(), 0.);
            SparseMatrix K_globale(mesh.nb_vertices());

            for ( int a = 0; a < mesh.nb_triangles(); ++a ) {
            	DenseMatrix Ke;
            	ElementMapping mapping(mesh, false, a);
            	Quadrature quadrature = Quadrature::get_quadrature(2);
            	ShapeFunctions shapefunctions(2, 1);
            	std::vector< double > Fe(shapefunctions.nb_functions());
            	assemble_elementary_matrix(mapping, shapefunctions, quadrature, unit_fct, Ke);
            	assemble_elementary_vector(mapping, shapefunctions, quadrature, unit_fct, Fe);
            	local_to_global_matrix(mesh, a, Ke, K_globale);
            	local_to_global_vector(mesh, true, a, Fe, F_globale);
            }
            
            std::vector< bool > attribute_dirichlet(2, false);
            attribute_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector < double > imposed_v(mesh.nb_vertices());
            
            for ( int b = 0; b<mesh.nb_vertices(); ++b) {
            	imposed_v[b] = zero_fct(mesh.get_vertex(b));
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribute_dirichlet, imposed_v, K_globale, F_globale);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_globale, F_globale, u);
            std::string export_name = "square_fine_source_dirichlet";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name +".bb");
            
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
        }
        
        void sinus_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            
            Mesh mesh;
            mesh.load(mesh_filename);
            
            std::vector<double> F_globale(mesh.nb_vertices(), 0.);
            SparseMatrix K_globale(mesh.nb_vertices());

            for ( int a = 0; a < mesh.nb_triangles(); ++a ) {
            	DenseMatrix Ke;
            	ElementMapping mapping(mesh, false, a);
            	Quadrature quadrature = Quadrature::get_quadrature(2);
            	ShapeFunctions shapefunctions(2, 1);
            	std::vector< double > Fe(shapefunctions.nb_functions());
            	assemble_elementary_matrix(mapping, shapefunctions, quadrature, unit_fct, Ke);
            	assemble_elementary_vector(mapping, shapefunctions, quadrature, sinus_fct, Fe);
            	local_to_global_matrix(mesh, a, Ke, K_globale);
            	local_to_global_vector(mesh, true, a, Fe, F_globale);
            }
            
            std::vector< bool > attribute_dirichlet(2, false);
            attribute_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector < double > imposed_v(mesh.nb_vertices());
            
            for ( int b = 0; b<mesh.nb_vertices(); ++b) {
            	imposed_v[b] = zero_fct(mesh.get_vertex(b));
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribute_dirichlet, imposed_v, K_globale, F_globale);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_globale, F_globale, u);
            std::string export_name = "square_sinus_bump_dirichlet";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name +".bb");
            
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
        }
        

    }

}
