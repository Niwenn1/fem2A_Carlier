#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature()
        {
        	Quadrature quad;
        	quad = quad.get_quadrature(2,false);
        	int i=0;
        	float somme = 0;
        	while (i<quad.nb_points()) {
        		somme = somme + quad.weight(i);
        		i = i+1;
        	}
        	std::cout << "somme : " << somme << std::endl;
        	return true;
        }

		bool test_ElementMapping()
		{
			Mesh carre;
			carre.load("data/square.mesh");
			// test sur le triangle 4 de square.mesh
			ElementMapping( carre, false, 4 );
			return true;
		}
		bool test_transform_ElementMaping() {
			Mesh carre;
			carre.load("data/square.mesh");
			vertex v;
			v.x = 0.2;
			v.y = 0.4;
			vertex u;
			ElementMapping U( carre, false, 4 );
			u = U.transform(v);
			std::cout << "x :" << u.x << "y : " << u.y <<std::endl;
			return true;
		}
		bool test_mat_et_jacobien() {
			Mesh carre;
			carre.load("data/square.mesh");
			vertex v;
			v.x = 0.2;
			v.y = 0.4;
			DenseMatrix D;
			ElementMapping EL( carre, false, 4 );
			D = EL.jacobian_matrix( v );
			D.print();
			double det = EL.jacobian (v);
			std::cout << "\n Jacobien : " << det << std::endl;
			return true;
		}
		
		bool test_shape_functions() {
			//ShapeFunctions( 3, 2 );
			ShapeFunctions( 2, 1 );
			//std::cout << ShapeFunctions( 2, 1 ).nb_functions() << std::endl;
			vertex v;
			v.x = 0.2;
			v.y = 0.4;
			std::cout << ShapeFunctions( 2, 1 ).evaluate( 2, v ) << std::endl;
			return true;
		}
		
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
		
		bool test_ass_elmt_matrix() {
			Mesh carre;
			carre.load("data/square.mesh");
			Quadrature quad;
        		quad = quad.get_quadrature(2,false);
			ShapeFunctions SF (2, 1);
			vertex v;
			v.x = 0.2;
			v.y = 0.4;
			DenseMatrix Ke;
			SparseMatrix K(carre.nb_vertices());
			ElementMapping EL( carre, false, 4 );
			assemble_elementary_matrix( EL, SF, quad, unit_fct, Ke);
			//Ke.print();
			// t=4 index du rectangle 4
			local_to_global_matrix( carre, 4, Ke, K);
			K.print();
			return true;
		}
		
		/*bool test_ass_elmt_vector() {
			Mesh carre;
			carre.load("data/square.mesh");
			Quadrature quad;
        		quad = quad.get_quadrature(2,false);
			ShapeFunctions SF (2, 1);
			vertex v;
			v.x = 0.2;
			v.y = 0.4;
			std::vector< double >& Fe;
			ElementMapping EL( carre, false, 4 );
			assemble_elementary_vector( EL, SF, quad, unit_fct, Fe);
			Fe.print();
			return true;
		*/
    	}
    }

