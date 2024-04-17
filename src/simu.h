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

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            Mesh M;
            M.load(mesh_filename);
            SparseMatrix K_glob(M.nb_vertices());
            std::vector<double> F_glob(M.nb_vertices(), 0.);
            for ( int t = 0; t < M.nb_triangles(); ++t ) {
            	ElementMapping my_map(M, false, t);
            	ShapeFunctions my_shpfct(2, 1);
            	Quadrature my_quad = Quadrature::get_quadrature(2);
            	DenseMatrix Ke;
            	assemble_elementary_matrix(my_map, my_shpfct, my_quad, unit_fct, Ke);
            	local_to_global_matrix(M, t, Ke, K_glob);
            }
            std::vector< bool > att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            M.set_attribute(unit_fct, 1, true);
            std::vector < double > imposed_values(M.nb_vertices());
            for ( int i = 0; i<M.nb_vertices(); ++i) {
            	imposed_values[i] = xy_fct(M.get_vertex(i));
            }
            apply_dirichlet_boundary_conditions(M, att_is_dirichlet, imposed_values, K_glob, F_glob);
            std::vector<double> u(M.nb_vertices());
            solve(K_glob, F_glob, u);
            std::string export_name = "pure_dirichlet";
            M.save(export_name+".mesh");
            save_solution(u, export_name +".bb");
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
        }
	bool test_pure_dirichlet_pb() {
			pure_dirichlet_pb( "data/square.mesh", false );
			return true;
		}
    }

}
