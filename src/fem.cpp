#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        if ( border==true ) {
        	for (int v = 0; v < 2; ++v ){
        		vertices_.push_back(M.get_edge_vertex(i,v));
        		//std::cout << "x : " << vertices_[v].x << " et y : " << vertices_[v].y << std::endl;
        	} 
        }
        else {
        	for ( int v = 0; v < 3; ++v ){
        		vertices_.push_back(M.get_triangle_vertex(i,v)) ;
        		//std::cout << "x : " << vertices_[v].x << " et y : " << vertices_[v].y << std::endl;
        	} 
        }
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        vertex r ;
        if ( border_ ) {
        	double xi = x_r.x ;
        	r.x = (1-xi)*vertices_[0].x + xi*vertices_[1].x ;
        	r.y = (1-xi)*vertices_[0].y + xi*vertices_[1].y ;
        }
        else {
        	double xi = x_r.x ;
        	double eta = x_r.y ;
        	r.x = (1-xi-eta)*vertices_[0].x + xi*vertices_[1].x + eta*vertices_[2].x ;
        	r.y = (1-xi-eta)*vertices_[0].y + xi*vertices_[1].y + eta*vertices_[2].y ;
		}
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        DenseMatrix J ;
        if ( border_ ) {
        	J.set_size(2, 1);
        	J.set(0, 0, -vertices_[0].x +vertices_[1].x);
        	J.set(1, 0, -vertices_[0].y +vertices_[1].y);
        }
        else {
        	J.set_size(2, 2);
        	J.set(0, 0, -vertices_[0].x +vertices_[1].x);
        	J.set(0, 1, -vertices_[0].x +vertices_[2].x);
        	J.set(1, 0, -vertices_[0].y +vertices_[1].y);
        	J.set(1, 1, -vertices_[0].y +vertices_[2].y);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        DenseMatrix J = jacobian_matrix(x_r);
        if ( border_ ) {
        	double JTrans_J = J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0);
        	return std::sqrt(JTrans_J);
        }
        else {
        	double JTrans_J  = J.det_2x2();
        	return JTrans_J;
        }
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        bool SF_construct = true;
        if ( dim != 1 && dim != 2 ) {
        	std::cout << "Implémentation en 1D et 2D" << std::endl;
        	SF_construct = false;
        }
        if ( order != 1 ) {
        	std::cout << "Fonctions de forme d'ordre 1 seulement" << std::endl;
        	SF_construct = false;
        }
        assert(SF_construct);
    }

    int ShapeFunctions::nb_functions() const
    {
        if ( dim_ == 1 ) {
        	return 2;
        }
        if ( dim_ == 2 ) {
        	return 3;
        }
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
    	if  ( i == 0 ) {
        	double xi = x_r.x ;
        	switch(i) {
        		case(0):
        			return 1-xi; break;
        		case(1):
        			return xi; break;
        	}
        }
        else {
        	double xi = x_r.x ; double eta = x_r.y;
        	switch(i) {
        		case(0):
        			return 1-xi-eta; break;
        		case(1):
        			return xi; break;
        		case(2):
        			return eta; break;
        	}
        }			
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        vec2 g ;
        if (dim_ == 1 ) {
        	switch(i) {
        		case 0:
        			g.x = -1.; break;
        		case 1:
        			g.x = 1.; break;
        	}
        	g.y = 0. ;
        }
        else {
        	switch(i) {
        		case 0:
        			g.x = -1.; g.y = -1.; break;
        		case 1:
        			g.x = 1.; g.y = 0.; break;
        		case 2:
        			g.x = 0.; g.y = 1.; break;
        	}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());
        for ( int i = 0; i < reference_functions.nb_functions(); ++i ) {
        	for ( int j = 0; j < reference_functions.nb_functions(); ++j ) {
        		Ke.set(i, j, 0.);
        		for (int k=0; k< quadrature.nb_points(); ++k ) {
        			vertex p_k = quadrature.point(k);
        			double w_k = quadrature.weight(k);
        			DenseMatrix inv_J = elt_mapping.jacobian_matrix(p_k).invert_2x2();
        			vec2 grad_i = inv_J.transpose().mult_2x2_2(reference_functions.evaluate_grad(i, p_k));
        			vec2 grad_j = inv_J.transpose().mult_2x2_2(reference_functions.evaluate_grad(j, p_k));
        			Ke.add(i, j, w_k * coefficient(elt_mapping.transform(p_k)) * dot(grad_i, grad_j) * elt_mapping.jacobian(p_k));			
        		}	
        	}
        }
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        for (int i=0; i<3; i++) {
        	for ( int j=0; j<3; j++) {
        		int a = M.get_triangle_vertex_index( t, i );
        		int b = M.get_triangle_vertex_index( t , j);
        		K.add (a,b, Ke.get( i , j ));
        	}
        }	
        	
        	

    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        for ( int i = 0; i < reference_functions.nb_functions(); ++i ) {
        	double s = 0;
        	for (int k=0; k< quadrature.nb_points(); ++k ) {
        		vertex p_k = quadrature.point(k);
        		double w_k = quadrature.weight(k);
        		s += w_k * source(elt_mapping.transform(p_k)) * reference_functions.evaluate(i, p_k) * elt_mapping.jacobian(p_k);			
        	}
        	Fe[i] = s;
        }
    }

    /*void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        for(int a = 0; a < reference_functions_1D.nb_functions(); ++a ){
        	double S = 0;
        	for (int b = 0; b < quadrature_1D.nb_points(); ++b ) {
        		vertex p_q = quadrature_1D.point(b);
        		double w_q = quadrature_1D.weight(b);
        		S = S + w_q * neumann(elt_mapping_1D.transform(p_q)) * reference_functions_1D.evaluate(a, p_q) * sqrt(dot(?));
        	}
        	Fe[a] = S;
        }   
    }*/

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        for (int y=0; y<Fe.size(); ++y) {
        	int a = M.get_triangle_vertex_index( i, y );
        	F[a] = Fe[y];
        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::vector< bool > vertices(values.size(), false);
        double p = 10000.;
        for (int ed = 0; ed < M.nb_edges(); ed++ ) {
        	int edge_attribute = M.get_edge_attribute(ed);
        	if( attribute_is_dirichlet[edge_attribute] ) {
        		for( int vert=0; vert < 2; vert++ ) {
        			int vertex_index = M.get_edge_vertex_index(ed, vert);
        			if ( !vertices[vertex_index] ) {
        				vertices[vertex_index] = true;
        				K.add(vertex_index, vertex_index, p);
        				F[vertex_index] += p*values[vertex_index];
        			}
        		}
        	}
        }
        
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
