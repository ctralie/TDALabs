/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

// This file contains a simple example that demonstrates the usage of the library interface

// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>
#include <phat/representations/bit_tree_pivot_column.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

int main( int argc, char** argv ) 
{
    std::cout << "We will build an ordered boundary matrix of this simplicial complex consisting of a single triangle: " << std::endl;
    std::cout << std::endl;
    std::cout << " 3" << std::endl;
    std::cout << " |\\" << std::endl;
    std::cout << " | \\" << std::endl;
    std::cout << " |  \\" << std::endl;
    std::cout << " |   \\ 4" << std::endl;
    std::cout << "5|    \\" << std::endl;
    std::cout << " |     \\" << std::endl;
    std::cout << " |  6   \\" << std::endl;
    std::cout << " |       \\" << std::endl;
    std::cout << " |________\\" << std::endl;
    std::cout << " 0    2    1" << std::endl;


    // first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix< phat::bit_tree_pivot_column > boundary_matrix;

    // set the number of columns (has to be 7 since we have 7 simplices)
    boundary_matrix.set_num_cols( 7 );
    
    // set the dimension of the cell that a column represents:
    boundary_matrix.set_dim( 0, 0 );
    boundary_matrix.set_dim( 1, 0 );
    boundary_matrix.set_dim( 2, 1 );
    boundary_matrix.set_dim( 3, 0 );
    boundary_matrix.set_dim( 4, 1 );
    boundary_matrix.set_dim( 5, 1 );
    boundary_matrix.set_dim( 6, 2 );

    // set the respective columns -- the columns entries have to be sorted
    std::vector< phat::index > temp_col;

    boundary_matrix.set_col( 0, temp_col );

    boundary_matrix.set_col( 1, temp_col );

    temp_col.push_back( 0 );
    temp_col.push_back( 1 );
    boundary_matrix.set_col( 2, temp_col );
    temp_col.clear();

    boundary_matrix.set_col( 3, temp_col );

    temp_col.push_back( 1 );
    temp_col.push_back( 3 );
    boundary_matrix.set_col( 4, temp_col );
    temp_col.clear();

    temp_col.push_back( 0 );
    temp_col.push_back( 3 );
    boundary_matrix.set_col( 5, temp_col );
    temp_col.clear();

    temp_col.push_back( 2 );
    temp_col.push_back( 4 );
    temp_col.push_back( 5 );
    boundary_matrix.set_col( 6, temp_col );
    temp_col.clear();

    // print some information of the boundary matrix:
    std::cout << std::endl;
    std::cout << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns: " << std::endl;
    for( phat::index col_idx = 0; col_idx < boundary_matrix.get_num_cols(); col_idx++ ) {
        std::cout << "Colum " << col_idx << " represents a cell of dimension " << (int)boundary_matrix.get_dim( col_idx ) << ". ";
        if( !boundary_matrix.is_empty( col_idx ) ) {
            std::vector< phat::index > temp_col;
            boundary_matrix.get_col( col_idx, temp_col ); 
            std::cout << "Its boundary consists of the cells";
            for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
                std::cout << " " << temp_col[ idx ];
        }
        std::cout << std::endl;
    }
    std::cout << "Overall, the boundary matrix has " << boundary_matrix.get_num_entries() << " entries." << std::endl;  
    

    // define the object to hold the resulting persistence pairs
    phat::persistence_pairs pairs;

    // choose an algorithm (choice affects performance) and compute the persistence pair
    // (modifies boundary_matrix)
    phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );
    
    // sort the persistence pairs by birth index 
    pairs.sort();

    // print the pairs:
    std::cout << std::endl;
    std::cout << "There are " << pairs.get_num_pairs() << " persistence pairs: " << std::endl;
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
        std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
}
