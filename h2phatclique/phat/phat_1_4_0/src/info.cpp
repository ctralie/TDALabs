/*  Copyright 2013 IST Austria
    Contributed by: Jan Reininghaus

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

#include <phat/boundary_matrix.h>

void print_help() {
    std::cerr << "Usage: " << "info " << "[options] input_filename_0 input_filename_1 ... input_filename_N" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--ascii   --  use ascii file format" << std::endl;
    std::cerr << "--binary  --  use binary file format (default)" << std::endl;
    std::cerr << "--help    --  prints this screen" << std::endl;
}

void print_help_and_exit() {
    print_help();
    exit( EXIT_FAILURE );
}

void parse_command_line( int argc, char** argv, bool& use_binary, std::vector< std::string >& input_filenames ) {

    if( argc < 2 ) print_help_and_exit();

    int number_of_options = 0;
    for( int idx = 1; idx < argc; idx++ ) {
        const std::string argument = argv[ idx ];
        if( argument.size() > 2 && argument[ 0 ] == '-' && argument[ 1 ] == '-' ) {
            if( argument == "--ascii" ) use_binary = false;
            else if( argument == "--binary" ) use_binary = true;
            else if( argument == "--help" ) print_help_and_exit();
            else print_help_and_exit();
        } else {
            input_filenames.push_back( argument );
        }
    }
}

int main( int argc, char** argv )
{
    bool use_binary = true; // interpret inputs as binary or ascii files
    std::vector< std::string > input_filenames; // name of file that contains the boundary matrix

    parse_command_line( argc, argv, use_binary, input_filenames );

    for( int idx_input = 0; idx_input < input_filenames.size(); idx_input++ ) {
        std::string input_filename = input_filenames[ idx_input ];
        phat::boundary_matrix<> matrix;
        bool read_successful = use_binary ? matrix.load_binary( input_filename ) : matrix.load_ascii( input_filename );
        if( !read_successful ) {
            std::cerr << std::endl << " Error opening file " << input_filename << std::endl;
            print_help_and_exit();
        }

        std::cout << input_filename << ":" << std::endl;
        std::cout << "\t" << "num_cols = " << matrix.get_num_cols() << std::endl;
        std::cout << "\t" << "max_dim = " << (int)matrix.get_max_dim() << std::endl;
        std::cout << "\t" << "max_col_entries = " << matrix.get_max_col_entries() << std::endl;
        std::cout << "\t" << "max_row_entries = " << matrix.get_max_row_entries() << std::endl;
        std::cout << "\t" << "total_nr_of_entries = " << matrix.get_num_entries() << std::endl;
    }
}
