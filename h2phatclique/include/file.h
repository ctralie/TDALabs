#pragma once

#ifndef _FILE
	#define _FILE
#endif

#include "ck.h"
#include "bmatrix.h"
#include "clique.h"

#include <phat/compute_persistence_pairs.h>

// New read file function
namespace file
{
ck::sparse* new_read_extdimacs_file(std::istream& dimacs, std::string& dimacs_err, bmatrix::cliquevector& cliquevector, bmatrix::indexvector& dimborders) {

    ck::sparse *gr = new ck::sparse; 

    cliques::index vertices = 0, edges = 0, cliques_added=0;
    cliques::weight_t weight = 0;
    cliques::index from, to;
    ck::edge tempedge;

    std::string line;
    while (getline(dimacs, line))
    {
     	std::istringstream iss(line);
        char ch;
        if (iss >> ch)
        {

            std::string format;

            switch(ch) {
                case 'c': break;
                case 'p':
                    if (vertices||edges)
                        { dimacs_err = "Corrupt file: More than one line starting with p!"; break; }
                    if (iss >> format >> vertices >> edges)
                    {
                     	if ("edge" != format)
                            { dimacs_err = "Corrupt file: Line beginning with p does not have 'edge' as its next word!"; }
                        else
                            {
                             	cliquevector.reserve(vertices+edges);
                                bmatrix::add_vertices_to_vector(vertices,cliquevector);
                                cliques_added += vertices;
                                dimborders.push_back(vertices);
                                dimborders.push_back(vertices+edges);
                                gr->n = vertices;
                                gr->e = edges;
                                gr->edges = std::vector<ck::edge>();
                                gr->edges.reserve(edges);
                            }
                    }
                    break;
                case 'e':
                    if ( edges-- )
                    {
                     	if (iss >> from >> to >> weight )
                            {
                             	// for each edge in the file, add it both the graph and$
                                if (to>from) {
                                        tempedge.s=from;
                                        tempedge.t=to;
                                }
                                else {
                                      	tempedge.s=to;
                                        tempedge.t=from;
				}
                                tempedge.weight=weight;
                                gr->edges.push_back(tempedge);

                                cliques::clique c({from,to},cliques_added,weight);
                                cliquevector.emplace_back(c);
                                ++cliques_added;
                            }
                        else
                            { dimacs_err = "Corrupt file: Line beginning with e does not contain two vertices and a weight!"; }
                    }
                    break;

                default:
                    dimacs_err = "Corrupt file: Line begins with unknown character (not c, p, or e)!";
            }
            if ( dimacs_err != "" )
                { break; }
        }
    }

    // Sort edges by weight and rewrite clique numbers accordindly
    __gnu_parallel::sort(cliquevector.begin()+vertices,cliquevector.end(),cliques::CliqueComparerByWeight());
   for ( cliques::index i = vertices; i<cliques_added; ++i )
        { cliquevector[i].set_cliquenum(i); }

   return gr;
}

void program_usage(char** argv)
{
            std::cout << "Program usage: " << argv[0] << " -i <input file> -m <max size> -o <pairs output file> [-b <boundary output file>]" << std::endl;
            exit(1);
}

bool file_exists(std::string& filename)
{
	if (std::ifstream(filename))
		{ return true; }
}

void parse_command_line( int argc, char** argv, std::string& infile, int& maxsize, std::string& outputfile, std::string& boundaryfile )
{

	// Check for correct number of command line arguments
	if ( !(argc == 7 | argc == 9) )
		{ program_usage(argv); }

	std::vector<std::string> allArgs(argv, argv + argc);

	for ( int i=1; i<argc; i=i+2 )
	{
		if ( allArgs[i] == "-i" )
		{
			infile = allArgs[i+1];
			if ( !std::ifstream(infile) )
			{
				std::cout << "Input file " << argv[i+1] << " does not exist!" << std::endl;
				program_usage(argv);
			}
			continue;
		}

		if ( allArgs[i] == "-m" )
		{
			try {
				maxsize = std::stoi(allArgs[i+1]);
			} catch (const std::exception &e) {
				std::cout << "Argument following -m must be an integer!" << std::endl;
				program_usage(argv);
			}
			continue;
		}

		// TODO: implement error checking for output files
		if ( allArgs[i] == "-o" )
			{ outputfile = allArgs[i+1]; continue; }

		if ( allArgs[i] == "-b" )
			{ boundaryfile = allArgs[i+1]; continue; }

		else
			{ program_usage(argv); }
	}
}

//void output_pairs(bmatrix::bmatrix& boundarymatrix, phat::persistence_pairs& pairs, std::string filename)
void output_pairs(std::vector<cliques::weight_t>& weightvector, bmatrix::indexvector dimborders, phat::persistence_pairs& pairs, std::string filename)
{
	std::ofstream output_stream( filename.c_str() );

	//bmatrix::indexvector dimborders = boundarymatrix.get_dimborders();
	bmatrix::indexvector dimbordersext = dimborders;
	dimbordersext.insert(dimbordersext.begin(),0);

	cliques::index birth,death,dim;
	cliques::weight_t birthweight,deathweight;

    for ( cliques::index i=0; i<pairs.get_num_pairs(); ++i )
    {
	birth = pairs.get_pair(i).first;
	death = pairs.get_pair(i).second;
	//birthweight = boundarymatrix[birth]->get_weight();
	birthweight = weightvector[birth];
	dim = bmatrix::cliquesize(dimbordersext,birth);

        if ( death != -1 )
//		{ deathweight = boundarymatrix[death]->get_weight(); }
		{ deathweight = weightvector[death]; }
	else
		{deathweight = -1; }

	// only output pairs with pers > 0
        if ( deathweight-birthweight != 0 )
		{ output_stream << dim << " " <<  birth << " " << death << " " << birthweight << " " << deathweight << std::endl; }
    }
}

} //end namespace
