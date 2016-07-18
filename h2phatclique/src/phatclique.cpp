#include <mytime.h>
#include <bmatrix.h>
#include <file.h>
#include <ck.h>

#ifndef _OPENMP
	#include <omp.h>
#endif // _OPENMP

#ifdef __APPLE__
	typedef unsigned int uint;
#endif

// We're not interested in clique sizes below 3 (triangles), as we already have them.
static int MinSize = 3;

//  This global will hold cliques of dim > 2
bmatrix::cliquearray highercliques;

// This global will hold the number of simplicies of each dimension
bmatrix::indexvector dimborders;

void find_inf_pers( phat::persistence_pairs& pairs, bmatrix::indexvector& dimborders )
{
	bmatrix::indexarray finpers(dimborders.size()-1,bmatrix::indexvector({}));
	bmatrix::indexarray infpers(dimborders.size()-1,bmatrix::indexvector({}));

	cliques::sz_t d;
	std::pair<cliques::index,cliques::index> p;
	bmatrix::indexvector dimbordersext = dimborders;
	dimbordersext.insert(dimbordersext.begin(),0);

	for (cliques::index i = 0; i<pairs.get_num_pairs(); ++i )
	{
		p = pairs.get_pair(i);
		d = bmatrix::cliquesize(dimbordersext,p.first);

		finpers[d].push_back(p.first);
		if ( d<(dimborders.size()-2) ) { finpers[d+1].push_back(p.second); }
	}

	for ( unsigned i=0; i<finpers.size(); ++i )
	{
		__gnu_parallel::sort(finpers[i].begin(),finpers[i].end());
//		#pragma omp parallel for schedule(guided) shared(infpers)
		for ( unsigned j=dimbordersext[i]; j<dimbordersext[i+1]; j++ )
		{
			// If a simplex is does not have finite persistence, it has infinite pers.
			// Note: this could take a very long time if there's an infinite pers.
			//       simplex with high simplex number, but the commented method below
			//       takes a ton of memory.
			if ( !(std::binary_search(finpers[i].begin(),finpers[i].end(),j)) )
				{ infpers[i].emplace_back(j); }
		}
/*
		bmatrix::indexvector range(dimbordersext[i+1]-dimbordersext[i]);
		std::iota(range.begin(),range.end(),dimbordersext[i]);

		infpers[i].reserve(range.size() - finpers[i].size());
		auto it=std::set_difference(range.begin(),range.end(),finpers[i].begin(),finpers[i].end(),infpers[i].begin());
		infpers[i].resize(it-infpers[i].begin());
*/
	}

	for ( auto i : infpers )
		{ for ( auto j : i) { pairs.append_pair(j,-1); }}
}


void run_phat(bmatrix::bmatrix& boundarymatrix,phat::persistence_pairs& pairs)
{
	phat::boundary_matrix< phat::bit_tree_pivot_column > phatboundarymatrix;

 	cliques::index bmatsize=boundarymatrix.size();
	phatboundarymatrix.set_num_cols( bmatsize );

	std::vector< phat::index > temp_col;
	std::vector< phat::index > zero_col = temp_col;

	cliques::index count=0;

	for ( auto c : boundarymatrix.get_boundaryvector() )
	{
		phatboundarymatrix.set_dim(count,(c->get_size())-1);
		if ( c->get_weight() == 0 ) { phatboundarymatrix.set_col(count,zero_col); }

		else {phatboundarymatrix.set_col(count,c->get_boundary()); }

		count++;

		// This clears the boundary so that the clique takes less memory....
		// It saves ~12% memory at a cost of adding 50% time to the PHAT compuation
		//cliques::boundary().swap(c->bdry);

	}

	bmatrix::indexvector dimborders = boundarymatrix.get_dimborders();
	boundarymatrix.destroyMe();

	// TODO: make the algorithm user-selectable
	// Is there any point in allowing dualization to be user-selectable?  Seems to always be faster.
	dualize( phatboundarymatrix );
	phat::compute_persistence_pairs< phat::chunk_reduction >( pairs, phatboundarymatrix );

	phatboundarymatrix = phat::boundary_matrix< phat::bit_tree_pivot_column >();

	phat::dualize_persistence_pairs( pairs, bmatsize );

	// Compute infinite persistence pairs
	find_inf_pers(pairs,dimborders);

	pairs.sort();
}

int main(int argc, char **argv)
{
	std::string infile, outputfile, boundaryfile;
	int maxsize;
	mytimer::timer timeme;

	file::parse_command_line(argc, argv, infile,maxsize,outputfile,boundaryfile);

	std::ifstream inputfile(infile);
	std::string dimacs_error;
	bmatrix::cliquevector cliquevector;
	bmatrix::bmatrix::bmatrixvector cliquepointers;

	omp_set_num_threads(omp_get_max_threads());
	// Use the below for debugging
	//omp_set_num_threads(1);

	// Read DIMACS file (with weights)
	timeme.start("Reading input file....");
	ck::sparse* gr = file::new_read_extdimacs_file(inputfile, dimacs_error, cliquevector, dimborders);

	if ( dimacs_error != "")
	{
		std::cerr << dimacs_error << "\n";
		return 1;
	}
	timeme.stop();

	// TODO: add vertices and edges to bmatrix here?

	timeme.start("Computing Clique Complex....");

	// Compute how many dim's above 2 we're working on
	int numhigherdim=maxsize-MinSize+1;

	// Initialize vector that holds higher dimensional cliques
	for ( int i=0; i<numhigherdim; ++i)
		{ bmatrix::cliquevector v; highercliques.emplace_back(v); }

	// ck initialization
	ck::mkgraph(gr);
	ck::kcore(gr);
	ck::mkspecial(gr);

	// run ck to count cliques
	ck::singlepass(gr,maxsize,highercliques,MinSize,dimborders,1);

	// reserve space for cliques of size > 2
	// TODO: can we insert cliques straight into cliquevector?
	//		This would save a ton of memory
	// 		I believe this is difficult at best, as random inserts into
	//		std::vector's are expensive.
	for (int i = 2; i<dimborders.size(); i++)
		{ highercliques[i-2].reserve(dimborders[i]-dimborders[i-1]); }

	// rerun ck to record cliques
	ck::singlepass(gr,maxsize,highercliques,MinSize,dimborders,0);

	ck::freesparse(gr);

	// Concatenate the higher dimensional vectors, releasing memory as we go
	for ( cliques::sz_t i=0; i<highercliques.size(); ++i )
	{
		cliquevector.reserve(cliquevector.size() + highercliques[i].size());
		cliquevector.insert(cliquevector.end(),highercliques[i].begin(),highercliques[i].end());
		// Release memory
		bmatrix::cliquevector().swap(highercliques[i]);
	}

	// release memory allocated to highercliques itself, as it's now in cliquevector.
	bmatrix::cliquearray().swap(highercliques);

	// Create vectors of pointers to cliques (and set clique numbers)
	cliquepointers.reserve(cliquevector.size());
	cliques::index counter = 0;
	for (auto it=cliquevector.begin(); it!=cliquevector.end(); ++it )
	{
		it->set_cliquenum(counter);
		cliquepointers.emplace_back(&(*it));
		counter++;
	}
	timeme.stop();

	// Initialize boundary matrix structure
	timeme.start("Transforming boundary matrix from vertices to faces....");
	bmatrix::bmatrix boundarymatrix(cliquepointers,dimborders);

	// Release memory
	bmatrix::bmatrix::bmatrixvector().swap(cliquepointers);

	// Rewrite matrix using faces instead of vertices, compute weights, and sort by them.
	boundarymatrix.transformboundarymatrix();

	// Release memory
	boundarymatrix.destroymap();
	timeme.stop();

	if ( !boundaryfile.empty() )
	{
		timeme.start("Outputting boundary matrix to ");
		std::cout << boundaryfile << "...." << std::flush;
		boundarymatrix.output_bmatrix(boundaryfile);
		timeme.stop();
	}

	timeme.start("Reducing boundary matrix using PHAT....");
	std::vector<cliques::weight_t> weightvector(boundarymatrix.size());
	for (cliques::index i = 0; i<boundarymatrix.size(); i++)
		{ weightvector[i] = boundarymatrix[i]->get_weight(); }

	phat::persistence_pairs pairs;
	run_phat(boundarymatrix, pairs);
	timeme.stop();

	// output pairs to file
	timeme.start("Outputting persistence pairs to ");
	std::cout<< outputfile << "...." << std::flush;
	file::output_pairs(weightvector, dimborders, pairs, outputfile);
	timeme.stop();

	return 0;
}
