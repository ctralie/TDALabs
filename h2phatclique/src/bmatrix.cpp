#include "clique.h"
#include "bmatrix.h"

namespace bmatrix 
{

// This merely adds the dim-0 simplicies to the cliquevector
void add_vertices_to_vector(cliques::index& vertices, cliquevector& cliquevector)
{
    for ( cliques::index i=0; i<vertices; ++i )
    {
        cliques::clique c({i},i);
        cliquevector.emplace_back(c);
    }
}

// This allows us to read the dimension of a simplex from the clique dimension list
cliques::sz_t cliquesize(indexvector& dimbordersext, cliques::index& simplex)
{
	cliques::sz_t d;
	for ( cliques::index j = dimbordersext.size() ; j>0; j-- )
	{
		if (simplex >= dimbordersext[j-1])
		{
			d = j-1;
			break;
		}
	}
	return d;
}

} //end namespace
