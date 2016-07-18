#pragma once
#ifndef _BMATRIX
	#define _BMATRIX
#endif

#include <clique.h>
#include <unordered_map>

#include <algorithm>

#ifndef _OPENMP
	#include <omp.h>
#endif // _OPENMP

#include <parallel/algorithm>
#include <fstream>

namespace bmatrix
{
typedef std::vector<cliques::clique> cliquevector;
typedef std::vector<cliquevector> cliquearray;
typedef std::vector<cliques::index> indexvector;
typedef std::vector<indexvector> indexarray;

class bmatrix
{
    public:
        typedef std::unordered_map<cliques::boundary, cliques::clique*, cliques::bHash> bmatrixmap;
        typedef std::vector<cliques::clique*> bmatrixvector;

    private:
        bmatrixmap      bmat;
        std::vector<cliques::index> dimborders;

    public:
        bmatrixvector   bvec;

    // Constructors and operators
    public:
        bmatrix():bmat({}),bvec({}) {};

        // If called with a vector of pointers to cliques, add all of the cliques
        explicit bmatrix(bmatrixvector& cliquematrix, std::vector<cliques::index>& dimborders):bmat({}),bvec(cliquematrix.begin(),cliquematrix.end())
            {
		set_dimborders(dimborders);
		// Reserve space for map for cliques of size maxsize-1
                bmat.reserve( *(dimborders.rbegin()+1) );
                createmap( bvec, bmat);
            }

	~bmatrix() {}

	void destroyMe() 
	{
		bmatrixmap().swap(bmat);
		bmatrixvector().swap(bvec);
	}

        // locate clique by a boundary, specified as a vector or brace-enclosed list
        template<class T>
        cliques::clique* operator[](T b)
            { return ( bmat[b] ); }

        // locate clique by its index
        cliques::clique* operator[](cliques::index i)
            { return ( bvec[i] ); }

        bmatrixvector& get_boundaryvector ()
            { return bvec; }

        void output_bmatrix(std::string filename)
        {
            std::ofstream output_stream( filename.c_str() );

            for ( auto c : bvec )
            {
                output_stream << (c->get_size())-1 << " ";
                if ( c->get_size() > 1 )
                {
                    for ( auto i : c->get_boundary() )
                        { output_stream << i << " "; }
                }
                output_stream << std::endl;
            }
        }

	void set_dimborders(std::vector<cliques::index> db)
	    { dimborders = db; }

	std::vector<cliques::index> get_dimborders()
	    { return dimborders; }

        cliques::index size()
            { return bvec.size(); }

        void destroymap()
            { bmatrixmap().swap (bmat); }

    // Methods
    private:

        void createmap(bmatrixvector& cliquematrix, bmatrixmap& bmat_priv)
        {
            // NOTE: this only works if our vectors are sorted by size!
            cliques::sz_t maxdim = cliquematrix.back()->get_size();
	    // We only need to create the map for dimensions up to maxsize-1
            for ( unsigned i=0; i<dimborders[maxdim-2]; ++i )
            {
                cliques::clique* c = cliquematrix[i];
		const cliques::boundary b = c->get_boundary();
		bmat_priv[b] = c;
	    }
	}

        void verticestofaces(cliques::clique*& c)
        {

            cliques::sz_t sz = c->get_size();

            if ( sz > 2 )
            {
                cliques::weight_t wt = 0;

                cliques::boundary vbdry = c->get_boundary();

                // facebdry will contain the faces of c
                cliques::boundary facebdry(sz);

                // temporary holder for faces of c listed as vertices
                cliques::boundary t(sz-1);

                // compute the boundary of c as faces, one face at a time
                for ( cliques::sz_t i=0; i<sz; ++i )
                {
                    // These for loops remove the i'th element from vbdry,
                    // resulting in a list of vertices of a face of c.
                    for ( cliques::sz_t j=0; j<i; ++j )
                            { t[j]=vbdry[j]; }
                    for ( cliques::sz_t j=i+1; j<sz; ++j )
                            { t[j-1]=vbdry[j]; }

                    // Look up the clique corresponding to the face
                    cliques::clique* face = bmat[t];

                    // Check if this face has greater weight.  If so, make it the new weight.
                    cliques::weight_t faceweight = face->get_weight();
                    if ( faceweight > wt )
                            { wt = faceweight; }

                    // Store the face's clique number in temp boundary vector
                    facebdry[sz-i-1] = face->get_cliquenum();
                }

                // Set properties for clique
                c->set_boundary(facebdry);
		c->set_weight(wt);
            }
        }

        // This gives an ordering on pointers to cliques by the weights of cliques pointed to
        struct CliquePointerComparerByWeight
        {
            bool operator() ( const cliques::clique* lhs, const cliques::clique* rhs ) const
                { return (lhs->get_weight()) < (rhs->get_weight()); }
        };

    public:

        // This function assumes that bvec is sorted by dimension.
        void transformboundarymatrix()
	{
		cliques::sz_t maxsize = bvec.back()->get_size();

		// start the transformation from triangles, as edges are already set

		cliques::index mincliquenum = dimborders[1];

	    	// Note: i in this next for loop is (clique size  - 2)
            	for (cliques::sz_t i = 1; i<(maxsize-1); ++i )
            	{
                	cliques::index maxcliquenum = dimborders[i+1];

                	// Rewrite cliques of size (i+2) with faces instead of vertices
                	// This also computes the weights
                	#pragma omp parallel for schedule(guided)
                	for ( auto c = bvec.begin()+mincliquenum; c < bvec.begin()+maxcliquenum; c++ )
                    		{ verticestofaces(*c); }

	                // Sort vector of pointers to cliques of size (i+2) by weight of cliques pointed to
        	        __gnu_parallel::sort(bvec.begin()+mincliquenum,bvec.begin()+maxcliquenum,CliquePointerComparerByWeight());

                	// Reset clique numbers for sorted cliques of size (i+2)
                	if ( i != (maxsize-2) )
                	{
	                	#pragma omp parallel for schedule(guided)
        	        	for ( cliques::index j=mincliquenum; j<maxcliquenum; ++j )
                	    		{ bvec[j]->set_cliquenum(j); }

                		mincliquenum = maxcliquenum;
			}
            	}
	}

};

void add_vertices_to_vector(cliques::index& vertices, std::vector<cliques::clique>& cliquevector);
cliques::sz_t cliquesize(std::vector<cliques::index>& dimbordersext, cliques::index& simplex);

} // end namespace
