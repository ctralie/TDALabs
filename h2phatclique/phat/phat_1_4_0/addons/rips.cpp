#include <topology/rips.h>
#include <topology/filtration.h>

#include <geometry/l2distance.h>
#include <geometry/distances.h>

#include <utilities/containers.h>           // for BackInsertFunctor

#include <vector>

#include <phat/boundary_matrix.h>

typedef PairwiseDistances<PointContainer, L2Distance> Pair_distances;
typedef Pair_distances::DistanceType Distance_type;
typedef Pair_distances::IndexType Vertex;

typedef Rips<Pair_distances> Generator;
typedef Generator::Simplex Smplx;
typedef Filtration<Smplx> Fltr;

int main(int argc, char* argv[])
{
    Dimension skeleton;
    Distance_type max_distance;
    std::string infilename;

    if(argc<2) {
        std::cerr << "Requires inputfile" << std::endl;
        std::exit(1);
    }
    infilename = argv[1];

    if(argc>=3) {
        skeleton=atoi(argv[2]);
        if(skeleton==0) {
            std::cerr << "# Command line argument 0 ignored" << std::endl;
            skeleton=std::numeric_limits<Dimension>::max();
        }
    } else {
        skeleton=std::numeric_limits<Dimension>::max();
    }

    if(argc>=4) {
        max_distance=atof(argv[3]);
    } else {
        max_distance=std::numeric_limits<Distance_type>::max();
    }

    PointContainer points;
    read_points(infilename, points);

    Pair_distances           distances(points);
    Generator               rips(distances);
    Generator::Evaluator    size(distances);
    Fltr                    f;
    
    // Generate 2-skeleton of the Rips complex for epsilon = 50
    rips.generate(skeleton, max_distance, make_push_back_functor(f));
    std::cerr << "# Generated complex of size: " << f.size() << std::endl;

    // Generate filtration with respect to distance and compute its persistence
    f.sort(Generator::Comparison(distances));


    std::map<Smplx,phat::index,Smplx::VertexComparison> simplex_map;
    phat::index size_of_simplex_map=0;
    
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    boundary_matrix.set_num_cols( f.size() );
    
    for(Fltr::Index it=f.begin();it!=f.end();it++) {
        phat::column boundary_indices;
        const Smplx& c = f.simplex(it);
        for(Smplx::BoundaryIterator bit = c.boundary_begin(); bit != c.boundary_end(); bit++) 
            boundary_indices.push_back( simplex_map[*bit] );
        std::sort(boundary_indices.begin(),boundary_indices.end());

        boundary_matrix.set_col( size_of_simplex_map, boundary_indices );

        phat::dimension dim_of_column = boundary_indices.size()==0 ? 0 : boundary_indices.size()-1;

        boundary_matrix.set_dim( size_of_simplex_map, dim_of_column );

        simplex_map[c] = size_of_simplex_map++;
    }

    boundary_matrix.save_binary( "rips.bin" );

    return 1;

}
