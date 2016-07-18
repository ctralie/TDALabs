#include <phat/boundary_matrix.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <fstream>
#include <list>
#include <cassert>
#include <sstream>

template<typename Index_>
struct Vertex_info_3 {

    typedef Index_ Index;

    Vertex_info_3() {
        index_=boost::none;
    }

    bool has_index() {
        return index_;
    }

    Index index() {
        CGAL_assertion(has_index());
        return index_.get();
    }

    void set_index(Index I) {
        index_=I;
    }

private: 
    boost::optional<Index> index_;
};


template<typename Index_>
struct Cell_info_3 {

    typedef Index_ Index;
    
    Cell_info_3() {
        for(std::size_t i=0; i<6; i++) {
            edge_index_[i] = boost::none;
        }
        for(std::size_t i=0; i<4; i++) {
            facet_index_[i] = boost::none;
        }
    }

    int edge_conv(int i, int j) {
        if(i>j) std::swap(i,j);
        if(i==0 && j==1) return 0;
        if(i==0 && j==2) return 1;
        if(i==0 && j==3) return 2;
        if(i==1 && j==2) return 3;
        if(i==1 && j==3) return 4;
        if(i==2 && j==3) return 5;
    }


    bool has_edge_index(int i, int j) {
        return edge_index_[edge_conv(i,j)];
    }

    Index edge_index(int i, int j) {
        CGAL_assertion(has_edge_index(i,j));
        int k = edge_conv(i,j);
        return edge_index_[k].get();
    }

    bool has_facet_index(int i) {
        CGAL_assertion(i>=0 && i<4);
        return facet_index_[i];
    }

    Index facet_index(int i) {
        CGAL_assertion(has_facet_index(i));
        return facet_index_[i].get();
    }

    void set_edge_index(int i, int j, Index I) {
        edge_index_[edge_conv(i,j)]=I;
    }

    void set_facet_index(int i, Index I) {
        facet_index_[i]=I;
    }

private:

    boost::optional<Index> edge_index_[6];
    boost::optional<Index> facet_index_[4];

};
    
typedef phat::dimension Dim;
typedef std::vector<Dim> Dim_container;
typedef phat::index Index;
typedef phat::column Column;
typedef std::vector<std::vector<Index> > Matrix;

typedef CGAL::Exact_predicates_exact_constructions_kernel Gt;
typedef Gt::FT FT;

typedef CGAL::Triangulation_cell_base_with_info_3<Cell_info_3<Index>,Gt> Cb;
typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_info_3<Index>,Gt> Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_data_structure_3<Vbh,Cb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Triangulation_hierarchy_3<Triangulation_3>       Triangulation_hierarchy;

typedef Gt::Point_3                                  Point;


typedef Triangulation_3::Vertex_handle Vertex_handle;
typedef Triangulation_3::Edge Edge;
typedef Triangulation_3::Facet Facet;
typedef Triangulation_3::Cell_handle Cell_handle;

typedef Triangulation_3::Cell_circulator Cell_circulator;

void set_index_of_edge(const Triangulation_3& T, const Edge& e, Index I) {

    Vertex_handle v1 = e.first->vertex(e.second);
    Vertex_handle v2 = e.first->vertex(e.third);

    Cell_circulator ch=T.incident_cells(e);
    Cell_circulator ch_start=ch;
    int count=0;
    do {
        ch->info().set_edge_index(ch->index(v1),ch->index(v2),I);
        ch++;
        count++;
    } while(ch!=ch_start);
    //std::cout << "Did " << count << " updates" << std::endl;
}

void set_index_of_facet(const Triangulation_3& T, const Facet& f, Index I) {

    f.first->info().set_facet_index(f.second,I);
    Facet mf = T.mirror_facet(f);
    mf.first->info().set_facet_index(mf.second,I);    
    
}

template<typename Triple>
struct Sort_triples {

    bool operator() (const Triple& a, const Triple& b) {
        if(a.first < b.first) return true;
        if(a.first > b.first) return false;
        return a.second < b.second;
    }

};

int main(int argc, char** argv)
{

    if(argc<2) {
        std::cerr << "Need input file" << std::endl;
        std::exit(1);
    }

    std::list<Point> lp;
    Point p;

    //read input
    std::ifstream is(argv[1]);
    std::string next_line;
    while( getline( is, next_line ) ) {
        if( next_line != "" && next_line[ 0 ] != '#' ) {
            //std::cerr << next_line;
            std::stringstream sstr(next_line);
            sstr >> p;
            lp.push_back(p);
        }
    }

    std::cerr << "Compute Delaunay triangulation..." << std::endl;
    Triangulation_hierarchy dt(lp.begin(),lp.end());
//    for (std::list<Point>::iterator it = lp.begin(); it != lp.end(); ++it) {
//        dt.insert(*it);
//    }/Users/uli/Downloads/neptune-raw.off/782_neptune-raw.pts.txt

    std::cerr << "Compute circumradii..." << std::endl;

    typedef CGAL::Triple<FT,int,CGAL::Object> Triple;

    std::vector<Triple > circumradii;

    for (Triangulation_hierarchy::Finite_vertices_iterator vertex = dt.finite_vertices_begin();
	vertex != dt.finite_vertices_end(); vertex++) {
        Vertex_handle vh(vertex);
        circumradii.push_back(CGAL::make_triple(FT(0),0,CGAL::make_object(vh)));
    }
    for (Triangulation_hierarchy::Finite_edges_iterator edge = dt.finite_edges_begin();
	edge != dt.finite_edges_end(); edge++) {
        Vertex_handle v1 = edge->first->vertex(edge->second);
        Vertex_handle v2 = edge->first->vertex(edge->third);
        circumradii.push_back(CGAL::make_triple(CGAL::squared_radius(v1->point(),v2->point()),1,CGAL::make_object(*edge)));
    }
    for (Triangulation_hierarchy::Finite_facets_iterator f = dt.finite_facets_begin();
	f != dt.finite_facets_end(); f++) {
        Vertex_handle v1 = f->first->vertex((f->second+1)%4);
        Vertex_handle v2 = f->first->vertex((f->second+2)%4);
        Vertex_handle v3 = f->first->vertex((f->second+3)%4);
        circumradii.push_back(CGAL::make_triple(CGAL::squared_radius(v1->point(),v2->point(),v3->point()),2,CGAL::make_object(*f)));
    }
    for (Triangulation_hierarchy::Finite_cells_iterator cit = dt.finite_cells_begin(); 
	cit != dt.finite_cells_end(); cit++) {
        Cell_handle ch(cit);
        Vertex_handle v1 = cit->vertex(0);
        Vertex_handle v2 = cit->vertex(1);
        Vertex_handle v3 = cit->vertex(2);
        Vertex_handle v4 = cit->vertex(3);
        circumradii.push_back(CGAL::make_triple(CGAL::squared_radius(v1->point(),v2->point(),v3->point(),v4->point()),3,CGAL::make_object(ch)));
    }
    
    std::cerr << "Sort circumradii..." << std::endl;

    std::sort(circumradii.begin(),circumradii.end(),Sort_triples<Triple>());

    std::cerr << "Filtration of size " << circumradii.size() << std::endl;

    phat::boundary_matrix< phat::vector_vector > boundary_matrix;

    Vertex_handle v;
    Edge e;
    Facet f;
    Cell_handle c;

    std::size_t filtration_index = 0;
    std::size_t filtration_size = circumradii.size();

    boundary_matrix.set_num_cols( filtration_size );

    Index curr_index = 0;

    for(std::vector<Triple>::const_iterator it = circumradii.begin();
        it != circumradii.end(); it++) {

        if(filtration_index % 100000 == 0) {
            std::cerr << filtration_index << " of " << filtration_size
                      << std::endl;
        }
        filtration_index++;

        const CGAL::Object& obj = it->third;

        Column col;

        if(CGAL::assign(v,obj)) {
            //std::cout << "Vertex " << it->second << std::endl;
            boundary_matrix.set_dim(curr_index, 0);
            v->info().set_index(curr_index);
            boundary_matrix.set_col(curr_index, col);
        }
        if(CGAL::assign(e,obj)) {
            //std::cout << "Edge " << it->second << std::endl;
            boundary_matrix.set_dim(curr_index, 1);
            Vertex_handle v1 = e.first->vertex(e.second);
            CGAL_assertion(v1->info().has_index());
            Vertex_handle v2 = e.first->vertex(e.third);
            CGAL_assertion(v2->info().has_index());
            Index i1 = v1->info().index();
            Index i2 = v2->info().index();
            if(i1>i2) {
                std::swap(v1,v2);
                std::swap(i1,i2);
            }
            col.push_back(i1);
            col.push_back(i2);
            boundary_matrix.set_col(curr_index, col);
            set_index_of_edge(dt,e,curr_index);
        }
        if(CGAL::assign(f,obj)) {
            //std::cout << "Facet " << it->second << std::endl;
            boundary_matrix.set_dim(curr_index, 2);

            Index i1= f.first->info().edge_index( (f.second+1)%4, (f.second+2)%4 );
            col.push_back(i1);
            Index i2= f.first->info().edge_index( (f.second+1)%4, (f.second+3)%4 );
            col.push_back(i2);
            Index i3= f.first->info().edge_index( (f.second+2)%4, (f.second+3)%4 );
            col.push_back(i3);
            std::sort(col.begin(),col.end());
            boundary_matrix.set_col(curr_index, col);
            set_index_of_facet(dt,f,curr_index);

        }
        if(CGAL::assign(c,obj)) {
            //std::cout << "Cell " << it->second << std::endl;
            boundary_matrix.set_dim(curr_index, 3);
            col.push_back(c->info().facet_index(0));
            col.push_back(c->info().facet_index(1));
            col.push_back(c->info().facet_index(2));
            col.push_back(c->info().facet_index(3));
            std::sort(col.begin(),col.end());
            boundary_matrix.set_col(curr_index, col);
        }
        curr_index++;
    }

    boundary_matrix.save_binary("alpha_filtration.bin");
    
    //phat::write(std::cout,M,dim_container);

    return 0;
}
