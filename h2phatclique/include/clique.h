#pragma once
#ifndef CLIQUE_H
	#define CLIQUE_H
#endif
#include <vector>
#include <cstdlib>
#include <iostream>

namespace cliques
{

// Type definitions
typedef int32_t index; // int64 to handle over a billion cliques
                        // change int32 to save memory
typedef uint8_t sz_t;   // size of column = (dim + 1)
typedef std::vector<index> boundary;
typedef float weight_t; // change to double if needed - float saves memory

class clique
{
private:
    index cliquenum;
    weight_t weight;  
public:
	boundary bdry;

    // Constructors and operators
public:
    clique(): bdry(std::vector<index>{}),cliquenum(0),weight(0) {};

    // Construct clique without weight
    clique( boundary& v, index cnum ):bdry(v),cliquenum(cnum),weight(0) {};

    // Construct clique without weight from array
    clique( index *v, index cnum ):bdry(v,v+sizeof(v)/sizeof(v[0])),cliquenum(cnum),weight(0) {};

    // Construct clique with weight
    clique( boundary& v, index cnum, weight_t wt ): bdry(v),cliquenum(cnum),weight(wt) {};

    // Construct clique with weight from brace-enclosed list
    clique( std::initializer_list<index> l, index cnum, weight_t wt ): bdry(l),cliquenum(cnum),weight(wt) {};

    // Construct clique without weight from brace-enclosed list
    clique( std::initializer_list<index> l, index cnum ): bdry(l),cliquenum(cnum),weight(0) {};

    // Destructor
    ~clique() {}

    // Copy Constructor
    clique( const clique& other ): bdry(other.bdry),cliquenum(other.cliquenum),weight(other.weight) {}

    // Copy assignment operator
    clique& operator=(const clique& other)
    {
        if ( *this != other )
        {
            bdry=other.bdry;
            cliquenum=other.cliquenum;
            weight=other.weight;
        }

        return *this;
    }

    // Move constructor
    clique(clique&& other): bdry(std::vector<index>{}),cliquenum(0),weight(0)
    {
        bdry=other.bdry;
        cliquenum=other.cliquenum;
        weight=other.weight;
    }

    // Move assignment operator
    clique& operator=(clique&& other)
    {
        if ( *this != other )
        {
            bdry={};
            cliquenum=0;
            weight=0;

            bdry=other.bdry;
            cliquenum=other.cliquenum;
            weight=other.weight;

            other.bdry={};
            other.cliquenum=0;
            other.weight=0;
        }

        return *this;
    }

    // Two cliques are equal if their boundaries are equal.
    bool operator==( const clique& other_clique ) const
	{ return ( bdry == other_clique.bdry ); }

    bool operator!=( const clique& other_clique ) const
	{ return !( *this == other_clique ); }

    // For sorting purposes, cliques are ordered by weight.
    bool operator<(const clique& other_clique) const
	{ return ( weight < other_clique.weight ); }

    bool operator>(const clique& other_clique) const
	{ return ( other_clique < *this ); }

    bool operator<=(const clique& other_clique) const
	{ return !( *this > other_clique ); }

    bool operator>=(const clique& other_clique) const
	{ return !( *this < other_clique ); }

    // Functions to get and set clique attributes
public:
    sz_t get_size() const
    { return bdry.size(); }

    void print_boundary() const
    {
        for ( const index i: bdry )
		{ std::cout << i << ' '; }
    }

    const boundary& get_boundary() const
	{ return bdry; }

    void set_boundary(boundary b)
	{ bdry = b; }

    weight_t get_weight() const
	{ return weight; }

    void set_weight(weight_t w)
	{ weight = w; }

    index get_cliquenum() const
	{ return cliquenum; }

    void set_cliquenum(index i)
	{ cliquenum = i; }
};

class bHash
{
public:
    // This is from Boost, but replicated to avoid a dependence on Boost
    std::size_t operator()( const boundary& bdry) const 
    {
	std::size_t seed = bdry.size();
        for(auto& i : bdry) 
        	{ seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2); }
	return seed;
    }

};

struct CliqueComparerByWeight
{
    bool operator() ( const cliques::clique& lhs, const cliques::clique& rhs ) const
        { return (lhs.get_weight()) < (rhs.get_weight()); }
};

}
