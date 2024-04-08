#ifndef ___JACKHEADER___
#define ___JACKHEADER___

#include <RcppArmadillo.h>
#include <boost/multiprecision/gmp.hpp>
#include "qspray.h"
#include "ratioOfQsprays.h"

using namespace QSPRAY;
using namespace RATIOOFQSPRAYS;

typedef std::vector<int>                    Partition;

// template <typename CoeffT>
// using Poly = std::unordered_map<Powers, CoeffT, vecHasher>;
// typedef Poly<int>  Zpoly;
// typedef Poly<gmpq> Qpoly;

class pairHasher {
public:
  size_t operator()(const std::pair<int, int>& ij) const {
    size_t seed = 0;
    seed ^= ij.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= ij.second + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
int _N(Partition, Partition);

template <typename T>
std::pair<T,T> _betaratio(Partition, Partition, int, T);

int weight(Partition);


template <typename numT>
numT ipow(numT, unsigned);


#endif
