#ifndef ___JACKHEADER___
#define ___JACKHEADER___

#include "qspray.h"
#include "ratioOfQsprays.h"

using namespace QSPRAY;
using namespace RATIOOFQSPRAYS;

typedef std::vector<int>                    Partition;
typedef std::vector<signed int>             Powers;
typedef boost::multiprecision::mpq_rational gmpq;
typedef boost::multiprecision::mpz_int      gmpi;

class vecHasher {
public:
  size_t operator()(const Powers& exponents) const {
    std::size_t seed = 0;
    for(auto& i : exponents) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

template <typename CoeffT>
using Poly = std::unordered_map<Powers, CoeffT, vecHasher>;
typedef Poly<int>  Zpoly;
typedef Poly<gmpq> Qpoly;


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

template <typename numT>
numT _betaratio(Partition, Partition, int, numT);
template <typename T>
RatioOfQsprays<T> _betaratio(Partition, Partition, int);

int weight(Partition);

template <typename numT>
numT ipow(numT, unsigned);

//std::string q2str(gmpq);


void simplifyPowers(Powers&);

template <typename CoeffT>
Poly<CoeffT> polyAdd(Poly<CoeffT>, Poly<CoeffT>);

Powers growPowers(Powers, int, int);

template <typename CoeffT>
Poly<CoeffT> polyMult(const Poly<CoeffT>, const Poly<CoeffT>);

template <typename CoeffT>
Poly<CoeffT> zeroPoly();

template <typename CoeffT>
Poly<CoeffT> unitPoly();

template <typename CoeffT>
Poly<CoeffT> constantPoly(CoeffT);

template <typename CoeffT>
Poly<CoeffT> lonePoly(const int);

template <typename CoeffT>
Poly<CoeffT> polyPow(const Poly<CoeffT>, int);



#endif
