#include "jack.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
void simplifyPowers(Powers& pows) {
  int n = pows.size();
  if(n == 0) {
    return;
  }
  n--;
  Powers::iterator it = pows.end();
  bool zero = pows[n] == 0;
  while(zero && n >= 0) {
    it--;
    n--;
    zero = pows[n] == 0;
  }
  if(n == -1) {
    pows = {};
  } else {
    pows.erase(it, pows.end());
  }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyAdd(Poly<CoeffT> P1, Poly<CoeffT> P2) {
  Powers pows;
  // Poly<CoeffT> P1copy; // USELESS
  // for(auto it = P1.begin(); it != P1.end(); ++it) {
  //   pows = it->first;
  //   P1copy[pows] = P1[pows];
  // }
  for(auto it = P2.begin(); it != P2.end(); ++it) {
    pows = it->first;
    P1[pows] += P2[pows];
    if(P1[pows] == 0) {
      P1.erase(pows);
    }
  }
  return P1;
}

template Zpoly polyAdd<int>(Zpoly, Zpoly);
template Qpoly polyAdd<gmpq>(Qpoly, Qpoly);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
Powers growPowers(Powers pows, int m, int n) {
  Powers gpows;
  gpows.reserve(n);
  for(int i = 0; i < m; i++) {
    gpows.emplace_back(pows[i]);
  }
  for(int i = m; i < n; i++) {
    gpows.emplace_back(0);
  }
  return gpows;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyMult(const Poly<CoeffT> P1, const Poly<CoeffT> P2) {

  Poly<CoeffT> Pout;
  Powers powssum;
  int i;

  for(auto it1 = P1.begin(); it1 != P1.end(); ++it1) {
    const CoeffT c1 = it1->second;
    if(c1 != 0) {
      Powers pows1 = it1->first;
      int n1 = pows1.size();
      for(auto it2 = P2.begin(); it2 != P2.end(); ++it2) {
        const CoeffT c2 = it2->second;
        if(c2 != 0) {
          Powers pows2 = it2->first;
          int n2 = pows2.size();
          powssum.clear();
          if(n1 < n2) {
            Powers gpows = growPowers(pows1, n1, n2);
            powssum.reserve(n2);
            for(i = 0; i < n2; i++) {
              powssum.emplace_back(gpows[i] + pows2[i]);
            }
          } else if(n1 > n2) {
            Powers gpows = growPowers(pows2, n2, n1);
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + gpows[i]);
            }
          } else {
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + pows2[i]);
            }
          }
          Pout[powssum] += c1 * c2;
        }
      }
    }
  }

  return Pout;
}

template Zpoly polyMult<int>(const Zpoly, const Zpoly);
template Qpoly polyMult<gmpq>(const Qpoly, const Qpoly);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> zeroPoly() {
  Poly<CoeffT> out;
  return out;
}

template Zpoly zeroPoly<int>();
template Qpoly zeroPoly<gmpq>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> unitPoly() {
  Poly<CoeffT> out;
  Powers pows(0);
  CoeffT one(1);
  out[pows] = one;
  return out;
}

template Zpoly unitPoly<int>();
template Qpoly unitPoly<gmpq>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> constantPoly(CoeffT a) {
  Poly<CoeffT> out;
  Powers pows(0);
  out[pows] = a;
  return out;
}

template Qpoly constantPoly<gmpq>(gmpq);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> lonePoly(const int n) {
  Poly<CoeffT> out;
  Powers pows(n, 0);
  pows[n-1] = 1;
  CoeffT one(1);
  out[pows] = one;
  return out;
}

template Zpoly lonePoly<int>(const int);
template Qpoly lonePoly<gmpq>(const int);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyPow(const Poly<CoeffT> P, int n) {
  Poly<CoeffT> out;
  if(n >= 1) {
    if(n == 1) {
      out = P;
    } else {
      out = unitPoly<CoeffT>();
      for(; n > 0; n--) {
        out = polyMult(P, out);
      }
    }
  } else {
    out = unitPoly<CoeffT>();
  }
  return out;
}

template Zpoly polyPow<int>(const Zpoly, int);
template Qpoly polyPow<gmpq>(const Qpoly, int);
