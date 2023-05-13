#include "jack.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename numT>
numT ipow(numT base, unsigned exp){
  numT result(1);
  while(exp) {
    if(exp & 1) {
      result *= base;
    }
    exp >>= 1;
    base *= base;
  }
  return result;
}

template gmpq ipow<gmpq>(gmpq, unsigned);
template double ipow<double>(double, unsigned);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
std::string q2str(gmpq r) {
  const gmpi numer = boost::multiprecision::numerator(r);
  const gmpi denom = boost::multiprecision::denominator(r);
  mpz_t p; mpz_init(p);
  mpz_set(p, numer.backend().data());
  mpz_t q; mpz_init(q);
  mpz_set(q, denom.backend().data());
  const size_t n = mpz_sizeinbase(p, 10) + 2;
  const size_t d = mpz_sizeinbase(q, 10) + 2;
  char* cnumer = new char[n];
  char* cdenom = new char[d];
  cnumer = mpz_get_str(cnumer, 10, p);
  cdenom = mpz_get_str(cdenom, 10, q);
  std::string snumer = cnumer;
  std::string sdenom = cdenom;
  delete[] cnumer; delete[] cdenom;
  mpz_clear(p); mpz_clear(q);
  return snumer + "/" + sdenom;
}
