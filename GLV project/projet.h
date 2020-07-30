#ifndef PROJET
#define PROJET

#include "operation.h"
#include <gmp.h>

void double_and_add(point* res, point* p , mpz_t k);
void ext_eucl(mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2, mpz_t m, mpz_t n);
void step_euc(mpz_t u1, mpz_t u2, mpz_t v1, mpz_t v2, mpz_t q, mpz_t r1, mpz_t r2);
void close_vect(mpz_t resx, mpz_t resy, mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2, mpz_t k);
void decomp_k(mpz_t k1, mpz_t k2, mpz_t k, mpz_t n, mpz_t lambda);
void simult_point_mul(point* r, mpz_t u, mpz_t v, point* p, point* q);


#endif // FINI
