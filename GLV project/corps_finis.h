#ifndef FINI
#define FINI

# include <gmp.h>

void add_fp(mpz_t res, mpz_t add1, mpz_t add2, mpz_t corps);

void sub_fp(mpz_t res, mpz_t val1, mpz_t val2, mpz_t corps);

void inv_fp(mpz_t res, mpz_t p, mpz_t corps);

void mul_fp(mpz_t res, mpz_t p ,mpz_t q, mpz_t corps);


#endif // FINI
