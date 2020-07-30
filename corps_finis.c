# include "corps_finis.h"
# include <gmp.h>
# include <stdio.h>

void add_fp(mpz_t res, mpz_t p, mpz_t q, mpz_t corps)   // res prend la valeur p+q [corps]
{
    mpz_add(res, p, q);
    mpz_mod(res, res, corps);
}

void sub_fp(mpz_t res, mpz_t val1, mpz_t val2, mpz_t corps)   // res prend la valeur p-q [corps]
{
    mpz_sub(res, val1, val2);
    mpz_mod(res, res, corps);
}

void inv_fp(mpz_t res, mpz_t p, mpz_t corps) // res prend la valeur p^-1 dans Fcorps
{
    mpz_t s,g;
    mpz_inits(s,g,NULL);
    mpz_gcdext(g,res,s,p,corps);
    mpz_mod(res, res, corps);
    mpz_clears(s,g,NULL);
}

void mul_fp(mpz_t res, mpz_t p ,mpz_t q, mpz_t corps)  // res prend la valeur p*q [corps]
{
    mpz_mul(res, p, q);
    mpz_mod(res, res, corps);
}
