# include "operation.h"
# include "corps_finis.h"
# include <gmp.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

void rand_prime(mpz_t res, int n)   // res devient un nombre premier entre 2^(n-1) et 2^(n+1)
{
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    gmp_randseed_ui(state,time.tv_nsec);

    mpz_rrandomb(res,state,n);
    mpz_nextprime(res,res);

    gmp_randclear(state);
}


void rand_val(mpz_t res, mpz_t corps, int n)    // res devient  un nombre aléatoire inférieur à 2^n
{
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    gmp_randseed_ui(state,time.tv_nsec);

    mpz_urandomb(res,state,n);
    mpz_mod(res,res,corps);

    gmp_randclear(state);
}

void phi1(point* res, point* p)  /* Phi1 correspond à la mult par un certain lambda (adapté à la courbe1) */
{
    mpz_t beta;
    mpz_init(beta);

    mpz_set_str(beta,"771473166210819779552257112796337671037538143582",10);
    mul_fp((*res).x,(*p).x,beta,(*p).corps);
    mpz_set((*res).y,(*p).y);
    (*res).z = 1;

    mpz_clear(beta);
}


void phi2(point* res, point* p)  /* Phi2 correspond à la mult par un certain lambda (adapté à la courbe2) */
{
    mpz_t beta;
    mpz_init(beta);

    mpz_set_str(beta,"100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",2);
    mul_fp((*res).x,(*p).x,beta,(*p).corps);
    mpz_set((*res).y,(*p).y);
    (*res).z = 1;

    mpz_clear(beta);
}
