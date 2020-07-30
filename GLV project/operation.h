#ifndef OPERATION
#define OPERATION

#include <gmp.h>
#include <stdio.h>

typedef struct point point;
struct point
{
    mpz_t x;
    mpz_t y;
    int z;         // vaut 0 ou 1 pour simplifier
    mpz_t corps;
    mpz_t a;          // Le a de l'équation y^2 = x^3 + ax + b
};

void init_point(point* p, mpz_t corps, mpz_t a);

void set_point_str(point* p, char* x, char* y, char* z);

void set_point_mpz(point* p, mpz_t x, mpz_t y, int z);

void set_point(point* p, point* q);

void clear_point(point* p);

int compare(point* p, point* q);

int compare_inv(point* p, point* q);

void print_point(point* p);

void fprint_point(FILE* fp, point*p);

void inv_point(point* res, point* p);

void add_point(point* res, point* p, point* q);

#endif // OPERATION
