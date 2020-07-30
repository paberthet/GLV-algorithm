#ifndef APPLI
#define APPLI

#include <gmp.h>
#include "operation.h"

void rand_prime(mpz_t res, int n);
void rand_val(mpz_t res, mpz_t corps, int n);
void phi1(point* res, point* p);
void phi2(point* res, point* p);

#endif // APPLI
