// On se place en carac différente de 2 sur une courbe d'équation y^2 = x^3 +a*x + b

# include "operation.h"
# include "corps_finis.h"
# include <gmp.h>
# include <stdio.h>
# include <stdlib.h>

void init_point(point* p, mpz_t corps, mpz_t a)         // initialise les arguments mpz_t et définit le corps et la courbe considérés
{
    mpz_inits((*p).x,(*p).y,(*p).corps,(*p).a, NULL);
    (*p).z = 0;
    mpz_set_si((*p).x,0);
    mpz_set_si((*p).y,1);
    mpz_set((*p).corps, corps);    // défini le corps sur lequel p est défini une fois pour toute
    mpz_set((*p).a,a);             // détermine la courbe où se trouve le point
}


void clear_point(point* p)
{
    mpz_clears((*p).x, (*p).y, (*p).corps, NULL);
}

void set_point_str(point* p, char* x, char* y, char* z)
{
    mpz_set_str((*p).x,x,10);    // première coordonnée du point : x
    mpz_set_str((*p).y,y,10);    // deuxième coordonnée du point : y
    (*p).z = atoi(z);       // 0 si p est le point à l'infini, 1 sinon
}


void set_point_mpz(point* p, mpz_t x, mpz_t y, int z)
{
    mpz_set((*p).x,x);
    mpz_set((*p).y,y);
    (*p).z = z;
}


void set_point(point* p, point* q) // copie q dans p
{
    mpz_set((*p).x,(*q).x);
    mpz_set((*p).y,(*q).y);
    (*p).z = (*q).z;
    mpz_set((*p).corps,(*q).corps);
    mpz_set((*p).a,(*q).a);
}

void print_point(point* p)      // affiche les coordonnées de p
{
    if ((*p).z == 0)
    {
        printf("infini\n");
    }
    else
    {
        gmp_printf("(p.x, p.y) = (%Zd, %Zd)\n", (*p).x, (*p).y);
    }
}


void fprint_point(FILE* fp, point*p)  // ecrit les coordonnées de p  dans le fichier fp
{
    if ((*p).z == 0)
    {
        fprintf(fp,"infini\n");
    }
    else
    {
        gmp_fprintf(fp,"(p.x, p.y) = (%Zd, %Zd)\n", (*p).x, (*p).y);
    }
}

int compare(point* p, point* q)   // renvoie 1 si p=q, 0 sinon
{
    if ((*p).z == (*q).z && (*p).z == 0)   // si p = q = 0
    {
        return(1);
    }
    else if ((*p).z != (*q).z)   // si l'un seulement vaut 0
    {
        return(0);
    }
    else if (mpz_cmp((*p).x,(*q).x) == 0 && mpz_cmp((*p).y,(*q).y) == 0)  // si p.x = q.x et p.y = q.y
    {
        return(1);
    }
    else
    {
        return(0);
    }
}


int compare_inv(point* p, point* q)   // renvoie 1 si p = -q, 0 sinon
{
    if ((*p).z == (*q).z && (*p).z == 0)   // si p et q valent O
    {
        return(1);
    }
    else if ((*p).z != (*q).z)    // si seulement l'un vaut 0
    {
        return(0);
    }
    else if (mpz_cmp((*p).x,(*q).x) == 0 && mpz_sgn((*p).y) == 0)   // si p = q sont des points spéciaux
    {
        return(1);
    }
    else if (mpz_cmp((*p).x,(*q).x) == 0 && mpz_cmp((*p).y,(*q).y) != 0)    // si p.x = q.x, p.y != q.y
    {
        return(1);
    }
    else
    {
        return(0);
    }
}



void inv_point(point* res, point* p)  // donne l'inverse de p dans res
{
    mpz_set((*res).x, (*p).x);
    mpz_neg((*res).y, (*p).y);
    mpz_mod((*res).y, (*res).y, (*p).corps);
    (*res).z = (*p).z;
}


void add_point(point* res, point* p, point* q)   // stocke p+q dans res (res peut être p ou q)
{
    mpz_t lambda,beta,resx,resy;
    mpz_inits(lambda,beta,resx,resy,NULL);

    if ((*p).z == 0)   // si p = 0
    {
        set_point(res,q);
    }
    else if ((*q).z == 0)   // si q = 0
    {
        set_point(res, p);
    }
    else if (compare_inv(p,q) == 1)   // si p = -q
    {
        (*res).z = 0;
        mpz_set_si((*res).x,0);
        mpz_set_si((*res).y,1);
    }
    else
    {
        if (compare(p,q) == 0) // si p != q, la pente est lambda = (y1-y2)/(x1-x2)
        {
            sub_fp(lambda,(*p).y,(*q).y,(*p).corps);
            sub_fp(beta,(*p).x,(*q).x,(*p).corps);
            inv_fp(beta,beta,(*p).corps);
            mul_fp(lambda,lambda,beta,(*p).corps);
        }
        else   // p = q points non spéciaux, la pente est lambda = (3x^2+a)/2y
        {
            mul_fp(lambda,(*p).x,(*p).x,(*p).corps);
            mpz_mul_si(lambda, lambda, 3);
            mpz_mod(lambda, lambda, (*p).corps);
            mpz_add(lambda, lambda, (*p).a);
            mpz_mod(lambda, lambda, (*p).corps);
            mpz_mul_si(beta, (*p).y, 2);
            mpz_mod(beta, beta, (*p).corps);
            inv_fp(beta,beta,(*p).corps);
            mul_fp(lambda, lambda, beta, (*p).corps);
        }
        mul_fp(resx, lambda, lambda, (*p).corps);
        sub_fp(resx, resx, (*p).x, (*p).corps);
        sub_fp(resx, resx, (*q).x, (*p).corps);
        sub_fp(resy, (*p).x, resx, (*p).corps);
        mul_fp(resy, lambda, resy, (*p).corps);
        sub_fp(resy, resy, (*p).y, (*p).corps);
        mpz_set((*res).x, resx);
        mpz_set((*res).y, resy);
        (*res).z = 1;
    }

    mpz_clears(lambda,beta,resx,resy,NULL);
}
