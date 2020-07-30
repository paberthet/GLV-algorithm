#include "operation.h"
#include "projet.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

void double_and_add(point* res, point* p , mpz_t k)
{
    int end = mpz_popcount(k);     // donne le nombre de bit à 1 dans la représentation de k
    int i = 0;                     // place du bit considéré dans k
	point temp;                  // le point qui va être doublé successivement
	init_point(&temp,(*p).corps, (*p).a);
	set_point(&temp,p);
	set_point_str(res,"0","1","0");

	while(end != 0)
    {
        if(mpz_tstbit(k,i) == 1)           // test si le ieme bit de k vaut 1
        {
            add_point(res,res,&temp);
            end -= 1;                      // on a trouvé un bit à 1 donc on décrémente end
        }
        add_point(&temp,&temp,&temp);
        i += 1;
    }

	clear_point(&temp);
}


void ext_eucl(mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2, mpz_t m, mpz_t n)    // algorithme d'euclide étendu adapté au problème
{
    mpz_t u1,u2,v1,v2,q,r1,r2,end,norm1,norm2,temp;
    mpz_inits(u1,u2,v1,v2,q,r1,r2,end,norm1,norm2,temp,NULL);
    mpz_set_si(u1,1);
    mpz_set_si(u2,0);
    mpz_set_si(v1,0);
    mpz_set_si(v2,1);
    mpz_set(r1,m);
    mpz_set(r2,n);
    mpz_sqrt(end,m);

    while(mpz_cmp(r2,end) >= 0)
    {
        step_euc(u1,u2,v1,v2,q,r1,r2);
    }
    mpz_set(x1,r1);
    mpz_set(x2,r2);
    mpz_neg(y1,v1);
    mpz_neg(y2,v2);
    // choix du vecteur le plus court entre (x1,y1) et le prochain (r2,v2)
    mpz_mul(temp,x1,x1);
    mpz_mul(norm1,y1,y1);
    mpz_add(norm1,norm1,temp);   // norm1 vaut le carré de la norme du vecteur (x1,y1) = (r1,-v1)
    step_euc(u1,u2,v1,v2,q,r1,r2);
    mpz_mul(temp,r2,r2);
    mpz_mul(norm2,v2,v2);
    mpz_add(norm2,norm2,temp);  // norm2 vaut le carré de la norme du nouveau vecteur (r2,v2)
    if(mpz_cmp(norm1,norm2)>0)
    {
        mpz_set(x1,r2);
        mpz_neg(y1,v2);
    }
// les vecteurs (x1,y1) et (x2,y2) sont ceux du lemme 1
    mpz_clears(u1,u2,v1,v2,q,r1,r2,end,norm1,norm2,temp,NULL);
}

// un pas de l'algorithme d'euclide en dessous
void step_euc(mpz_t u1, mpz_t u2, mpz_t v1, mpz_t v2, mpz_t q, mpz_t r1, mpz_t r2)  // 2 désigne l'étape précédente, 1 l'étape d'encore avant. q=quotient, r = reste
{
    mpz_t temp;
    mpz_init(temp);

    mpz_fdiv_qr(q,temp,r1,r2);              // On stocke dans q et temp les valeurs des nouveaux quotient et reste
    mpz_set(r1,r2);                        // On stocke dans r1 le reste précédent
    mpz_set(r2,temp);                      // On se retrouve de nouveau avec r2 le reste le plus récent
    mpz_mul(temp,q,u2);
    mpz_sub(temp,u1,temp);
    mpz_set(u1,u2);             // u1 vaut le u précédent
    mpz_set(u2,temp);          // u2 est la nouvelle valeur de u
    mpz_mul(temp,q,v2);
    mpz_sub(temp,v1,temp);
    mpz_set(v1,v2);
    mpz_set(v2,temp);     // même chose que pour u mais avec v

    mpz_clear(temp);
}

void close_vect(mpz_t resx, mpz_t resy, mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2, mpz_t k)   // en entrée, les vecteurs (x1,y1) et (x2,y2) du lemme1 et le k considéré
{
    mpz_t b1,b2,reste1,reste2;
    mpz_inits(b1,b2,reste1,reste2,NULL);

    mpz_mul(resy,x2,y1);
    mpz_submul(resy,x1,y2);       // resy vaut désormais x2*y1 - x1*y2
    mpz_mul(b1,k,y2);
    mpz_neg(b1,b1);
    mpz_mul(b2,k,y1);

    // on va arrondir à l'entier le plus proche b1/resy et b2/resy et les stocker dans b1 et b2
    mpz_cdiv_qr(resx,reste1,b1,resy);   // l'arrondi au-dessus dans resx
    mpz_abs(reste1,reste1);
    mpz_fdiv_qr(b1,reste2,b1,resy);    //l'arrondi en-dessous
    mpz_abs(reste2,reste2);
    if(mpz_cmp(reste1,reste2)<0)
    {
        mpz_set(b1,resx);
    }
    mpz_cdiv_qr(resx,reste1,b2,resy);   // l'arrondi au-dessus dans resx
    mpz_abs(reste1,reste1);
    mpz_fdiv_qr(b2,reste2,b2,resy);    //l'arrondi en-dessous
    mpz_abs(reste2,reste2);
    if(mpz_cmp(reste1,reste2)<0)
    {
        mpz_set(b2,resx);
    }
    // on construit le vecteur v du lemme2 du projet dans (resx,resy)
    mpz_mul(resx,b1,x1);
    mpz_addmul(resx,b2,x2);
    mpz_mul(resy,b1,y1);
    mpz_addmul(resy,b2,y2);
    // On construit le vecteur u du lemme 2 dans (resx,resy)
    mpz_sub(resx,k,resx);
    mpz_neg(resy,resy);

    mpz_clears(b1,b2,reste1,reste2,NULL);
}

void decomp_k(mpz_t k1, mpz_t k2, mpz_t k, mpz_t n, mpz_t lambda)   // décompose k sous la forme k1+lambda.k2[n]
{
    mpz_t x1,x2,y1,y2;
    mpz_inits(x1,x2,y1,y2,NULL);

    ext_eucl(x1,x2,y1,y2,n,lambda);  /* cette étape n'a pas besoin d'être recalculée à chaque nouvelle multiplication, on utilisera donc plutôt
    ext_eucl et close_vect séparément que via decomp_k */
    close_vect(k1,k2,x1,x2,y1,y2,k);

    mpz_clears(x1,x2,y1,y2,NULL);
}



void simult_point_mul(point* r, mpz_t u, mpz_t v, point* p, point* q)    // Calcule u.p + v.q dans r, u et v étant positifs. On a forcé w = 1 par rapport à l'algo du projet
{
    int i,j,ui,vi;
    point s,temp;
    init_point(&s,(*p).corps,(*p).a);
    init_point(&temp,(*p).corps,(*p).a);

    add_point(&s,p,q);  // Le précalcule
    // à partir d'ici, un double_and_add sur les 2 termes en même temps
    i = mpz_sizeinbase(u,2);
    j = mpz_sizeinbase(v,2);
    if(j<i)
    {
        j = i;
    }
    for(i = j-1; i >= 0; i--)
    {
        ui = mpz_tstbit(u,i);
        vi = mpz_tstbit(v,i);
        add_point(&temp,&temp,&temp);
        if((ui == 0) && (vi == 1))
        {
            add_point(&temp,&temp,q);
        }
        else if((ui == 1) && (vi == 0))
        {
            add_point(&temp,&temp,p);
        }
        else if((ui == 1) && (vi == 1))
        {
            add_point(&temp,&temp,&s);
        }
    }
    set_point(r,&temp);
    clear_point(&s);
    clear_point(&temp);
}
