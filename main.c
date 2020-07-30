#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "operation.h"
#include "corps_finis.h"
#include "projet.h"
#include "application.h"


int main(int argc, char** argv)
{
// VERIFICATION DE L'OUTPUT DE LA METHODE

/*

    int i;
    mpz_t zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2;
    mpz_inits(zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2,NULL);
    mpz_set_ui(zero,0);
    mpz_set_str(corps,"1461501637330902918203684832716283019655932313743",10);
    mpz_set_str(n,"1461501637330902918203687013445034429194588307251",10);
    mpz_set_str(lambda,"903860042511079968555273866340564498116022318806",10);
    point p,q,p2;
    init_point(&p,corps,zero);
    init_point(&q,corps,zero);
    init_point(&p2,corps,zero);

    set_point_str(&p,"1","2","1");
    set_point(&p2,&p);

    ext_eucl(x1,x2,y1,y2,n,lambda);   // ne dépend pas de k donc à calculer une seule fois

    for(i=0;i<10;i++)
    {
        rand_val(k,corps,160);
        printf("tour %d\n",i);

        double_and_add(&q,&p,k);
        print_point(&q);


        close_vect(k1,k2,x1,x2,y1,y2,k);
        phi1(&q,&p);
        if(mpz_sgn(k1)<0)
        {
            mpz_neg(k1,k1);
            inv_point(&p2,&p);
        }
        else
        {
            set_point(&p2,&p);
        }
        if(mpz_sgn(k2)<0)
        {
            mpz_neg(k2,k2);
            inv_point(&q,&q);
        }
        simult_point_mul(&q,k1,k2,&p2,&q);

        print_point(&q);
    }

    mpz_clears(zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2,NULL);
    clear_point(&p);
    clear_point(&q);
    clear_point(&p2);


/* COURBE 1, Y^2 = X^3 + 3 */


/*
    FILE* fp;
    fp = fopen("courbe1_1000.txt","a+");
    int i,j;
    float temps1,temps2;
    clock_t debut,fin;
    mpz_t zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2;
    mpz_inits(zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2,NULL);
    mpz_set_ui(zero,0);
    mpz_set_str(corps,"1461501637330902918203684832716283019655932313743",10);
    mpz_set_str(n,"1461501637330902918203687013445034429194588307251",10);
    mpz_set_str(lambda,"903860042511079968555273866340564498116022318806",10);
    point p,q,p2;
    init_point(&p,corps,zero);
    init_point(&q,corps,zero);
    init_point(&p2,corps,zero);

    set_point_str(&p,"1","2","1");
    set_point(&p2,&p);

    ext_eucl(x1,x2,y1,y2,n,lambda);   // ne dépend pas de k donc à calculer une seule fois

    for(i=0;i<1000;i++)
    {
        rand_val(k,corps,160);
        debut = clock();

        for(j=0;j<1000;j++)
        {
            double_and_add(&q,&p,k);
        }

        fin = clock();
        temps1 = (float)(fin - debut)/CLOCKS_PER_SEC;

        debut = clock();

        for(j=0;j<1000;j++)
        {
            close_vect(k1,k2,x1,x2,y1,y2,k);
            phi1(&q,&p);
            if(mpz_sgn(k1)<0)
            {
                mpz_neg(k1,k1);
                inv_point(&p2,&p);
            }
            else
            {
                set_point(&p2,&p);
            }
            if(mpz_sgn(k2)<0)
            {
                mpz_neg(k2,k2);
                inv_point(&q,&q);
            }
            simult_point_mul(&q,k1,k2,&p2,&q);
        }

        fin = clock();
        temps2 = (float)(fin - debut)/CLOCKS_PER_SEC;

        gmp_fprintf(fp,"%f\n",temps2/temps1);
    }

    mpz_clears(zero,corps,n,k,k1,k2,lambda,x1,x2,y1,y2,NULL);
    clear_point(&p);
    clear_point(&q);
    clear_point(&p2);
    fclose(fp);

    /* COURBE 2, Y^2 = X^3 - 2 */

/*
    FILE* fp;
    fp = fopen("courbe2_1000.txt","a+");
    int i,j;
    float temps1,temps2;
    clock_t debut,fin,deb,end;
    mpz_t zero,corps,n,k,k1,k2,lambda,st,x1,x2,y1,y2;
    mpz_inits(zero,corps,n,k,k1,k2,lambda,st,x1,x2,y1,y2,NULL);
    mpz_set_ui(st,63);
    mpz_set_ui(zero,0);
    mpz_set_str(corps,"1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011",2);
    mpz_set_str(n,"111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111",2);
    mpz_divexact_ui(n,n,63);
    mpz_set_str(lambda,"111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110",2);
    mpz_divexact_ui(lambda,lambda,3);
    point p,q,p2;
    init_point(&p,corps,zero);
    init_point(&q,corps,zero);
    init_point(&p2,corps,zero);

    deb = clock();

    set_point_str(&p,"3","5","1");
    double_and_add(&p,&p,st);  // p est désormais un point d'ordre n premier
    set_point(&p2,&p);

    ext_eucl(x1,x2,y1,y2,n,lambda); // ne dépend pas de k donc à calculer une seule fois

    for(i=0;i<1000;i++)
    {
        rand_val(k,corps,389);
        debut = clock();

        for(j=0;j<100;j++)
        {
            double_and_add(&q,&p,k);
        }

        fin = clock();
        temps1 = (float)(fin - debut)/CLOCKS_PER_SEC;

        debut = clock();

        for(j=0;j<100;j++)
        {
            close_vect(k1,k2,x1,x2,y1,y2,k);
            phi2(&q,&p);
            if(mpz_sgn(k1)<0)
            {
                mpz_neg(k1,k1);
                inv_point(&p2,&p);
            }
            else
            {
                set_point(&p2,&p);
            }
            if(mpz_sgn(k2)<0)
            {
                mpz_neg(k2,k2);
                inv_point(&q,&q);
            }
            simult_point_mul(&q,k1,k2,&p2,&q);
        }

        fin = clock();
        temps2 = (float)(fin - debut)/CLOCKS_PER_SEC;

        gmp_fprintf(fp,"%f\n",temps2/temps1);
    }


    // Libération de mémoire
    mpz_clears(zero,corps,n,k,k1,k2,lambda,st,x1,x2,y1,y2,NULL);
    clear_point(&p);
    clear_point(&q);
    clear_point(&p2);
    fclose(fp);
    end = clock();
    temps1 = (float)(end - deb)/CLOCKS_PER_SEC;
    printf("temps total = %f\n",temps1);


/* MOYENNE DES RATIOS DES DONNEES */

/*

    FILE* fp;
    fp = fopen("courbe2_1000.txt","r");

    float a = 0;
    char str[10];

    while(fgets(str,10,fp) != NULL)
    {
        a += atof(str);
    }
    a /= 1000;
    printf("%f\n",a);

    fclose(fp); // moyenne courbe 1: 0.594751, courbe 2: 0.571638, courbe 1 2eme round: 0,604026, courbe 2 2eme round: 0.565441  */


    return 0;
}
