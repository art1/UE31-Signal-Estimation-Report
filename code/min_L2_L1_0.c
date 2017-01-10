/*=============================================================================
 Interface Matlab/C de la fonction gauss_seidel
 >> [ x, J, NZ ] = Min_L2_L1(y,W,lambda[,n_it_max,seuil]);
 renvoie le vecteur x minimisant le critere 
        J(x) = 1/2||y-W*x||_2^2 + lambda*||x||_1
 avec J le vecteur donnant la valeur du critere a chaque iteration
 et NZ le vecteur donnant le nombre d'elements non nuls a chaque iteration
 
 Version specialisee pour le probleme d'analyse spectrale pour lequel la matrice 
 W'W est Toeplitz symetrique et W(:,k)'*W(:,k) = N
 Version sans penalisation pour le teme x_k0 avec k0 indice du milieu de x 
 Correspondant a la frequence nulle en analyse spectrale
 ==============================================================================*/

#include "mex.h"
#include "math.h"

/*=============================================================================
 declaration des strucures
 ==============================================================================*/

/* Structure pour les Matrices complexes */
struct Mat
{
   int        L;   /* Nombre de Lignes */
   int        C;   /* Nombre de Colonnes */
   double**   R;   /* Pointeur sur la partie r�elle */
   double**   I;   /* Pointeur sur la partie imaginaire */
   } ;

/* Structure pour les vecteurs complexes dim C */
struct vect
{
   int       C;   /* Nombre d'elements */
   double*   R;   /* Pointeur sur la partie r�elle */
   double*   I;   /* Pointeur sur la partie imaginaire */
} ;

/* Structure pour les variables complexes  */
struct var_comp
{
   double   R;   /* partie r�elle */
   double   I;   /* partie imaginaire */
} ;


/*=============================================================================
 Fonction qui calcule le carr� du module d'un nombre complexe
 >> mx = abs(x)^2;
 ==============================================================================*/
double abs2_C(struct var_comp x)
{
    return(x.R*x.R + x.I*x.I);
}

/*=============================================================================
 Fonction qui calcule le module d'un nombre complexe
 >> mx = abs(x);
 ==============================================================================*/
double abs_C(struct var_comp x)
{
    return(sqrt(x.R*x.R + x.I*x.I));
}

/*=============================================================================
 Fonction qui calcule le module du ki�me �l�ment d'un vecteur complexe
 >> mx = abs(x(k));
 ==============================================================================*/
double abs_VectCk(struct vect x,int k)
{
    return(sqrt(x.R[k]*x.R[k] + x.I[k]*x.I[k]));
}

/*=============================================================================
 Fonction qui calcule le carr� du module du ki�me �l�ment d'un vecteur complexe
 >> mx = abs(x(k))^2;
 ==============================================================================*/
double abs2_VectCk(struct vect x,int k)
{
    return(x.R[k]*x.R[k] + x.I[k]*x.I[k]);
}

/*=============================================================================
 Fonction qui calcule l'argument d'un nombre complexe
 >> argx = angle(x);
 ==============================================================================*/
double arg_C(struct var_comp x)
{
    return(atan2(x.I,x.R));
}
    
/*=============================================================================
 Fonction qui calcule le produit de la ki�me colonne de la matrice A par le 
 ki�me �l�ment du vecteur b
 >> c = A(:,k)*x(k)
 ==============================================================================*/
void Prod_MatColk_Vectk(struct vect c, struct Mat A, struct vect b,int k)
{
 int n;
 
 /* c est un vecteur colonne de meme nb d'elements que le nb de lignes de A */
 /* b est un vecteur colonne de meme nb d'elements que le nb de colonnes de A */
 
 c.C = A.L; 

     if ((A.I!=NULL) && (b.I!=NULL))
     /* 1er cas : A et b reels */
        for (n=0;n<A.L;n++) {
        c.R[n] = A.R[k][n]*b.R[k] - A.I[k][n]*b.I[k];
        c.I[n] = A.R[k][n]*b.I[k] + A.I[k][n]*b.R[k];
        }
     else /* 2�me cas : A reel et b complexe */
     if ((A.I==NULL) && (b.I!=NULL))
        for (n=0;n<A.L;n++) {
         c.R[n]=A.R[k][n]*b.R[k];
         c.I[n]=A.R[k][n]*b.I[k];
        }
     else /* 3�me cas : si A complexes et b reel */
     if ((A.I!=NULL) && (b.I==NULL))
        for (n=0;n<A.L;n++) {
         c.R[n]=A.R[k][n]*b.R[k];
         c.I[n]=A.I[k][n]*b.R[k];
        }
     else /* 4�me cas : si A et b reels */
     if ((A.I==NULL) && (b.I==NULL))
        for (n=0;n<A.L;n++) {
         c.R[n]=A.R[k][n]*b.R[k];
         c.I[n]=0;
        }
 
}

/*=============================================================================
 fonction qui calcule le produit de la ki�me colonne de la matrice A transpos�e  
 conjugu�e par le vecteur colonne b
 >> c = A(k,:)'*b;
 ==============================================================================*/
void Prod_MatColkConj_Vect(struct var_comp *c, struct Mat A, struct vect b,int k)
{
 int n;
  
 /* c est un scalaire */
 /* b est un vecteur colonne de meme nb d'elements que le nb de lignes de A */
 
     (*c).R = 0;
     (*c).I = 0;     
     
     if ((A.I!=NULL) && (b.I!=NULL))
     /* 1er cas : A et reels */
        for(n=0;n<A.L;n++) {
         (*c).R += A.R[k][n]*b.R[n] + A.I[k][n]*b.I[n];
         (*c).I += A.R[k][n]*b.I[n] - A.I[k][n]*b.R[n];
        }
     else /* 2�me cas : A reel et b complexe */
     if ((A.I==NULL) && (b.I!=NULL))
        for (n=0;n<A.L;n++) {
         (*c).R += A.R[k][n]*b.R[n];
         (*c).I += A.R[k][n]*b.I[n];
        }
     else /* 3�me cas: si A complexes et b reel */
     if ((A.I!=NULL) && (b.I==NULL))
        for (n=0;n<A.L;n++) {
         (*c).R += A.R[k][n]*b.R[n];
         (*c).I -= A.I[k][n]*b.R[n];
        }
     else /* 4�me cas: si A et b reels */
     if ((A.I==NULL) && (b.I==NULL))
        for (n=0;n<A.L;n++) {
         (*c).R += A.R[k][n]*b.R[n];
         (*c).I += 0.0;
        }
 
}


/*=============================================================================
 fonction qui calcule le produit de la matrice A transpos�e conjugu�e par le 
 vecteur b selon A et b reel ou complexe
 >> c = A'*b
 ==============================================================================*/
void Prod_MatConj_Vect(struct vect c, struct Mat A, struct vect b)
{
 int n,k;
     
 /* c est un vecteur colonne de m�me nb d'�lements que le nb de colonnes de A*/
 /* b est un vecteur colonne de meme nb d'elements que le nb de lignes de A */
 
     if ((A.I!=NULL) && (b.I!=NULL))
     /* 1er cas : A et b reels */
       for (k=0;k<A.C;k++) {
         /*initialisation du vecteur c */
         c.R[k] = 0;
         c.I[k] = 0;
         for (n=0;n<A.L;n++) {
             c.R[k] += A.R[k][n]*b.R[n] + A.I[k][n]*b.I[n];
             c.I[k] += A.R[k][n]*b.I[n] - A.I[k][n]*b.R[n];
         }
       }
     else
     /* 2�me cas : A reel et b complexe */
     if ((A.I==NULL) && (b.I!=NULL))
       for (k=0;k<A.C;k++) {
         /*initialisation du vecteur c */
         c.R[k] = 0;
         c.I[k] = 0;
         for (n=0;n<A.L;n++) {
             c.R[k] += A.R[k][n]*b.R[n];
             c.I[k] += A.R[k][n]*b.I[n];
         }
       }
     else
     /* 3�me cas: si A complexe et b reel */
     if ((A.I!=NULL) && (b.I==NULL))
       for (k=0;k<A.C;k++){
         /*initialisation du vecteur c */
         c.R[k] = 0;
         c.I[k] = 0;
         for (n=0;n<A.L;n++) {
             c.R[k] += A.R[k][n]*b.R[n];
             c.I[k] -= A.I[k][n]*b.R[n];
         }
       }
     else
     /* 4�me cas: si A et b reel */
     if ((A.I==NULL) && (b.I==NULL))
       for (k=0;k<A.C;k++) {
         /*initialisation du vecteur c */
         c.R[k] = 0;
         c.I[k] = 0;
         for (n=0;n<A.L;n++) {
             c.R[k] += A.R[k][n]*b.R[n];
             c.I[k] -= 0.0;
         }
       }
 
}
     
  
/*=============================================================================
   Algorithme de GaussSeidel 
   Gauss_Seidel(x,y,W,lambda,critere,nbNZ,n_it_max,seuil)
   Programme qui calcule le vecteur x, le critere et le nombre d'�l�ments non 
   nuls � chaque it�ration
/*=============================================================================*/
void Gauss_Seidel(struct vect x, struct vect y, struct Mat W, double lambda, 
                  double *critere, double *nbNZ, int n_it_max, double seuil)
{

/* Variables de l'algorithme */
struct vect e, r, Wx, Wdag_y, Wdag_W;
int n_it; 

/*variables utilis�es pour le calcul du crit�re */
double mod_xk, mod_rk, pre_terme, sec_terme;

/*variables utilis�es pour le principe de relaxation */
double  w = 1.5;
struct var_comp xk_prec,xk_relax;
double dj, mod2_xk_prec, mod2_xk_relax;

/*variable utilis�e pour le test d'arret */
int test;

/* Variables de stockage */
double module, argument;
struct var_comp condition;
int m,n,k;
int L; /* Ne balaye les parametres nuls qu'une fois sur L */

L = 100; /* W.C/50;*/
/*   mexPrintf("L = %d\n",L); */

/* Allocation de m�moire pour e, r, Wx et We */
   e.R  = (double *)malloc(W.L*sizeof(double));
   e.I  = (double *)malloc(W.L*sizeof(double));
   r.R  = (double *)malloc(W.C*sizeof(double));  
   r.I  = (double *)malloc(W.C*sizeof(double));
   Wx.R = (double *)malloc(W.L*sizeof(double));
   Wx.I = (double *)malloc(W.L*sizeof(double));
   Wdag_y.R = (double *)malloc(W.C*sizeof(double));
   Wdag_y.I = (double *)malloc(W.C*sizeof(double));
   Wdag_W.R = (double *)malloc(W.C*sizeof(double));
   Wdag_W.I = (double *)malloc(W.C*sizeof(double));

/* Calcul de Wdag_y = W'*y utilis� dans le test d'arr�t */
Prod_MatConj_Vect(Wdag_y,W,y);

/* Calcul de la premi�re colonne de la matrice Toeplitz Wdag_W = W'*W(:,1) 
   utilis�e dans le test d'arr�t */
for (n=0;n<W.L;n++) { /* Vecteur contenant la premi�re colonne de W  : e = W(:,1)*/
    e.R[n] = W.R[0][n];
    e.I[n] = W.I[0][n];
}
/* Produit de W conjugu� par la premi�re colonne de W : Wdag_W = W'*e */
Prod_MatConj_Vect(Wdag_W,W,e); 

/* Initialisation de nombre d'iterations */
   n_it=0;
/* Initialisation de valeur du test d'arret*/
   test=0;

/* Initialisation de l'erreur e = y - W*x = y pour x=0;*/
for (n=0;n<W.L;n++) {
    e.R[n] = y.R[n];
    if (y.I==NULL)
       e.I[n] = 0.0;
    else
       e.I[n] = y.I[n];
}
     
/* Boucle sur les it�rations */
while ((test!=W.C)&&(n_it<n_it_max)) {
                
     /****************************************************/
     /* Balayage pour la mise � jour des param�tres x(k) */
     /****************************************************/
     for (k=0;k<W.C;k++) {
     if ( (x.R[k]!=0) || (x.I[k]!=0) || ((double)(n_it/L)==n_it/((double)L)) ) {
         
         /* Mise � jour de l'erreur e = y - W*x sans x[k] */
         /* --------------------------------------------- */
         if ((x.R[k]!=0)||(x.I[k]!=0)) {
            Prod_MatColk_Vectk(Wx,W,x,k); /* Wx = W(:,k)*x(k) */
            for (n=0;n<W.L;n++) {         /* e = e + W(:,k)*x(k) */
                e.R[n] = e.R[n] + Wx.R[n];
                e.I[n] = e.I[n] + Wx.I[n];
            }
         }

         /* Mise � jour de x[k] suivant la valeur de W(k,:)'*e */
         /* -------------------------------------------------- */
         Prod_MatColkConj_Vect(&condition,W,e,k); /* condition = W(k,:)'*e */
         module = abs_C(condition);     
         argument = arg_C(condition);
      
         /* stokage de xk � l'iteration k-1 */
         xk_prec.R = x.R[k];
         xk_prec.I = x.I[k];

	 if (k==(W.C-1)/2) { /* Pour k = P+1 : frequence nulle, on ne penalise pas */
            x.R[k] = module*cos(argument)/W.L;
            x.I[k] = module*sin(argument)/W.L;
	 }
	 else { 
            /* calcul de xk selon le module(W(k,:)'e) */
            if (module < lambda) { /* x(k) = 0; */
               x.R[k]=0;
               x.I[k]=0; 
            } 
            else { /* x(k) = ( abs(W(k,:)'e) - lambda )/N*exp(j*arg(W(k,:)'e)) */
               x.R[k] = (module - lambda)*cos(argument)/W.L;
               x.I[k] = (module - lambda)*sin(argument)/W.L;
	    }
         }

        /* calcul de xk avec relaxation (param�tre w) */     
        /* -------------------------------------------------- */
        /* xk_relax = xk_prec +w(x(k) - xk_prec) */ 
        xk_relax.R = xk_prec.R + w*(x.R[k] - xk_prec.R);
        xk_relax.I = xk_prec.I + w*(x.I[k] - xk_prec.I);
      
        /* Difference du critere avec ou sans relaxation */
        mod2_xk_prec  = abs2_C(xk_prec);
        mod2_xk_relax = abs2_C(xk_relax);
        dj = W.L/2*(mod2_xk_relax - mod2_xk_prec) 
           - condition.R*(xk_relax.R - xk_prec.R)
           - condition.I*(xk_relax.I - xk_prec.I)
           + lambda*(sqrt(mod2_xk_relax) - sqrt(mod2_xk_prec));
      
        if (dj <0) { /* Si relaxation r�duit le crit�re, on relaxe */
          x.R[k] = xk_relax.R;
          x.I[k] = xk_relax.I;
        }

        /* Mise � jour de l'erreur e = y - W*x avec x(k) */
        /* --------------------------------------------- */
        if ((x.R[k]!=0)||(x.I[k]!=0)) {
          Prod_MatColk_Vectk(Wx,W,x,k); /* Wx = W(:,k)*x(k) */
          for (n=0;n<W.L;n++) {         /* e = e - W(:,k)*x(k) */
              e.R[n] = e.R[n]- Wx.R[n];
              e.I[n] = e.I[n]- Wx.I[n];
          }
        }

     }
     } /* Fin de balayage des param�tres x(k) */
     
     
     /*********************************************************/
     /* Calcul du vecteur r du test d'arr�t r = W'*W*x - W'y; */     
     /*********************************************************/
     if ((double)(n_it/L)==n_it/((double)L)) { /* Test d'arret une fois sur L */
     for (k=0;k<W.C;k++) {
         r.R[k] = - Wdag_y.R[k];
         r.I[k] = - Wdag_y.I[k];
     }
	 /* Calcul de W'*W*x pour les elements non nuls de x avec W'*W Toeplitz 
        dont on dispose de la premi�re colonne */
	 for (n=0;n<W.C;n++) 
         if ((x.R[n]!=0)||(x.I[n]!=0)) {
            for (k=0;k<n;k++) { /* W'W(k,n) = conj(W(0,n-k)) */
                r.R[k] += Wdag_W.R[n-k]*x.R[n] 
                        + Wdag_W.I[n-k]*x.I[n]; 
	            r.I[k] += Wdag_W.R[n-k]*x.I[n] 
                        - Wdag_W.I[n-k]*x.R[n]; 
            }
            for (k=n;k<W.C;k++) { /* W'W(k,n) = W(0,k-n) */
                r.R[k] += Wdag_W.R[k-n]*x.R[n] 
                        - Wdag_W.I[k-n]*x.I[n]; 
	            r.I[k] += Wdag_W.R[k-n]*x.I[n] 
                        + Wdag_W.I[k-n]*x.R[n]; 
            }
         } 
     
     /*****************************************************/
     /* Balayage des param�tres x(k) pour le test d'arret */
     /*****************************************************/
     test = 0;
     /* Second terme du critere : somme des modules des x(k) */
     sec_terme = 0;
     for (k=0;k<W.C;k++) {           
        mod_xk = abs_VectCk(x,k); /* mod_xk = abs(x(k)) */
        /* Test d'arret */
        if (k==(W.C-1)/2) { /* Pour k = P+1 : frequence nulle, on ne penalise pas */
             mod_rk = abs_VectCk(r,k); /* mod_rk = abs(r(k)) */
             if (mod_rk<=  seuil)
                test++;
             if (nbNZ!=NULL) {   
                if (mod_xk>seuil)
                   nbNZ[n_it] = nbNZ[n_it] + 1;
             }
	} 
	else {
           if (mod_xk<=seuil) {/* abs(rk) < lambda + seuil si x(k)=0 */
               mod_rk = abs_VectCk(r,k); /* mod_rk = abs(r(k)) */
               if (mod_rk<= lambda + seuil)
                  test++;
           }
           else { /* abs( rk + lambda*x(k)/abs(xk)) < seuil sinon */
               condition.R = r.R[k] + lambda*x.R[k]/mod_xk;
               condition.I = r.I[k] + lambda*x.I[k]/mod_xk;
               module = abs_C(condition);     
               if (module<seuil)
                  test++;
               /* Mise � jour du second terme du crit�re */
               sec_terme += mod_xk;
               if (nbNZ!=NULL) 
                  nbNZ[n_it] = nbNZ[n_it] + 1;
           }
	}
     }
          
     /***********************************/
     /* Calcul et stockage du crit�re J */
     /***********************************/
     if (critere!=NULL) {
        /* Calcul du premier terme du crit�re : norme de l'erreur */
        pre_terme = 0;
        for (n=0;n<W.L;n++) 
            pre_terme += abs2_VectCk(e,n); 
            critere[n_it] = pre_terme/2 + lambda*sec_terme;
     mexPrintf("%d : J = %f; test = %d\n",n_it,critere[n_it],test);
     fflush(NULL); 
     }        
     }
     n_it +=1;
          
 } /*fin de boucle while sur les it�rations */
   
if (n_it<n_it_max)
   mexPrintf("min_L2_L1_0 converged in %d iterations\n",n_it);
else
   mexPrintf("min_L2_L1_0 did not converged in %d iterations\n",n_it);
   

/* Lib�ration de la m�moire */     
   free(e.R);
   free(e.I);
   free(r.R);
   free(r.I);
   free(Wx.R);
   free(Wx.I);
   free(Wdag_y.R);
   free(Wdag_y.I);
   free(Wdag_W.R);
   free(Wdag_W.I);
  
} /* fin du programme */
     
/*=============================================================================
 Interface Matlab/C de la fonction gauss_seidel
 >> [ x, J, NZ ] = Min_L2_L1(y,W,lambda[,n_it_max,seuil]);
 renvoie le vecteur x minimisant le crit�re 
        J(x) = 1/2||y-W*x||_2^2 + lambda*||x||_1
 avec J le vecteur donnant la valeur du crit�re � chaque it�ration
 et NZ le vecteur donnant le nombre d'�l�ments non nuls � chaque it�ration
 
 ==============================================================================*/
      
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  struct Mat W;
  struct vect x, y; 
  double *J, *NZ;
  double lambda, seuil, *tmp;
  int    n_it_max, mrows, ncols;

  /* Verification du nombre d'arguments */
  if ((nrhs<3) || (nrhs > 5)) {
     mexPrintf("min_L2_L1 : three to five input required:\n");
     mexPrintf("  [x, J, NZ] = min_L2_L1(y,W,lambda,[n_it_max,thres]);\n");
     mexPrintf("  Minimise criterion J(x) = 1/2||y-W*x||_2^2 + lambda*||x||_1\n");
     mexErrMsgTxt("  Error using min_L2_L1");
  } 
  if ((nlhs<1) || (nlhs > 3)) {
     mexPrintf("min_L2_L1 : one to three output requires:\n");
     mexPrintf("  [x, J, NZ] = min_L2_L1(y,W,lambda,[n_it_max,thres]);\n");
     mexPrintf("  Minimise criterion J(x) = 1/2||y-W*x||_2^2 + lambda*||x||_1\n");
     mexErrMsgTxt("  Error using min_L2_L1");
  }



  /* R�cup�ration de l'entree y */
  y.C = mxGetM(prhs[0]);  /* N nb lignes de y */
  y.R = (double *)mxGetPr(prhs[0]);
  y.I = (double *)mxGetPi(prhs[0]);
  /*printf("y  = [");
  for (ncols=0;ncols<y.C;ncols++) {
      if (y.I!=NULL)
          printf(" %.4f + %.4fi ;",y.R[ncols],y.I[ncols]);
      else
          printf(" %.4f ;",y.R[ncols]);
   }  printf("]\n");*/
  
  
  /* R�cup�ration de l'entree W */
  W.C = mxGetN(prhs[1]);  /* nb colonnes de W */
  W.L = mxGetM(prhs[1]);  /* nb lignes de W */
  tmp = (double *)mxGetPr(prhs[1]); 
  /* Stockage sous forme default:matrice complexe */
  W.R = (double **)malloc(W.C*sizeof(double *));  
  for (ncols=0;ncols<W.C;ncols++)
      W.R[ncols] = &tmp[ncols*W.L];
  tmp = (double *)mxGetPi(prhs[1]);
  if (tmp!=NULL) {
     W.I = (double **)malloc(W.C*sizeof(double *));  
     for (ncols=0;ncols<W.C;ncols++)
         W.I[ncols] = &tmp[ncols*W.L];
  } else W.I=NULL;
  /*printf("W  = ");
  for (mrows=0;mrows<W.L;mrows++) {
      printf("\t[ ");
      for (ncols=0;ncols<W.C;ncols++) {
      if (W.I!=NULL)
          printf(" %.4f + %.4fi ;",W.R[ncols][mrows],W.I[ncols][mrows]);
      else
          printf(" %.4f ;",W.R[ncols][mrows]);
      } printf("]\n");
  }*/
      
  
 /* Verification de l'entree lambda.*/
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
     !(mrows == 1 && ncols == 1)) {
     mexPrintf("Input lambda must be a noncomplex scalar:\n");
     mexErrMsgTxt("  [x, J, NZ] = Gauss_seidel(y,W,lambda,[n_it_max,thres]);");
  } else { lambda = (double)mxGetScalar(prhs[2]);
  }
  /* printf("lambda  = %e\n",lambda); */
  
  /* Verification de l'entree n_it_max */
  if (nrhs == 3) {
     n_it_max = (int)100000;
  } else {
    mrows = mxGetM(prhs[3]);
    ncols = mxGetN(prhs[3]);
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
       !(mrows == 1 && ncols == 1)) {
       mexPrintf("Input n_it_max must be a noncomplex scalar:\n");
       mexErrMsgTxt("  [x, J, NZ] = Gauss_seidel(y,W,lambda,[n_it_max,thres]);");
    } else { n_it_max = (int)mxGetScalar(prhs[3]);
    }
  }
  /* printf("nb it max  = %d\n",n_it_max); */

 /* Verification de l'entree seuil.*/
  if ((nrhs == 3)||(nrhs == 4)) {
     seuil = (double)1e-10;
  } else {
    mrows = mxGetM(prhs[4]);
    ncols = mxGetN(prhs[4]);
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
       !(mrows == 1 && ncols == 1)) {
       mexPrintf("Input thres seuil must be a noncomplex scalar:\n");
       mexErrMsgTxt("  [x, J, NZ] = Gauss_seidel(y,W,lambda,[n_it_max,thres]);");
    } else { seuil = (double)mxGetScalar(prhs[4]);
    }
  }
  /*printf("thres  = %e\n",thres);*/

  /* Creation du vecteur pour l'argument de sortie x */
  /* x est un vecteur colonne de meme nb d'elements que le nb de colonnes de W */
  x.C = W.C; 
  /*  if ((W.I==NULL) && (y.I==NULL)) 
     plhs[0] = mxCreateDoubleMatrix(x.C,1, mxREAL);
  else */
  plhs[0] = mxCreateDoubleMatrix(x.C,1, mxCOMPLEX);
  x.R = mxGetPr(plhs[0]); x.I = mxGetPi(plhs[0]);

  /* Creation du vecteur pour l'argument de sortie J */
  if (nlhs>1) {
     plhs[1] = mxCreateDoubleMatrix(n_it_max,1, mxREAL);
     J = mxGetPr(plhs[1]); 
  } else J = NULL;   

  /* Creation du vecteur pour l'argument de sortie NZ */
  if (nlhs>2) {
     plhs[2] = mxCreateDoubleMatrix(n_it_max,1, mxREAL);
     NZ = mxGetPr(plhs[2]); 
  } else NZ = NULL;   

  /* Appel a la fonction. */
    Gauss_Seidel(x,y,W,lambda,J,NZ,n_it_max,seuil); 
      
  free(W.R);
  free(W.I);
} /* Fin MexFunction */

