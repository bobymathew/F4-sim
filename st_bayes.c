#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include "ranlib.h"
float rank(float *A,int a1,float *B);
void inverse(float *a,int li,float *CI);
void inv_A(float **pr,int in,float *Ai);
void sqrtm(float **pp,int li,float *H);

extern void spotrf_(char *uplo, int *n, float *a, int *lda, int *info);
extern void ssyev_( char *jobz, char *uplo, int *n, float *a, int *lda,
        float *w, float *work, int *lwork, int *info );
extern void spotri_(char *uplo, int *n, float *a, int *lda,int *info);
extern void sgetri_(int *n, float *a, int *lda, int *ipiv,float *work, int *lwork, int *info);
extern void sgetrf_(int *m, int *n, float *a, int *lda,int *ipiv, int *info);
/*extern void sgesv_(int *n, int *nrhs, float *a, int *lda, 
	int *ipiv, float *b, int *ldb, int *info);*/

struct data{
int records,order,fix,ord_fix;
float **pr,*Y,*C,*null_space,*CI,*A_in,*A_sqr,*ran_ord,*a_star1,*ran_reco,*diag;
float *sa,*se,*Wt,*Y_Wt,va,ve,k_rat,*Z_1,*z_star,*Wf,*LH,*C_Wf,*t,*t_tr,*A_tr,*Y_1,*x_sol;

};

int  main()
{
struct data base,*pbase;
pbase=&base;
int i,j,*sir,*dam,*p,tt,y=1,va=1,ve=1,l=0,*ipiv,info,bb=1;
float *in,k;
/*Matrices for W*/
int *X1,**Z,no_colu_record=4,bval=0,store=6,iteration=10,th=5,niter;
float **data,**W,**E,*W_1,*Cw,*E_1,an=0.0,d=1.0,**bvalue,temp=0.0,bsum=0.0,*bmean,*bmedian;

pbase->order=100,pbase->records=100;
pbase->fix=1;
pbase->va=1.0;pbase->ve=1.0;pbase->k_rat=0.0;
pbase->ord_fix=pbase->order+pbase->fix;
FILE *fpx,*fpd,*fpy;

/*pedgree data*/
sir= malloc((pbase->order+1) * sizeof (int));
dam= malloc((pbase->order+1) * sizeof (int));
p= malloc((pbase->order+1) * sizeof (int));
in= malloc((pbase->order+1) * sizeof (float));

time_t now;
time(&now);
printf("It’s %s", ctime(&now));

/*Matrices for W*/

X1=malloc((pbase->records) * sizeof (int));
pbase->Y=malloc((pbase->records) * sizeof (float));
pbase->Y_1=malloc((pbase->records) * sizeof (float));

bmean=(float *)malloc(pbase->ord_fix* sizeof (float));
bmedian=(float *)malloc(pbase->ord_fix* sizeof (float));

data= (float**)calloc(pbase->records,sizeof(float *));
     for (i = 0; i < pbase->records; i++) {
    data[i] = (float *)calloc(no_colu_record ,sizeof(float));  }
Z= (int**)calloc((pbase->records),sizeof(int *));
     for (i = 0; i < pbase->records; i++) {
    Z[i] = (int *)calloc(pbase->order,sizeof(int));  }

bvalue= (float**)calloc(pbase->ord_fix,sizeof(float *));
     for (i = 0; i < pbase->ord_fix; i++) {
    bvalue[i] = (float *)calloc(store,sizeof(float));  }

pbase->Z_1=(float *)malloc(pbase->records*pbase->order* sizeof (float));

W= (float**)calloc((pbase->records),sizeof(float *));
     for (i = 0; i < (pbase->records); i++) {
     W[i] = (float *)calloc(pbase->ord_fix,sizeof(float));  }

W_1=(float *)malloc(pbase->records*pbase->ord_fix* sizeof (float));

pbase->Wf=(float *)malloc(pbase->ord_fix*1* sizeof (float));

E= (float**)calloc((pbase->ord_fix),sizeof(float *));
     for (i = 0; i < pbase->ord_fix; i++) {
    E[i] = (float *)calloc(pbase->ord_fix,sizeof(float));  }

E_1=(float *)malloc(pbase->ord_fix*pbase->ord_fix* sizeof (float));

Cw=(float *)malloc(pbase->ord_fix*pbase->ord_fix* sizeof (float));

pbase->C=(float *)malloc(pbase->ord_fix*pbase->ord_fix* sizeof (float));
pbase->x_sol=(float *)malloc(1*pbase->ord_fix* sizeof (float));

pbase->CI=(float *)malloc(pbase->ord_fix*pbase->ord_fix* sizeof (float));

pbase->null_space=(float *)malloc(pbase->ord_fix*pbase->ord_fix* sizeof (float));

pbase->pr= (float**)calloc((pbase->order+1),sizeof(float *));
     for (i = 0; i < pbase->order+1; i++) {
    pbase->pr[i] = (float *)calloc((pbase->order+1) ,sizeof(float));  }

pbase->A_in = (float *)malloc(pbase->order*pbase->order*sizeof(float));
pbase->A_sqr = (float *)malloc(pbase->order*pbase->order*sizeof(float));

pbase->ran_ord = (float *)malloc(pbase->order*1*sizeof(float));
pbase->ran_reco = (float *)malloc(pbase->records*1*sizeof(float));

pbase->a_star1= (float *)malloc(pbase->order*1*sizeof(float));

pbase->z_star=(float *)malloc(pbase->records*1*sizeof(float));

pbase->LH=(float *)malloc(pbase->ord_fix*1* sizeof (float));

pbase->C_Wf=(float *)malloc(pbase->ord_fix*1* sizeof (float));

pbase->t=(float *)malloc(pbase->ord_fix*1* sizeof (float));

pbase->t_tr=(float *)malloc(pbase->order*1* sizeof (float));

pbase->A_tr=(float *)malloc(pbase->order*1* sizeof (float));

pbase->diag=(float *)malloc(pbase->ord_fix*sizeof (float));

pbase->sa = (float *)malloc(1*1*sizeof(float));
pbase->se = (float *)malloc(1*1*sizeof(float));

pbase->Wt=(float *)malloc(pbase->records*1* sizeof (float));
pbase->Y_Wt=malloc(pbase->records*1* sizeof (float));

if((fpx=fopen("pedgree_test.txt","r+"))==NULL){printf("Cant open file");}
if((fpd=fopen("data_test.txt","r+"))==NULL){printf("Cant open file");}

fpy=fopen("bvalue_5569.txt","w+");
for(i=0;i<pbase->order+1;i++)
	{
 		fscanf(fpx,"%d %d %d %f",&p[i],&sir[i],&dam[i],&in[i]);
	}
printf("%f\n",sir[2]);

/*calculation the A matrix*/		
	for(i=1;i< pbase->order+1;i++)
	{

 		for(j=i;j< pbase->order+1;j++)
		 {
 	 	if(i!=j)
                  {
                     if(( sir[j]&& dam[j])==0)
                     {
                      k=0.0;
                      pbase->pr[i][j]=k+pbase->pr[i][j];
                      pbase->pr[j][i]=pbase->pr[i][j];
                     }
		     if(sir[j]==0) 
			{
				if(dam[j]!=0)
				{ 
				k=0.0;
				k=0.5*pbase->pr[i][dam[j]];
				pbase->pr[i][j]=k+pbase->pr[i][j];
				pbase->pr[i][j]=pbase->pr[j][i];
				}
			}
		 if(sir[j]!=0)
			{
			if(dam[j]==0)
				{
				k=0.0;
				k=0.5*(pbase->pr[i][sir[j]]);
				pbase->pr[i][j]=k+pbase->pr[i][j];
				pbase->pr[j][i]=pbase->pr[j][i];
				}
			}
                     if((dam[j] && sir[j])!=0)
                     {
                        k=0.5*(pbase->pr[i][sir[j]]+pbase->pr[i][dam[j]]);
                        pbase->pr[i][j]=k+pbase->pr[i][j];
                        pbase->pr[j][i]=pbase->pr[i][j];
                      }
		
		}
	else
	    	{
		pbase->pr[i][j]=1+in[i];
			
		}		
		}
        }
free(sir);free(dam);free(in);free(p);
/*Creating W matrix*/

for(i=0;i<pbase->records;i++)
        {
                for(j=0;j<no_colu_record;j++)
                {
                fscanf(fpd,"%f",&data[i][j]);
		}
        }
printf("%f\n",data[2][2]);
for(i=0;i<pbase->records;i++)
        {
	X1[i]=1;
	pbase->Y[i]=data[i][no_colu_record-1];
	}

for(i=0;i<pbase->records;i++)
        {
                for(j=1;j<pbase->order+1;j++)
                {
                if(j==data[i][0]){Z[i][j-1]=1;}
                else{Z[i][j-1]=0;}
		}
        }
for(i=0;i<pbase->records;i++)
{
	tt=0;
	W[i][tt]=X1[i];
		 
	 for(j=0;j<pbase->order;j++)
        {
                tt++;
		W[i][tt]=Z[i][j];
		pbase->Z_1[i*pbase->order+j]=(float)Z[i][j];
        }
}
l=0;

/*one dimension W*/
for(i=0;i<pbase->records;i++)
{	
	for(j=0;j<pbase->ord_fix;j++)
        {
               	W_1[i*pbase->ord_fix+j]=W[i][j];
	}
}
cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, pbase->ord_fix, pbase->ord_fix,
                 pbase->records, 1.0f, W_1,pbase->ord_fix, W_1,pbase->ord_fix, 0.0f,
                 Cw, pbase->ord_fix);/* W*W'*/

		free(X1);
		for (i = 0; i < pbase->records; i++)
	 		{
        		free(data[i]);
	 		free(Z[i]);
	 		free(W[i]);
	 		}

/*E matrix*/
for(i=0;i<pbase->ord_fix;i++)
{                for(j=0;j<pbase->ord_fix;j++)
                 {
                      E[i][j]=0.0;
                 }
}
for(i=0;i<pbase->fix;i++)
{
                for(j=0;j<pbase->ord_fix;j++)
                 {
                        E[i][j]=0.0;
                 }
}
for(i=pbase->fix;i<pbase->ord_fix;i++)
{
tt=pbase->fix;
		for(j=1;j< pbase->order+1;j++)
		{
		E[i][tt]=pbase->pr[y][j];
		tt++;
		}
               y++;
}
/*one dimension E*/
for(i=0;i<pbase->ord_fix;i++)
{	
	for(j=0;j<pbase->ord_fix;j++)
        {
               	E_1[i*pbase->ord_fix+j]=E[i][j];
	}
}

		for(i=0;i<pbase->ord_fix;i++)
		{
		free(E[i]);
		}
	
for(i=0;i<pbase->ord_fix*pbase->ord_fix;i++)
{	
	pbase->C[i]=Cw[i]+E_1[i];
}



/*rank(pbase->C,pbase->ord_fix,pbase->null_space);*/

inv_A(pbase->pr,pbase->order,pbase->A_in);/*inv(A)*/
sqrtm(pbase->pr,pbase->order,pbase->A_sqr);/*sqrtm(A)*/

/* LOOP*/

/*for(niter=1;niter<=10;niter++)
{*/
for(i=0;i<pbase->order;i++)
	{
 		pbase->ran_ord[i]=gennor(an,d);
	}

for(i=0;i<pbase->records;i++)
	{
 		pbase->ran_reco[i]=gennor(an,d);
	}

/* sqrtm(A)*randn(n,1) */
cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, pbase->order, 1,
                 pbase->order, 1.0f, pbase->A_sqr, pbase->order, pbase->ran_ord, pbase->order, 0.0f,pbase->a_star1,1);

for(i=0;i<pbase->order;i++)
	{
	pbase->a_star1[i]=sqrt(pbase->va)*pbase->a_star1[i];/* (sqrt(va)*a_star1')*/
	}

for(i=0;i<pbase->records;i++)
	{
	pbase->ran_reco[i]=sqrt(pbase->ve)*pbase->ran_reco[i];/* sqrt(ve)*randn(records,1) */
	}
/* z_star=Z*a_star' */
cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, pbase->records, 1,
                 pbase->order, 1.0f, pbase->Z_1, pbase->order, pbase->a_star1, pbase->order, 0.0f,pbase->z_star, 1);


for(i=0;i<pbase->records;i++)
	{
	pbase->z_star[i]=pbase->z_star[i]+pbase->ran_reco[i];/* z_star=Z*a_star'+sqrt(ve)*randn(records,1); */
	}

for(i=0;i<pbase->records;i++)
	{
	pbase->Y_1[i]=pbase->Y[i]-pbase->z_star[i];/* y-z_star; */
	}
/* Wf=W'*(y-z_star); */
cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, pbase->ord_fix, 1,
                 pbase->records, 1.0f, W_1,pbase->ord_fix, pbase->Y_1,1, 0.0f,
                pbase->Wf, 1);
for(i=0;i<pbase->ord_fix;i++)
{
 printf("t=%f\n",pbase->Wf[i]);
}

for(i=0;i<pbase->fix;i++)/* LH=[zeros(1,n_fix),a_star]' */
	{
pbase->LH[i]=0.0;
	}
for(i=pbase->fix;i<pbase->ord_fix;i++)
                 {
		pbase->LH[i]=pbase->a_star1[i-pbase->fix];
		 }
ipiv = (int *) malloc (pbase->ord_fix * sizeof(int));

/*sgesv_(&pbase->ord_fix,&bb, pbase->C, &pbase->ord_fix, 
	ipiv, pbase->Wf, &pbase->ord_fix, &info);
printf("%d\n",info);
/*inverse(pbase->C,pbase->ord_fix,pbase->CI);/*inv(C)*/
/*(inv(C)*Wf) */
/*cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, pbase->ord_fix, 1,
                 pbase->ord_fix, 1.0f, pbase->CI, pbase->ord_fix, pbase->Wf, pbase->ord_fix, 0.0f,pbase->C_Wf,1);*/
/*for(i=0;i<pbase->ord_fix;i++)
{
 fprintf(fpy,"t=%f\n",pbase->Wf[i]);
}*/
for(i=0;i<pbase->ord_fix;i++)
                 {
		pbase->t[i]=pbase->LH[i]+pbase->LH[i];/* t=LH+(pinv(C)*Wf); */
		 }


for(i=pbase->fix;i<pbase->ord_fix;i++)/*t((n_fix+1):(n+n_fix))*/
	{
	pbase->t_tr[i-pbase->fix]=pbase->t[i];
	}

/* t((n_fix+1):(n+n_fix))'*inv(A) */
cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, 1,pbase->order ,
                 pbase->order, 1.0f, pbase->t_tr, 1, pbase->A_in, pbase->order, 0.0f,pbase->A_tr,pbase->order);
/* t((n_fix+1):(n+n_fix))'*inv(A)*t((n_fix+1):(n+n_fix)) */
cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1,1 ,
                 pbase->order, 1.0f, pbase->A_tr, pbase->order, pbase->t_tr, pbase->order, 0.0f,pbase->sa,1);
/* W*t */
cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,pbase->records ,1 ,
                 pbase->ord_fix, 1.0f, W_1, pbase->ord_fix, pbase->t, 1, 0.0f,pbase->Wt,1);

for(i=0;i<pbase->records;i++)/* y-W*t */
	{
	pbase->Y_Wt[i]=pbase->Y[i]-pbase->Wt[i];
	}

/*(y-W*t)'*(y-W*t)*/
cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1,1 ,
                 pbase->records, 1.0f, pbase->Y_Wt, pbase->records, pbase->Y_Wt, pbase->records, 0.0f,pbase->se,1);

pbase->va=(pbase->sa[0])/genchi(pbase->order);
pbase->ve=(pbase->se[0])/genchi(pbase->records);

pbase->k_rat=pbase->ve/pbase->va;


for(i=0;i<pbase->fix;i++)/* LH=[zeros(1,n_fix),a_star]' */
	{
pbase->diag[i]=0.0;
	}
for(i=pbase->fix;i<pbase->ord_fix;i++)
                 {
		pbase->diag[i]=pbase->k_rat;
		
		 }

for(i=0;i<pbase->ord_fix;i++)
                 {
for(j=0;j<pbase->ord_fix;j++)
                 {
		if(i==j){pbase->C[i*pbase->ord_fix+j]=Cw[i*pbase->ord_fix+j]+pbase->diag[i];}
		else{pbase->C[i*pbase->ord_fix+j]=Cw[i*pbase->ord_fix+j];}
		 }
		}
/*if(niter>=th && niter%1==0)
		{
		for(j=0;j<pbase->ord_fix;j++)
			{
			bvalue[j][bval]=pbase->t[j];
			}
		bval++;
		}
}*/
free(pbase->Y);free(pbase->Z_1);free(W_1);free(pbase->Wf);free(E_1);
free(Cw);free(pbase->C);free(pbase->CI);free(pbase->null_space);free(pbase->A_in);
free(pbase->A_sqr);free(pbase->ran_ord);free(pbase->ran_reco);free(pbase->a_star1);free(pbase->z_star);
free(pbase->LH);free(pbase->C_Wf);free(pbase->t);free(pbase->t_tr);free(pbase->A_tr);
free(pbase->diag);free(pbase->sa);free(pbase->se);free(pbase->Wt);free(pbase->Y_Wt);
for (i = 0; i < pbase->order+1; i++)
	 		{
        		free(pbase->pr[i]);
			}
/*for(i=0;i<pbase->ord_fix;i++)
	{
	bsum=0.0;
                 for(j=0;j<store;j++)
                 {
                        bsum=bsum+bvalue[i][j];
						
	         }bmean[i]=bsum/store;
		
	}
	for(i=0;i<pbase->ord_fix;i++)
	{
	      for(j=0;j<store;j++)
		for(l=j+1;l<store;l++)
		{
		if(bvalue[i][j]>bvalue[i][l])
		{
		temp=bvalue[i][l];
		bvalue[i][l]=bvalue[i][j];
		bvalue[i][j]=temp;
		}
		}
		bmedian[i]=(bvalue[i][store/2]+bvalue[i][(store/2)-1])/2;
		
	}*

for(i=0;i<pbase->ord_fix;i++)
{
 fprintf(fpy,"mn=%f\n",bmean[i]);
}
for(i=0;i<pbase->ord_fix;i++)
{
 fprintf(fpy,"md=%f\n",bmedian[i]);
}
/*
free(bmedian);free(bmean);

for(i=0;i<pbase->ord_fix;i++)
	{
	free(bvalue[i]);
	}*/

fclose(fpx);fclose(fpd);fclose(fpy);
time_t new;
time(&new);
printf("It’s %s", ctime(&new));
return 0;
}
/* functions*/

void inverse(float *a,int li,float *CI)
{

    	int info,i,n=li;
	for(i=0;i<n*n;i++)
	{
    	CI[i]=a[i];
	}
	spotrf_("L", &n, CI, &n,&info );
	printf("InverseA1=%d\n",info);
        spotri_("L", &n, CI, &n,&info );
	printf("InverseA2=%d\n",info);
	
}

void inv_A(float **pr,int in,float *B)
{
    float *a;
	int info,i,j,q=0,n=in;
		        
        a = (float *)malloc(n*n*sizeof(float));
	
	for (i = 1; i <n+1; i++){
				for (j = 1; j <n+1; j++){
							a[q]=pr[i][j];
							++q;
							}
				}
	spotrf_("L", &n, a, &n,&info );
	printf("InverseA1=%d\n",info);
        spotri_("L", &n, a, &n,&info );
	printf("InverseA2=%d\n",info);
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
		B[i*n+j]=B[j*n+i]=a[j*n+i];
		}
	}

  	free(a);
}
void sqrtm(float **prr,int uu,float *H)
{
float *a,*work, *w,*re;
        int lwork,info=1,i,j,q=0,n=uu;
	
	
        re= (float*)calloc(n*n,sizeof(float *));
	a = (float *)malloc(n*n*sizeof(float));
        w = (float *)malloc(n*sizeof(float));
	for (i = 1; i < n+1; i++){
				for (j = 1; j < n+1; j++){
							a[q]=prr[i][j];
							++q;
							}
				}
	work = (float *)malloc(1*sizeof(float));
        lwork=-1;
        ssyev_( "V", "U", &n, a, &n, w, work, &lwork, &info );
	printf("sqrtA1=%d\n",info);
        lwork= work[0];
	free(work);
        work = (float *)malloc(lwork*sizeof(float));
	ssyev_( "V", "U", &n, a, &n, w, work, &lwork, &info );
	printf("sqrtA2=%d\n",info);
		for(i=0;i<n*n;i++)
		{
		re[i]=a[i]*sqrt(w[i/n]);
		}
cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n,
                 n, 1.0f, re, n, a, n, 0.0f,H, n);


	free(re);
        free(work);
        free(w);
	free(a);
}
/*float rank(float *A,int a1,float *B)

{
int i,rank=0,N=a1;
FILE *fpr;
fpr=fopen("ZE_1.txt","w+");

int LWORK,INFO,temp,j,k=0,tem2,l=0,*ord,rank_def;
float *s,*wk,*uu,*vt,*vR,*null,null_space;

ord=malloc(N* sizeof (int));
s=malloc(N* sizeof (float));
uu=malloc(N* sizeof (float));
vt=malloc(N*N*sizeof (float));
vR=malloc(N*N* sizeof (float));
wk = malloc(1*sizeof(float));


for(i=0;i<N;i++)
{
	ord[i]=i;
}

LWORK = -1;

sgeev_( "V", "V", &N,A, &N, s, uu,vt, &N, vR, &N,wk,&LWORK,&INFO);

LWORK=wk[0];

free(wk);
wk = malloc(LWORK*sizeof(float));

sgeev_( "V", "V", &N,A, &N, s, uu,vt, &N, vR, &N,wk,&LWORK,&INFO);
for (i = 0; i < N; ++i)
     {
	if(fabs(s[i]) > 0.000001) {rank++;}
	fprintf(fpr,"%f\n",s[i]);
      }

rank_def=N-rank;
printf("rank=%d",rank);
null=malloc(N*rank_def* sizeof (float));
for (i = 0; i < (N- 1); ++i)
     {
	
          for (j = 0; j < N- 1 - i; ++j )
          {
               if (s[j] < s[j+1])
               {
                    temp = s[j+1];
                    s[j+1] = s[j];
                    s[j] = temp;
		    tem2 = ord[j+1];
                    ord[j+1] = ord[j];
                    ord[j] = tem2;
				
               }
		
          }
	
     }

for (i = 0; i < rank_def; ++i)
	{
for (j = 0; j < N; ++j)
	{	
vR[l]=vt[i*N+j];
	}
	}
for (i = 0; i < N*rank_def; ++i)
	{
null[i]=vR[i];
	}
cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N,
                 rank, 1.0f, null, N, null, N, 0.0f,B, N);


free(ord);free(s);free(uu);free(vt);free(vR);
}*/


