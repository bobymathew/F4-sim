#include "ranlib.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
void shuffle( int *array, int n,int k);
void insert(int *a,int *b, int l);
float sum_base(float *TT,int kk,int oo);
float crossing(int SIR1,int DAM1,int first1);
float second_crossing(int rr,int ll);
float calculate1(int r);
float calculate3(int r);
float calculate2(int r);
float third_crossing(int third,int first);
float fourth_crossing(int third1,int first1);
int *base_cros,n=5000,*first_zero,*base_zero,*second_zero,*homo,suff,suff2=50,*third_zero,*fourth_zero,base=20;
float **matrx,*bpr,*first_loc,*base_loc,*second_loc,*homo_loc,*third_loc,*homo_loc2,*fourth_loc,**fourth_sum,**ped;
int main()
{
int size_y=11,t=200,dum[n],h=0,kk=0,f1=0,f2=0,f3=0,f4=0,si=0,da=0,par,pp=0,par2,par3,t1=0,t2=0,t3=0,t4=0,g1=5,g2=500,g3=5000,g4=50000;
int i,j,final=15790,sir[5],dam[5],l=0,a[55570],first=0,generation=1,first_gene=0,second_gene=0,third_gene=0,fourth_gene=0,first_gene_second=0,second_value=1,in,gg=0,in_l,track=0,pd=0,f_2=0,f_3=0;
float an=0.0,d=1.0,sum=0.0,loc_eff[6];
FILE *fp,*fz,*sz,*fl,*sl;

bpr= malloc(n * sizeof (int));
homo_loc=malloc(n * sizeof (float));
homo_loc2=malloc(n * sizeof (float));

matrx= (float**)calloc(final,sizeof(float*));
	 for (i = 0; i < final; i++)
	 {
	matrx[i] = (float*)calloc(size_y ,sizeof(float));
	 }

ped= (float**)calloc(final,sizeof(float*));
	 for (i = 0; i < final; i++)
	 {
	ped[i] = (float*)calloc(5 ,sizeof(float));
	 }
base_loc=(float *)malloc(base*n* sizeof (float));
base_cros=(int *)malloc(base*n* sizeof (int));
base_zero=(int *)malloc(base*n* sizeof (int));
first_zero=(int *)malloc(g1*n* sizeof (int));

first_loc=(float *)malloc(g1*n* sizeof (float));

second_zero=(int *)malloc(g2*n* sizeof (int));

second_loc=(float *)malloc(g2*n* sizeof (float));

third_loc=(float *)malloc(g3*n* sizeof (float));

third_zero=(int *)malloc(g3*n* sizeof (int));

fourth_loc=(float *)malloc(g4*n* sizeof (float));

fourth_zero=(int *)malloc(g4*n* sizeof (int));

fourth_sum=(float**)calloc(50000,sizeof(float*));
	 for (i = 0; i < 50000; i++)
	 {
	 fourth_sum[i] = (float*)calloc(3,sizeof(float));
	 }

fp=fopen("f4.txt","w+");
/*fz=fopen("ped.txt","w+");
fl=fopen("secon_loc.txt","w+");
sz=fopen("fourth_zero.txt","w+");
sl=fopen("fourth_loc.txt","w+");*/
for (i=0;i<5;i++)
       {
	sir[i]=rand()%18+1;
	dam[i]=rand()%18+1;
	/*printf("%d\t%d\n",sir[i],dam[i]);*/
	}
for (i=0;i<6;i++)
       {
	loc_eff[i]=gennor(0.0,1.0);
	/*printf("%d\t%d\n",sir[i],dam[i]);*/
	}
for(i=0; i<n;i++)
		 {
				 dum[i]=i;	 
			bpr[i]=	gennor(an,d);	
		}	
for(i=0; i<n;i++)
		 {
				 
			homo_loc[i]=gennor(0.0,1.0);	
			homo_loc2[i]=gennor(0.0,1.0);
			/*fprintf(sl,"%f\t%f\n",homo_loc[i],homo_loc2[i]);*/
		}	

for(i=0;i<base;i++)
				 {
				 shuffle(dum,n,t);
				insert(base_cros,dum,i);
                 t+=300;
			}

/*for(i=0;i<base;i++)
		
				  {  
					fprintf(fp,"\n\n");
					for(j=0; j<n;j++)
                    			{	
					fprintf(fp,"%d\t",base_cros[i*n+j]);			
					}
				}*/


for(i=0;i<base;i++)
		
				  {  
					for(j=0; j<n;j++)
                    			{
						if(base_cros[i*n+j]%2==0){base_loc[i*n+j]=bpr[j];base_zero[i*n+j]=1;}
						else{base_loc[i*n+j]=-1*bpr[j];base_zero[i*n+j]=0;}
					}
				}


for(i=0;i<=24;i++)
{
if(i<base){a[l]=0;l++;}
if(i>=base)
{
	for(f1=0;f1<10;f1++)
	{
	a[l]=1;l++;
		{
		for(f2=0;f2<10;f2++)
			{
			a[l]=2;l++;
				for(f3=0;f3<10;f3++)
				{
				a[l]=3;l++;
					for(f4=0;f4<10;f4++)
						{
						a[l]=4;l++;
						}
				}
			}
		}
	}
}
}
for(i=0;i<55570;i++)
{
 	
	if(a[i]==0)	
	{
	ped[pd][0]=a[i];ped[pd][1]=pd;ped[pd][2]=0;ped[pd][3]=0;ped[pd][4]=0.99;pd++;
	for(in=0;in<6;in++)
	{
	for(in_l=0;in_l<2;in_l++)
	{
	matrx[gg][0]=a[i];
	matrx[gg][1]=i;
	matrx[gg][2]=0;
	matrx[gg][3]=0;
	matrx[gg][4]=sum_base(base_loc,n,i);
	matrx[gg][5]=0.99;
	matrx[gg][6]=in+1;
	matrx[gg][7]=loc_eff[in];
	matrx[gg][8]=in_l+1;
	matrx[gg][9]=gennor(an,d);
	matrx[gg][10]=matrx[gg][4]+matrx[gg][7]+matrx[gg][9];
	gg++;
	}
	}
	}
	if(a[i]==1)
	{
	matrx[gg][0]=a[i];
	matrx[gg][1]=i;
	matrx[gg][2]=sir[si];
	matrx[gg][3]=dam[da];
	matrx[gg][4]=crossing(sir[si],dam[da],first_gene);
	matrx[gg][5]=0.0;
	matrx[gg][6]=1;
	matrx[gg][7]=loc_eff[0];
	matrx[gg][8]=1;
	matrx[gg][9]=gennor(an,d);
	matrx[gg][10]=matrx[gg][4]+matrx[gg][7]+matrx[gg][9];
	par=matrx[gg][1];
	
	gg++;
	generation++;
	if(generation>10)
		{
		si++;da++;
		generation=1;
		first_gene++;
		}
	if(second_value>10)
		{
		second_value=1;
		first_gene_second++;
		}
	second_value++;
	}
	if(a[i]==2)
	{
	matrx[gg][0]=a[i];
	matrx[gg][1]=i;
	matrx[gg][2]=par;
	matrx[gg][3]=par;
	matrx[gg][4]=second_crossing(second_gene,first_gene_second);
	matrx[gg][5]=0.5;
	matrx[gg][6]=2;
	matrx[gg][7]=loc_eff[1];
	matrx[gg][8]=1;
	matrx[gg][9]=gennor(an,d);
	matrx[gg][10]=matrx[gg][4]+matrx[gg][7]+matrx[gg][9];
	par2=matrx[gg][1];
	gg++;
	second_gene++;
	}
	if(a[i]==3)
	{
	matrx[gg][0]=a[i];
	matrx[gg][1]=i;
	matrx[gg][2]=par2;
	matrx[gg][3]=par2;
	matrx[gg][4]=third_crossing(third_gene,second_gene-1);
	/*printf("%d\t%d\n",third_gene,second_gene-1);*/
	matrx[gg][5]=0.75;
	matrx[gg][6]=3;
	matrx[gg][7]=loc_eff[2];
	matrx[gg][8]=1;
	matrx[gg][9]=gennor(an,d);
	matrx[gg][10]=matrx[gg][4]+matrx[gg][7]+matrx[gg][9];
	par3=matrx[gg][1];
	gg++;
	third_gene++;
	}
if(a[i]==4)
	{
	
	fourth_sum[fourth_gene][0]=fourth_crossing(fourth_gene,third_gene-1);
	fourth_sum[fourth_gene][1]=gennor(an,d);
	fourth_sum[fourth_gene][2]=fourth_sum[fourth_gene][0]+fourth_sum[fourth_gene][1]+loc_eff[3];
	track++;
	if(track==5)
		{
		track=0;
		matrx[gg][0]=a[i];
		matrx[gg][1]=i;
		matrx[gg][2]=par3;
		matrx[gg][3]=par3;
		matrx[gg][4]=calculate1(fourth_gene);
		matrx[gg][5]=0.87;
		matrx[gg][6]=4;
		matrx[gg][7]=loc_eff[3];
		matrx[gg][8]=1;
		matrx[gg][9]=calculate2(fourth_gene);
		matrx[gg][10]=calculate3(fourth_gene);
		gg++;
		}
	fourth_gene++;
	}
}
printf("%f\n%f\n%f\n%f\n%f",fourth_sum[49999][0],fourth_sum[49998][0],fourth_sum[49997][0],fourth_sum[49996][0],fourth_sum[49995][0]);
for(i=0;i<final;i++)
	{
if(matrx[i][0]==1){
			ped[pd][0]=matrx[i][0];ped[pd][1]=pd;ped[pd][2]=matrx[i][2];ped[pd][3]=matrx[i][3];ped[pd][4]=0.0;t2=ped[pd][1];pd++;f_2++;
			if(f_2==10){f_2=0;t2=ped[pd-1][1];}
		}	
if(matrx[i][0]==2){
		ped[pd][0]=matrx[i][0];ped[pd][1]=pd;ped[pd][2]=t2;ped[pd][3]=t2;ped[pd][4]=0.5;t3=ped[pd][1];pd++;f_3++;
			if(f_3==10){f_3=0;t3=ped[pd-1][1];}
		}
if(matrx[i][0]==3){ped[pd][0]=matrx[i][0];ped[pd][1]=pd;ped[pd][2]=t3;ped[pd][3]=t3;ped[pd][4]=0.75;pd++;t4=1;}
if(matrx[i][0]==4){if(t4==1)
			{ped[pd][0]=matrx[i][0];ped[pd][1]=pd;ped[pd][2]=ped[pd-1][1];ped[pd][3]=ped[pd-1][1];ped[pd][4]=0.87;pd++;t4=0;}
			else{ped[pd][0]=matrx[i][0];ped[pd][1]=pd;ped[pd][2]=ped[pd-2][1];ped[pd][3]=ped[pd-2][1];ped[pd][4]=0.87;pd++;}
			}		
}
pd=20;
fprintf(fp,"Gene\tLine\tgenotpic value\tIn Br.\tloc\tloc effect\trep\treplication\tphenotype\n");
for(i=0;i<final;i++)
	{
	if((int)matrx[i][0]==0){
	fprintf(fp,"%d\t%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\n",(int)matrx[i][0],(int)matrx[i][1],matrx[i][4],matrx[i][5],(int)matrx[i][6],matrx[i][7],(int)matrx[i][8],matrx[i][9],matrx[i][10]);
				}
else		{
fprintf(fp,"%d\t%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\n",(int)matrx[i][0],(int)ped[pd][1],matrx[i][4],matrx[i][5],(int)matrx[i][6],matrx[i][7],(int)matrx[i][8],matrx[i][9],matrx[i][10]);
	pd++;
		}
	}

/*for(i=0;i<15570;i++)
	{
	fprintf(fp,"%d\t%d\t%d\t%d\t%f\n",(int)ped[i][0],(int)ped[i][1],(int)ped[i][2],(int)ped[i][3],ped[i][4]);
	}
/*for(i=0;i<15790;i++)
	{
	fprintf(fp,"%d\t%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\n",(int)matrx[i][0],(int)matrx[i][1],matrx[i][4],matrx[i][5],(int)matrx[i][6],matrx[i][7],(int)matrx[i][8],matrx[i][9],matrx[i][10]);
	}
/*for (i=0;i<625;i++)
	{fprintf(fz,"\n");
	for (j=0;j<n;j++)
	{
	fprintf(fz,"%d",third_zero[i][j]);
	}
	}
for (i=0;i<3125;i++)
	{fprintf(sz,"\n");
	for (j=0;j<n;j++)
	{
	fprintf(sz,"%d",fourth_zero[i][j]);
	}
	}*/
/*for (i=0;i<3125;i++)
	{fprintf(fl,"\n");
	for (j=0;j<3;j++)
	{
	fprintf(fl,"%f\t",fourth_sum[i][j]);
	}
	}/*
for (i=0;i<3125;i++)
	{fprintf(sl,"\n");
	for (j=0;j<n;j++)
	{
	fprintf(sl,"%f\t",fourth_loc[i][j]);
	}
	}*/
return 0;
}

/*functions*/
void shuffle(int *dum,int n,int k)
{
 
 srand(k);
int i,j,t;
    if (n > 1) 
    {
       	for (i = 0; i < n - 1; i++)
         {
	  j = rand()%n;
	  t = dum[j];
	  dum[j] = dum[i];
	  dum[i] = t;
      	 }
    }
}
 void insert(int *base_cros,int *dum,int l)
      {
      	int k;
		for(k=0; k<n;k++)
				 {
					base_cros[l*n+k]=dum[k];

      	         }
      }

float sum_base(float *base_loc,int kk,int oo)

{
float hh=0.0;
int i;
for(i=0;i<kk;i++)
	{
	hh=hh+base_loc[oo*n+i];
	}
return hh;
}

float crossing(int SIR,int DAM,int first) 
       {
float result=0.0;
int count=0;
int j,h,kk;
h=SIR;
kk=DAM;		
		for(j=0; j<n;j++)
		
                    			{
						if(base_zero[h*n+j]==1) 
						{
						if(base_zero[kk*n+j]==0) 
						{first_zero[first*n+j]=2;count++;first_loc[first*n+j]=bpr[j];}
						if(base_zero[kk*n+j]==1) 
						{first_zero[first*n+j]=1;first_loc[first*n+j]=bpr[j];}
						}
						if(base_zero[h*n+j]==0) 
						{
						if(base_zero[kk*n+j]==0) 
						{first_zero[first*n+j]=0;first_loc[first*n+j]=-1*bpr[j];}
						if(base_zero[kk*n+j]==1) 
						{first_zero[first*n+j]=2;count++;first_loc[first*n+j]=bpr[j];}
						}
					}

for(j=0; j<n;j++)
	{
	result+=first_loc[first*n+j];
	}
return result;
}

float second_crossing(int second,int first) 
 {
float result=0.0;
int two=0,j,AA,Aa,aa,i,l,k,u,in=0;

for(j=0; j<n;j++)
		
                    			{
				if(first_zero[first*n+j]==2){two++;}
					}
homo=(int*) calloc(two,sizeof(int));

if(two%2==0)
	{
	Aa=two/2;
	if(Aa%2==0)
	{
	AA=aa=Aa/2;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
	if(Aa%2!=0)
	{
	u=Aa-1;
	AA=aa=Aa/2;
	AA=AA+1;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
	}

if(two%2!=0)
{
	u=two+1;
	Aa=u/2;
	AA=Aa/2;
	aa=AA-1;
	l=Aa+AA;
	k=AA+aa+Aa;
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
}
shuffle(homo,two,suff);
suff+=100;
		for(j=0; j<n;j++)
		
                    			{
						if(first_zero[first*n+j]==1) 
						{
						second_zero[second*n+j]=1; 
						second_loc[second*n+j]=first_loc[first*n+j];
						}
						if(first_zero[first*n+j]==0) 
						{
						second_loc[second*n+j]=0;
						second_loc[second*n+j]=first_loc[first*n+j];
						}
						if(first_zero[first*n+j]==2) 
						{
						second_zero[second*n+j]=homo[in];
						if(homo[in]==1)
						{					
						second_loc[second*n+j]=1*homo_loc[in];
						}
						if(homo[in]==2) 
						{second_loc[second*n+j]=1*homo_loc[in];
						}
						if(homo[in]==0) 
						{second_loc[second*n+j]=-1*homo_loc[in];
						}
						in++;
						}
					}
for(j=0;j<n;j++)
	{
	result+=second_loc[second*n+j];
	}
return result;
}


float third_crossing(int third,int first) 
 {
float result=0.0;
int two=0,j,AA,Aa,aa,i,l,k,u,in=0;

for(j=0; j<n;j++)
		
                    			{
				if(second_zero[first*n+j]==2){two++;}
					}
homo=(int*) calloc(two,sizeof(int));

if(two%2==0)
	{
	Aa=two/2;
	if(Aa%2==0)
	{
	AA=aa=Aa/2;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
if(Aa%2!=0)
	{
	u=Aa-1;
	AA=aa=Aa/2;
	AA=AA+1;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
	}

if(two%2!=0)
{
	u=two+1;
	Aa=u/2;
	AA=Aa/2;
	aa=AA-1;
	l=Aa+AA;
	k=AA+aa+Aa;
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
}
shuffle(homo,two,suff2);
suff2+=50;
		for(j=0; j<n;j++)
		
                    			{
						if(second_zero[first*n+j]==1) 
						{
						third_zero[third*n+j]=1; 
						third_loc[third*n+j]=second_loc[first*n+j];
						}
						if(second_zero[first*n+j]==0) 
						{
						third_zero[third*n+j]=0;
						third_loc[third*n+j]=second_loc[first*n+j];
						}
						if(second_zero[first*n+j]==2) 
						{
						third_zero[third*n+j]=homo[in];
						if(homo[in]==1)
						{					
						third_loc[third*n+j]=1*homo_loc2[in];
						}
						if(homo[in]==2) 
						{third_loc[third*n+j]=1*homo_loc2[in];
						}
						if(homo[in]==0) 
						{third_loc[third*n+j]=-1*homo_loc2[in];
						}
						in++;
						}
					}
for(j=0;j<n;j++)
	{
	result+=third_loc[third*n+j];
	}
return result;
}
float fourth_crossing(int third1,int first1) 
 {
float result=0.0;
int two=0,j,AA,Aa,aa,i,l,k,u,in=0;

for(j=0; j<n;j++)
		
                    			{
				if(third_zero[first1*n+j]==2){two++;}
					}

homo=(int*) calloc(two,sizeof(int));

if(two%2==0)
	{
	Aa=two/2;
	if(Aa%2==0)
	{
	AA=aa=Aa/2;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
if(Aa%2!=0)
	{
	u=Aa-1;
	AA=aa=Aa/2;
	AA=AA+1;
	l=Aa+AA;
	k=AA+aa+Aa;
	}
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
	}

if(two%2!=0)
{
	u=two+1;
	Aa=u/2;
	AA=Aa/2;
	aa=AA-1;
	l=Aa+AA;
	k=AA+aa+Aa;
for(i=0;i<Aa;i++){homo[i]=2;}
for(i=Aa;i<l;i++){homo[i]=1;}
for(i=l;i<two;i++){homo[i]=0;}
}
shuffle(homo,two,suff2);
suff2+=100;
		for(j=0; j<n;j++)
		
                    			{
						if(third_zero[first1*n+j]==1) 
						{
						fourth_zero[third1*n+j]=1; 
						fourth_loc[third1*n+j]=third_loc[first1*n+j];
						}
						if(third_zero[first1*n+j]==0) 
						{
						fourth_zero[third1*n+j]=0;
						fourth_loc[third1*n+j]=third_loc[first1*n+j];
						}
						if(third_zero[first1*n+j]==2) 
						{
						fourth_zero[third1*n+j]=homo[in];
						if(homo[in]==1)
						{					
						fourth_loc[third1*n+j]=1*homo_loc2[in];
						}
						if(homo[in]==2) 
						{fourth_loc[third1*n+j]=1*homo_loc2[in];
						}
						if(homo[in]==0) 
						{fourth_loc[third1*n+j]=-1*homo_loc2[in];
						}
						in++;
						}
					}
for(j=0;j<n;j++)
	{
	result+=fourth_loc[third1*n+j];
	}
return result;
}
float calculate1(int r)
{

float sum=0.0,res=0.0;
sum=fourth_sum[r][0]+fourth_sum[r-1][0]+fourth_sum[r-2][0]+fourth_sum[r-3][0]+fourth_sum[r-4][0];
res=sum/5;
return res;
}
float calculate2(int r)
{
float sum=0.0,res=0.0;
sum=fourth_sum[r][1]+fourth_sum[r-1][1]+fourth_sum[r-2][1]+fourth_sum[r-3][1]+fourth_sum[r-4][1];
res=sum/5;
return res;
}
float calculate3(int r)
{
float sum=0.0,res=0.0;
sum=fourth_sum[r][2]+fourth_sum[r-1][2]+fourth_sum[r-2][2]+fourth_sum[r-3][2]+fourth_sum[r-4][2];
res=sum/5;
return res;
}
