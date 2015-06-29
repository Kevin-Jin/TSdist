
#include <R.h>
#include<math.h> 


//Function that calculates the DTW cost matrix with a temporal condition.
void dtw(int *tamx, int *tamy, int *sigma, double *costMatrix, double *distMatrix, int *pathMatrix, int *xGaps, int *yGaps){

	int i,j,max,siggma;
	siggma=*sigma+1;

	//Fill position (0,0) of cost matrix.
	costMatrix[0]=distMatrix[0];

	//Fill line i=0.
	for(j=1;j<siggma;j++){
	costMatrix[j]=distMatrix[j]+costMatrix[j-1];}

	//Fill column j=0 until i=sigma+1
	for(i=1;i<siggma;i++){costMatrix[i*(*tamy)]=distMatrix[i*(*tamy)]+costMatrix[(i-1)*(*tamy)];}

	//Fill all lines until i=sigma+1
	for(i=1;i<siggma;i++){
	max=fmin(i+*sigma+1,*tamy);
		for(j=1;j<max;j++){
		  costMatrix[i*(*tamy)+j]=distMatrix[i*(*tamy)+j]+fmin(fmin(costMatrix[(i-1)*(*tamy)+j],costMatrix[i*(*tamy)+(j-1)]),costMatrix[(i-1)*(*tamy)+(j-1)]);}
	}
	
	//Fill the rest of the matrix
	for(i=(siggma);i<(*tamx);i++){
	  max=fmin(i+*sigma+1,*tamy);
		for(j=(i-*sigma);j<max;j++){
		  costMatrix[i*(*tamy)+j]=distMatrix[i*(*tamy)+j]+fmin(fmin(costMatrix[(i-1)*(*tamy)+j],costMatrix[i*(*tamy)+(j-1)]),costMatrix[(i-1)*(*tamy)+(j-1)]);}
	}

}

//Function that calculates the DTW cost matrix with no temporal condition.
void dtwnw(int *tamx, int *tamy, double *costMatrix, double *distMatrix, int *pathMatrix, int *xGaps, int *yGaps){

	int i,j;

	//Fill position (0,0) of cost matrix.
	costMatrix[0]=distMatrix[0];;

	//Fill the edges of the matrix.
	for(j=1;j<(*tamy);j++){
		costMatrix[j]=distMatrix[j]+costMatrix[j-1];
		}

	for(i=1;i<(*tamx);i++){
		costMatrix[i*(*tamy)]=distMatrix[i*(*tamy)]+costMatrix[(i-1)*(*tamy)];
		}


	//Fill the rest of the matrix.
	for(i=1;i<(*tamx);i++){
		for(j=1;j<(*tamy);j++){
		costMatrix[i*(*tamy)+j]=distMatrix[i*(*tamy)+j]+fmin(fmin(costMatrix[(i-1)*(*tamy)+j],costMatrix[i*(*tamy)+(j-1)]),costMatrix[(i-1)*(*tamy)+(j-1)]);
		}
	}
}







