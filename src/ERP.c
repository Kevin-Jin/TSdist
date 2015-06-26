#include <R.h>
#include<math.h> 

#include "tsdist.h"
//costMatrix is offset by 1, so x[i - 1] and y[j - 1] are the current values and costMatrix[i - 1, j] is actually the prior cost
#define Y_GAP_ERP_COST(gapPenalty, x, costMatrix, lengthY, i, j) (fabs((gapPenalty) - (x)[i - 1]) + ELT(costMatrix, lengthY, i - 1, j))
#define X_GAP_ERP_COST(gapPenalty, y, costMatrix, lengthY, i, j) (fabs((gapPenalty) - (y)[j - 1]) + ELT(costMatrix, lengthY, i, j - 1))
#define NO_GAP_ERP_COST(distMatrix, costMatrix, lengthY, i, j) (ELT(distMatrix, lengthY - 1, i - 1, j - 1) + ELT(costMatrix, lengthY, i - 1, j - 1))
#define OPTIMIZE_ERP_COST(gapPenalty, x, y, distMatrix, costMatrix, xGaps, yGaps, lengthY, i, j, yGapCost, xGapCost, noGapCost) \
	yGapCost = Y_GAP_ERP_COST(gapPenalty, x, costMatrix, lengthY, i, j); \
	xGapCost = X_GAP_ERP_COST(gapPenalty, y, costMatrix, lengthY, i, j); \
	noGapCost = NO_GAP_ERP_COST(distMatrix, costMatrix, lengthY, i, j); \
	if (yGapCost <= xGapCost && yGapCost <= noGapCost) { \
		ELT(costMatrix, lengthY, i, j) = yGapCost; \
		ELT(pathMatrix, lengthY, i, j) = PATH_Y_GAP; \
		ELT(xGaps, tam2, i, j) = ELT(xGaps, tam2, i - 1, j); \
		ELT(yGaps, tam2, i, j) = ELT(yGaps, tam2, i - 1, j) + 1; \
	} else if (xGapCost <= yGapCost && xGapCost <= noGapCost) { \
		ELT(costMatrix, lengthY, i, j) = xGapCost; \
		ELT(pathMatrix, lengthY, i, j) = PATH_X_GAP; \
		ELT(xGaps, tam2, i, j) = ELT(xGaps, tam2, i, j - 1) + 1; \
		ELT(yGaps, tam2, i, j) = ELT(yGaps, tam2, i, j - 1); \
	} else if (noGapCost <= yGapCost && noGapCost <= xGapCost) { \
		ELT(costMatrix, lengthY, i, j) = noGapCost; \
		ELT(pathMatrix, lengthY, i, j) = PATH_NO_GAP; \
		ELT(xGaps, tam2, i, j) = ELT(xGaps, tam2, i - 1, j - 1); \
		ELT(yGaps, tam2, i, j) = ELT(yGaps, tam2, i - 1, j - 1); \
	}

//Function that calculates the cost matrix of the ERP distance.
void erp(double *x, double *y, int *tamx, int *tamy, int *sigma, double *costMatrix, double *distMatrix, int *pathMatrix, int *xGaps, int *yGaps, double *g){

	int i,j,tam1,tam2, siggma, max;
	tam1=*tamx+1;
	tam2=*tamy+1;
	siggma=*sigma+2;
	double yGapCost, xGapCost, noGapCost;
	double gapPenalty = *g;

	//The (0,0) position of the matrix is filled
	ELT(costMatrix, tam2, 0, 0) = 0.0;
	ELT(xGaps, tam2, 0, 0) = 0;
	ELT(yGaps, tam2, 0, 0) = 0;

	//The edges of the matrix are filled.
	for(i=1;i<siggma;i++){
		yGapCost = Y_GAP_ERP_COST(gapPenalty, x, costMatrix, tam2, i, 0);
		ELT(costMatrix, tam2, i, 0) = yGapCost;
		ELT(pathMatrix, tam2, i, 0) = PATH_Y_GAP;
		ELT(yGaps, tam2, i, 0) = ELT(yGaps, tam2, i - 1, 0) + 1;
	}

	for(j=1;j<siggma;j++){
		xGapCost = X_GAP_ERP_COST(gapPenalty, y, costMatrix, tam2, 0, j);
		ELT(costMatrix, tam2, 0, j) = xGapCost;
		ELT(pathMatrix, tam2, 0, j) = PATH_X_GAP;
		ELT(xGaps, tam2, 0, j) = ELT(xGaps, tam2, 0, j - 1) + 1;
	}

	//Fill all lines until i=sigma+2
	for(i=1;i<siggma;i++){
		max=MIN(i+siggma - 1,tam2);
		for(j=1;j<max;j++){
			OPTIMIZE_ERP_COST(gapPenalty, x, y, distMatrix, costMatrix, xGaps, yGaps, tam2, i, j, yGapCost, xGapCost, noGapCost);
		}
	}

	//Fill the rest of the matrix
	for(i=siggma;i<tam1;i++){
		max=MIN(i+siggma - 1,tam2);
		for(j=i-(siggma - 2);j<max;j++){
			OPTIMIZE_ERP_COST(gapPenalty, x, y, distMatrix, costMatrix, xGaps, yGaps, tam2, i, j, yGapCost, xGapCost, noGapCost);
		}
	}
}

//Function that calculates the cost matrix of the ERP distance, without temporal constraints.
void erpnw(double *x, double *y, int *tamx, int *tamy, double *costMatrix, double *distMatrix, int *pathMatrix, int *xGaps, int *yGaps, double *g){

	int i,j,tam1,tam2;
	tam1=*tamx+1;
	tam2=*tamy+1;
	double yGapCost, xGapCost, noGapCost;
	double gapPenalty = *g;
	

	//The (0,0) position of the matrix is filled
	ELT(costMatrix, tam2, 0, 0) = 0.0;
	ELT(xGaps, tam2, 0, 0) = 0;
	ELT(yGaps, tam2, 0, 0) = 0;

	//The edges of the matrix are filled.
	for(i=1;i<tam1;i++){
		yGapCost = Y_GAP_ERP_COST(gapPenalty, x, costMatrix, tam2, i, 0);
		ELT(costMatrix, tam2, i, 0) = yGapCost;
		ELT(pathMatrix, tam2, i, 0) = PATH_Y_GAP;
		ELT(yGaps, tam2, i, 0) = ELT(yGaps, tam2, i - 1, 0) + 1;
	}

	for(j=1;j<tam2;j++){
		xGapCost = X_GAP_ERP_COST(gapPenalty, y, costMatrix, tam2, 0, j);
		ELT(costMatrix, tam2, 0, j) = xGapCost;
		ELT(pathMatrix, tam2, 0, j) = PATH_X_GAP;
		ELT(xGaps, tam2, 0, j) = ELT(xGaps, tam2, 0, j - 1) + 1;
	}

	//The rest of the matrix is filled.
	for(i=1;i<tam1;i++){
		for(j=1;j<tam2;j++){
			OPTIMIZE_ERP_COST(gapPenalty, x, y, distMatrix, costMatrix, xGaps, yGaps, tam2, i, j, yGapCost, xGapCost, noGapCost);
		}
	}
}






