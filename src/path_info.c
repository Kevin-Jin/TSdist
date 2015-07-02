#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "tsdist.h"

SEXP ts_xLeadOverY(SEXP tamx, SEXP tamy, SEXP pathMatrixVec, SEXP xGapsVec, SEXP yGapsVec) {
	int i, j, k, tam1, tam2, pathLen, *pathMatrix, *xGaps, *yGaps;
	SEXP fullPathDimNames, fullPathRowNames, fullPathColNames, fullPath, xLeads, yLeads, result, resultNames;
	
	pathMatrix = INTEGER(pathMatrixVec);
	xGaps = INTEGER(xGapsVec);
	yGaps = INTEGER(yGapsVec);
	tam1 = *INTEGER(tamx) + 1;
	tam2 = *INTEGER(tamy) + 1;
	i = tam1 - 1;
	j = tam2 - 1;
	
	if ((tam2 - 1) - ELT(xGaps, tam2, i, j) != (tam1 - 1) - ELT(yGaps, tam2, i, j)) {
		error("reconstruct_path: gaps did not converge. two possible PATH_NO_GAP counts are %d and %d", (tam2 - 1) - ELT(xGaps, tam2, i, j) - 1, (tam1 - 1) - ELT(yGaps, tam2, i, j) - 1);
		return allocVector(VECSXP, 0);
	}
	pathLen = ((tam2 - 1) - ELT(xGaps, tam2, i, j))	//number of PATH_NO_GAPs on optimal path
		+ ELT(xGaps, tam2, i, j)					//number of PATH_X_GAPS on optimal path
		+ ELT(yGaps, tam2, i, j)					//number of PATH_Y_GAPS on optimal path
	;
	PROTECT(fullPathDimNames = allocVector(VECSXP, 2));
	PROTECT(fullPathRowNames = allocVector(STRSXP, 0));
	PROTECT(fullPathColNames = allocVector(STRSXP, 2));
	PROTECT(fullPath = allocMatrix(INTSXP, pathLen, 2));
	PROTECT(xLeads = allocVector(INTSXP, tam1));
	PROTECT(yLeads = allocVector(INTSXP, tam2));
	PROTECT(result = allocVector(VECSXP, 3));
	PROTECT(resultNames = allocVector(STRSXP, 3));
	for (k = tam1 - 1; k >= 0; --k)
		INTEGER(xLeads)[k] = NA_INTEGER;
	for (k = tam2 - 1; k >= 0; --k)
		INTEGER(yLeads)[k] = NA_INTEGER;
	
	INTEGER(xLeads)[i] = INTEGER(yLeads)[j] = 0;
	k = pathLen;
	while (i >= 0 && j >= 0 && k > 0) {
		k--;
		ELT(INTEGER(fullPath), pathLen, 0, k) = i;
		ELT(INTEGER(fullPath), pathLen, 1, k) = j;
		if (INTEGER(xLeads)[i] == NA_INTEGER || INTEGER(yLeads)[j] == NA_INTEGER) {
			error("reconstruct_path: assertion failed. xLeads[i] or yLeads[j] is NA");
			goto end;
		}
		switch (ELT(pathMatrix, tam2, i, j)) {
			case PATH_X_GAP: //stay still on X, move on on Y
				j--;
				INTEGER(xLeads)[i]++;
				INTEGER(yLeads)[j] = 0;
				break;
			case PATH_Y_GAP: //stay still on Y, move on on X
				i--;
				INTEGER(xLeads)[i] = 0;
				INTEGER(yLeads)[j]++;
				break;
			case PATH_NO_GAP: //move on on X and Y
				i--;
				j--;
				INTEGER(xLeads)[i] = INTEGER(yLeads)[j] = 0;
				break;
			default:
				error("reconstruct_path: bad path value. pathMatrix[%d, %d] == %d", i, j, ELT(pathMatrix, tam2, i, j));
				goto end;
		}
	}
	if (i != 0 || j != 0 || k != 0) {
		error("reconstruct_path: gaps did not converge. end %d short at pathMatrix[%d, %d]", k, i, j);
		goto end;
	}
	
end:
	SET_STRING_ELT(fullPathColNames, 0, mkChar("x"));
	SET_STRING_ELT(fullPathColNames, 1, mkChar("y"));
	SET_VECTOR_ELT(fullPathDimNames, 0, fullPathRowNames);
	SET_VECTOR_ELT(fullPathDimNames, 1, fullPathColNames);
	setAttrib(fullPath, R_DimNamesSymbol, fullPathDimNames);
	SET_VECTOR_ELT(result, 0, fullPath);
	SET_VECTOR_ELT(result, 1, xLeads);
	SET_VECTOR_ELT(result, 2, yLeads);
	SET_STRING_ELT(resultNames, 0, mkChar("full.path"));
	SET_STRING_ELT(resultNames, 1, mkChar("x.leads"));
	SET_STRING_ELT(resultNames, 2, mkChar("y.leads"));
	setAttrib(result, R_NamesSymbol, resultNames);
	UNPROTECT(8);
	return result;
}






