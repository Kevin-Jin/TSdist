#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "tsdist.h"

SEXP ts_xLeadOverY(SEXP tamx, SEXP tamy, SEXP pathMatrixVec, SEXP xGapsVec, SEXP yGapsVec) {
	int i, j, k, tam1, tam2, pathLen, *pathMatrix, *xGaps, *yGaps;
	SEXP xLeadOverY, fullPath, xLeads, yLeads, result, resultNames;
	
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
	PROTECT(xLeadOverY = allocVector(INTSXP, pathLen));
	PROTECT(fullPath = allocMatrix(INTSXP, 2, pathLen));
	PROTECT(xLeads = allocVector(INTSXP, tam1));
	PROTECT(yLeads = allocVector(INTSXP, tam2));
	PROTECT(result = allocVector(VECSXP, 4));
	PROTECT(resultNames = allocVector(STRSXP, 4));
	for (k = pathLen - 1; k >= 0; --k)
		INTEGER(xLeadOverY)[k] = NA_INTEGER;
	for (k = tam1 - 1; k >= 0; --k)
		INTEGER(xLeads)[k] = NA_INTEGER;
	for (k = tam2 - 1; k >= 0; --k)
		INTEGER(yLeads)[k] = NA_INTEGER;
	
	INTEGER(xLeads)[i] = INTEGER(yLeads)[j] = 0;
	while (i >= 0 && j >= 0 && pathLen > 0) {
		pathLen--;
		INTEGER(xLeadOverY)[pathLen] = ELT(xGaps, tam2, i, j) - ELT(yGaps, tam2, i, j);
		ELT(INTEGER(fullPath), 2, pathLen, 0) = i;
		ELT(INTEGER(fullPath), 2, pathLen, 1) = j;
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
	if (i != 0 || j != 0 || pathLen != 0) {
		error("reconstruct_path: gaps did not converge. end %d short at pathMatrix[%d, %d]", pathLen, i, j);
		goto end;
	}
	
end:
	SET_VECTOR_ELT(result, 0, xLeadOverY);
	SET_VECTOR_ELT(result, 1, fullPath);
	SET_VECTOR_ELT(result, 2, xLeads);
	SET_VECTOR_ELT(result, 3, yLeads);
	SET_STRING_ELT(resultNames, 0, mkChar("x.cumulative.lead"));
	SET_STRING_ELT(resultNames, 1, mkChar("full.path"));
	SET_STRING_ELT(resultNames, 2, mkChar("x.leads"));
	SET_STRING_ELT(resultNames, 3, mkChar("y.leads"));
	setAttrib(result, R_NamesSymbol, resultNames);
	UNPROTECT(6);
	return result;
}






