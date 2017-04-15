#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Function declarations */
int sizeLinearSystem(FILE*);
void backSubstitution(float**, int, int*, float**);
int findNextPivot(float**, int, float*, int*, int);
int gaussianElimination(float**, int, int*, int, int);
int solveLinearSystem(float**, int, float**, int);

int main(int argc, char **argv) {
	/*
		Description:
			Requires a path to a file which contains a linear system of equations in matrix format, parses the matrix then attepts to solve the system using both Gaussian Elimation and Gaussian Elimation with Scaled Partial Pivoting. Displays the resulting solution to stdout.
		Inputs:
			argv[1] (string): the filepath/filename where the linear system of equations matrix is stored
		Outputs:
			Displays the resulting solutions for Gaussian Elimation and Scaled Partial Pivoting to stdout
		Error Codes:
			returns 0 upon success
			-1 : indicates error opening file
			-2 : indicates file was empty
			-3 : indicates error with subfunction which will be described in stdout
			-4 : indicates that no filepath was provided or too many inputs were provided
	*/

	/* Declare variables */
	int i, j, linsize;
	char _filepath[512];
	float temp;
	float **linSystemGaus, **linSystemSPP;
	float *rootsGaus, *rootsSPP;
	FILE *inFile;

	/* Get filename from command line */
	// Check # of inputs
	if(argc != 2) {																												
		printf("\tError not enough arguments: provide path to file containing linear system.\n");
		return -4;
	}
	// Copy filepath from cmd
	strcpy(_filepath, argv[1]);																									

	/* Open file */
	if((inFile = fopen(_filepath, "r")) == NULL) {
		printf("\tError opening %s. Exiting to system.\n", _filepath);
	}

	/* Determine size of linear system */
	if((linsize = sizeLinearSystem(inFile)) < 1) {
		printf("\tError file does not contain any information. Exiting to system.\n");
		return -2;
	}

	/* Allocate Memory */
	linSystemGaus = (float**) malloc((linsize-1) * sizeof(float*));
	linSystemSPP = (float**) malloc((linsize-1) * sizeof(float*));
	for(i = 0; i < linsize-1; i++) {
		linSystemGaus[i] = (float*) malloc(linsize * sizeof(float));
		linSystemSPP[i] = (float*) malloc(linsize * sizeof(float));
	}

	/* Parse Linear System */
	for(i = 0; i < linsize-1; i++) {
		for(j = 0; j < linsize; j++) {
			fscanf(inFile, "%f", &temp);
			linSystemGaus[i][j] = temp;
			linSystemSPP[i][j] = temp;
		}
	}

	/* Solve Linear System */
	solveLinearSystem(linSystemGaus, linsize, &rootsGaus, 0);
	solveLinearSystem(linSystemSPP, linsize, &rootsSPP, 1);

	/* Display Results */
	printf("\tGaussian Elimation:      ");
	for(i = 0; i < linsize-1; i++) {
		printf("x%d = %f\t", i, rootsGaus[i]);
	}
	printf("\n");

	printf("\tScaled Partial Pivoting: ");
	for(i = 0; i < linsize-1; i++) {
		printf("x%d = %f\t", i, rootsSPP[i]);
	}
	printf("\n");

	/* Clean-up and Return */
	for(i = 0; i < linsize-1; i++) {
		free(linSystemGaus[i]);
		free(linSystemSPP[i]);
	}
	free(linSystemGaus);
	free(linSystemSPP);
	free(rootsGaus);
	free(rootsSPP);
	return 0;
}

int solveLinearSystem(float **_linearSystem, int _linearsize, float **_solutions, int _spp) {
 	/*
 		Description: Will solve the system of linear equations specified using either scaled-partial-pivoting or simple gaussian elimination
 		Inputs: 
			_linearSystem: the 2D array containing the linear system of equations to be solved
			_linearsize: the number of equations contained in the linearsystem
			_solutions: a 2D array pointer to the solutions of the system
			_spp: if 1 the scaled-partial-pivoting is used if 0 basic gaussian elimination is used
 		Outputs:
			Solutions to the system are placed in _solutions
 	*/

 	/* Declare variables */
	int i, j, curr_pivot, indexcurr_pivot;
	int *pivotList, *pivotIndices;
	float *scalingFactors;

 	/* Allocate Memory */
	scalingFactors = (float*) malloc((_linearsize-1) * sizeof(float));
	pivotList = (int*) malloc((_linearsize-1) * sizeof(int));
	pivotIndices = (int*) malloc((_linearsize-1) * sizeof(int));

	/* Set pivots and scaling factors */
	for(i = 0; i < _linearsize-1; i++) {
		// set initial pivots in order
		pivotList[i] = i;																									
		scalingFactors[i] = _linearSystem[i][0];		
		for(j = 1; j < _linearsize; j++) {
			// find the largest factor in each row
			scalingFactors[i] = scalingFactors[i] > _linearSystem[i][j] ? scalingFactors[i] : _linearSystem[i][j];
		}	
	}	

 	/* Perform elimation */
	for(i = 0; i < _linearsize-1; i++) {
		if(_spp) indexcurr_pivot = findNextPivot(_linearSystem, _linearsize, scalingFactors, pivotList, i);
		else indexcurr_pivot = i;
		gaussianElimination(_linearSystem, _linearsize, pivotList, i, indexcurr_pivot);
		
	}

	/* Determine the indices of the pivots in ascending order */
	for(i = 0; i < _linearsize-1; i++) {
		pivotIndices[pivotList[i]] = i;
	}

 	/* Find solutions */
 	backSubstitution(_linearSystem, _linearsize, pivotIndices, _solutions);

 	/* Clean-up and return */
	free(scalingFactors);
	free(pivotList);
	free(pivotIndices);
 	return 0;
}

int findNextPivot(float **_linSys, int _linsize, float *_scaleFactors, int *_pivotOrder, int _currpivot) {
	/* 
		Description: Will find the index of the next pivot equation
		Inputs:
			_linSys: the 2D matrix of the linearsystem
			_linsize: the number of equations in the linearsystem
			_scaleFactors: the scaling factors of the array
			_pivotOrder: the current ordering of the matrix
			_currPivot: the current pivot number of the matrix
		Outputs: returns the index of the next pivot equation
	*/
	
	/* Declare variables */
	int i, indexbestratio, indexcurrpivot;
	float bestratio=-1;

	/* Determine curr_pivot */
	for(i = 0; i < _linsize-1; i++) {
		if(_pivotOrder[i] >= _currpivot) {
			if(bestratio < fabsf(_linSys[i][_currpivot]) / _scaleFactors[i]) {
				indexbestratio = i;
				bestratio = fabsf(_linSys[i][_currpivot]) / _scaleFactors[i];
			}
			if(_pivotOrder[i] == _currpivot) {
				indexcurrpivot = i;
			}
		}
	}

	/* Swap current pivot to location of greatest ratio */
	_pivotOrder[indexcurrpivot] = _pivotOrder[indexbestratio];
	_pivotOrder[indexbestratio] = _currpivot;

	/* Clean-up and return */
	return indexbestratio;
}


int gaussianElimination(float **_linSys, int _linsize, int *_pivotOrder, int _currpivot, int _indexcurrpivot) {
	/*
		Description:
			Takes a matrix which describes a linear system, a set of scaling factors, a set of pivotOrders and a current pivot which defines the column/row (currcol=currrow=_currpivot-1)at which to evaluate and then uses eliminates all values at that column with row indices greater than the current row. The algorithm can be described through the following steps:

			1. Determine the largest ratio of coefficent/scalefactor for all _pivotOrders greater than or equal to currpivot
			2. The row which has the largest ratio is swapped with the row which is currently the _currpivot and becomes the _currpivot
			3. For every row with _pivotOrder greater than _currpivot eliminate the values in the column currpivot-1

		Inputs:
			_linSys (float**): a n*(n+1) matrix which describes a linear system
			_linsize (int): the n+1 value of the matrix
			_scaleFactors (float*): the scaling ratios for each row
			_pivotOrder (int*): a list containing the order of the rows
			_currpivot (int): the current row value which is being evaluated

		Outputs:
			Modifies _linSys and _pivotOrder to reflect the new state of the matrix as it approaches row-eschelon form.
		Error Codes:
			-1: Function exits to system if any of the pointers passed are NULL
			-2: Function exits to system if _currpivot or _linsize are < 1
	*/

	/* Declare variables */
	int i, j, indexcurrpivot;
	float factorDiv;

	/* Check inputs */
	if(_linSys == NULL) {
		printf("\tError in gaussianElimation. _linSys was NULL. Exiting to system.\n");
		exit(-1);
	}
	if(_pivotOrder == NULL) {
		printf("\tError in gaussianElimation. _pivotOrder was NULL. Exiting to system.\n");
		exit(-1);
	}
	if(_currpivot < 0) {
		printf("\tError in gaussianElimation. _currpivot < 1. Exiting to system.\n");
		exit(-2);
	}
	if(_linsize < 1) {
		printf("\tError in gaussianElimation. _linsize < 1. Exiting to system.\n");
		exit(-2);
	}


	/* Cancel out pivots greater than curr_pivot at currcolumn */
	for(i = 0; i < _linsize; i++) {
		if(_pivotOrder[i] > _currpivot) {
			factorDiv = _linSys[i][_currpivot] / _linSys[_indexcurrpivot][_currpivot];
			_linSys[i][_currpivot] = 0;
			for(j = _currpivot; j < _linsize; j++) {
				_linSys[i][j] -= factorDiv * _linSys[_indexcurrpivot][j];
			}
		}
	}

	/* Clean-up and return */
	return 0;
}

void backSubstitution(float **_rowEschelon, int _sizesys, int *_pivotIndices, float **_sols) {
	/*
		Description:
			Finds the solution to a system of linear equations which is in the row eschelon form. The algorithm works as follows:
				1. i = n; i -> 0 start at final row and move towards first, solve for x_n first since row
						eschelon form ensures that the solution is trival.
				2. find the right hand side of equation by looping from j=i, j->n subtracting coefficients*x_j
						from the known equate where we have already solved for x_j in previous iterations. Basically 6x_1+5x_2+3x_3=2 becomes 	
						6x_1=2-5x_2-3x_3
				3. divide resulting right hand side of equation by coefficient of current x_j; i.e.
						6x_1=2-5x_2-3x_3 -> x_1=(2-5x_2-3x_3)/6 store this value as the solution for x_j.
				4. repeat steps 2-4 until i = 0, now all the solutions are known.

		Inputs:
			_rowEschelon (float**): a n*n+1 matrix which describes a linear system in row-eschelon form
			_sizesys (int): the n+1 value or the number of columns in the linear system
			_pivotIndices (int*): the indices at which the pivots may be found shall be ordered by increasing pivot order
			_sols (float*): a pointer to where the solutions of the system will be stored. This will be
											allocated by the program to contain _sizesys-1 elements
		Outputs:
			Places the solutions to the linear system at the location of _sols.
		Errors:
			if _rowEschelon is a nullpointer then program will exit(-1) and display an error
			if _sizesys is not a positive integer greater than 1 then program will exit(-2)
	*/

	/* Declare variables */
	int i, j, pivindex;
	float rheq;

	/* Check inputs */
	if(_rowEschelon == NULL) {
		printf("\tError in backSubstitution(), _rowEschelon was NULL. Exiting to System.\n");
		exit(-1);
	}
	if(_sizesys < 1) {
		printf("\tError in backSubstitution(), _sizesys was < 1. Exiting to System.\n");
		exit(-2);
	}

	/* Allocate memory */
	*_sols = (float*) malloc((_sizesys-1) * sizeof(float));

	for(i = _sizesys-2; i >= 0; i--) {
		// get index of pivot i
		pivindex = _pivotIndices[i];
		// set right hand side of equation equal to the constants
		rheq = _rowEschelon[pivindex][_sizesys-1];
		for(j=i+1; j < _sizesys-1; j++) {
			// subtract values of known variables from right hand side of equation
			rheq -= _rowEschelon[pivindex][j]*_sols[0][j];
		}
		// divide right hand side of equation by coefficient of solution i
		_sols[0][i] = rheq/_rowEschelon[pivindex][i];
	}

	// Clean-up and return
	return;
}

int sizeLinearSystem(FILE *_fp) {
	/*
		Description:
			Reads in the first line of the file as a string then parses floats from it until the end of the string is reached. This will determine how many columns the matrix has, by the property of the linear system being fully defined where there are n equations and n unknowns there must be columns-1 rows in the matrix as well.
		Inputs:
			_fp (FILE): the file from which the linear system will be read the number of characters per row shall not exceed 512 or the program will terminate with error.
		Outputs:
			Returns the number of columns in the matrix.
		Error Codes:
			-1 if the *_fp is NULL error stating so will be printed to stdout.
			-2 if the file is empty and therefore describes no linear system. No comment on this is displayed it is left to the calling function to determine how to handle this case.
			-3 if the number of characters contained per row in the file exceeds 512.
	*/

	/* Declare variables */
	char *buffRow;
	size_t buffsize=512;
	int count=0, bytesread=0, bytestot=0;
	float temp;

	/* Allocate memory */
	buffRow = (char*) malloc(buffsize * sizeof(char));

	/* Determine number of columns */
	getline(&buffRow, &buffsize, _fp);																					
	// Read first line of file
	while(sscanf(buffRow+bytestot, "%f%n", &temp, &bytesread) != EOF && count < 10) {
		// Count number of floats
		count++;	
		// Count number of bytes read 																															
		bytestot += bytesread;
	}

	/* Clean-up and return */
	// Return to start of file
	rewind(_fp);																																
	return count;

}
