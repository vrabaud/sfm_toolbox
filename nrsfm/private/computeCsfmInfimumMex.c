#include "mex.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*
mkoctfile --mex -v  -DCXX=g++-4.1 -DCC=g++-4.1 -DLD_CXX=g++-4.1 -o computeCsfmInfimum.mex computeCsfmInfimum.c;
load('../../W.mat');
[ Sim2, q ] = computeCsfmInfimum(W);

load('../data/shark/jawSource.mat'); tic; [ Sim2, q ] = computeCsfmInfimum(W); toc


 mex computeAnimSimilarity.c -lmwlapack -lmwblas -DCXX=g++-4.2 -DCC=g++-4.2 -DLD=g++-4.2

 load('shark.mat'); Sim2 = computeAnimSimilarity(reshape(permute(animGT.W,[3,1,2]),size(W,3),[]));*/

#ifdef MATLAB_MEX_FILE
#include "blas.h"
#else
double ddot_(int *, double *, int *, double *, int *);
double dgemm_(char*, char*, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
double dnrm2_(int *, double *, int *);
#endif

const double NUM_TOL = 1e-8;
const double NUM_TOL_SQ = 1e-8*1e-8;

void printMatrix(double *array, int m, int n) {
	int i, j;

	for (j = 0; j < m; ++j) {
		for (i = 0; i < n; ++i)
			fprintf(stderr,"%f ", array[j + i * m]);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
}

inline void normalizeQuaternion(double *quaternion) {
	/* normalize */
	double norm = sqrt(quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);
	int k;
	
	for( k = 0; k < 4; ++k )
		quaternion[k] /= norm;
}

inline void cleanQuaternion(double *quaternion) {
	if ((quaternion[0] * quaternion[0] + quaternion[3] * quaternion[3]) < NUM_TOL_SQ) {
		quaternion[0] = 1e-5;
		quaternion[3] = -1e-5;
	}
	if ((quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2]) < NUM_TOL_SQ) {
		quaternion[1] = 1e-5;
		quaternion[2] = -1e-5;
	}
	/* normalize */
	normalizeQuaternion(quaternion);
}

inline void copyQuaternion(double *quaternionIn, double *quaternionOut, char doClean) {
	int k;
	
	for (k = 0; k < 4; ++k)
		quaternionOut[k] = quaternionIn[k];

	if (doClean)
		cleanQuaternion(quaternionOut);
}

inline double normSq(int n, double *A) {
	int one = 1;

	return ddot_(&n, A, &one, A, &one);

	/* 	double res = dnrm2_(&n, A, &one);
	return res*res; */
}

/* Compute A*B', A is of size m *p and B n * p */
inline double matrixMultiply(double *A, double *B, double *C, int m, int p, int n) {
	double one = 1.0, zero = 0.0;
	char *chn = "N";
	char *chnb = "C";

	dgemm_(chn, chnb, &m, &n, &p, &one, A, &m, B, &n, &zero, C, &m);
}

/* Compute trace(A*B'), A and B are of size m*n */
inline double matrixTraceProduct4By4(double *A, double *B) {
	int elemNbr = 4 * 4, one = 1;

	return ddot_(&elemNbr, A, &one, B, &one);
}

inline void _copywPrimeWPrimeTransposeBlock(double *wPrimeWPrimeTranspose, int jBlock, int iBlock, int iFrame, int jFrame, double *wPairwiseProd, int nFrame) {
	int i, j;
	for (i = 0; i < 2; ++i)
		for(j = 0; j < 2; ++j)
			wPrimeWPrimeTranspose[2*4*iBlock + 2*jBlock + 4*i + j] = wPairwiseProd[(2*iFrame+i)*2*nFrame + 2*jFrame + j];
}

inline void updateW(double *wPairwiseProd, int iFrame, int jFrame, int nFrame, double *wPrimeWPrimeTranspose) {
	int i, j;
	/* update wPrimeWPrimeTranspose */
	_copywPrimeWPrimeTransposeBlock(wPrimeWPrimeTranspose, 0, 0, iFrame, iFrame, wPairwiseProd, nFrame);
	_copywPrimeWPrimeTransposeBlock(wPrimeWPrimeTranspose, 0, 1, iFrame, jFrame, wPairwiseProd, nFrame);
	_copywPrimeWPrimeTransposeBlock(wPrimeWPrimeTranspose, 1, 0, jFrame, iFrame, wPairwiseProd, nFrame);
	_copywPrimeWPrimeTransposeBlock(wPrimeWPrimeTranspose, 1, 1, jFrame, jFrame, wPairwiseProd, nFrame);
}

inline double reconstrPairVal(double *quaternion, double *wPrimeWPrimeTranspose, double *RTransposeR) {
	double a = quaternion[0], b = quaternion[1], c = quaternion[2], d =
			quaternion[3];

	double  t19, t18, t17, t16, t9, t15, t13, t10, t14, t8, t7, t6, t5, t4, t3, t2, t1, t;
	t19 = b*d;                                                                      
	t18 = a*b;                                                                         
	t17 = c*d;                                                                         
	t16 = 1/(d*d+a*a)/(c*c+b*b);                                                       
	t9 = t18+t17;                                                                      
	t15 = t9*t16;                                                                      
	t13 = a*c;                                                                         
	t10 = t13+t19;                                                                     
	t14 = t10*t16;
	t8 = -t13+t19;
	t7 = -t17+t18;
	t6 = t7*t14;
	t5 = t8*t15;
	t4 = t7*t15;
	t3 = t8*t14;
	t2 = t9*t14;
	t1 = t8*t7*t16;
	RTransposeR[0] = t7*t7*t16;
	RTransposeR[1] = t6;
	RTransposeR[2] = -t4;
	RTransposeR[3] = t1;
	RTransposeR[4] = t6;
	RTransposeR[5] = t10*t10*t16;
	RTransposeR[6] = -t2;
	RTransposeR[7] = t3;
	RTransposeR[8] = -t4;
	RTransposeR[9] = -t2;
	RTransposeR[10] = t9*t9*t16;
	RTransposeR[11] = -t5;
	RTransposeR[12] = t1;
	RTransposeR[13] = t3;
	RTransposeR[14] = -t5;
	RTransposeR[15] = t8*t8*t16;
	
	return matrixTraceProduct4By4(RTransposeR, wPrimeWPrimeTranspose);
}

inline double reconstrPairValGrad(double *quaternion, double *wPrimeWPrimeTranspose, double *gradient, double *RTransposeR, double *dRTransposeR) {
	double a = quaternion[0], b = quaternion[1], c = quaternion[2], d =
			quaternion[3];
	int i, j, k;

	double  t86, t85, t25, t28, t21, t19, t26, t27, t22, t20, t84, t17, t83, t24, t16, t82, t23, t13, t81, t18, t80, t14, t79, t11, t78, t15, t7, t77, t12, t9, t76, t8, t75, t74, t73, t72, t71, t70, t10, t69, t68, t67, t66, t65, t64, t63, t62, t61, t60, t59, t58, t57, t56, t55, t54, t53, t52, t51, t50, t49, t48, t47, t46, t45, t44, t43, t42, t41, t40, t39, t38, t37, t36, t35, t34, t33, t32, t31, t30, t29, t6, t5, t4, t3, t2, t1, t;     
	t86 = c*d;                                                                                             
	t85 = b*d;                                                                                                
	t25 = a*a;                                                                                                
	t28 = d*d;                                                                                                
	t21 = t28+t25;                                                                                            
	t19 = 1/t21;                                                                                              
	t26 = b*b;                                                                                                
	t27 = c*c;                                                                                                
	t22 = t27+t26;                                                                                            
	t20 = 1/t22;                                                                                              
	t84 = t19*t20;                                                                                            
	t17 = 1/(t22*t22);                                                                                        
	t83 = t17*t19;                                                                                            
	t24 = a*b;                                                                                                
	t16 = -t86+t24;                                                                                           
	t82 = t16*t19;                                                                                            
	t23 = a*c;                                                                                                
	t13 = t23-t85;                                                                                            
	t81 = t20*t13;                                                                                            
	t18 = 1/(t21*t21);                                                                                        
	t80 = t18*t20;                                                                                            
	t14 = t25*t26-t27*t28;                                                                                    
	t79 = t20*t14;                                                                                            
	t11 = t25*t27-t26*t28;                                                                                    
	t78 = t11*t20;                                                                                            
	t15 = t23+t85;                                                                                            
	t7 = t15*t15;                                                                                             
	t77 = t7*t83;                                                                                             
	t12 = t24+t86;                                                                                            
	t9 = t12*t12;                                                                                             
	t76 = t9*t83;                                                                                             
	t8 = t16*t16;                                                                                             
	t75 = t8*t83;                                                                                             
	t74 = t7*t80;                                                                                             
	t73 = t9*t80;                                                                                             
	t72 = t8*t80;                                                                                             
	t71 = t18*t78;                                                                                            
	t70 = t11*t83;                                                                                            
	t10 = t13*t13;                                                                                            
	t69 = t10*t83;                                                                                            
	t68 = t12*t83;                                                                                            
	t67 = t14*t83;                                                                                            
	t66 = t18*t79;                                                                                            
	t65 = t10*t80;                                                                                            
	t64 = t17*t82;                                                                                            
	t63 = t19*t79;                                                                                            
	t62 = t15*t84;                                                                                            
	t61 = t19*t78;                                                                                            
	t60 = t12*t81;                                                                                            
	t59 = t16*t80;                                                                                            
	t58 = c*t70;                                                                                              
	t57 = b*t67;                                                                                              
	t56 = d*t71;                                                                                              
	t55 = d*t66;                                                                                              
	t54 = a*t71;                                                                                              
	t53 = a*t66;                                                                                              
	t52 = t15*t68;                                                                                            
	t51 = t15*t64;                                                                                            
	t50 = t13*t64;                                                                                            
	t49 = t13*t68;                                                                                            
	t48 = t18*t60;                                                                                            
	t47 = t12*t15*t80;                                                                                        
	t46 = t13*t59;                                                                                            
	t45 = t15*t59;                                                                                            
	t44 = t81*t82;                                                                                            
	t43 = t12*t62;                                                                                            
	t42 = c*t51;                                                                                              
	t41 = c*t52;                                                                                              
	t40 = c*t50;                                                                                              
	t39 = c*t49;                                                                                              
	t38 = b*t51;                                                                                              
	t37 = b*t52;                                                                                              
	t36 = b*t50;                                                                                              
	t35 = b*t49;                                                                                              
	t34 = d*t45;                                                                                              
	t33 = d*t47;                                                                                              
	t32 = d*t48;                                                                                              
	t31 = a*t45;                                                                                              
	t30 = a*t46;                                                                                              
	t29 = a*t48;                                                                                              
	t6 = c*t67;                                                                                               
	t5 = b*t70;                                                                                               
	t4 = t16*t62;                                                                                             
	t3 = t19*t60;                                                                                             
	t2 = d*t46;                                                                                               
	t1 = a*t47;                                                                                               
	RTransposeR[0] = t8*t84;                                                                                           
	RTransposeR[1] = t4;                                                                                               
	RTransposeR[2] = -t63;                                                                                             
	RTransposeR[3] = -t44;                                                                                             
	RTransposeR[4] = t4;                                                                                               
	RTransposeR[5] = t7*t84;                                                                                           
	RTransposeR[6] = -t43;                                                                                             
	RTransposeR[7] = -t61;                                                                                             
	RTransposeR[8] = -t63;                                                                                             
	RTransposeR[9] = -t43;                                                                                             
	RTransposeR[10] = t9*t84;                                                                                          
	RTransposeR[11] = t3;                                                                                              
	RTransposeR[12] = -t44;                                                                                            
	RTransposeR[13] = -t61;                                                                                            
	RTransposeR[14] = t3;                                                                                              
	RTransposeR[15] = t10*t84;                                                                                         
	dRTransposeR[0] = t34;                                                                                             
	dRTransposeR[1] = d*t74;                                                                                           
	dRTransposeR[2] = -t33;                                                                                            
	dRTransposeR[3] = -t56;                                                                                            
	dRTransposeR[4] = -d*t72;                                                                                          
	dRTransposeR[5] = -t34;                                                                                            
	dRTransposeR[6] = t55;                                                                                             
	dRTransposeR[7] = t2;                                                                                              
	dRTransposeR[8] = t2;                                                                                              
	dRTransposeR[9] = t56;                                                                                             
	dRTransposeR[10] = -t32;                                                                                           
	dRTransposeR[11] = -d*t65;                                                                                         
	dRTransposeR[12] = -t55;                                                                                           
	dRTransposeR[13] = -t33;                                                                                           
	dRTransposeR[14] = d*t73;                                                                                          
	dRTransposeR[15] = t32;                                                                                            
	dRTransposeR[16] = t42;                                                                                            
	dRTransposeR[17] = c*t77;                                                                                          
	dRTransposeR[18] = -t41;                                                                                           
	dRTransposeR[19] = -t58;                                                                                           
	dRTransposeR[20] = -c*t75;                                                                                         
	dRTransposeR[21] = -t42;                                                                                           
	dRTransposeR[22] = t6;                                                                                             
	dRTransposeR[23] = t40;                                                                                            
	dRTransposeR[24] = -t40;                                                                                           
	dRTransposeR[25] = -t58;                                                                                           
	dRTransposeR[26] = t39;                                                                                            
	dRTransposeR[27] = c*t69;                                                                                          
	dRTransposeR[28] = t6;                                                                                             
	dRTransposeR[29] = t41;                                                                                            
	dRTransposeR[30] = -c*t76;                                                                                         
	dRTransposeR[31] = -t39;                                                                                           
	dRTransposeR[32] = -t38;                                                                                           
	dRTransposeR[33] = -b*t77;                                                                                         
	dRTransposeR[34] = t37;                                                                                            
	dRTransposeR[35] = t5;                                                                                             
	dRTransposeR[36] = b*t75;                                                                                          
	dRTransposeR[37] = t38;                                                                                            
	dRTransposeR[38] = -t57;                                                                                           
	dRTransposeR[39] = -t36;                                                                                           
	dRTransposeR[40] = t36;                                                                                            
	dRTransposeR[41] = t5;                                                                                             
	dRTransposeR[42] = -t35;                                                                                           
	dRTransposeR[43] = -b*t69;                                                                                         
	dRTransposeR[44] = -t57;                                                                                           
	dRTransposeR[45] = -t37;                                                                                           
	dRTransposeR[46] = b*t76;                                                                                          
	dRTransposeR[47] = t35;                                                                                            
	dRTransposeR[48] = -t31;                                                                                           
	dRTransposeR[49] = -a*t74;                                                                                         
	dRTransposeR[50] = t1;                                                                                             
	dRTransposeR[51] = t54;                                                                                            
	dRTransposeR[52] = a*t72;                                                                                          
	dRTransposeR[53] = t31;                                                                                            
	dRTransposeR[54] = -t53;                                                                                           
	dRTransposeR[55] = -t30;                                                                                           
	dRTransposeR[56] = -t30;                                                                                           
	dRTransposeR[57] = -t54;                                                                                           
	dRTransposeR[58] = t29;                                                                                            
	dRTransposeR[59] = a*t65;                                                                                          
	dRTransposeR[60] = t53;                                                                                            
	dRTransposeR[61] = t1;                                                                                             
	dRTransposeR[62] = -a*t73;                                                                                         
	dRTransposeR[63] = -t29;
	
	/* Compute the gradient */
	for (k = 0; k < 4; ++k)
		gradient[k] = 2*matrixTraceProduct4By4(dRTransposeR + 4 * 4 * k, wPrimeWPrimeTranspose);

	/* Return the error */
	return matrixTraceProduct4By4(RTransposeR, wPrimeWPrimeTranspose);
}

void gradientDescent(double *wPrimeWPrimeTranspose, int iFrame, int jFrame, double * wPairwiseProd, double *errArray, double *quaternionArray, char *isDone, int nFrame, int nPoint, double *RTransposeR, double *dRTransposeR, double *quaternionTmp, double *gradient) {
	int i, ii, j, k;

	/* update the wPrimeWPrimeTranspose matrix */
	updateW(wPairwiseProd, iFrame, jFrame, nFrame, wPrimeWPrimeTranspose);

	/* Figure out the best guess to start from */
	double errMin, errTmp;

	double *quaternionMin = quaternionArray + 4 * (iFrame * nFrame + jFrame);

	/* Go over all the neighbors around the current location to find the best one */
	if (isDone[jFrame * nFrame + iFrame])
		errMin = errArray[jFrame * nFrame + iFrame];
	else
		errMin = -1;
	for (i = iFrame - 1; i <= iFrame + 1; ++i)
		for (j = jFrame - 1; j <= jFrame + 1; ++j)
			if ((i >= 0) && (i < nFrame) && (j >= 0) && (j < nFrame)
					&& isDone[j * nFrame + i]) {
				/* and check which one already gives the smallest value */
				copyQuaternion(quaternionArray + 4 * (i * nFrame + j), quaternionTmp, 1);
			
				errTmp = reconstrPairVal(quaternionTmp, wPrimeWPrimeTranspose,RTransposeR);

				/* Keep the best pick */
				if ((errTmp < errMin) || (errMin < 0)) {
					errMin = errTmp;
					copyQuaternion(quaternionTmp, quaternionMin, 0);
				}
			}

	/* Start from that optimal value and perform gradient descent */
	int nItr;
	double lambda=1.0;

	/*fprintf(stderr,"Start error %f\n", errMin);*/
	for (nItr = 0; nItr < 40; ++nItr) {
		errMin = reconstrPairValGrad(quaternionMin, wPrimeWPrimeTranspose, gradient, RTransposeR, dRTransposeR);

		/* stop if the gradient is too small or the error too */
		double normGradSq = normSq(4,gradient);
		if ((normGradSq < NUM_TOL_SQ) || (errMin < NUM_TOL ))
			break;

		/* Perform line search */
		int nItrLine;
		double errMinOld=errMin;
		for (nItrLine = 0; nItrLine < 20; ++nItrLine) {
			/* get the value for the current quaternion */
			for (k = 0; k < 4; ++k)
				quaternionTmp[k] = quaternionMin[k] - lambda * gradient[k];
			cleanQuaternion(quaternionTmp);
			errTmp = reconstrPairVal(quaternionTmp, wPrimeWPrimeTranspose,RTransposeR);

			/* check if the error is better than what we had before */
			if (errTmp < errMin) {
				errMin = errTmp;
				copyQuaternion(quaternionTmp, quaternionMin, 0);
				break;
			} else
				lambda /= 2;
		}

		/* stop if the error does not change much, < 0.1 percent */
		if ((errMinOld-errMin)<0.001*errMinOld)
			break;
	}

	/*fprintf(stderr,"%i ", nItr);*/
	/* update if nothing before or if the current result is better than before */
	if ((!isDone[iFrame * nFrame + jFrame]) || ((isDone[iFrame * nFrame + jFrame]) && (errMin <errArray[iFrame * nFrame + jFrame]))) {
		errArray[iFrame * nFrame + jFrame] = errMin;
		errArray[jFrame * nFrame + iFrame] = errMin;
		normalizeQuaternion(quaternionMin);
		copyQuaternion(quaternionMin, quaternionArray + 4 * (jFrame * nFrame + iFrame),0);
		isDone[iFrame * nFrame + jFrame] = 1;
		isDone[jFrame * nFrame + iFrame] = 1;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
	/* compute the infimum in CSFM, provided that the W
	(which is 2 x nPoint x nFrame) is centered */
	double *wOri, *q, *res;
	int i, ii, j, k, l, m, n, p, iFrame, jFrame, nFrame, nPoint;
	mxArray *pIsDone;

	/* Retrieve the useful data */
	wOri= mxGetPr(prhs[0]);
	const mwSize *dim_array = mxGetDimensions(prhs[0]);
	nPoint = dim_array[1];
	nFrame = dim_array[2];

	/* Create the output data */
	plhs[0] = mxCreateDoubleMatrix(nFrame, nFrame, mxREAL);
	double *errArray = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(4, nFrame * nFrame, mxREAL);
	double *quaternionArray = mxGetPr(plhs[1]);

	pIsDone = mxCreateLogicalMatrix(nFrame, nFrame);
	char *isDone = (char*) mxGetPr(pIsDone);
	
	for (iFrame = 0, k = 0; iFrame < nFrame; ++iFrame)
		for (jFrame = 0; jFrame < nFrame; ++jFrame, ++k) {
			errArray[k] = 0;
			if (iFrame == jFrame) {
				/* set the rotation to identity */
				isDone[k] = 1;
				quaternionArray[4 * k + 0] = 1;
				for (i = 1; i < 4; ++i)
					quaternionArray[4 * k + i] = 0;
			} else
				isDone[k] = 0;
		}

	/* create the pairwise W multiplications between frames */
	mxArray *pWPairwiseProd = mxCreateDoubleMatrix(2*nFrame, 2*nFrame, mxREAL);
	double *wPairwiseProd = mxGetPr(pWPairwiseProd);
	mxArray *pTmp2By2 = mxCreateDoubleMatrix(2, 2, mxREAL);
	double *tmp2By2 = mxGetPr(pTmp2By2);
	
	for (i = 0; i < nFrame; ++i )
		for (j = 0; j <= i; ++j ) {
			matrixMultiply(wOri + i * 2 * nPoint, wOri + j * 2 * nPoint, tmp2By2, 2, nPoint, 2);
			for( k = 0; k < 2; ++k )
				for( l = 0; l < 2; ++l )
					wPairwiseProd[(2*i+k)*2*nFrame + 2*j + l] = tmp2By2[l+2*k];
		}
	/* symmetrize the matrix */
	for (i = 0; i < 2*nFrame; ++i )
		for (j = i+1; j < 2*nFrame; ++j )
			wPairwiseProd[i*2*nFrame+j] = wPairwiseProd[j*2*nFrame+i];

	/* create some cache some temporary matrices */
	mxArray *pRTransposeR = mxCreateDoubleMatrix(4, 4, mxREAL);
	double *RTransposeR = mxGetPr(pRTransposeR);
	mxArray *pDRTransposeR = mxCreateDoubleMatrix(4, 4*4, mxREAL);
	double *dRTransposeR = mxGetPr(pDRTransposeR);
	mxArray *pQuaternionTmp = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *quaternionTmp = mxGetPr(pQuaternionTmp);
	mxArray *pGradient = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *gradient = mxGetPr(pGradient);
	mxArray *pWPrimeWPrimeTranspose = mxCreateDoubleMatrix(4, 4, mxREAL);
	double *wPrimeWPrimeTranspose = mxGetPr(pWPrimeWPrimeTranspose);

	/* do everything else */
	for (iFrame = 0; iFrame < nFrame; ++iFrame) {
		for (jFrame = iFrame+1; jFrame < nFrame; ++jFrame)
			gradientDescent(wPrimeWPrimeTranspose, iFrame, jFrame, wPairwiseProd, errArray, quaternionArray, isDone, nFrame, nPoint, RTransposeR, dRTransposeR, quaternionTmp, gradient);
		++iFrame;
		for (jFrame = nFrame - 1; jFrame > iFrame; --jFrame)
			gradientDescent(wPrimeWPrimeTranspose, iFrame, jFrame, wPairwiseProd, errArray, quaternionArray, isDone, nFrame, nPoint, RTransposeR, dRTransposeR, quaternionTmp, gradient);
	}
	/* do everything else in a different order */
	for (jFrame = nFrame-1; jFrame > 0; --jFrame) {
		for (iFrame = jFrame+1; iFrame < nFrame; ++iFrame)
			gradientDescent(wPrimeWPrimeTranspose, iFrame, jFrame, wPairwiseProd, errArray, quaternionArray, isDone, nFrame, nPoint, RTransposeR, dRTransposeR, quaternionTmp, gradient);

		--jFrame;
		for (iFrame = nFrame - 1; iFrame > jFrame; --iFrame)
			gradientDescent(wPrimeWPrimeTranspose, iFrame, jFrame, wPairwiseProd, errArray, quaternionArray, isDone, nFrame, nPoint, RTransposeR, dRTransposeR, quaternionTmp, gradient);
	}

	/* Free memory */
	mxDestroyArray(pIsDone);
	mxDestroyArray(pWPrimeWPrimeTranspose);
	mxDestroyArray(pTmp2By2);
	mxDestroyArray(pWPairwiseProd);
 	mxDestroyArray(pRTransposeR);
 	mxDestroyArray(pDRTransposeR);
	mxDestroyArray(pQuaternionTmp);
	mxDestroyArray(pGradient);
}
