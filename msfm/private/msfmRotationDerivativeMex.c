#include "mex.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
	/* compute the infimum in CSFM, provided that the W
	(which is 2 x nPoint x nFrame) is centered */
	double a, b, c, d;
	double  t25, t28, t62, t53, t26, t27, t1, t61, t60, t59, t58, t57, t56, t55, t54, t20, t24, t51, t19, t22, t50, t49, t48, t17, t23, t47, t46, t13, t14, t45, t15, t21, t44, t43, t42, t41, t40, t39, t16, t38, t18, t37, t36, t35, t34, t33, t12, t11, t10, t5, t4, t3;                                          
	int n, nFrame;

	/* Retrieve the useful data */
	double *quaternion = mxGetPr(prhs[0]);
	nFrame = mxGetN(prhs[0]);

	/* Create the output data */
	mwSize dim_array[4] = {2, 3, 4, nFrame};
	plhs[0] = mxCreateNumericArray(4, dim_array, mxDOUBLE_CLASS, mxREAL);
	double *dR = mxGetPr(plhs[0]);
	
	/* fill all the rotations */
	for( n = 0; n < nFrame; ++n) {
		a = quaternion[0];
		b = quaternion[1];
		c = quaternion[2];
		d = quaternion[3];
		
		t25 = b*b;                                       
		t28 = c*c;                                                
		t62 = (t25-t28)*d;                                        
		t53 = -t25-t28;                                           
		t26 = a*a;                                                
		t27 = d*d;                                                
		t1 = 1/(pow(t26+t27-t53,2.0));                            
		t61 = 4.0*t1;                                             
		t60 = 2.0*t1;                                             
		t59 = -2.0*t1;                                            
		t58 = -4.0*t1;                                            
		t57 = a*d;                                                
		t56 = c*b;                                                
		t55 = c*t25;                                              
		t54 = b*t28;                                              
		t20 = d*t26;                                              
		t24 = d*t27;                                              
		t51 = t20+t24;                                            
		t19 = b*t26;                                              
		t22 = b*t25;                                              
		t50 = t19+t22;                                            
		t49 = c*t57;                                              
		t48 = d*t56;                                              
		t17 = c*t27;                                              
		t23 = c*t28;                                              
		t47 = -t23-t17;                                           
		t46 = b*t57;                                              
		t13 = a*t28;                                              
		t14 = a*t27;                                              
		t45 = t14-t13;                                            
		t15 = a*t25;                                              
		t21 = a*t26;                                              
		t44 = t15+t21;                                            
		t43 = a*t56;                                              
		t42 = (t25+t27)*t61;                                      
		t41 = (t28+t27)*t61;                                      
		t40 = (t26+t28)*t58;
		t39 = (t26+t25)*t58;
		t16 = b*t27;
		t38 = t16-t54;
		t18 = c*t26;
		t37 = t18-t55;
		t36 = t44-t45;
		t35 = t38+t50;
		t34 = t37-t47;
		t33 = t20+t53*d-t24;
		t12 = -2.0*t43;
		t11 = -2.0*t49;
		t10 = -2.0*t48;
		t5 = 2.0*t43;
		t4 = 2.0*t46;
		t3 = 2.0*t48;
		dR[0] = a*t41;
		dR[1] = (t5+t33)*t59;
		dR[2] = (t12+t33)*t60;
		dR[3] = a*t42;
		dR[4] = (t4+t37+t47)*t59;
		dR[5] = (t19-t22-t54-t16+t11)*t60;
		dR[6] = b*t41;
		dR[7] = (-2.0*t46+t34)*t60;
		dR[8] = (t4+t34)*t60;
		dR[9] = b*t40;
		dR[10] = (t12-t62+t51)*t60;
		dR[11] = (t21-t15+t13+t14+t3)*t59;
		dR[12] = c*t39;
		dR[13] = (t11+t35)*t60;
		dR[14] = (2.0*t49+t35)*t60;
		dR[15] = c*t42;
		dR[16] = (t10+t44+t45)*t60;
		dR[17] = (t5+t62+t51)*t60;
		dR[18] = d*t39;
		dR[19] = (t10+t36)*t60;
		dR[20] = (t3+t36)*t59;
		dR[21] = d*t40;
		dR[22] = (t11-t38+t50)*t60;
		dR[23] = (t18+t55+t23-t17+t4)*t60;
		
		dR += 24;
		quaternion += 4;
	}
}
