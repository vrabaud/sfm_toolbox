#include "mex.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef MATLAB_MEX_FILE
#include "blas.h"
#if !defined(dgemm)
#define dgemm dgemm_
#define ddot ddot_
#endif
#else
/* if we are in Octave */
#if defined(_WIN32) || defined(_WIN64) || defined(__hpux)
#define FORTRAN_WRAPPER(x) x
#else
#define FORTRAN_WRAPPER(x) x ## _
#endif
typedef int mwSignedIndex;

#define ddot FORTRAN_WRAPPER(ddot)
double ddot(mwSignedIndex *, double *, mwSignedIndex *, double *,
			mwSignedIndex *);

#define dgemm FORTRAN_WRAPPER(dgemm)
double dgemm(char*, char*, mwSignedIndex *, mwSignedIndex *,
			 mwSignedIndex *, double *, double *, mwSignedIndex *,
			 double *, mwSignedIndex *, double *, double *,
			 mwSignedIndex *);
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

__inline void normalizeQuaternion(double *quaternion) {
	/* normalize */
	double norm = sqrt(quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);
	int k;
	
	for( k = 0; k < 4; ++k )
		quaternion[k] /= norm;
}

__inline void copyQuaternion(double *quaternionIn, double *quaternionOut) {
	int k;
	
	for (k = 0; k < 4; ++k)
		quaternionOut[k] = quaternionIn[k];
}

__inline double normSq(mwSignedIndex n, double *A) {
	mwSignedIndex one = 1;
	return ddot(&n, A, &one, A, &one);
	/* 	double res = dnrm2_(&n, A, &one);
	return res*res; */
}

/* Compute A*B', A is of size m *p and B n * p */
__inline double matrixMultiplyABt(double *A, double *B, double *C,
								  mwSignedIndex m, mwSignedIndex p,
								  mwSignedIndex n) {
	double one = 1.0, zero = 0.0;
	char chn[] = "N";
	char chnb[] = "C";

	dgemm(chn, chnb, &m, &n, &p, &one, A, &m, B, &n, &zero, C, &m);
}

/* Compute trace(A*B'), A and B are of size 2*2 */
__inline double matrixTraceProductABt(double *A, double *B, mwSignedIndex elemNbr) {
  double res = 0;
  double *A_end = A+elemNbr;
  for(;A<A_end;++A,++B)
    res += (*A)*(*B);
  return res;
  /* The following crashes on my new Matlab for some reason ... */
	/*mwSignedIndex one = 1;
	return ddot(&elemNbr, A, &one, B, &one);*/
}

__inline double reconstrPairVal(double *quaternion, double wWtTrace, double *sSt, double *wSt, double *r, double *rtR) {
	double a = quaternion[0], b = quaternion[1], c = quaternion[2], d = quaternion[3];
	
	double  t16, t18, t20, t15, t17, t19, t7, t5, t23, t10, t22, t21, t14, t12, t11, t9, t8, t6, t4, t3, t2, t1;
	t16 = c*c;                                                                             
	t18 = a*a;                                                                                
	t20 = t16+t18;                                                                            
	t15 = d*d;                                                                                
	t17 = b*b;                                                                                
	t19 = t15+t17;                                                                            
	t7 = t19+t20;                                                                             
	t5 = 1/t7;                                                                                
	t23 = 2.0*t5;                                                                             
	t10 = c*d-a*b;                                                                            
	t22 = 2.0*t10;                                                                            
	t21 = a*d;                                                                                
	t14 = b*c;                                                                                
	t12 = t21+t14;                                                                            
	t11 = t14-t21;
	t9 = a*c+b*d;
	t8 = t18+t17-t16-t15;
	t6 = -t19+t20;
	t4 = 1/(t7*t7);
	t3 = (2.0*t8*t9+4.0*t12*t10)*t4;
	t2 = (4.0*t11*t9+t6*t22)*t4;
	t1 = 2.0*(t8*t11+t12*t6)*t4;
	r[0] = t8*t5;
	r[1] = t12*t23;
	r[2] = t11*t23;
	r[3] = t6*t5;
	r[4] = t9*t23;
	r[5] = t5*t22;
	rtR[0] = (t8*t8+4.0*t12*t12)*t4;
	rtR[1] = t1;
	rtR[2] = t3;
	rtR[3] = t1;
	rtR[4] = (4.0*t11*t11+t6*t6)*t4;
	rtR[5] = t2;
	rtR[6] = t3;
	rtR[7] = t2;
	rtR[8] = 4.0*(t9*t9+t10*t10)*t4;
	
	return (wWtTrace - 2*matrixTraceProductABt(r, wSt, 6) + matrixTraceProductABt(rtR, sSt, 9));
}

__inline double reconstrPairValGrad(double *quaternion, double *gradient, double wWtTrace, double *sSt, double *wSt, double *r, double *dR, double *rtR, double *rtDR ) {
	double a = quaternion[0], b = quaternion[1], c = quaternion[2], d =
			quaternion[3];
	int k;
	
	double  t111, t107, t221, t176, t110, t72, t113, t106, t156, t253, t112, t117, t222, t259, t258, t109, t227, t66, t84, t76, t47, t130, t115, t87, t101, t88, t257, t121, t124, t255, t220, t73, t254, t252, t200, t199, t171, t39, t251, t172, t194, t62, t250, t196, t173, t191, t64, t249, t79, t179, t41, t197, t248, t192, t57, t184, t165, t247, t246, t189, t170, t190, t56, t245, t74, t244, t204, t243, t82, t158, t78, t242, t206, t241, t85, t163, t240, t27, t28, t7, t5, t238, t6, t4, t237, t236, t235, t234, t233, t232, t231, t230, t218, t48, t229, t228, t226, t83, t225, t49, t8, t224, t9, t223, t71, t77, t75, t10, t219, t217, t216, t215, t214, t213, t212, t211, t210, t80, t209, t198, t17, t208, t193, t169, t13, t207, t202, t81, t201, t26, t195, t29, t188, t187, t186, t185, t183, t182, t180, t178, t177, t175, t174, t168, t167, t166, t164, t162, t161, t153, t152, t151, t150, t149, t148, t147, t146, t145, t144, t143, t142, t95, t96, t138, t137, t93, t89, t136, t135, t69, t94, t70, t92, t134, t68, t133, t132, t67, t131, t129, t128, t60, t54, t52, t51, t46, t45, t44, t43, t25, t24, t23, t22, t21, t19, t16, t14, t3, t2, t1;                                                                             
	t111 = b*b;                                                                            
	t107 = b*t111;                                                                            
	t221 = t107*c;                                                                            
	t176 = a*t221;                                                                            
	t110 = d*d;                                                                               
	t72 = a*t110;                                                                             
	t113 = a*a;                                                                               
	t106 = a*t113;                                                                            
	t156 = -t106-t72;                                                                         
	t253 = t156*c;                                                                            
	t112 = c*c;                                                                               
	t117 = c*t112;                                                                            
	t222 = t117*a;                                                                            
	t259 = (-t222-t253)*b-t176;                                                               
	t258 = b*t253;                                                                            
	t109 = d*t110;                                                                            
	t227 = b*c;                                                                               
	t66 = t109*t227;                                                                          
	t84 = c*d;                                                                                
	t76 = b*t113;                                                                             
	t47 = t76*t84;                                                                            
	t130 = (t117*b+t221)*d;                                                                   
	t115 = t111*t111;                                                                         
	t87 = a*t115;                                                                             
	t101 = t112*t112;                                                                         
	t88 = a*t101;                                                                             
	t257 = 2.0*t130+6.0*t47-2.0*t66-t88+t87;                                                  
	t121 = t110*t110;                                                                         
	t124 = t113*t113;                                                                         
	t255 = t124+t121;                                                                         
	t220 = t109*a;                                                                            
	t73 = a*t111;                                                                             
	t254 = d*t222+(-t220+(-t106+t73)*d)*c;                                                    
	t252 = (6.0*t72-2.0*t106)*c;                                                              
	t200 = t113*t107;                                                                         
	t199 = t112*t113;                                                                         
	t171 = b*t199;                                                                            
	t39 = 2.0*t171;                                                                           
	t251 = 2.0*t200+t39;                                                                      
	t172 = d*t199;                                                                            
	t194 = t109*t112;                                                                         
	t62 = 2.0*t194;                                                                           
	t250 = t62+2.0*t172;                                                                      
	t196 = t113*t111;                                                                         
	t173 = d*t196;                                                                            
	t191 = t109*t111;                                                                         
	t64 = 2.0*t191;                                                                           
	t249 = t64+2.0*t173;                                                                      
	t79 = c*t113;                                                                             
	t179 = 2.0*t79;                                                                           
	t41 = t111*t179;                                                                          
	t197 = t113*t117;                                                                         
	t248 = t41+2.0*t197;                                                                      
	t192 = t107*t110;                                                                         
	t57 = 2.0*t192;                                                                           
	t184 = t110*t112;                                                                         
	t165 = 2.0*t184;                                                                          
	t247 = t57+b*t165;                                                                        
	t246 = t47+t66;                                                                           
	t189 = t110*t111;                                                                         
	t170 = a*t189;                                                                            
	t190 = t106*t111;                                                                         
	t56 = 2.0*t190;                                                                           
	t245 = 2.0*t170+t56;                                                                      
	t74 = b*t112;                                                                             
	t244 = -t107-t74;                                                                         
	t204 = (t115+t101)*b;                                                                     
	t243 = t204+2.0*t107*t112;                                                                
	t82 = d*t113;                                                                             
	t158 = t109+t82;                                                                          
	t78 = c*t111;                                                                             
	t242 = -t117-t78;                                                                         
	t206 = c*t115+t117*t112;                                                                  
	t241 = 2.0*t117*t111+t206;                                                                
	t85 = a*d;                                                                                
	t163 = 4.0*t85;                                                                           
	t240 = t242*t163;                                                                         
	t27 = t113+t112;                                                                          
	t28 = t111+t110;                                                                          
	t7 = t28+t27;                                                                             
	t5 = 1/(t7*t7);                                                                           
	t238 = 2.0*t5;                                                                            
	t6 = 1/t7;                                                                                
	t4 = t6*t5;                                                                               
	t237 = -4.0*t4;                                                                           
	t236 = 2.0*t4;                                                                            
	t235 = -4.0*t5;                                                                           
	t234 = -2.0*t5;                                                                           
	t233 = 4.0*t5;                                                                            
	t232 = 2.0*t6;                                                                            
	t231 = 4.0*t4;                                                                            
	t230 = -2.0*t4;                                                                           
	t218 = 2.0*t84;                                                                           
	t48 = t110+t113;                                                                          
	t229 = 4.0*t48;                                                                           
	t228 = -2.0*t107;                                                                         
	t226 = a*c;                                                                               
	t83 = b*d;                                                                                
	t225 = a*b;                                                                               
	t49 = t111+t112;                                                                          
	t8 = -t49+t48;                                                                            
	t224 = t4*t8;                                                                             
	t9 = t225+t84;                                                                            
	t223 = t4*t9;                                                                             
	t71 = a*t112;                                                                             
	t77 = c*t110;                                                                             
	t75 = b*t110;                                                                             
	t10 = t83-t226;                                                                           
	t219 = t10*t4;                                                                            
	t217 = 2.0*t225;                                                                          
	t216 = t77-t79;                                                                           
	t215 = 4.0*t246;                                                                          
	t214 = -2.0*t226;                                                                         
	t213 = -2.0*t84;                                                                          
	t212 = -2.0*t225;                                                                         
	t211 = -t76+t75;                                                                          
	t210 = -t73+t72;                                                                          
	t80 = d*t112;                                                                             
	t209 = t80-t82;                                                                           
	t198 = t117*t110;                                                                         
	t17 = -2.0*c*t189;                                                                        
	t208 = -2.0*t198+t17;                                                                     
	t193 = t106*t112;                                                                         
	t169 = a*t184;                                                                            
	t13 = -2.0*t169;                                                                          
	t207 = t13-2.0*t193;                                                                      
	t202 = t71+t106;                                                                          
	t81 = d*t111;                                                                             
	t201 = t109+t81;                                                                          
	t26 = t113+t111;                                                                          
	t195 = -t111+t112;                                                                        
	t29 = t112+t110;                                                                          
	t188 = t8*t234;                                                                           
	t187 = t112*t111;                                                                         
	t186 = -t113+t110;                                                                        
	t185 = -4.0*t223;                                                                         
	t183 = t27*t235;                                                                          
	t182 = t28*t233;                                                                          
	t180 = t29*t233;                                                                          
	t178 = 2.0*t222;                                                                          
	t177 = t26*t235;                                                                          
	t175 = 4.0*t219;                                                                          
	t174 = 2.0*t83+t214;                                                                      
	t168 = 4.0*t196;                                                                          
	t167 = 2.0*t106*t110+t255*a;                                                              
	t166 = t218+t217;                                                                         
	t164 = 2.0*t109*t113+t255*d;                                                              
	t162 = -4.0*t49*t224;                                                                     
	t161 = t224*t229;                                                                         
	t153 = t77+t79+t117-t78;                                                                  
	t152 = t244+t211;                                                                         
	t151 = -t202+t210;                                                                        
	t150 = -t71+t73-t156;                                                                     
	t149 = t201+t209;                                                                         
	t148 = t202+t210;                                                                         
	t147 = -t244+t211;                                                                        
	t146 = t75+t76-t74+t107;                                                                  
	t145 = t242+t216;                                                                         
	t144 = t80-t81+t158;                                                                      
	t143 = t201-t209;                                                                         
	t142 = -t242+t216;                                                                        
	t95 = t101*d;                                                                             
	t96 = t115*d;                                                                             
	t138 = -2.0*t176-t95+t96;                                                                 
	t137 = t244*t85;                                                                          
	t93 = t121*b;                                                                             
	t89 = t124*b;                                                                             
	t136 = d*t178+t109*t214+t93-6.0*t73*t84+t106*t213-t89;                                    
	t135 = 2.0*a*t187+t87+t88+t167;                                                           
	t69 = t106*t83;                                                                           
	t94 = c*t121;                                                                             
	t70 = b*t220;                                                                             
	t92 = t124*c;                                                                             
	t134 = 2.0*t69+6.0*t71*t83+t94+t85*t228-t92+2.0*t70;                                      
	t68 = t113*t75;                                                                           
	t133 = 2.0*t68+t93+t89+t243;                                                              
	t132 = t110*t179+t94+t92+t241;                                                            
	t67 = d*t187;                                                                             
	t131 = 2.0*t67+t95+t96+t164;                                                              
	t129 = t69+t137+t70;                                                                      
	t128 = -t130+t246;                                                                        
	t60 = 2.0*t198;                                                                           
	t54 = 2.0*t193;                                                                           
	t52 = -2.0*t197;                                                                          
	t51 = -2.0*t200;                                                                          
	t46 = c*t217;                                                                             
	t45 = d*t217;                                                                             
	t44 = a*t218;                                                                             
	t43 = b*t218;                                                                             
	t25 = c*t212;                                                                             
	t24 = d*t212;                                                                             
	t23 = a*t213;                                                                             
	t22 = b*t213;                                                                             
	t21 = -2.0*t173;                                                                          
	t19 = -2.0*t172;                                                                          
	t16 = -2.0*t170;                                                                          
	t14 = -2.0*b*t184;                                                                        
	t3 = t10*t9*t235;                                                                         
	t2 = t10*t188;                                                                            
	t1 = t9*t188;                                                                             
	r[0] = -(t186+t195)*t6;                                                                   
	r[1] = (t85+t227)*t232;                                                                   
	r[2] = -2.0*(-t227+t85)*t6;                                                               
	r[3] = -(t186-t195)*t6;                                                                   
	r[4] = (t226+t83)*t232;                                                                   
	r[5] = (t84-t225)*t232;                                                                   
	dR[0] = a*t180;                                                                           
	dR[1] = (t25+t149)*t238;                                                                  
	dR[2] = (t46+t149)*t234;                                                                  
	dR[3] = a*t182;                                                                           
	dR[4] = (t24+t142)*t238;                                                                  
	dR[5] = (t44+t147)*t234;                                                                  
	dR[6] = b*t180;                                                                           
	dR[7] = (t24+t153)*t238;                                                                  
	dR[8] = (t45+t153)*t238;                                                                  
	dR[9] = b*t183;                                                                           
	dR[10] = (t25+t144)*t238;                                                                 
	dR[11] = (t43+t148)*t234;                                                                 
	dR[12] = c*t177;                                                                          
	dR[13] = (t23+t146)*t238;                                                                 
	dR[14] = (t44+t146)*t238;                                                                 
	dR[15] = c*t182;                                                                          
	dR[16] = (t22+t150)*t238;                                                                 
	dR[17] = (t46+t143)*t238;                                                                 
	dR[18] = d*t177;                                                                          
	dR[19] = (t43+t151)*t234;                                                                 
	dR[20] = (t22+t151)*t238;                                                                 
	dR[21] = d*t183;                                                                          
	dR[22] = (t44+t152)*t234;                                                                 
	dR[23] = (t24+t145)*t234;                                                                 
	rtR[0] = (t7+t174)*(t7-t174)*t5;                                                          
	rtR[1] = t3;                                                                              
	rtR[2] = t2;                                                                              
	rtR[3] = t3;                                                                              
	rtR[4] = (t7+t166)*(t7-t166)*t5;                                                          
	rtR[5] = t1;                                                                              
	rtR[6] = t2;                                                                              
	rtR[7] = t1;                                                                              
	rtR[8] = t49*t5*t229;                                                                     
	rtDR[0] = (t45+t142)*t175;                                                                
	rtDR[1] = (t21+4.0*t172+t64+(-2.0*t222-t252)*b+t138+t164)*t230;                           
	rtDR[2] = (t94+(t117+t26*c)*t110+t129+t248)*t231;                                         
	rtDR[3] = (t19+t62+d*t168+(t178+t252)*b-t138+t164)*t236;                                  
	rtDR[4] = (t23+t147)*t185;                                                                
	rtDR[5] = (t68+t93-t244*t110+t251+t254)*t237;                                             
	rtDR[6] = (t52+t60-t244*t163+2.0*t186*t78+t132)*t230;                                     
	rtDR[7] = (-2.0*t171+t51+t240+t133+t247)*t236;                                            
	rtDR[8] = a*t162;                                                                         
	rtDR[9] = -4.0*(t46+t144)*t219;                                                           
	rtDR[10] = (t41+t52+(-2.0*t117-4.0*t77)*t111+t134-t206)*t230;                             
	rtDR[11] = (t95+t67+t158*t112+t249-t259)*t231;                                            
	rtDR[12] = (t17+t60+c*t168+t134+t241)*t230;                                               
	rtDR[13] = (t22+t148)*t185;
	rtDR[14] = (t88+(t106+t28*a)*t112+t128+t245)*t231;
	rtDR[15] = (-2.0*t191+t21-4.0*t258+t131+t250)*t230;
	rtDR[16] = (t16+a*t165-2.0*t190+t54+t135-t215)*t230;
	rtDR[17] = b*t161;
	rtDR[18] = (t43+t150)*t175;
	rtDR[19] = (4.0*t171+t57+t14+t136+t243)*t230;
	rtDR[20] = (-t87+(-t106-t29*a)*t111+t128+t207)*t231;
	rtDR[21] = (t51+t39+(t228-4.0*t75)*t112+t136-t204)*t230;
	rtDR[22] = (t25+t143)*t185;
	rtDR[23] = (t96+t67+t158*t111+t250+t259)*t231;
	rtDR[24] = (t135+t207+t215+t245)*t236;
	rtDR[25] = (t19-2.0*t194+4.0*t258+t131+t249)*t230;
	rtDR[26] = c*t161;
	rtDR[27] = (t23+t152)*t175;
	rtDR[28] = (t54+4.0*t170+t13+t167-t257)*t236;
	rtDR[29] = (t68+t89-t244*t113+t247-t254)*t237;
	rtDR[30] = (t56+t16+4.0*t169+t167+t257)*t230;
	rtDR[31] = 4.0*(t45+t145)*t223;
	rtDR[32] = (-t92+(-t117-t28*c)*t113+t129+t208)*t231;
	rtDR[33] = (-2.0*t192+t14-t240+t133+t251)*t236;
	rtDR[34] = (4.0*t137+t132+t208+t248)*t236;
	rtDR[35] = d*t162;
	
	/* Compute the gradient */
	for (k = 0; k < 4; ++k, rtDR += 9, dR += 6)
		gradient[k] = 2*matrixTraceProductABt(sSt, rtDR, 9) - 2*matrixTraceProductABt(dR, wSt, 6);

	/* Return the error */
	return (wWtTrace - 2*matrixTraceProductABt(r, wSt, 6) + matrixTraceProductABt(rtR, sSt, 9));
}

void gradientDescent(double wWtTrace, double *sSt, double *wSt, double *errArray, double *quaternionArray, int iFrame, int iFrameIni, double *r, double *dR, double *rtR, double *rtDR, double *quaternionTmp, double *gradient) {
	int k;

	/* Figure out the best guess to start from */
	double errMin, errTmp, errMinOld;
	double *quaternionMin = quaternionArray + 4 * iFrame;
	int nItr;
	double lambda=1.0;
	int nItrLine;
	double normGradSq;

	copyQuaternion(quaternionArray + 4 * iFrameIni, quaternionMin);
	
	/* Start from that optimal value and perform gradient descent */

	for (nItr = 0; nItr < 40; ++nItr) {
		errMin = reconstrPairValGrad(quaternionMin, gradient, wWtTrace, sSt, wSt, r, dR, rtR, rtDR );

    /* stop if the gradient is too small or the error too */
		normGradSq = normSq(4,gradient);

		if ((normGradSq < NUM_TOL_SQ) || (errMin < NUM_TOL ))
			break;

		/* Perform line search */
		errMinOld=errMin;
		for (nItrLine = 0; nItrLine < 20; ++nItrLine) {
			/* get the value for the current quaternion */
			for (k = 0; k < 4; ++k)
				quaternionTmp[k] = quaternionMin[k] - lambda * gradient[k];
			normalizeQuaternion(quaternionTmp);
			errTmp = reconstrPairVal(quaternionTmp, wWtTrace, sSt, wSt, r, rtR);

			/* check if the error is better than what we had before */
			if (errTmp < errMin) {
				errMin = errTmp;
				copyQuaternion(quaternionTmp, quaternionMin);
				break;
			} else
				lambda /= 2;
		}

		/* stop if the error does not change much, < 0.1 percent */
		if ((errMinOld-errMin)<0.001*errMinOld)
			break;
	}

	/* update if nothing before or if the current result is better than before */
	if ((errMin < errArray[iFrame]) || (errArray[iFrame]<0))
		errArray[iFrame] = errMin;
/*   fprintf( stdout, "Hello %f\n", errArray[0] );*/
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
	/* refines the exterior orientation given
	 * w (which is 2 x nPoint x nFrame)
	 * s (which is 3 x nPoint x nFrame)
	 * q (which is 4 x nFrame) an initial estimate of the quaternion
	 */
	double *wOri, *q, *res, *w, *s, *quaternionIn;
	int i, j, k, iFrame, nFrame, nPoint;
  const mwSize *dim_array;
  double *quaternionOut;
  double *wWtTrace;
  double *errArray;
  mxArray *pWWtTrace;
  mxArray *pSSt;
  double *sSt;
  mxArray *pWSt;
  double *wSt;
	mxArray *pR, *pDR, *pRtR, *pRtDR, *pQuaternionTmp, *pGradient;
  double *r, *dR, *rtR, *rtDR, *quaternionTmp, *gradient;

  /* Retrieve the useful data */
	w = mxGetPr(prhs[0]);
	dim_array = mxGetDimensions(prhs[0]);
	nPoint = dim_array[1];
	if (mxGetNumberOfDimensions(prhs[0])==2)
		nFrame = 1;
	else
		nFrame = dim_array[2];
	s = mxGetPr(prhs[1]);
	quaternionIn = mxGetPr(prhs[2]);

	/* Create the output data */
	plhs[0] = mxCreateDoubleMatrix(4, nFrame, mxREAL);
	quaternionOut = mxGetPr(plhs[0]);
	for( i = 0, k = 0; i < nFrame; ++i )
		for( j = 0; j < 4; ++j, ++k )
			quaternionOut[k] = quaternionIn[k];

	/* fill the errors with -1 */
	plhs[1] = mxCreateDoubleMatrix(1, nFrame, mxREAL);
	errArray = mxGetPr(plhs[1]);
	for( i = 0; i < nFrame; ++i )
		errArray[i] = -1.0;

	/* create the pairwise ww',ss', sw' multiplications between frames */
	pWWtTrace = mxCreateDoubleMatrix(1, nFrame, mxREAL);
	wWtTrace = mxGetPr(pWWtTrace);
	pSSt = mxCreateDoubleMatrix(9, nFrame, mxREAL);
	sSt = mxGetPr(pSSt);
	pWSt = mxCreateDoubleMatrix(6, nFrame, mxREAL);
	wSt = mxGetPr(pWSt);
	
	for ( i = 0; i < nFrame; ++i, w += 2 * nPoint, s += 3 * nPoint, sSt += 9, wSt += 6 ) {
		wWtTrace[i] = matrixTraceProductABt(w, w, 2*nPoint);
		matrixMultiplyABt(s, s, sSt, 3, nPoint, 3);
		matrixMultiplyABt(w, s, wSt, 2, nPoint, 3);
	}
	
	/* create some cache some temporary matrices */
	pR = mxCreateDoubleMatrix(2, 3, mxREAL);
	r = mxGetPr(pR);
	pDR = mxCreateDoubleMatrix(6, 4, mxREAL);
	dR = mxGetPr(pDR);
	pRtR = mxCreateDoubleMatrix(3, 3, mxREAL);
	rtR = mxGetPr(pRtR);
	pRtDR = mxCreateDoubleMatrix(9, 4, mxREAL);
	rtDR = mxGetPr(pRtDR);
	pQuaternionTmp = mxCreateDoubleMatrix(4, 1, mxREAL);
	quaternionTmp = mxGetPr(pQuaternionTmp);
	pGradient = mxCreateDoubleMatrix(4, 1, mxREAL);
	gradient = mxGetPr(pGradient);

	for (iFrame = 0; iFrame < nFrame; ++iFrame, sSt += 9, wSt += 6)
		gradientDescent(wWtTrace[iFrame], sSt, wSt, errArray, quaternionOut, iFrame, iFrame, r, dR, rtR, rtDR, quaternionTmp, gradient);

  /* do everything else a few times */
	for( i = 0; i < 5; ++i ) {
		for (iFrame = 1; iFrame < nFrame; ++iFrame)
			gradientDescent(wWtTrace[iFrame], sSt+9*iFrame, wSt+6*iFrame, errArray, quaternionOut, iFrame, iFrame-1, r, dR, rtR, rtDR, quaternionTmp, gradient);

    /* do everything else in reverse order */
		for (iFrame = nFrame-2; iFrame >= 0; --iFrame)
			gradientDescent(wWtTrace[iFrame], sSt+9*iFrame, wSt+6*iFrame, errArray, quaternionOut, iFrame, iFrame+1, r, dR, rtR, rtDR, quaternionTmp, gradient);
	}

	/* Free memory */
	mxDestroyArray(pWWtTrace);
	mxDestroyArray(pSSt);
	mxDestroyArray(pWSt);
 	mxDestroyArray(pR);
	mxDestroyArray(pDR);
	mxDestroyArray(pRtR);
	mxDestroyArray(pRtDR);
	mxDestroyArray(pQuaternionTmp);
	mxDestroyArray(pGradient);
}
