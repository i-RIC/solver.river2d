#pragma options cpluscmt	
#include "Fe_PS.h"


extern struct control  N ;
extern double 	UW ;

int	get_shape(el_type,x,psf)
int	el_type ;
struct gausspts		*x ;
struct shapefuncs	*psf ;

{
	int	err ;
	switch (el_type) {
	
		case 210 :
			err = shap_2t1(x,psf) ;
			break ;
		
		case 211 :
			err = shap_2s1(x,psf) ;
			break ;
		
		case 216 :
			err = shap_2s1(x,psf) ;
			break ;
		
		case 220 :
			err = shap_2t2(x,psf) ;
			break ;
		
		case 221 :
			err = shap_2s2(x,psf) ;
			break ;
		
		case 226 :
			err = shap_2s2(x,psf) ;
			break ;
			
		case 222 :
			err = shap_2l2(x,psf) ;
			break ;
		
		case 227 :
			err = shap_2l2(x,psf) ;
			break ;
			
		case 228 :
			err = shqu6_2l2(x,psf) ;
			break ;
			
		case 229 :
			err = shqu9_2l2(x,psf) ;
			break ;
			
		case 224 :
			err = shqu8_2l2(x,psf) ;
			break ;
			
		case 219 :
			err = shqu_2l1(x,psf) ;
			break ;
			
		case 231 :
			err = shap_2s3(x,psf) ;
			break ;
		
		case 111 :
			err = shap_1s1(x,psf) ;
			break ;
			
		case 118 :
			err = SU1D_1s1(x,psf) ;
			break ;

		case 119 :
			err = SU1D_1s1(x,psf) ;
			break ;

		case 121 :
			err = shap_1s2(x,psf) ;
			break ;
		
		case 131 :
			err = shap_1s3(x,psf) ;
			break ;
		
		case 128 :
			err = shcu3r_1s2(x,psf) ;
			break ;
		
		case 129 :
			err = shcu3_1s2(x,psf) ;
			break ;
  		case 139 :
			err = shapcu_1s3(x,psf) ;
			break ;
		
		case 138 :
			err = shapcur_1s3(x,psf) ;
			break ;
/*		
		case 149 :
			err = shapcu_1s4(x,psf) ;
			break ;
		
		case 159 :
			err = shapcu_1s5(x,psf) ;
			break ;
		
		case 169 :
			err = shapcu_1s6(x,psf) ;
			break ;
		
		case 179 :
			err = shapcu_1s7(x,psf) ;
			break ;
		
		case 178 :
			err = shcu_1s7r(x,psf) ;
			break ;
*/			
		case 1 :
			err = shap_pt(x,psf) ;
			break ;

		case 230 :
			err = shcu2_2x1(x,psf) ;
			break ;

		case 2311 :
			err = shcu3_2x2(x,psf) ;
			break ;
		
		case 232 :
			err = shcu4_2x2(x,psf) ;
			break ;
		
		case 233 :
			err = shcu3_3x1(x,psf) ;
			break ;
		
		case 2341 :
			err = shcu4_3x2(x,psf) ;
			break ;
		
		case 2342 :
			err = sh4_3x2b(x,psf) ;
			break ;
				
		case 2351 :
			err = shcu5_3x2(x,psf) ;
			break ;

		case 2353 :
			err = sh5_3x2b(x,psf) ;
			break ;
				
		case 2352 :
			err = shcu5_3x3(x,psf) ;
			break ;
		
		case 2361 :
			err = shcu6_3x2(x,psf) ;
			break ;
		
		case 2363 :
			err = sh6_3x2b(x,psf) ;
			break ;
				
		case 2362 :
			err = shcu6_3x3(x,psf) ;
			break ;
		
		case 2371 :
			err = shcu7_3x3(x,psf) ;
			break ;

		case 2372 :
			err = sh7_3x3b(x,psf) ;
			break ;
				
		case 238 :
			err = shcu8_3x3(x,psf) ;
			break ;
	
		case 239 :
			err = shcu9_3x3(x,psf) ;
			break ;
				
		case 235 :
			err = hermite(x,psf) ;
			break ;
	}
	return(err) ;
}


int	nsf(el_type)
int	el_type ;

{
	struct gausspts		x ;
	struct shapefuncs	psf ;
	int	err ;

	x.x = 100.0 ;
	x.y = 100.0 ;
	x.z = 100.0 ;
	err = get_shape(el_type,&x,&psf) ;
	return(err) ;
}


int	nbounds(el_type)
int	el_type ;

{
	int	err ;
	
	err = 0 ;
	switch (el_type) {
	
		case 210 :
			err = 3 ;
			break ;
		
		case 211 :
			err = 4 ;
			break ;
		
		case 216 :
			err = 4 ;
			break ;
		
		case 220 :
			err = 6 ;
			break ;
		
		case 221 :
			err = 8 ;
			break ;
		
		case 222 :
			err = 8 ;
			break ;
		
		case 226 :
			err = 8 ;
			break ;
		
		case 227 :
			err = 8 ;
			break ;
		
		case 224 :
			err = 4 ;
			break ;
		
		case 228 :
			err = 4 ;
			break ;
		
		case 229 :
			err = 4 ;
			break ;
		
		case 239 :
			err = 4 ;
			break ;
		
		case 230 :
			err = 4 ;
			break ;
		
		case 2311 :
			err = 4 ;
			break ;
		
		case 232 :
			err = 4 ;
			break ;
		
		case 233 :
			err = 4 ;
			break ;
		
		case 2341 :
			err = 4 ;
			break ;
		
		case 2342 :
			err = 4 ;
			break ;
				
		case 2351 :
			err = 4 ;
			break ;
				
		case 2353 :
			err = 4 ;
			break ;
		
		case 2352 :
			err = 4 ;
			break ;
		
		case 2361 :
			err = 4 ;
			break ;

		case 2363 :
			err = 4 ;
			break ;
				
		case 2362 :
			err = 4 ;
			break ;
		
		case 2371 :
			err = 4 ;
			break ;
		
		case 2372 :
			err = 4 ;
			break ;
			
		case 238 :
			err = 4 ;
			break ;
	
		case 231 :
			err = 12 ;
			break ;
		
		case 111 :
			err = 2 ;
			break ;
		
		case 118 :
			err = 2 ;
			break ;
		
		case 119 :
			err = 2 ;
			break ;
		
		case 128 :
			err = 2 ;
			break ;
		
		case 129 :
			err = 2 ;
			break ;
		
		case 139 :
			err = 2 ;
			break ;
		
		case 138 :
			err = 2 ;
			break ;
		
		case 121 :
			err = 2 ;
			break ;
		
		case 131 :
			err = 2 ;
			break ;
		
		case 1 :
			err = 1 ;
			break ;
		
		case 235 :
			err = 4 ;
			break ;
	}
	return(err) ;
}

int		get_nlcs(el_type,k,r,s) 
int		k, el_type ;
double	*r, *s ;

{
	double	athird = 0.33333333333 ;
			
	switch (el_type) {
	
		case 210 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
			case 1 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
			}
			return(3) ;
		
		case 211 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
			}
			return(4) ;
		
		case 216 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
			}
			return(4) ;
			
		case 220 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 1 :
					*r = 0.5 ;
					*s = 0.5 ;
					break ;
				case 2 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 0.5 ;
					break ;
				case 4 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
				case 5 :
					*r = 0.5 ;
					*s = 0.0 ;
					break ;
			}
			return(6) ;
	
		case 221 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 2 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 7 :
					*r = 0.0 ;
					*s = -1.0 ;
					break ;
			}
			return(8) ;
	
		case 226 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 2 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 7 :
					*r = 0.0 ;
					*s = -1.0 ;
					break ;
			}
			return(8) ;
	
		case 222 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 2 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 7 :
					*r = 0.0 ;
					*s = -1.0 ;
					break ;
				case 8 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
			}
			return(9) ;
	
		case 229 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 7 :
					*r = 0.0 ;
					*s = -1.0 ;
					break ;
				case 8 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
			}
			return(9) ;
	
		case 228 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
			}
			return(9) ;
	
		case 227 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = 0.0 ;
					break ;
				case 2 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 4 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -1.0 ;
					*s = 0.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 7 :
					*r = 0.0 ;
					*s = -1.0 ;
					break ;
				case 8 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
			}
			return(9) ;
	
		case 231 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;
				case 1 :
					*r = 1.0 ;
					*s = -athird ;
					break ;
				case 2 :
					*r = 1.0 ;
					*s = athird ;
					break ;
				case 3 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 4 :
					*r = athird ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -athird ;
					*s = 1.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 7 :
					*r = -1.0 ;
					*s = athird ;
					break ;
				case 8 :
					*r = -1.0 ;
					*s = -athird ;
					break ;
				case 9 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 10 :
					*r = -athird ;
					*s = -1.0 ;
					break ;
				case 11 :
					*r = athird ;
					*s = -1.0 ;
					break ;
			}
			return(8) ;
	
		case 111 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 1 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
			}
			return(2) ;
	
		
		case 121 :
			switch (k) {
				case 0 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
				case 1 :
					*r = 0.0 ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
			}
			return(3) ;
	
		case 131 :
			switch (k) {
				case 0 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;
				case 1 :
					*r = -athird ;
					*s = 1.0 ;
					break ;
				case 2 :
					*r = athird ;
					*s = 1.0 ;
					break ;
				case 3 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;
			}
			return(4) ;
			
			
		case 1 :
			switch (k) {
				case 0 :
					*r = 0.0 ;
					*s = 0.0 ;
					break ;
			}
			return(1) ;
			
		case 239 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -5.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 10 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 11 :
					*r = -5.0 ;
					*s = -5.0 ;
					break ;	
				case 12 :
					*r = -3.0 ;
					*s = -5.0 ;
					break ;	
				case 13 :
					*r = -1.0 ;
					*s = -5.0 ;
					break ;	
				case 14 :
					*r = 1.0 ;
					*s = -5.0 ;
					break ;	
				case 15 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (16) ;
		case 238 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -5.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 10 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 11 :
					*r = -3.0 ;
					*s = -5.0 ;
					break ;	
				case 12 :
					*r = -1.0 ;
					*s = -5.0 ;
					break ;	
				case 13 :
					*r = 1.0 ;
					*s = -5.0 ;
					break ;	
				case 14 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (15) ;
		case 2371 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -5.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 10 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 11 :
					*r = -1.0 ;
					*s = -5.0 ;
					break ;	
				case 12 :
					*r = 1.0 ;
					*s = -5.0 ;
					break ;	
				case 13 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (14) ;	
				
		case 2372 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*s = -3.0 ;
					*r = 1.0 ;
					break ;	
				case 5 :
					*s = -5.0 ;
					*r = 1.0 ;
					break ;	
				case 6 :
					*s = -5.0 ;
					*r = -1.0 ;
					break ;	
				case 7 :
					*s = -3.0 ;
					*r = -1.0 ;
					break ;	
				case 8 :
					*s = -5.0 ;
					*r = -3.0 ;
					break ;	
				case 9 :
					*s = -3.0 ;
					*r = -3.0 ;
					break ;	
				case 10 :
					*s = -1.0 ;
					*r = -3.0 ;
					break ;	
				case 11 :
					*s = -1.0 ;
					*r = -5.0 ;
					break ;	
				case 12 :
					*s = 1.0 ;
					*r = -5.0 ;
					break ;	
				case 13 :
					*s = 1.0 ;
					*r = -3.0 ;
					break ;	
				}
				return (14) ;		
		case 2362 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 10 :
					*r = -1.0 ;
					*s = -5.0 ;
					break ;	
				case 11 :
					*r = 1.0 ;
					*s = -5.0 ;
					break ;	
				case 12 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (13) ;
		case 2361 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -5.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 10 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;		
				case 11 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (12) ;
				
		case 2363 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*s = -3.0 ;
					*r = 1.0 ;
					break ;	
				case 5 :
					*s = -5.0 ;
					*r = 1.0 ;
					break ;	
				case 6 :
					*s = -5.0 ;
					*r = -1.0 ;
					break ;	
				case 7 :
					*s = -3.0 ;
					*r = -1.0 ;
					break ;	
				case 8 :
					*s = -5.0 ;
					*r = -3.0 ;
					break ;	
				case 9 :
					*s = -3.0 ;
					*r = -3.0 ;
					break ;	
				case 10 :
					*s = -1.0 ;
					*r = -3.0 ;
					break ;		
				case 11 :
					*s = 1.0 ;
					*r = -3.0 ;
					break ;	
				}
				return (12) ;		
		case 2351 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;
				case 8 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;
				case 10 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (11) ;
				
		case 2353 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*s = -3.0 ;
					*r = 1.0 ;
					break ;	
				case 5 :
					*s = -5.0 ;
					*r = 1.0 ;
					break ;	
				case 6 :
					*s = -5.0 ;
					*r = -1.0 ;
					break ;	
				case 7 :
					*s = -3.0 ;
					*r = -1.0 ;
					break ;
				case 8 :
					*s = -3.0 ;
					*r = -3.0 ;
					break ;	
				case 9 :
					*s = -1.0 ;
					*r = -3.0 ;
					break ;
				case 10 :
					*s = 1.0 ;
					*r = -3.0 ;
					break ;	
				}
				return (11) ;
		
		case 2352 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = -1.0 ;
					*s = -5.0 ;
					break ;	
				case 10 :
					*r = 1.0 ;
					*s = -5.0 ;
					break ;	
				case 11 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (12) ;
		case 2341 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 8 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;	
				case 9 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;	
				}
				return (10) ;		
		case 2342 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*s = -3.0 ;
					*r = 1.0 ;
					break ;	
				case 5 :
					*s = -5.0 ;
					*r = 1.0 ;
					break ;	
				case 6 :
					*s = -5.0 ;
					*r = -1.0 ;
					break ;	
				case 7 :
					*s = -3.0 ;
					*r = -1.0 ;
					break ;	
				case 8 :
					*s = -1.0 ;
					*r = -3.0 ;
					break ;	
				case 9 :
					*s = 1.0 ;
					*r = -3.0 ;
					break ;	
				}
				return (10) ;		
						
		case 233 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;	
				case 5 :
					*r = -5.0 ;
					*s = 1.0 ;
					break ;	
				case 6 :
					*r = -5.0 ;
					*s = -1.0 ;
					break ;	
				case 7 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;
				}
				return (8) ;
				
		case 232 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				case 6 :
					*r = -3.0 ;
					*s = -3.0 ;
					break ;	
				case 7 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;
				case 8 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;
				}
				return (9) ;
				
		case 2311 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;
				case 6 :
					*r = -1.0 ;
					*s = -3.0 ;
					break ;
				case 7 :
					*r = 1.0 ;
					*s = -3.0 ;
					break ;
				}
				return (8) ;
				
		case 230 :
			switch (k) {
				case 3 :
					*r = -1.0 ;
					*s = -1.0 ;
					break ;
				case 0 :
					*r = 1.0 ;
					*s = -1.0 ;
					break ;	
				case 1 :
					*r = 1.0 ;
					*s = 1.0 ;
					break ;	
				case 2 :
					*r = -1.0 ;
					*s = 1.0 ;
					break ;	
				case 4 :
					*r = -3.0 ;
					*s = 1.0 ;
					break ;
				case 5 :
					*r = -3.0 ;
					*s = -1.0 ;
					break ;	
				}
				return (6) ;
		}
		return(0) ;
}


int	shap_2t1(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t ;
	
	r = x->x ;
	s = x->y ;
	t = 1.0 - r - s ;
	if( (r<-0.0001) || (r>1.0001) || (s<-0.00001) || (s>1.00001) || (t<-0.00001) ) return(3) ;
	
	p->dof = 3 ;
	
	p->f[0] = r ;
	p->f[1] = s ;
	p->f[2] = t ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[0] =  1.0 ;
	p->dfdr[1] =  0.0 ;
	p->dfdr[2] = -1.0 ;
	
	p->dfds[0] =  0.0 ;
	p->dfds[1] =  1.0 ;
	p->dfds[2] = -1.0 ;
	
	return(0) ;
}


int	shap_2s1(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) {
//printf("r=%g\ts=%g\n",r,s);
	return(4) ;
	}
	
	t = 0.25 ;
	mr = t * (1.0 - r) ;
	ms = 1.0 - s ;
	pr = t * (1.0 + r) ;
	ps = 1.0 + s ;
	
	p->dof = 4 ;
	
	p->f[0] = pr * ms ;
	p->f[1] = pr * ps ;
	p->f[2] = mr * ps ;
	p->f[3] = mr * ms ;
	
	if(x->w < 0.0)
		return(0) ;
	ms *= t ;
	ps *= t ;
	
	p->dfdr[0] =  ms ;
	p->dfdr[1] =  ps ;
	p->dfdr[2] = -ps ;
	p->dfdr[3] = -ms ;
	
	p->dfds[0] = -pr ;
	p->dfds[1] =  pr ;
	p->dfds[2] =  mr ;
	p->dfds[3] = -mr ;
	
	return(0) ;
}
int	shap_2t2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t ;
	
	r = x->x ;
	s = x->y ;
	t = 1.0 - r - s ;
	if( (r<-0.00001) || (r>1.00001) || (s<-0.00001) || (s>1.00001) || (t<-0.00001) ) return(6) ;
	
	p->dof = 6 ;
	
	p->f[0] = r * (2.0*r - 1.0) ;
	p->f[1] = 4.0 * r * s ;
	p->f[2] = s * (2.0*s - 1.0) ;
	p->f[3] = 4.0 * s * t ;
	p->f[4] = t * (2.0*t - 1.0) ;
	p->f[5] = 4.0 * r * t ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[0] =  4.0 * r - 1.0 ;
	p->dfdr[1] =  4.0 * s ;
	p->dfdr[2] =  0.0 ;
	p->dfdr[3] = -4.0 * s ;
	p->dfdr[4] = -4.0 * t + 1.0 ;
	p->dfdr[5] =  4.0 * (t - r) ;
	
	p->dfds[0] =  0.0 ;
	p->dfds[1] =  4.0 * r ;
	p->dfds[2] =  4.0 * s - 1.0 ;
	p->dfds[3] =  4.0 * (t - s) ;
	p->dfds[4] = -4.0 * t + 1.0 ;
	p->dfds[5] = -4.0 * r ;
	
	return(0) ;
}

int	shap_2s2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t,tr,ts,tprps,tmrps,tprms,tmrms ;
	double	trms,trps,tsmr,tspr,tps,tms,tpr,tmr ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(8) ;
	
	t = 0.25 ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	tprms = t * pr * ms ;
	tprps = t * pr * ps ;
	tmrps = t * mr * ps ;
	tmrms = t * mr * ms ;
	
	p->dof = 8 ;
	
	p->f[0] = tprms * (r - ps) ;
	p->f[2] = tprps * (r - ms) ;
	p->f[4] = tmrps * (s - pr) ;
	p->f[6] = tmrms * (-s - pr) ;
	tprms *= 2.0 ;
	tmrps *= 2.0 ;
	p->f[1] = tprms * ps ;
	p->f[3] = tmrps * pr ;
	p->f[5] = tmrps * ms ;
	p->f[7] = tprms * mr ;
	
	if(x->w < 0.0)
		return(0) ;
		
	ts = 2.0 * s ;
	tr = 2.0 * r ;
	trms = tr - s ;
	trps = tr + s ;
	tsmr = ts - r ;
	tspr = ts + r ;
	tps = t * ps ;
	tms = t * ms ;
	tpr = t * pr ;
	tmr = t * mr ;
	
	p->dfdr[0] = tms * (trms) ;
	p->dfdr[2] = tps * (trps) ;
	p->dfdr[4] = tps * (trms) ;
	p->dfdr[6] = tms * (trps) ;
	p->dfdr[1] =  0.5 * ms * ps ;
	p->dfdr[3] = -r * ps ;
	p->dfdr[5] = -(p->dfdr[1]) ;
	p->dfdr[7] = -r * ms ;
	
	p->dfds[0] = tpr * (tsmr) ;
	p->dfds[2] = tpr * (tspr) ;
	p->dfds[4] = tmr * (tsmr) ;
	p->dfds[6] = tmr * (tspr) ;
	p->dfds[1] = -s * pr ;
	p->dfds[3] = 0.5 * pr * mr ;
	p->dfds[5] = -s * mr ;
	p->dfds[7] = - (p->dfds[3]) ;
	
	return(0) ;
}

int	shap_2l2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,t,tr,ts,pmr,pms ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(9) ;
	
	t = 0.25 ;
	tr = 2.0 * r ;
	ts = 2.0 * s ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	pmr = pr * mr ;
	pms = ps * ms ;
	
	p->dof = 9 ;
	
	p->f[0] = -t * pr * ms * r * s ;
	p->f[2] = t * pr * ps * r * s ;
	p->f[4] = -t * r * mr * s * ps ;
	p->f[6] = t * r * mr * s * ms ;
	p->f[1] = 0.5 * r * pr * pms ;
	p->f[3] = 0.5 * pmr * s * ps ;
	p->f[5] = -0.5 * r * mr * pms ;
	p->f[7] = -0.5 * pmr * s * ms ;
	p->f[8] = pmr * pms ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[0] = -t * s * ms * (tr + 1.0) ;
	p->dfdr[2] = t * s * ps * (1.0 + tr) ;
	p->dfdr[4] = -t * s * ps * (1.0 - tr) ;
	p->dfdr[6] = t * s * ms * (1.0 - tr) ;
	p->dfdr[1] =  0.5 * pms * (1.0 + tr) ;
	p->dfdr[3] = -r * s * ps ;
	p->dfdr[5] = -0.5 * (1.0 - tr) * pms ;
	p->dfdr[7] =  r * s * ms ;
	p->dfdr[8] = -2.0 * r * pms ;
	
	p->dfds[0] = -t * r * pr * (1.0 - ts) ;
	p->dfds[2] = t * r * pr * (1.0 + ts) ;
	p->dfds[4] = -t * r * mr * (1.0 + ts) ;
	p->dfds[6] = t * r * mr * (1.0 - ts) ;
	p->dfds[1] = -s * r * pr ;
	p->dfds[3] = 0.5 * pmr * (1.0 + ts) ;
	p->dfds[5] = s * r * mr ;
	p->dfds[7] = -0.5 * pmr * (1.0 - ts) ;
	p->dfds[8] = -2.0 * s * pmr ;
	
	return(0) ;
}
int	shap_2s3(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,s,mr,ms,pr,ps,m3r,p3r,m3s,p3s,r2,mr2,s2,ms2,frs,t,nt ;
	double	temp, mr2p3r,mr2m3r,ms2p3s,ms2m3s ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(12) ;
	
	t = 0.03125 ;
	nt = 0.28125 ;
	p3r = 1.0 + 3.0 * r ;
	m3r = 1.0 - 3.0 * r ;
	p3s = 1.0 + 3.0 * s ;
	m3s = 1.0 - 3.0 * s ;
	mr = 1.0 - r ;
	ms = 1.0 - s ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	r2 = r * r ;
	s2 = s * s ;
	mr2 = nt * (1.0 - r2) ;
	ms2 = nt * (1.0 - s2) ;
	frs = t * (-10 + 9 * (r2 + s2)) ;
	ms2m3s = ms2 * m3s ;
	ms2p3s = ms2 * p3s ;
	mr2m3r = mr2 * m3r ;
	mr2p3r = mr2 * p3r ;
	
	p->dof = 12 ;
	
	p->f[0] = pr * ms * frs ;
	p->f[3] = pr * ps * frs ;
	p->f[6] = mr * ps * frs ;
	p->f[9] = mr * ms * frs ;
	p->f[1] = pr * ms2m3s ;
	p->f[2] = pr * ms2p3s ;
	p->f[4] = ps * mr2p3r ;
	p->f[5] = ps * mr2m3r ;
	p->f[7] = mr * ms2p3s ;
	p->f[8] = mr * ms2m3s ;
	p->f[10] = ms * mr2m3r ;
	p->f[11] = ms * mr2p3r ;
	
	if(x->w < 0.0)
		return(0) ;
	
	nt *= 2.0 ;
	r *= nt ;
	s *= nt ;
	m3s *= -s ;
	p3s *= -s ;
	m3r *= -r ;
	p3r *= -r ;
	ms2 *= 3.0 ;
	mr2 *= 3.0 ;
	
	temp = frs + pr * r ;
	p->dfdr[0] = ms * temp ;
	p->dfdr[3] = ps * temp ;
	temp = mr * r - frs ;
	p->dfdr[6] = ps * temp ;
	p->dfdr[9] = ms * temp ;
	
	p->dfdr[1] =  ms2m3s ;
	p->dfdr[2] =  ms2p3s ;
	p->dfdr[4] =  ps * (p3r + mr2) ;
	p->dfdr[5] =  ps * (m3r - mr2) ;
	p->dfdr[7] = -ms2p3s ;
	p->dfdr[8] = -ms2m3s ;
	p->dfdr[10] =  ms * (m3r - mr2) ;
	p->dfdr[11] =  ms * (p3r + mr2) ;
	
	temp = frs + ps * s ;
	p->dfds[3] = pr * temp ;
	p->dfds[6] = mr * temp ;
	temp = ms * s - frs ;
	p->dfds[0] = pr * temp ;
	p->dfds[9] = mr * temp ;
	
	p->dfds[1] =  pr * (m3s - ms2) ;
	p->dfds[2] =  pr * (p3s + ms2) ;
	p->dfds[4] =  mr2p3r ;
	p->dfds[5] =  mr2m3r ;
	p->dfds[7] =  mr * (p3s + ms2) ;
	p->dfds[8] =  mr * (m3s - ms2) ;
	p->dfds[10] =  -mr2m3r ;
	p->dfds[11] =  -mr2p3r ;

	return(0) ;
}
int	shap_1s1(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,mr,pr,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(2) ;
	
	t = 0.5 ;
	mr = 1.0 - r ;
	pr = 1.0 + r ;
	
	p->dof = 2 ;
	
	p->f[0] =  t * mr ;
	p->f[1] =  t * pr ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -0.5 ;
	p->dfdr[1] =  0.5 ;
	
	return(0) ;
}
int	SU1D_1s1(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,mr,pr,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(2) ;
	
	t = 0.5 ;
	mr = 1.0 - r ;
	pr = 1.0 + r ;
	
	p->dof = 2 ;
	
/*	Bubble Functions
	p->f[0] =  - UW * mr * pr ;
	p->f[1] =    UW * mr * pr;
*/

/*Upstream Weighting Functions*/
	p->f[0] = -0.5 * UW ;
	p->f[1] =  0.5 * UW ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] =  2 * UW * r ;
	p->dfdr[1] = -2 * UW * r;
	
	return(0) ;
}
int	shap_1s2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,mr,pr,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(3) ;
	
	t = 0.5 ;
	mr = 1.0 - r ;
	pr = 1.0 + r ;
	
	p->dof = 3 ;
	
	p->f[0] = -t * r * mr ;
	p->f[1] = pr * mr ;
	p->f[2] = t * pr * r ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -t * (1.0 - 2.0*r) ;
	p->dfdr[1] = -2.0 * r  ;
	p->dfdr[2] = t * (1.0 + 2.0*r) ;
	
	return(0) ;
}
int	shap_1s3(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,rm,rp,r3m,r3p,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(4) ;
	
	t = 9.0/144.0  ;
	rm = r - 1.0 ;
	rp = 1.0 + r ;
	r3m = 3.0 * r - 1.0 ;
	r3p = 3.0 * r + 1.0 ;
	
	p->dof = 4 ;
	
	p->f[0] = -t * rm * r3p * r3m ;
	p->f[1] = 9.0 * t * rp * rm * r3m ;
	p->f[2] = -9.0 * t * rp * rm * r3p ;
	p->f[3] = t * rp * r3p * r3m ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -t * ( r3p * r3m + 3.0*rm*r3m + 3.0*rm*r3p) ;
	p->dfdr[1] = 9.0 * t * ( 3.0*rp * rm + rm*r3m + rm*r3p)  ;
	p->dfdr[2] = -9.0 * t * ( 3.0*rp * rm + rm*r3p + rp*r3p) ;
	p->dfdr[3] = t * ( r3p * r3m + 3.0*rm*r3m + 3.0*rm*r3p) ;
	
	return(0) ;
}
int	shap_pt(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(1) ;
	
	p->dof = 1 ;
	p->f[0] =  1 ;
	p->dfdr[0] = 0.0 ;
	
	return(0) ;
}

int	shcu3_1s2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r, mr, pr, rp3, t, tt ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(3) ;
	
	t = 0.125 ;
	tt = 0.25 ;
	mr = 1.0 - r ;
	pr = 1.0 + r ;
	rp3 = 3.0 + r ;
	
	p->dof = 3 ;
	
	p->f[0] = tt * rp3 * mr ;
	p->f[1] = t * rp3 * pr ;
	p->f[2] = -t * pr * mr ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -0.5 * pr ;
	p->dfdr[1] = tt *(2.0+r)  ;
	p->dfdr[2] = tt * r ;
	
	p->d2fdr[0] = -0.5 ;
	p->d2fdr[1] = tt  ;
	p->d2fdr[2] = tt ;
	
	return(0) ;
}

int	shcu3r_1s2(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,mr,pr,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(3) ;
	
	t = 0.125 ;
	mr = 1.0 - r ;
	pr = 1.0 + r ;
	
	p->dof = 3 ;
	
	p->f[0] = t *(3.0- r) * mr ;
	p->f[1] = 2*t* (3.0- r) * pr ;
	p->f[2] = -t * pr * mr ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -2.0*t * (2.0-r) ;
	p->dfdr[1] = 4.0 * t*mr  ;
	p->dfdr[2] = 2.0*t * r ;
	
	return(0) ;
}

int	shapcu_1s3(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,rm,rp,rp3,rp5,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(4) ;
	
	t = 1.0/48.0  ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	rp3 = r + 3.0 ;
	rp5 = r + 5.0 ;
	
	p->dof = 4 ;
	
	p->f[0] =  3.0*t * rm * rp3 * rp5 ;
	p->f[1] =  t * rp * rp3 * rp5 ;
	p->f[2] = -3.0*t * rp * rm * rp5 ;
	p->f[3] = t * rp * rm * rp3 ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = -3.0 * t * ( 7.0 + (14.0 + 3.0*r)*r) ;
	p->dfdr[1] = t * ( 23.0 + r * (18.0 + 3.0*r))  ;
	p->dfdr[2] = -3.0 * t * ( 1.0 - r*(10.0 +3.0*r)) ;
	p->dfdr[3] = t * ( 1.0 -3.0*r*(2.0 + r)) ;
	
	return(0) ;
}

int	shapcur_1s3(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
	double	r,rm,rp,rp3,rp5,t ;
	
	r = x->x ;
	if( (r<-1.0) || (r>1.0) ) return(4) ;
	
	t = 1.0/48.0  ;
	rm = 1.0 - r ;
	rp = 1.0 + r ;
	rp3 = 3.0 - r ;
	rp5 = 5.0 - r;
	
	p->dof = 4 ;
	
	p->f[0] = t * rm * rp3 * rp5 ;
	p->f[1] =  3.0 * t * rp * rp3 * rp5 ;
	p->f[2] = -3.0*t * rp * rm * rp5 ;
	p->f[3] = t * rp * rm * rp3 ;
	
	if(x->w < 0.0)
		return(0) ;

	p->dfdr[0] = t * ( -23.0 + (6.0 - r)*r*3.0) ;
	p->dfdr[1] = 3.0 * t * (7.0 - r * (14.0 - 3.0*r))  ;
	p->dfdr[2] = 3.0 * t * ( 1.0 + r*(10.0 -3.0*r)) ;
	p->dfdr[3] = t * ( -1.0 -3.0*r*(2.0 - r)) ;
	
	return(0) ;
}


int	hermite(x,p)
struct gausspts		*x ;
struct shapefuncs	*p ;

{
/*
	double	r,s,mr,ms,pr,ps,t,tr,ts ;
	
	r = x->x ;
	s = x->y ;
	if( (r<-1.00001) || (r>1.00001) || (s<-1.00001) || (s>1.00001) ) return(4) ;
	
	t = 1.0/16.0 ;
	tr = 2.0 * r ;
	ts = 2.0 * s ;
	mr = r - 1.0 ;
	ms = s - 1.0 ;
	pr = 1.0 + r ;
	ps = 1.0 + s ;
	
	p->dof = 12 ;
	
	p->f[0] = t * pr * pr * (r - 2.0) * ms * ms * (-s - 2.0) ;
	p->f[3] = t * pr * pr * (r - 2.0) * ps * ps * (s - 2.0) ;
	p->f[6] = t * mr * mr * (-r - 2.0) * ps * ps * (s - 2.0) ;
	p->f[9] = t * mr * mr * (-r - 2.0) * ms * ms * (-s - 2.0) ;
	p->f[1] = -t * pr * pr * (r - 1.0) * ms * ms * (-s - 2.0) ;
	p->f[4] = -t * pr * pr * (r - 1.0) * ps * ps * (s - 2.0) ;
	p->f[7] = t * mr * mr * (-r - 1.0) * ps * ps * (s - 2.0) ;
	p->f[10] = t * mr * mr * (-r - 1.0) * ms * ms * (-s - 2.0) ;
	p->f[2] = t * pr * pr * (r - 2.0) * ms * ms * (-s - 1.0) ;
	p->f[5] = -t * pr * pr * (r - 2.0) * ps * ps * (s - 1.0) ;
	p->f[8] = -t * mr * mr * (-r - 2.0) * ps * ps * (s - 1.0) ;
	p->f[11] = t * mr * mr * (-r - 2.0) * ms * ms * (-s - 1.0) ;
	
	if(x->w < 0.0)
		return(0) ;
	
	p->dfdr[0] = t * pr * (3*r - 3.0) * ms * ms * (-s - 2.0) ;
	p->dfdr[3] = t * pr * (3*r - 3.0) * ps * ps * (s - 2.0) ;
	p->dfdr[6] = t * mr * (-3*r - 3.0) * ps * ps * (s - 2.0) ;
	p->dfdr[9] = t * mr * (-3*r - 3.0) * ms * ms * (-s - 2.0) ;
	p->dfdr[1] = -t * pr * (3*r - 1.0) * ms * ms * (-s - 2.0) ;
	p->dfdr[4] = -t * pr * (3*r - 1.0) * ps * ps * (s - 2.0) ;
	p->dfdr[7] = t * mr * (-3*r - 1.0) * ps * ps * (s - 2.0) ;
	p->dfdr[10] = t * mr * (-3*r - 1.0) * ms * ms * (-s - 2.0) ;
	p->dfdr[2] = t * pr * (3*r - 3.0) * ms * ms * (-s - 1.0) ;
	p->dfdr[5] = -t * pr * (3*r - 3.0) * ps * ps * (s - 1.0) ;
	p->dfdr[8] = -t * mr * (-3*r - 3.0) * ps * ps * (s - 1.0) ;
	p->dfdr[11] = t * mr * (-3*r - 3.0) * ms * ms * (-s - 1.0) ;
	
	p->dfds[0] = t * pr * pr * (r - 2.0) * ms * (-3*s - 3.0) ;
	p->dfds[3] = t * pr * pr * (r - 2.0) * ps * (3*s - 3.0) ;
	p->dfds[6] = t * mr * mr * (-r - 2.0) * ps * (3*s - 3.0) ;
	p->dfds[9] = t * mr * mr * (-r - 2.0) * ms * (-3*s - 3.0) ;
	p->dfds[1] = -t * pr * pr * (r - 1.0) * ms * (-3*s - 3.0) ;
	p->dfds[4] = -t * pr * pr * (r - 1.0) * ps * (3*s - 3.0) ;
	p->dfds[7] = t * mr * mr * (-r - 1.0) * ps * (3*s - 3.0) ;
	p->dfds[10] = t * mr * mr * (-r - 1.0) * ms * (-3*s - 3.0) ;
	p->dfds[2] = t * pr * pr * (r - 2.0) * ms * (-3*s - 1.0) ;
	p->dfds[5] = -t * pr * pr * (r - 2.0) * ps * (3*s - 1.0) ;
	p->dfds[8] = -t * mr * mr * (-r - 2.0) * ps * (3*s - 1.0) ;
	p->dfds[11] = t * mr * mr * (-r - 2.0) * ms * (-3*s - 1.0) ;
	
	*/
	
	return(0) ;
}
