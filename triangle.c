#include <stdio.h>
#include <math.h>

#define CONSTRAINT_SOLVER_MAX_ITERATIONS 1000
#define CONSTRAINT_SLOPE_INITIAL_GUESS 1.0
#define CONSTRAINT_ERROR_TRANSFER_RATIO 0.1
#define CONSTRAINT_ERROR_THRESHOLD 0.0000000001

#define RADIANS_TO_DEGREES (180 / M_PI)
#define FIFTY_DEGREES_AS_RADIANS (50 * M_PI / 180)

typedef long double real;

typedef struct {
	real x,y;
} vec2;
vec2 vec2Subtract(const vec2 a, const vec2 b){
	const vec2 r = {
		.x = a.x - b.x,
		.y = a.y - b.y
	};
	return r;
}
real vec2Dot(const vec2 a, const vec2 b){
	return a.x*b.x + a.y*b.y;
}
real vec2Magnitude(const vec2 v){
	return sqrtl(v.x*v.x + v.y*v.y);
}
real vec2Angle(const vec2 a, const vec2 b, const vec2 c){
	const vec2 ab = vec2Subtract(b, a);
	const vec2 cb = vec2Subtract(b, c);
	return acosl(vec2Dot(ab, cb) / (vec2Magnitude(ab)*vec2Magnitude(cb)));
}

real ConstraintTest(const real m){

	// Calculates the current value of <BKL.

	// BK is some line from the point (1, 0) to K
	// with slope m, given by the following equation:
	// y = m(k - K.x).
	//
	// Where K lies on the circle:
	// y^2 = 2kK.x - K.x^2.
	//
	// And k = |CB| = |BK| = 1.
	//
	// Thus, K has the following form:
	// K.x = 1 - 1/sqrt(m^2 + 1)
	// K.y = m/sqrt(m^2 + 1)
	//
	// From this, we can calculate the slope of CK:
	// K.y / K.x = (sqrt(m^2 + 1) + 1)/m
	//
	// Which has magnitude sqrt(2.0*K.x).
	const real msqrt = sqrtl(m*m + 1.0);
	const vec2 K = {
		.x = 1.0 - 1.0/sqrtl(m*m + 1.0),
		.y = m/sqrtl(m*m + 1.0)
	};
	const real slopeK = (msqrt+1.0)/m;
	const real j = sqrtl(2.0*K.x);

	// The point A has y coordinate 0.5, as the
	// triangle is an isosceles. The intersection
	// of the infinite line CK with y = 1/2 gives
	// the x coordinate, and the magnitude of the
	// line segment CA - j gives the length i.
	const vec2 A = {
		.x = 0.5,
		.y = slopeK*0.5
	};
	const real magnitudeA = vec2Magnitude(A);
	const real i = magnitudeA - j;

	// From here, we can find point L, which is
	// i from B towards A:
	// L = B + BA/|A| * i
	const vec2 L = {
		.x = 1.0 - i*0.5/magnitudeA,
		.y = i*A.y/magnitudeA
	};

	// Finally, we can calculate the angle BKL.
	// The constraint on this angle is 50 degrees.
	const vec2 B = {
		.x = 1.0,
		.y = 0.0
	};
	return vec2Angle(B, K, L);

}

real alpha(const real m){

	// Calculate the final value for alpha
	// based on the slope of line segment BK.

	// First, calculate j, the magnitude of
	// the line segment CK.
	const real j = sqrtl(2.0 - 2.0/sqrtl(m*m + 1.0));

	// Using the lengths of the triangle BCK,
	// calculate the angle <BCK.
	// <BCK = arccos(j/2)
	const real angleBCK = acosl(j*0.5) * RADIANS_TO_DEGREES;

	// The sum of a 3-sided polygon's angles
	// add up to 180 degrees. Therefore:
	// alpha = 180 - 2(<BCK)
	return 180.0 - 2.0*angleBCK;

}

int main(int argc, char **argv){

	// The following can all be graphed on Desmos
	// for a visual representation of the algorithm.
	//
	// Equations:
	//
	// K Circle : x^2 + y^2 = 2x
	// BK : y = m(1 - x)
	// CA : y = ((sqrt(m^2 + 1) + 1)/m)x
	// BA : y = ((sqrt(m^2 + 1) + 1)/m)(1 - x)
	//
	// Coordinates of each point on the triangle
	// in terms of the slope of BK:
	//
	// B = (1, 0)
	// C = (0, 0)
	// K = (1 - 1/sqrt(m^2 + 1), m/sqrt(m^2 + 1))
	// A = (1/2, (sqrt(m^2 + 1) + 1)/2m)
	// L = B + BA/|A| * i
	//
	// Side lengths:
	//
	// One dash     : i = |A| - j
	// Two dashes   : j = sqrt(2K.x)
	// Three dashes : k = 1

	unsigned int i = CONSTRAINT_SOLVER_MAX_ITERATIONS;
	real m = CONSTRAINT_SLOPE_INITIAL_GUESS;
	real error;

	// Iteratively solve the constraint.
	while(i > 0){
		const real theta = ConstraintTest(m);
		error = theta - FIFTY_DEGREES_AS_RADIANS;
		if(fabsl(error) <= CONSTRAINT_ERROR_THRESHOLD){
			// The error is within some threshold
			// that we're happy with. Exit the loop.
			break;
		}
		// Modify the slope according to feedback
		// from the error. As it turns out, simply
		// adding the product of the error and some
		// transfer ratio gets the slope closer to
		// the final solution.
		m += error * CONSTRAINT_ERROR_TRANSFER_RATIO;
		--i;
	}

	// Using our current triangle state, we
	// can attempt to solve for alpha.
	printf("Total iterations = %d\n", CONSTRAINT_SOLVER_MAX_ITERATIONS-i);
	printf("Slope of Line Segment BK = %.20Lf\n", m);
	printf("Angle BKL Error = %.20Lf\n", error);
	printf("Alpha = %.20Lf\n", alpha(m));
	return 0;

}