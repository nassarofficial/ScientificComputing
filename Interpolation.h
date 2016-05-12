#ifndef INTERPOLATION_H
#define INTERPOLATION_H


class Interpolation
{
    public:
        /**
        @param x[]  all points
        @param y[]  corresponding f(x)
        @param a    point to be interpolated
        @param n size of the input array
        @return interpolated point
        */
        double NewtInt(double x[], double y[], double a, int n);
        /**
        @param x[]  all points
        @param y[]  corresponding f(x)
        @param a    point to be interpolated
        @param n size of the input array
        @return interpolated point
        */
        double Spline(double x[], double y[], double a, int n);
    protected:

    private:
    	/**
    	@param x[] array of interval limits
    	@param a   point to be interpolated
    	@param n   number of elements in array x[]
    	
    	Finds the interval of the interpolated point to calculate the interpolant
    	*/
    	double findingInterval(double x[], double a, int n);
    	/**
    	@param a[] array to be printed
    	@param n   size of the array
    	*/
		void printArray(double a[], int n);
};

#endif // INTERPOLATION_H
