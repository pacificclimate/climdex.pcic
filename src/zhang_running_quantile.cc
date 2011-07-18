#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>
#include <math.h>
#include <limits>
#include <algorithm>

extern "C" {
  // Implements R's type=8
  double c_quantile(double* data, const int n, const double quantile) {
    if(n < 2) {
      return(R_NaReal);
    }

    // Constants for quantiles. Can be modified if needed.
    const double a = 1/3;
    const double b = 1/3;
    double* data_end = data + n;

    const double fuzz = 4 * std::numeric_limits<double>::epsilon();
    const double nppm = a + quantile * (n + 1 - a - b) - 1;
    const int j = floor(nppm + fuzz);
    // Variance from R: Should probably be <= not < here.
    const double h = (fabs(nppm - (double)j) <= fuzz) ? 0 : nppm - (double)j;
    double* right_elem = std::max(data, std::min(data + j + 1, data_end - 1));;
    double* left_elem = std::max(data, std::min(data + j, data_end - 1));;

    if(h == 1) {
      std::nth_element(data, right_elem, data_end);
      return(*right_elem);
    } else {
      std::nth_element(data, left_elem, data_end);
      if(h == 0) {
	return(*left_elem);
      } else {
	std::nth_element(data, right_elem, data_end);
	return((1 - h) * (*left_elem) + h * (*right_elem));
      }
    }
  }

  // Expects data in date sequence
  void running_quantile_windowed_365day(const double* data, double* quantiles, const int* n, const double* q, const int* data_length, const int* num_quantiles) {
    const int window = floor(*n / 2);
    const int true_data_length = *data_length - 2 * window;
    const int days_per_year = 365;
    const int num_years = ceil(true_data_length / days_per_year);
    const int repeated_data_size = num_years * (days_per_year + 2 * window);
    
    // This data will be 2 dimensional, row major, with the major dimension being the day of the year
    double* buf = new double[num_years * (*n)];
    
    for(int day = 0; day < days_per_year; ++day) {
      int count = 0;
      // Fill buffer with data with NA removed
      for(int winday = day; winday < day + *n; ++winday) {
	for(int yr = 0; yr < num_years; ++yr) {
	  const double d = data[yr * days_per_year + winday];
	  if(d != R_NaReal)
	    buf[count++] = d;
	}
      }
      // Quantiles on said buffer
      for(int q_idx = 0; q_idx < *num_quantiles; ++q_idx) {
	quantiles[*num_quantiles * day + q_idx] = c_quantile(buf, count, q[q_idx]);
      }
    }
    delete[] buf;
  }

  void R_init_mylib(DllInfo* info) {
    R_CMethodDef cMethods[] = { 
      {"running_quantile_windowed_365day", (DL_FUNC) &running_quantile_windowed_365day, 6, (R_NativePrimitiveArgType[]) { REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, }, (R_NativeArgStyle[]) { R_ARG_IN, R_ARG_OUT, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN } }, 
      { NULL, NULL, 0 } 
    };
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    /*R_RegisterCCallable("pcicspatial", "get_coverage", get_coverage);*/
  }
  
  void R_unload_mylib(DllInfo* info) {
  }
}
