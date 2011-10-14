#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>
#include <math.h>
#include <limits>
#include <algorithm>

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
    // No guarantee that 2nd nth_element call will preserve order such that the pointer used by the 1st call still points to the same thing; so store the result before calling nth_element again.
    std::nth_element(data, left_elem, data_end);
    const double left = *left_elem;
    if(h == 0) {
      return(left);
    } else {
      std::nth_element(data, right_elem, data_end);
      const double right = *right_elem;
      return((1 - h) * left + h * right);
    }
  }
}

extern "C" {
  // Expects data in date sequence
  //void running_quantile_windowed_365day(const double* data, double* quantiles, const int* n, const double* q, const int* data_length, const int* num_quantiles) {
  SEXP running_quantile_windowed(SEXP data, SEXP n, SEXP q, SEXP dpy) {
    PROTECT(n = coerceVector(n, INTSXP));
    PROTECT(dpy = coerceVector(dpy, INTSXP));
    
    const int win_size = INTEGER(n)[0];
    
    // Implicit pseudo-floor (nearest number to zero) by integer division
    const int win_border = win_size / 2;
    const int true_data_length = length(data) - 2 * win_border;
    const int days_per_year = INTEGER(dpy)[0];
    const int num_years = (int)ceil((double)true_data_length / (double)days_per_year);
    const int repeated_data_size = num_years * (days_per_year + 2 * win_border);
    const int nq = length(q);
    const double* q_ptr = REAL(q);
    const double* data_ptr = REAL(data);

    SEXP quantiles = allocVector(REALSXP, nq * days_per_year);
    double* quantiles_ptr = REAL(quantiles);

    // This data will be 2 dimensional, row major, with the major dimension being the day of the year
    double* buf = new double[num_years * win_size];

    // Comment for preservation of sanity...
    // The input data does not stay on day 'day'; it starts on day 'day' - win_border. This gets around the problem of an off-by-two error.
    for(int day = 0; day < days_per_year; ++day) {
      int count = 0;
      // Fill buffer with data with NA removed
      for(int yr = 0; yr < num_years; ++yr) {
	const double* yroff_ptr = &data_ptr[yr * days_per_year];
	for(int winday = day; winday < day + win_size; ++winday) {
	  const double d = yroff_ptr[winday];
	  if(!ISNA(d))
	    buf[count++] = d;
	}
      }
      // Quantiles on said buffer
      for(int q_idx = 0; q_idx < nq; ++q_idx) {
	quantiles_ptr[nq * day + q_idx] = c_quantile(buf, count, q_ptr[q_idx]);
      }
    }
    delete[] buf;

    UNPROTECT(2);
    return(quantiles);
  }

  void R_init_mylib(DllInfo* info) {
    R_CallMethodDef callMethods[] = { 
      {"running_quantile_windowed", (DL_FUNC) &running_quantile_windowed, 4 },
      { NULL, NULL, 0 } 
    };
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    /*R_RegisterCCallable("pcicspatial", "get_coverage", get_coverage);*/
  }
  
  void R_unload_mylib(DllInfo* info) {
  }
}
