#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

static const double min_data_fraction = 0.1;

// Implements R's type=8
double c_quantile(double* data, const int n, const double quantile, const bool sorted=false) {
  if(n == 0 || quantile < 0 || quantile > 1) {
    return(R_NaReal);
  }

  if(n == 1)
    return(data[0]);
  
  // Constants for quantiles. Can be modified if needed.
  static const double a = 1.0/3.0;
  static const double b = 1.0/3.0;
  double* const data_end = data + n;
  
  static const double fuzz = 4 * numeric_limits<double>::epsilon();
  const double nppm = a + quantile * (n + 1 - a - b) - 1;
  const double j = floor(nppm + fuzz);
  // Variance from R: Should probably be <= not < here.
  const double h = (fabs(nppm - j) <= fuzz) ? 0 : nppm - j;
  double* right_elem = max(data, min(data + (int)j + 1, data_end - 1));
  double* left_elem = max(data, min(data + (int)j, data_end - 1));
  
  if(h == 1) {
    if(!sorted)
      nth_element(data, right_elem, data_end);
    return(*right_elem);
  } else {
    // No guarantee that 2nd nth_element call will preserve order such that the pointer used by the 1st call still points to the same thing; so store the result before calling nth_element again.
    if(!sorted)
      nth_element(data, left_elem, data_end);
    const double left = *left_elem;
    if(h == 0) {
      return(left);
    } else {
      if(!sorted)
	nth_element(data, right_elem, data_end);
      const double right = *right_elem;
      return((1 - h) * left + h * right);
    }
  }
}

class DatYrTuple {
public:
  DatYrTuple() { dat = nan(""); yr = -1; }
  DatYrTuple(double dat, int yr, int day): dat(dat), yr(yr), day(day) { }
  bool operator < (const DatYrTuple& d) const { return dat < d.dat; }
  double dat;
  int yr;
  int day;
};

class IdxDupflagPair {
public:
  IdxDupflagPair() { idx = -1; dup = false; }
  bool operator < (const IdxDupflagPair& d) const { return idx < d.idx; }
  int idx;
  bool dup;
};

class IdxDayPair {
public:
  IdxDayPair() { idx = -1; day = -1; }
  IdxDayPair(int idx, int day): idx(idx), day(day) {}
  int idx;
  int day;
};

class ClimdexBootstrapper {
public:
  const int win_size, nyr, dpy;
  const double* dat;
  const int* notna_map;
  const int half_win;

  ClimdexBootstrapper(const double* dat, const int* notna_map, const int win_size, const int nyr, const int dpy): win_size(win_size), nyr(nyr), dpy(dpy), dat(dat), notna_map(notna_map), half_win(win_size / 2) { }
  
  // Extracts an n-day window into the data, removing NAs and generating a 2-tuple.
  vector<DatYrTuple> extract_window_with_year(const int day) {
    const int min_day = day - half_win;
    const int max_day = day + half_win;
    vector<DatYrTuple> out_dat;
    out_dat.reserve(win_size * nyr);

    for(int yr = 0; yr < nyr; ++yr) {
      const int yr_base = yr * dpy;
      for(int day_idx = min_day; day_idx <= max_day; ++day_idx) {
	const int wday = (day_idx + dpy) % dpy;
	const int idx = wday + yr_base;
	if(notna_map[idx])
	  out_dat.push_back(DatYrTuple(dat[idx], yr, wday));
      }
    }

    return out_dat;
  }

  vector<vector<IdxDayPair> > create_yrs_index(const vector<DatYrTuple>& sorted_in) {
    const int max_elems = sorted_in.size() / nyr;
    vector<IdxDayPair> temp(max_elems);
    temp.resize(0);
    vector<vector<IdxDayPair> > yidx(nyr, temp);
    int idx = 0;
    for(vector<DatYrTuple >::const_iterator i = sorted_in.begin(); i != sorted_in.end(); ++i, ++idx) {
      const int yr = (*i).yr;
      const int day = (*i).day;
      // Should raise an exception or something here.
      if(yr >= 0 && yr < nyr)
	yidx[yr].push_back(IdxDayPair(idx, day));
    }
    return yidx;
  }

  vector<IdxDupflagPair> get_index_tuples(const vector<vector<IdxDayPair> >& yrs_index, const int rm_year, const int dup_year, const int day) {
    vector<IdxDupflagPair> df_pairs(yrs_index[rm_year].size() + yrs_index[dup_year].size() + half_win + 1);
  
    int df_idx = 0;
    for(vector<IdxDayPair>::const_iterator i = yrs_index[rm_year].begin(); i != yrs_index[rm_year].end(); ++i, ++df_idx)
      df_pairs[df_idx].idx = (*i).idx;

    for(vector<IdxDayPair>::const_iterator i = yrs_index[dup_year].begin(); i != yrs_index[dup_year].end(); ++i, ++df_idx) {
      df_pairs[df_idx].dup = true;
      df_pairs[df_idx].idx = (*i).idx;
    }

    df_pairs.resize(df_idx);
    sort(df_pairs.begin(), df_pairs.end());
    return df_pairs;
  }

  void replace_data_year(const vector<double>& in_data, vector<double>& out, const vector<vector<IdxDayPair> >& yrs_index, const int rm_year, const int dup_year, const int day) {
    const vector<IdxDupflagPair>& df_pairs = get_index_tuples(yrs_index, rm_year, dup_year, day);
    int numdup = 0, numnotdup = 0;
    for(vector<IdxDupflagPair>::const_iterator i = df_pairs.begin(); i != df_pairs.end(); ++i) { numdup += (int)((*i).dup); numnotdup += (int)(!(*i).dup); }
    out.resize(in_data.size() + numdup - numnotdup);

    int last_in_idx = -1, next_out_idx = 0;
    for(vector<IdxDupflagPair>::const_iterator i = df_pairs.begin(); i != df_pairs.end(); ++i) {
      const bool isdup = (*i).dup;
      const int idx = (*i).idx;
      const int nonspecial_block_length = idx - last_in_idx - 1;

      if(isdup) {
	copy(&in_data[last_in_idx + 1], &in_data[idx + 1], &out[next_out_idx]);
	out[next_out_idx + nonspecial_block_length + 1] = in_data[idx];
	next_out_idx += 2;
      } else {
	copy(&in_data[last_in_idx + 1], &in_data[idx], &out[next_out_idx]);
      }
      next_out_idx += nonspecial_block_length;
      last_in_idx = idx;
    }

    if((unsigned int)(last_in_idx + 1) < in_data.size()) {
      copy(&in_data[last_in_idx + 1], &in_data[in_data.size()], &out[next_out_idx]);
    }
  }

  vector<double> get_data_only(const vector<DatYrTuple>& in_dat) {
    vector<double> out_dat(in_dat.size());
    const unsigned int dat_size = in_dat.size();
    for(unsigned int i = 0; i < dat_size; ++i)
      out_dat[i] = in_dat[i].dat;
    return out_dat;
  }
};

RcppExport SEXP c_quantile2(SEXP data_, SEXP quantile_) {
  const NumericVector q(quantile_);
  const NumericVector data(data_);
  const int n = data.size();
  const int nq = q.size();
  NumericVector res(nq);
  
  for(int i = 0; i < nq; ++i)
    res[i] = c_quantile(const_cast<NumericVector&>(data).begin(), n, q[i]);
  
  return res;
}

RcppExport SEXP running_quantile_windowed_bootstrap(SEXP data_, SEXP n_, SEXP q_, SEXP dpy_, SEXP min_fraction_) {
  const int win_size = as<int>(n_);
  const int days_per_year = as<int>(dpy_);
  const double min_fraction = as<double>(min_fraction_);
  const NumericVector q(q_);
  const NumericVector data(data_);
  const int nq = q.size();
  const int data_length = data.size();
  
  // Implicit pseudo-floor (nearest number to zero) by integer division
  const int num_years = (int)ceil((double)data_length / (double)days_per_year);
  const int slab_size = win_size * num_years;
  
  const int day_mul = 1;
  const int rmyr_mul = days_per_year;
  const int dupyr_mul = days_per_year * num_years;
  const int q_mul = days_per_year * num_years * (num_years - 1);

  NumericVector quantiles(nq * days_per_year * num_years * (num_years - 1));
  LogicalVector notna_map = !is_na(data);
  ClimdexBootstrapper bs(data.begin(), notna_map.begin(), win_size, num_years, days_per_year);
  
  vector<double> rep_dat(slab_size);
  for(int day = 0; day < days_per_year; ++day) {
    vector<DatYrTuple> win_tuples(bs.extract_window_with_year(day));
    std::sort(win_tuples.begin(), win_tuples.end());
    const vector<vector<IdxDayPair> >& yrs_index = bs.create_yrs_index(win_tuples);
    const vector<double>& win_dat = bs.get_data_only(win_tuples);

    const int off_day = day * day_mul;
    for(int rm_year = 0; rm_year < num_years; ++rm_year) {
      const int off_day_rmyr = off_day + rm_year * rmyr_mul;
      int dup_idx = 0;
      for(int dup_year = 0; dup_year < num_years; ++dup_year) {
	if(dup_year != rm_year) {
	  bs.replace_data_year(win_dat, rep_dat, yrs_index, rm_year, dup_year, day);
	  const int off_day_rmyr_dupyr = off_day_rmyr + dup_idx * dupyr_mul;
	  if(((double)rep_dat.size() / (double)slab_size) < min_fraction)
	    for(int q_idx = 0; q_idx < nq; ++q_idx)
	      quantiles[off_day_rmyr_dupyr + q_idx * q_mul] = R_NaReal;
	  else
	    for(int q_idx = 0; q_idx < nq; ++q_idx)
	      quantiles[off_day_rmyr_dupyr + q_idx * q_mul] = c_quantile(&rep_dat[0], rep_dat.size(), q[q_idx], true);
	  ++dup_idx;
	}
      }
    }
  }

  return quantiles;
}

// Expects data in date sequence
//void running_quantile_windowed_365day(const double* data, double* quantiles, const int* n, const double* q, const int* data_length, const int* num_quantiles) {
RcppExport SEXP running_quantile_windowed(SEXP data_, SEXP n_, SEXP q_, SEXP dpy_, SEXP min_fraction_) {
  const int win_size = as<int>(n_);
  const int dpy = as<int>(dpy_);
  const double min_fraction = as<double>(min_fraction_);
  const NumericVector q(q_);
  const NumericVector data(data_);
  const int nq = q.size();
  const int data_length = data.size();
  
  // Implicit pseudo-floor (nearest number to zero) by integer division
  const int half_win = win_size / 2;
  const int num_years = (int)ceil((double)data_length / (double)dpy);
  const int slab_size = win_size * num_years;

  NumericVector quantiles(nq * dpy);
  
  // This data will be 2 dimensional, row major, with the major dimension being the day of the year
  double* buf = new double[slab_size];
  
  LogicalVector notna_map = !is_na(data);
  
  for(int day = 0; day < dpy; ++day) {
    // Fill buffer with data with NA removed
    int count = 0;
    const int min_day = day - half_win;
    const int max_day = day + half_win;
    for(int yr = 0; yr < num_years; ++yr) {
      const int yr_base = yr * dpy;
      for(int day_idx = min_day; day_idx <= max_day; ++day_idx) {
	const int wday = (day_idx + dpy) % dpy;
	const int idx = wday + yr_base;
	if(notna_map[idx])
	  buf[count++] = data[idx];
      }
    }

    // Quantiles on said buffer
    const int nq_day = nq * day;
    if(((double)count / (double)slab_size) < min_fraction)
      for(int q_idx = 0; q_idx < nq; ++q_idx)
	quantiles[nq_day + q_idx] = R_NaReal;
    else
      for(int q_idx = 0; q_idx < nq; ++q_idx)
	quantiles[nq_day + q_idx] = c_quantile(buf, count, q[q_idx]);
  }
  delete[] buf;
  
  return quantiles;
}
