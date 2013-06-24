#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// Implements R's type=8
double c_quantile(double* data, const int n, const double quantile, bool sorted=false) {
  if(n < 2 || quantile < 0 || quantile > 1) {
    return(R_NaReal);
  }
  
  // Constants for quantiles. Can be modified if needed.
  const double a = 1.0/3.0;
  const double b = 1.0/3.0;
  double* data_end = data + n;
  
  const double fuzz = 4 * numeric_limits<double>::epsilon();
  const double nppm = a + quantile * (n + 1 - a - b) - 1;
  const int j = (int)floor(nppm + fuzz);
  // Variance from R: Should probably be <= not < here.
  const double h = (fabs(nppm - (double)j) <= fuzz) ? 0 : nppm - (double)j;
  double* right_elem = max(data, min(data + j + 1, data_end - 1));
  double* left_elem = max(data, min(data + j, data_end - 1));
  
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

// NOTE: Takes data with floor(win_size / 2) elements attached to beginning and end.
// Extracts an n-day window into the data, removing NAs and generating a 2-tuple.
vector<DatYrTuple> extract_window_with_year(const double* dat, const int* notna_map, const int day, const int win_size, const int nyr, const int dpy) {
  const int min_day = day;
  const int max_day = day + win_size;
  const int half_win = win_size / 2;
  vector<DatYrTuple> out_dat;
  out_dat.reserve(win_size * nyr + half_win);

  for(int yr = 0; yr < nyr; ++yr) {
    const int yr_base = yr * dpy;
    for(int day_idx = min_day; day_idx < max_day; ++day_idx) {
      const int idx = day_idx + yr_base;
      if(notna_map[idx]) {
	// This is kind of ugly. But it works.
	const int actual_yr = (int)floor((float)(idx - half_win) / (float)dpy);
	out_dat.push_back(DatYrTuple(dat[idx], actual_yr, (day_idx - half_win) % dpy));
      }
    }
  }

  // In the case of data at the edge, we need to copy additional values for replacement purposes
  // and omit them as needed later.
  if(day < half_win) {
    const int min_idx = nyr * dpy + day;
    const int max_idx_plusone = nyr * dpy + half_win;
    for(int idx = min_idx; idx < max_idx_plusone; ++idx) {
      if(notna_map[idx]) {
	const int actual_yr = (int)floor((float)(idx - half_win) / (float)dpy);
	out_dat.push_back(DatYrTuple(dat[idx], actual_yr, (idx - half_win) % dpy));
      }
    }
    // For cases where the day will wrap around the end, if we're going to duplicate the 1st year,
    // we will need the first 1 or 2 days sorted into the array too. Because we can only deal with sorted data later,
    // that data has to be added regardless and removed as needed.
  } else if(day > dpy - half_win - 1) {
    const int min_idx = (dpy - 1) - day + half_win;
    const int max_idx_plusone = half_win + half_win;
    for(int idx = min_idx; idx < max_idx_plusone; ++idx) {
      if(notna_map[idx]) {
	const int actual_yr = (int)floor((float)(idx - half_win) / (float)dpy);
	out_dat.push_back(DatYrTuple(dat[idx], actual_yr, (idx - half_win) % dpy));
      }
    }
  }

  return out_dat;
}

RcppExport SEXP extract_win_yr(SEXP data_, SEXP notna_map_, SEXP day_, SEXP win_size_, SEXP nyr_, SEXP dpy_) {
  const NumericVector data(data_);
  const LogicalVector notna_map(notna_map_);
  const int day = as<int>(day_);
  const int nyr = as<int>(nyr_);
  const int win_size = as<int>(win_size_);
  const int dpy = as<int>(dpy_);

  const vector<DatYrTuple>& dyt = extract_window_with_year(data.begin(), notna_map.begin(), day, win_size, nyr, dpy);

  NumericVector dat(dyt.size());
  IntegerVector yr(dyt.size());
  IntegerVector day_no(dyt.size());
  List ret;

  int idx = 0;
  for(vector<DatYrTuple>::const_iterator i = dyt.begin(); i != dyt.end(); ++i, ++idx) {
    dat[idx] = (*i).dat;
    yr[idx] = (*i).yr;
    day_no[idx] = (*i).day;
  }
  
  ret["dat"] = dat;
  ret["yr"] = yr;
  ret["day"] = yr;

  return ret;
}

vector<vector<IdxDayPair> > create_yrs_index(const vector<DatYrTuple>& sorted_in, const int nyr) {
  const int max_elems = sorted_in.size() / nyr;
  vector<IdxDayPair> temp(max_elems);
  temp.resize(0);
  vector<vector<IdxDayPair> > yidx(nyr, temp);
  int idx = 0;
  for(vector<DatYrTuple >::const_iterator i = sorted_in.begin(); i != sorted_in.end(); ++i, ++idx) {
    const int yr = (*i).yr;
    const int day = (*i).day;
    if(yr >= 0 && yr < nyr)
      yidx[yr].push_back(IdxDayPair(idx, day));
  }
  return yidx;
}

/*
RcppExport SEXP create_yrs_idx(SEXP sorted_in_, SEXP nyr_) {
  List sorted_in(sorted_in_);
  const int nyr = as<int>(nyr_);
  const NumericVector& dat = sorted_in["dat"];
  const IntegerVector& yr = sorted_in["yr"];
  const IntegerVector& day_no = sorted_in["day"];
  const int datsize = dat.size();
  vector<DatYrTuple> d(datsize);

  for(int i = 0; i < datsize; ++i) {
    d[i].dat = dat[i];
    d[i].yr = yr[i];
    d[i].day = day_no[i];
  }

  return wrap(create_yrs_index(d, nyr));
}
*/

vector<IdxDupflagPair> get_index_tuples(const vector<vector<IdxDayPair> >& yrs_index, const int rm_year, const int dup_year, const int nyr, const int day, const int win_size, const int dpy) {
  const int half_win = win_size / 2;
  const int day_check_lower = half_win;
  const int day_check_upper = dpy - half_win - 1;
  const int min_day = 0;
  const int max_day = dpy - 1;
  vector<IdxDupflagPair> id(yrs_index[rm_year].size() + yrs_index[dup_year].size() + half_win + 1);
  
  const bool need_rm_first_year_extras = day > day_check_upper && rm_year != 0 && dup_year != 0;
  const bool need_rm_last_year_extras = day < day_check_lower && rm_year != (nyr - 1) && dup_year != (nyr - 1);
  const bool need_dup_omit = (day > day_check_upper && (dup_year == 0 || rm_year == 0)) || (day < day_check_lower && (dup_year == (nyr - 1) || rm_year == (nyr - 1)));
  int id_idx = 0;
  
  for(vector<IdxDayPair>::const_iterator i = yrs_index[rm_year].begin(); i != yrs_index[rm_year].end(); ++i, ++id_idx)
    id[id_idx].idx = (*i).idx;

  // This runs if we need to ensure that only stuff within the OK range is copied (the other stuff is simply preserved. Viva le WTF!
  if(need_dup_omit) {
    const int ok_range_min = max(min_day, day - half_win);
    const int ok_range_max = min(max_day, day + half_win);
    for(vector<IdxDayPair>::const_iterator i = yrs_index[dup_year].begin(); i != yrs_index[dup_year].end(); ++i) {
      const int cur_day = (*i).day;
      if(cur_day <= ok_range_max && cur_day >= ok_range_min) {
	id[id_idx].dup = true;
	id[id_idx].idx = (*i).idx;
	++id_idx;
      }
    }
  } else {
    for(vector<IdxDayPair>::const_iterator i = yrs_index[dup_year].begin(); i != yrs_index[dup_year].end(); ++i, ++id_idx) {
      id[id_idx].dup = true;
      id[id_idx].idx = (*i).idx;
    }
  }

  // This clears out the duplicate values we throw in at the edges so that we have enough data to work with if we don't need them.
  if(need_rm_first_year_extras) {
    for(vector<IdxDayPair>::const_iterator i = yrs_index[0].begin(); i != yrs_index[0].end(); ++i) {
      const int cur_day = (*i).day;
      if(cur_day < day_check_lower) {
	id[id_idx].idx = (*i).idx;
	++id_idx;
      }
    }
  }
  if(need_rm_last_year_extras) {
    for(vector<IdxDayPair>::const_iterator i = yrs_index[nyr - 1].begin(); i != yrs_index[nyr - 1].end(); ++i) {
      const int cur_day = (*i).day;
      if(cur_day > day_check_upper) {
	id[id_idx].idx = (*i).idx;
	++id_idx;
      }
    }
  }

  id.resize(id_idx);
  sort(id.begin(), id.end());
  return id;
}

/*
RcppExport SEXP get_idx_pairs(SEXP yrs_index_, SEXP rm_year_, SEXP dup_year_) {
  const int rm_year = as<int>(rm_year_);
  const int dup_year = as<int>(dup_year_);
  List yrs_index_list(yrs_index_);
  vector<vector<int> > yrs_index(yrs_index_list.size());

  int yrid = 0;
  for(List::iterator i = yrs_index_list.begin(); i != yrs_index_list.end(); ++i, ++yrid) {
    const IntegerVector& d = (*i);
    yrs_index[yrid].resize(d.size());
    std::copy(d.begin(), d.end(), yrs_index[yrid].begin());
  }

  const vector<IdxDupflagPair>& idxt = get_index_tuples(yrs_index, rm_year, dup_year);
  LogicalVector dup(idxt.size());
  IntegerVector idx(idxt.size());
  List ret;

  for(int i = 0; i < dup.size(); ++i) {
    dup[i] = idxt[i].dup;
    idx[i] = idxt[i].idx;
  }

  ret["dup"] = dup;
  ret["idx"] = idx;

  return ret;
}
*/

vector<double> replace_data_year(const vector<double>& in_data, vector<double>& out, const vector<vector<IdxDayPair> >& yrs_index, const int rm_year, const int dup_year, const int nyr, const int day, const int win_size, const int dpy) {
  const vector<IdxDupflagPair> dfp = get_index_tuples(yrs_index, rm_year, dup_year, nyr, day, win_size, dpy);
  int numdup = 0, numnotdup = 0;
  for(vector<IdxDupflagPair>::const_iterator i = dfp.begin(); i != dfp.end(); ++i) { numdup += (int)((*i).dup); numnotdup += (int)(!(*i).dup); }
  out.resize(in_data.size() + numdup - numnotdup);

  int last_in_idx = -1, next_out_idx = 0;
  for(vector<IdxDupflagPair>::const_iterator i = dfp.begin(); i != dfp.end(); ++i) {
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

  if(last_in_idx + 1 < in_data.size()) {
    copy(&in_data[last_in_idx + 1], &in_data[in_data.size()], &out[next_out_idx]);
  }

  return out;
}

/*
RcppExport SEXP replace_data_year_R(SEXP in_data_, SEXP yrs_index_, SEXP rm_year_, SEXP dup_year_, SEXP win_size_) {
  const int rm_year = as<int>(rm_year_);
  const int dup_year = as<int>(dup_year_);
  const int win_size = as<int>(win_size_);
  const vector<double> in_data = as<vector<double> >(in_data_);

  List yrs_index_list(yrs_index_);
  vector<vector<int> > yrs_index(yrs_index_list.size());

  int yrid = 0;
  for(List::iterator i = yrs_index_list.begin(); i != yrs_index_list.end(); ++i, ++yrid) {
    const IntegerVector& d = (*i);
    yrs_index[yrid].resize(d.size());
    std::copy(d.begin(), d.end(), yrs_index[yrid].begin());
  }
  vector<double> out(yrs_index.size() * win_size);
  
  replace_data_year(in_data, out, yrs_index, rm_year, dup_year);

  return wrap(out);
}
*/

vector<double> get_data_only(const vector<DatYrTuple>& in_dat) {
  vector<double> out_dat(in_dat.size());
  const int dat_size = in_dat.size();
  for(int i = 0; i < in_dat.size(); ++i)
    out_dat[i] = in_dat[i].dat;
  return out_dat;
}

RcppExport SEXP c_quantile2(SEXP data_, SEXP quantile_) {
  const NumericVector q(quantile_);
  const NumericVector data(data_);
  const int n = data.size();
  const int nq = q.size();
  NumericVector res(nq);
  
  for(int i = 0; i < nq; ++i)
    res[i] = c_quantile(data.begin(), n, q[i]);
  
  return res;
}

RcppExport SEXP running_quantile_windowed_bootstrap(SEXP data_, SEXP n_, SEXP q_, SEXP dpy_) {
  const int win_size = as<int>(n_);
  const int days_per_year = as<int>(dpy_);
  const NumericVector q(q_);
  const NumericVector data(data_);
  const int nq = q.size();
  const int data_length = data.size();
  
  // Implicit pseudo-floor (nearest number to zero) by integer division
  const int win_border = win_size / 2;
  const int true_data_length = data_length - 2 * win_border;
  const int num_years = (int)ceil((double)true_data_length / (double)days_per_year);
  
  const int day_mul = 1;
  const int rmyr_mul = days_per_year;
  const int dupyr_mul = days_per_year * num_years;
  const int q_mul = days_per_year * num_years * (num_years - 1);

  NumericVector quantiles(nq * days_per_year * num_years * (num_years - 1));
  LogicalVector notna_map = !is_na(data);
  
  // Comment for preservation of sanity...
  // The input data does not stay on day 'day'; it starts on day 'day' - win_border. This gets around the problem of an off-by-two error.
  vector<double> rep_dat(win_size * num_years);
  for(int day = 0; day < days_per_year; ++day) {
    vector<DatYrTuple> win_tuples(extract_window_with_year(data.begin(), notna_map.begin(), day, win_size, num_years, days_per_year));
    std::sort(win_tuples.begin(), win_tuples.end());
    const vector<vector<IdxDayPair> >& yrs_index = create_yrs_index(win_tuples, num_years);
    const vector<double>& win_dat = get_data_only(win_tuples);

    const int off_day = day * day_mul;
    for(int rm_year = 0; rm_year < num_years; ++rm_year) {
      const int off_day_rmyr = off_day + rm_year * rmyr_mul;
      int dup_idx = 0;
      for(int dup_year = 0; dup_year < num_years; ++dup_year) {
	if(dup_year != rm_year) {
	  replace_data_year(win_dat, rep_dat, yrs_index, rm_year, dup_year, num_years, day, win_size, days_per_year);
	  const int off_day_rmyr_dupyr = off_day_rmyr + dup_idx * dupyr_mul;
	  const int nq_day = nq * day;
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
RcppExport SEXP running_quantile_windowed(SEXP data_, SEXP n_, SEXP q_, SEXP dpy_) {
  const int win_size = as<int>(n_);
  const int days_per_year = as<int>(dpy_);
  const NumericVector q(q_);
  const NumericVector data(data_);
  const int nq = q.size();
  const int data_length = data.size();
  
  // Implicit pseudo-floor (nearest number to zero) by integer division
  const int win_border = win_size / 2;
  const int true_data_length = data_length - 2 * win_border;
  const int num_years = (int)ceil((double)true_data_length / (double)days_per_year);
  
  NumericVector quantiles(nq * days_per_year);
  
  // This data will be 2 dimensional, row major, with the major dimension being the day of the year
  double* buf = new double[num_years * win_size];
  
  LogicalVector notna_map = !is_na(data);
  
  // Comment for preservation of sanity...
  // The input data does not stay on day 'day'; it starts on day 'day' - win_border. This gets around the problem of an off-by-two error.
  for(int day = 0; day < days_per_year; ++day) {
    int count = 0;
    // Fill buffer with data with NA removed
    for(int yr = 0; yr < num_years; ++yr) {
      const int ydayoff = yr * days_per_year + day;
      const int* ydayoff_notna = &notna_map.begin()[ydayoff];
      const double* ydayoff_data = &data.begin()[ydayoff];
      int winday = win_size;
      while(winday--)
	if(ydayoff_notna[winday])
	  buf[count++] = ydayoff_data[winday];
    }
    // Quantiles on said buffer
    const int nq_day = nq * day;
    for(int q_idx = 0; q_idx < nq; ++q_idx)
      quantiles[nq_day + q_idx] = c_quantile(buf, count, q[q_idx]);
  }
  delete[] buf;
  
  return quantiles;
}
