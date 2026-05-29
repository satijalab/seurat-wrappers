// Zero-copy CSC sparse × dense operations for BANKSY lazy operator.
// All functions read dgCMatrix slots (@i, @p, @x) directly.
// OpenMP parallelizes within each matmul across columns.

#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstring>
#include <vector>

using namespace Rcpp;

// ── Helper: CSC column-subset matmul ────────────────────────────────────────
// Compute result[nrow x k] += mat[, col_idx] %*% rhs[n_cols x k]
// Accumulates into pre-zeroed result buffer. No allocation.
static void csc_matmul_accum(const int* mi, const int* mp, const double* mx,
                             int nrow, const int* cidx, int n_cols,
                             const double* rhs, int k,
                             double* result) {
    #pragma omp parallel
    {
        std::vector<double> local(nrow * k, 0.0);
        #pragma omp for schedule(dynamic, 64)
        for (int c = 0; c < n_cols; c++) {
            int j = cidx[c] - 1;
            int p_start = mp[j];
            int p_end = mp[j + 1];
            for (int p = p_start; p < p_end; p++) {
                int row = mi[p];
                double val = mx[p];
                for (int kk = 0; kk < k; kk++) {
                    local[row + kk * nrow] += val * rhs[c + kk * n_cols];
                }
            }
        }
        #pragma omp critical
        for (int i = 0; i < nrow * k; i++) result[i] += local[i];
    }
}

// ── Helper: CSC column-subset crossprod ─────────────────────────────────────
// Compute result[n_cols x k] = t(mat[, col_idx]) %*% lhs[nrow x k]
// Each output row is independent — no race condition.
static void csc_crossprod_accum(const int* mi, const int* mp, const double* mx,
                                int nrow, const int* cidx, int n_cols,
                                const double* lhs, int k,
                                double* result) {
    #pragma omp parallel for schedule(dynamic, 64)
    for (int c = 0; c < n_cols; c++) {
        int j = cidx[c] - 1;
        int p_start = mp[j];
        int p_end = mp[j + 1];
        for (int p = p_start; p < p_end; p++) {
            int row = mi[p];
            double val = mx[p];
            for (int kk = 0; kk < k; kk++) {
                result[c + kk * n_cols] += val * lhs[row + kk * nrow];
            }
        }
    }
}

// ── Helper: full-matrix CSC matmul (no column subset) ───────────────────────
// result[nrow x k] += mat %*% rhs[ncol x k]
static void csc_full_matmul_accum(const int* mi, const int* mp, const double* mx,
                                  int nrow, int ncol,
                                  const double* rhs, int k,
                                  double* result) {
    #pragma omp parallel
    {
        std::vector<double> local(nrow * k, 0.0);
        #pragma omp for schedule(dynamic, 64)
        for (int c = 0; c < ncol; c++) {
            int p_start = mp[c];
            int p_end = mp[c + 1];
            for (int p = p_start; p < p_end; p++) {
                int row = mi[p];
                double val = mx[p];
                for (int kk = 0; kk < k; kk++) {
                    local[row + kk * nrow] += val * rhs[c + kk * ncol];
                }
            }
        }
        #pragma omp critical
        for (int i = 0; i < nrow * k; i++) result[i] += local[i];
    }
}

// ── Forward: all groups in one C++ call ─────────────────────────────────────
// Computes the full forward pass for irlba: result = rbind(own, h0)
//
// Takes per-group gcm CSC slots. For each group gr:
//   x_g = x[cid, ]           (cid maps group-local → global cell index)
//   Wx  = W_gr %*% x_g
//   own += gcm_gr %*% x_g    (direct per-group gcm access)
//   h0  += gcm_gr %*% Wx
//
// If split_scale: center/scale own and h0 per group before accumulating.
// If !split_scale: center/scale globally after the loop.
// Excess corrections are applied in R (sparse, negligible cost).
//
// [[Rcpp::export]]
NumericMatrix banksy_forward_cpp(
        List gcm_i_list, List gcm_p_list, List gcm_x_list,
        int gcm_nrow, int n_cells,
        List W_i_list, List W_p_list, List W_x_list, IntegerVector W_ncol_vec,
        List group_idx_list,
        NumericMatrix x,
        bool split_scale,
        NumericMatrix mu_own_mat, NumericMatrix sd_own_mat,
        NumericMatrix mu_h0_mat, NumericMatrix sd_h0_mat,
        NumericVector lam,
        Nullable<LogicalVector> valid_own_,
        Nullable<LogicalVector> valid_h0_) {

    int ng = gcm_nrow;
    int k = x.ncol();
    int n_groups = group_idx_list.size();

    std::vector<double> own(ng * k, 0.0);
    std::vector<double> h0(ng * k, 0.0);

    int max_ng = 0;
    for (int gr = 0; gr < n_groups; gr++) {
        int sz = ((IntegerVector)group_idx_list[gr]).size();
        if (sz > max_ng) max_ng = sz;
    }
    std::vector<double> x_g(max_ng * k);
    std::vector<double> Wx(max_ng * k);
    std::vector<double> own_g(ng * k);
    std::vector<double> h0_g(ng * k);

    for (int gr = 0; gr < n_groups; gr++) {
        IntegerVector cid_r = group_idx_list[gr];
        int n_g = cid_r.size();
        const int* cid = cid_r.begin();

        // Per-group gcm slots
        IntegerVector gmi_r = gcm_i_list[gr];
        IntegerVector gmp_r = gcm_p_list[gr];
        NumericVector gmx_r = gcm_x_list[gr];
        const int* gmi = gmi_r.begin();
        const int* gmp = gmp_r.begin();
        const double* gmx = gmx_r.begin();

        IntegerVector wi = W_i_list[gr];
        IntegerVector wp = W_p_list[gr];
        NumericVector wx = W_x_list[gr];
        int w_ncol = W_ncol_vec[gr];

        // x_g = x[cid, ] (cid maps group-local → global)
        for (int kk = 0; kk < k; kk++)
            for (int c = 0; c < n_g; c++)
                x_g[c + kk * n_g] = x[cid[c] - 1 + kk * n_cells];

        // Wx = W_gr %*% x_g
        std::memset(Wx.data(), 0, n_g * k * sizeof(double));
        csc_full_matmul_accum(wi.begin(), wp.begin(), wx.begin(),
                              n_g, w_ncol, x_g.data(), k, Wx.data());

        // own_g = gcm_gr %*% x_g (direct per-group access)
        std::memset(own_g.data(), 0, ng * k * sizeof(double));
        csc_full_matmul_accum(gmi, gmp, gmx, ng, n_g,
                              x_g.data(), k, own_g.data());

        // h0_g = gcm_gr %*% Wx
        std::memset(h0_g.data(), 0, ng * k * sizeof(double));
        csc_full_matmul_accum(gmi, gmp, gmx, ng, n_g,
                              Wx.data(), k, h0_g.data());

        if (split_scale) {
            std::vector<double> cs(k, 0.0);
            for (int kk = 0; kk < k; kk++)
                for (int c = 0; c < n_g; c++)
                    cs[kk] += x_g[c + kk * n_g];

            const double* mu_o = &mu_own_mat[gr * ng];
            const double* sd_o = &sd_own_mat[gr * ng];
            const double* mu_h = &mu_h0_mat[gr * ng];
            const double* sd_h = &sd_h0_mat[gr * ng];

            for (int r = 0; r < ng; r++) {
                double sdo = sd_o[r];
                double sdh = sd_h[r];
                for (int kk = 0; kk < k; kk++) {
                    own_g[r + kk * ng] =
                        (own_g[r + kk * ng] - mu_o[r] * cs[kk]) / sdo;
                    h0_g[r + kk * ng] =
                        (h0_g[r + kk * ng] - mu_h[r] * cs[kk]) / sdh;
                }
            }
        }

        for (int i = 0; i < ng * k; i++) {
            own[i] += own_g[i];
            h0[i] += h0_g[i];
        }
    }

    double lam_own = lam[0];
    double lam_h0 = lam[1];

    if (split_scale) {
        if (valid_own_.isNotNull()) {
            LogicalVector valid_own(valid_own_);
            for (int r = 0; r < ng; r++)
                if (!valid_own[r])
                    for (int kk = 0; kk < k; kk++)
                        own[r + kk * ng] = 0.0;
        }
        if (valid_h0_.isNotNull()) {
            LogicalVector valid_h0(valid_h0_);
            for (int r = 0; r < ng; r++)
                if (!valid_h0[r])
                    for (int kk = 0; kk < k; kk++)
                        h0[r + kk * ng] = 0.0;
        }
        for (int i = 0; i < ng * k; i++) {
            own[i] *= lam_own;
            h0[i] *= lam_h0;
        }
    }

    NumericMatrix result(2 * ng, k);
    for (int kk = 0; kk < k; kk++) {
        std::memcpy(&result[kk * 2 * ng], &own[kk * ng], ng * sizeof(double));
        std::memcpy(&result[kk * 2 * ng + ng], &h0[kk * ng], ng * sizeof(double));
    }
    return result;
}

// ── Adjoint: all groups in one C++ call ─────────────────────────────────────
// Computes the full adjoint pass: result[n_cells x k]
// Takes per-group gcm CSC slots. cid used for output indexing only.
//
// [[Rcpp::export]]
NumericMatrix banksy_adjoint_cpp(
        List gcm_i_list, List gcm_p_list, List gcm_x_list,
        int gcm_nrow, int n_cells,
        List W_i_list, List W_p_list, List W_x_list, IntegerVector W_ncol_vec,
        List group_idx_list,
        NumericMatrix xo_s, NumericMatrix xh_s,
        NumericVector adj_o, NumericVector adj_h,
        bool split_scale,
        NumericMatrix mu_own_mat, NumericMatrix sd_own_mat,
        NumericMatrix mu_h0_mat, NumericMatrix sd_h0_mat,
        NumericVector lam) {

    int ng = gcm_nrow;
    int nc = n_cells;
    int k = xo_s.ncol();
    int n_groups = group_idx_list.size();
    double lam_own = lam[0];
    double lam_h0 = lam[1];

    NumericMatrix result(nc, k);

    int max_ng = 0;
    for (int gr = 0; gr < n_groups; gr++) {
        int sz = ((IntegerVector)group_idx_list[gr]).size();
        if (sz > max_ng) max_ng = sz;
    }
    std::vector<double> xo_local(ng * k);
    std::vector<double> xh_local(ng * k);
    std::vector<double> adj_o_local(k);
    std::vector<double> adj_h_local(k);
    std::vector<double> own_raw(max_ng * k);
    std::vector<double> ht_raw(max_ng * k);
    std::vector<double> ht(max_ng * k);

    for (int gr = 0; gr < n_groups; gr++) {
        IntegerVector cid_r = group_idx_list[gr];
        int n_g = cid_r.size();
        const int* cid = cid_r.begin();

        // Per-group gcm slots
        IntegerVector gmi_r = gcm_i_list[gr];
        IntegerVector gmp_r = gcm_p_list[gr];
        NumericVector gmx_r = gcm_x_list[gr];
        const int* gmi = gmi_r.begin();
        const int* gmp = gmp_r.begin();
        const double* gmx = gmx_r.begin();

        IntegerVector wi = W_i_list[gr];
        IntegerVector wp = W_p_list[gr];
        NumericVector wx = W_x_list[gr];
        int w_ncol = W_ncol_vec[gr];

        const double* xo_ptr;
        const double* xh_ptr;
        const double* adj_o_ptr;
        const double* adj_h_ptr;

        if (split_scale) {
            const double* mu_o = &mu_own_mat[gr * ng];
            const double* sd_o = &sd_own_mat[gr * ng];
            const double* mu_h = &mu_h0_mat[gr * ng];
            const double* sd_h = &sd_h0_mat[gr * ng];

            std::fill(adj_o_local.begin(), adj_o_local.end(), 0.0);
            std::fill(adj_h_local.begin(), adj_h_local.end(), 0.0);
            for (int kk = 0; kk < k; kk++) {
                for (int r = 0; r < ng; r++) {
                    double xo_val = xo_s[r + kk * ng] / sd_o[r];
                    double xh_val = xh_s[r + kk * ng] / sd_h[r];
                    xo_local[r + kk * ng] = xo_val;
                    xh_local[r + kk * ng] = xh_val;
                    adj_o_local[kk] += mu_o[r] * xo_val;
                    adj_h_local[kk] += mu_h[r] * xh_val;
                }
            }
            xo_ptr = xo_local.data();
            xh_ptr = xh_local.data();
            adj_o_ptr = adj_o_local.data();
            adj_h_ptr = adj_h_local.data();
        } else {
            xo_ptr = &xo_s[0];
            xh_ptr = &xh_s[0];
            adj_o_ptr = adj_o.begin();
            adj_h_ptr = adj_h.begin();
        }

        std::memset(own_raw.data(), 0, n_g * k * sizeof(double));
        std::memset(ht_raw.data(), 0, n_g * k * sizeof(double));

        // Fused crossprod: direct per-group column access
        #pragma omp parallel for schedule(dynamic, 64)
        for (int c = 0; c < n_g; c++) {
            int p_start = gmp[c];
            int p_end = gmp[c + 1];
            for (int p = p_start; p < p_end; p++) {
                int row = gmi[p];
                double val = gmx[p];
                for (int kk = 0; kk < k; kk++) {
                    own_raw[c + kk * n_g] += val * xo_ptr[row + kk * ng];
                    ht_raw[c + kk * n_g] += val * xh_ptr[row + kk * ng];
                }
            }
            for (int kk = 0; kk < k; kk++)
                own_raw[c + kk * n_g] -= adj_o_ptr[kk];
        }

        // ht = t(W_gr) %*% ht_raw - adj_h
        std::memset(ht.data(), 0, n_g * k * sizeof(double));
        #pragma omp parallel for schedule(dynamic, 64)
        for (int c = 0; c < w_ncol; c++) {
            int p_start = wp[c];
            int p_end = wp[c + 1];
            for (int p = p_start; p < p_end; p++) {
                int row = wi[p];
                double val = wx[p];
                for (int kk = 0; kk < k; kk++)
                    ht[c + kk * n_g] += val * ht_raw[row + kk * n_g];
            }
        }
        for (int c = 0; c < n_g; c++)
            for (int kk = 0; kk < k; kk++)
                ht[c + kk * n_g] -= adj_h_ptr[kk];

        // r[cid, ] = lam_own * own + lam_h0 * ht (cid for output)
        for (int c = 0; c < n_g; c++) {
            int ci = cid[c] - 1;
            for (int kk = 0; kk < k; kk++)
                result[ci + kk * nc] =
                    lam_own * own_raw[c + kk * n_g] +
                    lam_h0 * ht[c + kk * n_g];
        }
    }

    return result;
}

// ── H0 scaling: per-group stats without materializing gcm_gr %*% W ─────────
// Takes per-group gcm CSC slots. Column-wise sweep: for each output column j,
//   product_col = sum over W_gr's nnz in col j: gcm_gr[, w_row] * w_val
// then accumulate row-wise mu (sum), ss (sum of squares), max.
// Memory: O(n_genes) per thread. No intermediate matrix allocation.

// [[Rcpp::export]]
List h0_group_stats_cpp(
        List gcm_i_list, List gcm_p_list, List gcm_x_list,
        int gcm_nrow,
        List W_i_list, List W_p_list, List W_x_list,
        List group_idx_list) {

    int n_genes = gcm_nrow;
    int n_groups = group_idx_list.size();

    NumericMatrix mu_out(n_genes, n_groups);
    NumericMatrix ss_out(n_genes, n_groups);
    NumericMatrix max_out(n_genes, n_groups);
    std::fill(max_out.begin(), max_out.end(), -INFINITY);

    for (int gr = 0; gr < n_groups; gr++) {
        IntegerVector gmi_r = gcm_i_list[gr];
        IntegerVector gmp_r = gcm_p_list[gr];
        NumericVector gmx_r = gcm_x_list[gr];
        const int* gmi = gmi_r.begin();
        const int* gmp = gmp_r.begin();
        const double* gmx = gmx_r.begin();

        int n_g = ((IntegerVector)group_idx_list[gr]).size();

        IntegerVector wi = W_i_list[gr];
        IntegerVector wp = W_p_list[gr];
        NumericVector wx = W_x_list[gr];

        double* mu_gr = &mu_out[gr * n_genes];
        double* ss_gr = &ss_out[gr * n_genes];
        double* mx_gr = &max_out[gr * n_genes];

        #pragma omp parallel
        {
            std::vector<double> col_buf(n_genes, 0.0);
            std::vector<double> lmu(n_genes, 0.0);
            std::vector<double> lss(n_genes, 0.0);
            std::vector<double> lmx(n_genes, -INFINITY);

            #pragma omp for schedule(dynamic, 256)
            for (int j = 0; j < n_g; j++) {
                int wp_s = wp[j], wp_e = wp[j + 1];
                for (int p = wp_s; p < wp_e; p++) {
                    int w_row = wi[p];
                    // w_row is a local column index in per-group gcm (0-based)
                    int gp_s = gmp[w_row], gp_e = gmp[w_row + 1];
                    for (int gp = gp_s; gp < gp_e; gp++)
                        col_buf[gmi[gp]] += gmx[gp] * wx[p];
                }
                for (int r = 0; r < n_genes; r++) {
                    double v = col_buf[r];
                    lmu[r] += v;
                    lss[r] += v * v;
                    if (v > lmx[r]) lmx[r] = v;
                    col_buf[r] = 0.0;
                }
            }

            #pragma omp critical
            for (int r = 0; r < n_genes; r++) {
                mu_gr[r] += lmu[r];
                ss_gr[r] += lss[r];
                if (lmx[r] > mx_gr[r]) mx_gr[r] = lmx[r];
            }
        }

        double inv_ng = 1.0 / n_g;
        for (int r = 0; r < n_genes; r++)
            mu_gr[r] *= inv_ng;
    }

    return List::create(
        Named("mu") = mu_out,
        Named("ss") = ss_out,
        Named("max_val") = max_out);
}

// ── H0 scaling: per-group excess (clipping) without materializing product ──
// Same column-wise sweep as stats. Takes per-group gcm CSC slots.
// Records entries where z > scale_max for clip genes.

// [[Rcpp::export]]
List h0_group_excess_cpp(
        List gcm_i_list, List gcm_p_list, List gcm_x_list,
        int gcm_nrow,
        List W_i_list, List W_p_list, List W_x_list,
        List group_idx_list,
        IntegerVector clip_genes_r,
        bool split_scale,
        NumericMatrix mu_mat,
        NumericMatrix sd_mat,
        double scale_max) {

    int n_genes = gcm_nrow;
    int n_groups = group_idx_list.size();
    int n_clip = clip_genes_r.size();

    std::vector<int> clip_0(n_clip);
    for (int i = 0; i < n_clip; i++)
        clip_0[i] = clip_genes_r[i] - 1;

    std::vector<int> all_ei, all_ej;
    std::vector<double> all_ex;

    for (int gr = 0; gr < n_groups; gr++) {
        IntegerVector cid_r = group_idx_list[gr];
        int n_g = cid_r.size();
        const int* cid = cid_r.begin();

        IntegerVector gmi_r = gcm_i_list[gr];
        IntegerVector gmp_r = gcm_p_list[gr];
        NumericVector gmx_r = gcm_x_list[gr];
        const int* gmi = gmi_r.begin();
        const int* gmp = gmp_r.begin();
        const double* gmx = gmx_r.begin();

        IntegerVector wi = W_i_list[gr];
        IntegerVector wp = W_p_list[gr];
        NumericVector wx = W_x_list[gr];

        const double* mu_ptr = split_scale ? &mu_mat[gr * n_genes] : &mu_mat[0];
        const double* sd_ptr = split_scale ? &sd_mat[gr * n_genes] : &sd_mat[0];

        #pragma omp parallel
        {
            std::vector<double> col_buf(n_genes, 0.0);
            std::vector<int> local_i, local_j;
            std::vector<double> local_x;

            #pragma omp for schedule(dynamic, 256)
            for (int j = 0; j < n_g; j++) {
                int wp_s = wp[j], wp_e = wp[j + 1];
                for (int p = wp_s; p < wp_e; p++) {
                    int w_row = wi[p];
                    int gp_s = gmp[w_row], gp_e = gmp[w_row + 1];
                    for (int gp = gp_s; gp < gp_e; gp++)
                        col_buf[gmi[gp]] += gmx[gp] * wx[p];
                }

                int cell_1b = cid[j];  // cid still needed for output
                for (int ci = 0; ci < n_clip; ci++) {
                    int r = clip_0[ci];
                    double z = (col_buf[r] - mu_ptr[r]) / sd_ptr[r];
                    if (z > scale_max) {
                        local_i.push_back(r + 1);
                        local_j.push_back(cell_1b);
                        local_x.push_back(z - scale_max);
                    }
                }

                for (int r = 0; r < n_genes; r++)
                    col_buf[r] = 0.0;
            }

            #pragma omp critical
            {
                all_ei.insert(all_ei.end(), local_i.begin(), local_i.end());
                all_ej.insert(all_ej.end(), local_j.begin(), local_j.end());
                all_ex.insert(all_ex.end(), local_x.begin(), local_x.end());
            }
        }
    }

    return List::create(
        Named("i") = wrap(all_ei),
        Named("j") = wrap(all_ej),
        Named("x") = wrap(all_ex));
}

// ── Own expression excess (clipping) via C++ ───────────────────────────────
// For each nonzero in per-group gcm, compute z = (x - mu) / sd.
// If z > scale_max, record (row+1, global_cell, z - scale_max).
// Eliminates ~70GB of R temporary vectors at 10M cells.

// [[Rcpp::export]]
List own_excess_cpp(
        List gcm_i_list, List gcm_p_list, List gcm_x_list,
        List group_idx_list,
        bool split_scale,
        NumericMatrix mu_mat,
        NumericMatrix sd_mat,
        NumericVector mu_vec,
        NumericVector sd_vec,
        Nullable<LogicalVector> valid_,
        double scale_max,
        int n_genes, int n_cells) {

    int n_groups = group_idx_list.size();

    std::vector<int> all_ei, all_ej;
    std::vector<double> all_ex;

    bool has_valid = valid_.isNotNull();
    const int* valid_ptr = NULL;
    LogicalVector valid_lv;
    if (has_valid) {
        valid_lv = LogicalVector(valid_);
        valid_ptr = valid_lv.begin();
    }

    for (int gr = 0; gr < n_groups; gr++) {
        IntegerVector cid_r = group_idx_list[gr];
        int n_g = cid_r.size();
        const int* cid = cid_r.begin();

        IntegerVector gi = gcm_i_list[gr];
        IntegerVector gp = gcm_p_list[gr];
        NumericVector gx = gcm_x_list[gr];

        const double* mu_ptr;
        const double* sd_ptr;
        if (split_scale) {
            mu_ptr = &mu_mat[gr * n_genes];
            sd_ptr = &sd_mat[gr * n_genes];
        } else {
            mu_ptr = mu_vec.begin();
            sd_ptr = sd_vec.begin();
        }

        // Iterate over columns of per-group gcm
        for (int c = 0; c < n_g; c++) {
            int cell_1b = cid[c];  // global cell index (1-based)
            int p_start = gp[c];
            int p_end = gp[c + 1];
            for (int p = p_start; p < p_end; p++) {
                int row = gi[p];  // 0-based gene index
                if (has_valid && !valid_ptr[row]) continue;
                double val = gx[p];
                double m = mu_ptr[row];
                double s = sd_ptr[row];
                if (val > m + scale_max * s) {
                    double z = (val - m) / s;
                    all_ei.push_back(row + 1);
                    all_ej.push_back(cell_1b);
                    all_ex.push_back(z - scale_max);
                }
            }
        }
    }

    return List::create(
        Named("i") = wrap(all_ei),
        Named("j") = wrap(all_ej),
        Named("x") = wrap(all_ex));
}

// ── Lanczos Q buffer management ─────────────────────────────────────────────
// Q is a flat NumericVector (n * work), column-major. Managed in C++ to avoid
// R's copy-on-modify. Peak memory = Q + one column. No spike at extraction.

// [[Rcpp::export]]
void lanczos_set_col(NumericVector Q, int j, NumericVector v, int n) {
    std::memcpy(&Q[j * (long long)n], v.begin(), n * sizeof(double));
}

// [[Rcpp::export]]
NumericVector lanczos_get_col(NumericVector Q, int j, int n) {
    NumericVector v(n);
    std::memcpy(v.begin(), &Q[j * (long long)n], n * sizeof(double));
    return v;
}

// [[Rcpp::export]]
NumericVector lanczos_reorth(NumericVector Q, int ncols, NumericVector v, int n) {
    // v = v - Q[,0:ncols-1] * (Q[,0:ncols-1]^T * v)
    // Two-pass CGS for numerical stability
    for (int pass = 0; pass < 2; pass++) {
        for (int c = 0; c < ncols; c++) {
            const double* qc = &Q[c * (long long)n];
            double dot = 0.0;
            #pragma omp parallel for reduction(+:dot) schedule(static)
            for (int i = 0; i < n; i++)
                dot += qc[i] * v[i];
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; i++)
                v[i] -= dot * qc[i];
        }
    }
    return v;
}

// [[Rcpp::export]]
double lanczos_norm(NumericVector v) {
    double s = 0.0;
    int n = v.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (int i = 0; i < n; i++)
        s += v[i] * v[i];
    return std::sqrt(s);
}

// [[Rcpp::export]]
NumericMatrix lanczos_extract(NumericVector Q, NumericMatrix vectors,
                              int n, int work, int npcs) {
    // V[n x npcs] = Q[n x work] %*% vectors[work x npcs]
    NumericMatrix V(n, npcs);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int col = 0; col < npcs; col++) {
            double sum = 0.0;
            for (int k = 0; k < work; k++)
                sum += Q[k * (long long)n + i] * vectors(k, col);
            V(i, col) = sum;
        }
    }
    return V;
}

