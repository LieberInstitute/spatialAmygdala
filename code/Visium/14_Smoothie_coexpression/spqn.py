"""
spqn.py

Pure Python/NumPy implementation of Spatial Quantile Normalization (SpQN)
for correcting mean-correlation bias in co-expression matrices.

Based on: Wang, Hicks & Hansen (2022) PLoS Comp Biol.
Original R implementation: https://github.com/hansenlab/spqn

Usage:
    from spqn import normalize_correlation
    cor_adj = normalize_correlation(cor_mat, ave_expr)
"""

import numpy as np


def _get_grps(ngene, ngrp, size_grp):
    """Build overlapping outer bins. Returns list of index arrays."""
    if size_grp >= ngene:
        return [np.arange(ngene)]

    d = (ngene - size_grp) / (ngrp - 1)
    grps = []
    for k in range(ngrp):
        start = int(round(k * d))
        end = min(start + size_grp, ngene)
        grps.append(np.arange(start, end))
    return grps


def _get_grps_inner(grps, ngene):
    """
    Build disjoint inner bins from overlapping outer bins.
    Follows the R implementation logic exactly.
    Returns list of index arrays.
    """
    ngrp = len(grps)
    if ngrp == 1:
        return [np.arange(ngene)]

    size_bin = len(grps[0])
    d = grps[1][0] - grps[0][0]  # step between outer bin starts

    inner = []
    # First inner bin
    length_inner_1 = int(round(size_bin / 2 + d / 2))
    inner.append(np.arange(0, length_inner_1))

    # Middle inner bins
    for k in range(1, ngrp - 1):
        prev_tail = inner[-1][-1]
        d_k = grps[k + 1][0] - grps[k][0]
        width = int(round(d_k))
        inner.append(np.arange(prev_tail + 1, prev_tail + 1 + width))

    # Last inner bin
    prev_tail = inner[-1][-1]
    inner.append(np.arange(prev_tail + 1, ngene))

    return inner


def normalize_correlation(cor_mat, ave_expr, ngrp=20, size_grp=400, ref_grp=None):
    """
    Spatial quantile normalization of a correlation matrix.

    Parameters
    ----------
    cor_mat : np.ndarray, shape (G, G)
        Symmetric correlation matrix.
    ave_expr : np.ndarray, shape (G,)
        Mean expression per gene.
    ngrp : int
        Number of overlapping bins (default 20).
    size_grp : int
        Size of each overlapping bin in genes (default 400).
    ref_grp : int or None
        Which bin to use as reference (0-indexed). Default ngrp-1.

    Returns
    -------
    cor_adj : np.ndarray, shape (G, G)
        SpQN-corrected correlation matrix.
    """
    if ref_grp is None:
        ref_grp = ngrp - 1

    ngene = cor_mat.shape[0]
    assert cor_mat.shape == (ngene, ngene)
    assert len(ave_expr) == ngene

    if size_grp >= ngene:
        size_grp = ngene
        ngrp = 1

    # Sort genes by expression
    idx_sort = np.argsort(ave_expr)
    cor_sorted = cor_mat[np.ix_(idx_sort, idx_sort)].copy()

    # Build bins
    outer = _get_grps(ngene, ngrp, size_grp)
    inner = _get_grps_inner(outer, ngene)

    # Reference distribution: upper triangle of reference block
    ref_block = cor_sorted[np.ix_(outer[ref_grp], outer[ref_grp])]
    ref_upper = ref_block[np.triu_indices_from(ref_block, k=1)]
    ref_sorted_vals = np.sort(ref_upper)
    n_ref = len(ref_sorted_vals)

    # Output: start with copy, overwrite inner blocks
    cor_adj = np.full_like(cor_sorted, np.nan)

    for i in range(ngrp):
        for j in range(i, ngrp):
            # --- Rank within the OUTER block ---
            outer_block = cor_sorted[np.ix_(outer[i], outer[j])]
            n_block = outer_block.size

            flat = outer_block.ravel()
            rank_flat = np.empty(len(flat), dtype=np.float64)
            rank_flat[np.argsort(flat)] = np.arange(len(flat), dtype=np.float64)
            rank_2d = rank_flat.reshape(outer_block.shape)

            # Number of diagonal elements (where outer[i] and outer[j] overlap)
            n_diag = len(np.intersect1d(outer[i], outer[j]))
            n_scale = max(n_block - n_diag - 1, 1)

            # --- Extract inner block ranks ---
            # Find where inner[i] genes sit within outer[i], and inner[j] within outer[j]
            # outer[i] is a sorted array of gene indices; inner[i] is a subset
            outer_i_set = set(outer[i].tolist())
            outer_j_set = set(outer[j].tolist())

            # Which inner genes are in the outer bin
            inner_i_in = np.array([g for g in inner[i] if g in outer_i_set])
            inner_j_in = np.array([g for g in inner[j] if g in outer_j_set])

            if len(inner_i_in) == 0 or len(inner_j_in) == 0:
                continue

            # Positions of those genes within the outer bin arrays
            outer_i_list = outer[i].tolist()
            outer_j_list = outer[j].tolist()
            pos_i = np.array([outer_i_list.index(g) for g in inner_i_in])
            pos_j = np.array([outer_j_list.index(g) for g in inner_j_in])

            inner_ranks = rank_2d[np.ix_(pos_i, pos_j)]

            # --- Scale to [0, 1] then to reference ---
            # R code: (rank - 1) / n_scale, then 1 + scaled * (n_ref - 1)
            scaled = (inner_ranks - 1) / n_scale
            scaled = np.clip(scaled, 0.0, 1.0)
            ref_idx = 1.0 + scaled * (n_ref - 1)
            ref_idx = np.clip(ref_idx, 1.0, n_ref)

            # --- Interpolate from reference distribution ---
            # Convert to 0-indexed for array access
            idx0 = ref_idx - 1.0
            idx_low = np.floor(idx0).astype(int)
            idx_low = np.clip(idx_low, 0, n_ref - 1)
            idx_high = np.minimum(idx_low + 1, n_ref - 1)
            frac = idx0 - idx_low
            adjusted = (1 - frac) * ref_sorted_vals[idx_low] + frac * ref_sorted_vals[idx_high]

            # --- Write into output at the inner gene positions ---
            cor_adj[np.ix_(inner_i_in, inner_j_in)] = adjusted
            if i != j:
                cor_adj[np.ix_(inner_j_in, inner_i_in)] = adjusted.T

    # Diagonal = 1
    np.fill_diagonal(cor_adj, 1.0)

    # Unsort back to original gene order
    idx_unsort = np.argsort(idx_sort)
    cor_adj = cor_adj[np.ix_(idx_unsort, idx_unsort)]

    return cor_adj
