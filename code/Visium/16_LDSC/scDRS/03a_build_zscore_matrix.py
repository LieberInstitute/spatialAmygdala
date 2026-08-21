#!/usr/bin/env python
"""
03a_build_zscore_matrix.py -- collect the per-trait MAGMA .genes.out files into
the single gene x trait Z-score matrix that `scdrs munge-gs` consumes.

    python 03a_build_zscore_matrix.py <magma_out_dir> <traits.tsv> <out.tsv>

INPUT   $MAGMA_OUT/<trait>.genes.out, whitespace-delimited, one row per gene:
            GENE  CHR  START  STOP  NSNPS  NPARAM  N  ZSTAT  P
        GENE is the gene SYMBOL because 01_prep_magma_ref.sh built the
        gene-loc file with the symbol in column 1.

OUTPUT  a tab-delimited matrix
            GENE   ADHD   Alcohol   ...   T2D
        with the MAGMA ZSTAT in each cell. This is exactly the
        `--zscore-file` format for scdrs munge-gs, which then takes the top
        --n-max genes per trait and writes a .gs file.

WHY Z AND NOT P
    munge-gs can take either. Z keeps the sign and the ranking resolution at
    the top of the distribution, and it is what the scDRS authors used to
    build their published gene sets, so results stay comparable to the paper.

Traits with no .genes.out yet are reported and skipped rather than failing,
so this can be run against a partially finished MAGMA array.
"""
import os
import sys
import pandas as pd


def main(magma_dir, traits_tsv, out_tsv):
    traits = pd.read_csv(traits_tsv, sep="\t")
    series, missing = {}, []

    for trait in traits["TRAIT"]:
        path = os.path.join(magma_dir, f"{trait}.genes.out")
        if not os.path.exists(path):
            missing.append(trait)
            continue
        g = pd.read_csv(path, sep=r"\s+", engine="c")
        if "ZSTAT" not in g.columns or "GENE" not in g.columns:
            raise SystemExit(f"{path}: expected GENE and ZSTAT, got {list(g.columns)}")
        # a symbol should appear once; if MAGMA emitted duplicates keep the strongest
        g = g.sort_values("ZSTAT", ascending=False).drop_duplicates("GENE", keep="first")
        series[trait] = g.set_index("GENE")["ZSTAT"]

    if not series:
        raise SystemExit("no .genes.out files found -- run 02_magma_array.sh first")

    mat = pd.DataFrame(series)
    mat.index.name = "GENE"
    mat.to_csv(out_tsv, sep="\t", float_format="%.6g")

    print(f"traits written : {mat.shape[1]}")
    print(f"genes (union)  : {mat.shape[0]}")
    print(f"genes complete : {int(mat.notna().all(axis=1).sum())}")
    if missing:
        print(f"MISSING ({len(missing)}): {', '.join(missing)}")
    print("\nper-trait non-missing genes and max Z:")
    for t in mat.columns:
        col = mat[t].dropna()
        print(f"  {t:24s} n={len(col):6d}  maxZ={col.max():7.3f}")
    print(f"\n-> {out_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
