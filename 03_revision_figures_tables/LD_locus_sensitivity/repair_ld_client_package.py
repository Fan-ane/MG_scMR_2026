#!/usr/bin/env python3
"""Repair and validate the interim LD client-delivery package.

The source analysis is an Ensembl/1000G Phase 3 EUR LD sensitivity analysis,
not a PLINK clumping run.  This script reconstructs the two membership tables,
adds client-facing documentation, verifies every expected row count and
rebuilds the ZIP only after all files have been closed and validated.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile


EXPECTED_ROWS = {
    "01_tables/01_locked_39_gene_cell_pairs.csv": 39,
    "01_tables/02_locked_33_lead_variants.csv": 33,
    "01_tables/03_pairwise_LD_55_with_sensitivity_annotation.csv": 55,
    "01_tables/04_locus_membership_observed_support.csv": 33,
    "01_tables/05_locus_summary_observed_support.csv": 13,
    "01_tables/06_gene_cell_mapping_observed_support.csv": 39,
    "01_tables/07_locus_membership_forced_merge.csv": 33,
    "01_tables/08_locus_summary_forced_merge.csv": 13,
    "01_tables/09_gene_cell_mapping_forced_merge.csv": 39,
    "01_tables/10_MHC_boundary_LD_summary.csv": 3,
    "01_tables/11_unresolved_pair_EUR_subpopulation_audit.csv": 6,
    "01_tables/12_LD_locus_count_sensitivity_summary.csv": 2,
}


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def md5(path: Path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pairwise_lookup(rows):
    out = {}
    for row in rows:
        try:
            out[frozenset((row["rsid_1"], row["rsid_2"]))] = float(row["r2"])
        except (KeyError, TypeError, ValueError):
            pass
    return out


def reconstruct_membership(root: Path, summary_name: str, output_name: str):
    table = root / "01_tables"
    variants = read_csv(table / "02_locked_33_lead_variants.csv")
    summaries = read_csv(table / summary_name)
    ld_lookup = pairwise_lookup(
        read_csv(table / "03_pairwise_LD_55_with_sensitivity_annotation.csv")
    )
    variant_by_rsid = {x["rsid"]: x for x in variants}
    rows = []
    for locus in summaries:
        index = locus["index_rsid"]
        members = [x for x in locus["lead_variants"].split(";") if x]
        for member in members:
            source = variant_by_rsid[member]
            r2 = 1.0 if member == index else ld_lookup.get(frozenset((index, member)))
            if r2 is None:
                raise RuntimeError(f"Missing index-member LD: {index} / {member}")
            rows.append({
                "ld_locus": locus["ld_locus"],
                "index_rsid": index,
                "member_rsid": member,
                "chr_grch37": source["chr_grch37"],
                "pos_grch37": source["pos_grch37"],
                "r2_to_index": f"{r2:.6f}" if member != index else "1",
            })
    if len(rows) != 33 or len({x["member_rsid"] for x in rows}) != 33:
        raise RuntimeError("Reconstructed membership did not contain 33 unique variants")
    write_csv(
        table / output_name,
        rows,
        ["ld_locus", "index_rsid", "member_rsid", "chr_grch37",
         "pos_grch37", "r2_to_index"],
    )


def write_texts(root: Path):
    text_dir = root / "03_text"
    code_dir = root / "04_code"
    code_dir.mkdir(exist_ok=True)

    (root / "README_先读.md").write_text(
        """# MG LD 位点收敛敏感性分析（客户核对版）

## 结论

- 39 个显著基因–细胞类型配对对应 33 个 lead variants。
- 33 个变异在观测支持情景和最保守强制合并情景下均收敛为 **13 个 LD-based loci**。
- 55 个同染色体且距离不超过 1 Mb 的变异对中，54 个获得 pooled-EUR LD；唯一未解析配对为 rs3893044–rs12453217。
- 该未解析配对的最高可用欧洲亚群估计来自 GBR（r²=0.055765），而且两个变异已通过 rs7359623 归入同一位点 LD03，所以强制合并不会改变 13 个位点的结论。
- 这是 **1000G Phase 3 EUR/Ensembl GRCh37 LD-based sensitivity analysis**，不是 PLINK clumping；稿件及回复中不得写成 PLINK 结果。

## MHC 边界敏感性

- rs13191296 与扩展 MHC 内变异存在强 LD（最大 r²=0.965696）。
- rs34991172 与扩展 MHC 内变异存在较弱但超过 0.1 的 LD（最大 r²=0.153467）。
- rs116418977 在评估窗口内没有直接 MHC 配对估计，但其与 rs34991172 的 r²=0.639146。
- 建议保留坐标定义的 post-MHC 分析为主分析，并补充同时排除上述 3 个边界邻近变异的保守敏感性计数；不要宣称 MHC 影响已被完全消除。

## 文件导航

- `01_tables/05_locus_summary_observed_support.csv`：13 个位点摘要。
- `01_tables/06_gene_cell_mapping_observed_support.csv`：39 个配对到位点的映射。
- `01_tables/04_locus_membership_observed_support.csv`：33 个 lead variants 的位点成员表。
- `01_tables/10_MHC_boundary_LD_summary.csv`：3 个边界变异的 MHC LD 摘要。
- `02_figures/Main_LD_sensitivity_summary.*`：主汇总图（PDF/PNG）。
- `03_text/R1_comment3_response_interim.txt`：审稿回复草稿。
- `PACKAGE_QA_REPORT.csv`：逐文件行数、大小与 MD5 验收结果。
""",
        encoding="utf-8",
    )

    (text_dir / "RESULTS_interim_LD_sensitivity.txt").write_text(
        """Interim LD-based locus results

The 39 globally significant post-MHC gene–cell-type pairs corresponded to 33 distinct lead variants. Of 55 same-chromosome lead-variant pairs separated by no more than 1 Mb, 54 had pairwise LD estimates in pooled 1000 Genomes Phase 3 Europeans. The remaining pair, rs3893044–rs12453217, had no pooled-EUR estimate; its maximum available European-subpopulation estimate was observed in GBR (r2=0.055765). Both variants were already assigned to LD03 through the index variant rs7359623 (r2=0.424284 and 0.122050, respectively). Consequently, both the observed-support scenario and a conservative scenario that forced the unresolved pair to be dependent yielded 13 LD-based loci.

For the three variants immediately outside the coordinate-defined extended-MHC boundary, rs13191296 showed strong LD with an MHC variant (maximum r2=0.965696), rs34991172 showed weaker LD exceeding 0.1 (maximum r2=0.153467), and no direct MHC partner was observed for rs116418977 within the evaluated 500-kb window. However, rs116418977 was correlated with rs34991172 (r2=0.639146). The coordinate-based post-MHC analysis was therefore retained as the primary analysis, accompanied by a conservative boundary sensitivity analysis excluding all three variants.
""",
        encoding="utf-8",
    )

    (text_dir / "R1_comment3_response_interim.txt").write_text(
        """Response to Reviewer 1, Comment 3

We thank the reviewer for requesting locus-level reporting. Without reselecting the prespecified cis instruments, we consolidated the 33 lead variants representing 39 globally significant post-MHC gene–cell-type pairs using 1000 Genomes Phase 3 European pairwise LD estimates obtained through the Ensembl GRCh37 REST API (r2 threshold 0.10; 1-Mb window). Fifty-four of 55 required comparisons were evaluable in pooled EUR data. The only unresolved pair, rs3893044–rs12453217, had a maximum available European-subpopulation estimate of r2=0.055765 in GBR. Because both variants were already assigned to the same locus through rs7359623, both an observed-support analysis and a conservative forced-dependency analysis yielded 13 LD-based loci. Complete variant and gene–cell-type membership is provided in the Supplementary Tables. We describe this result as an LD-based sensitivity analysis and not as PLINK clumping; confirmation using PLINK and a complete ancestry-matched GRCh37 reference should be performed if that reference becomes available.
""",
        encoding="utf-8",
    )

    (text_dir / "给客户的说明.txt").write_text(
        """给客户的核对说明

1. 本次可以报告“33 个 lead variants 在两种敏感性情景下均收敛为 13 个 LD-based loci”。不要写“13–13 个”，也不要写“PLINK clumping 得到 13 个”。
2. 唯一缺少 pooled-EUR 配对估计的是 rs3893044–rs12453217；最高可用欧洲亚群估计为 GBR r²=0.055765。但这两个变异本来已通过 rs7359623 归入 LD03，所以最保守强制合并也不改变位点数。
3. MHC 边界结果应作为透明的敏感性分析：rs13191296 强 LD，rs34991172 较弱但 r²>0.1，rs116418977 无直接估计但与 rs34991172 的 r²=0.639146。
4. 手稿主线建议仍以坐标定义的 post-MHC 分析为主；把边界三变异排除结果放补充材料并在局限性中说明，不要把文章改写成方法学文章，也不要声称 MHC 生物学不重要。
5. 若后续获得完整 1000G EUR GRCh37 基因型，可再用 PLINK 对 13 个位点数进行确认；在此之前应使用“LD-based locus consolidation sensitivity analysis”。
""",
        encoding="utf-8",
    )

    (root / "FIGURE_LEGEND_LD.md").write_text(
        """## Main LD sensitivity figure

**LD-based locus consolidation and MHC-boundary sensitivity.** (A) Counts of globally significant gene–cell-type pairs, distinct lead variants, and loci under the observed-support and conservative forced-dependency scenarios. (B) Number of gene–cell-type pairs assigned to each of the 13 LD-based loci. (C) GRCh37 positions of the three variants immediately outside the coordinate-defined extended-MHC boundary; the shaded region denotes the extended MHC. (D) Maximum pairwise r² with variants inside the extended MHC and the LD estimate between rs116418977 and rs34991172. Dashed and dotted lines denote r²=0.10 and r²=0.50. The analysis used 1000 Genomes Phase 3 EUR estimates from the Ensembl GRCh37 REST API and is not a PLINK clumping result.

## Supplementary pairwise LD figures

**Pairwise LD among lead variants located on the same chromosome and within 1 Mb.** The 55 required comparisons are split across two panels for readability. The vertical dashed line denotes r²=0.10. Orange marks the single unresolved pooled-EUR comparison represented by the maximum available European-subpopulation estimate (GBR r²=0.055765); blue denotes pooled-EUR estimates.
""",
        encoding="utf-8",
    )

    source_script = Path(__file__).resolve()
    shutil.copy2(source_script, code_dir / source_script.name)
    finalizer = source_script.parent / "delivery" / "ld_interim_client" / \
        "MG_finalize_LD_sensitivity_client_package.py"
    if finalizer.exists():
        shutil.copy2(finalizer, code_dir / finalizer.name)


def validate_and_package(root: Path):
    qa = []
    for relative, expected in EXPECTED_ROWS.items():
        path = root / relative
        observed = len(read_csv(path)) if path.exists() and path.stat().st_size else 0
        ok = path.exists() and path.stat().st_size > 0 and observed == expected
        qa.append({
            "relative_path": relative,
            "expected_data_rows": expected,
            "observed_data_rows": observed,
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "status": "PASS" if ok else "FAIL",
        })
    write_csv(root / "PACKAGE_QA_REPORT.csv", qa)
    failed = [x for x in qa if x["status"] != "PASS"]
    if failed:
        raise RuntimeError(f"Package QA failed: {failed}")

    manifest_path = root / "FILE_MANIFEST.csv"
    zip_path = root / "MG_LD_interim_client_delivery_CORRECTED.zip"
    files = sorted(
        p for p in root.rglob("*")
        if p.is_file() and p not in {manifest_path, zip_path}
    )
    manifest = [
        {"relative_path": str(p.relative_to(root)),
         "size_bytes": p.stat().st_size, "md5": md5(p)}
        for p in files
    ]
    write_csv(manifest_path, manifest)
    files.append(manifest_path)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            archive.write(path, path.relative_to(root))

    # Validate the archive itself, including the two previously empty tables.
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        for relative, expected in EXPECTED_ROWS.items():
            if relative not in names:
                raise RuntimeError(f"ZIP missing {relative}")
            payload = archive.read(relative).decode("utf-8-sig")
            rows = list(csv.DictReader(payload.splitlines()))
            if len(rows) != expected:
                raise RuntimeError(
                    f"ZIP row-count failure for {relative}: {len(rows)} != {expected}"
                )
    return zip_path


def main(source: str, destination: str):
    src = Path(source).resolve()
    dst = Path(destination).resolve()
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

    reconstruct_membership(
        dst, "05_locus_summary_observed_support.csv",
        "04_locus_membership_observed_support.csv",
    )
    reconstruct_membership(
        dst, "08_locus_summary_forced_merge.csv",
        "07_locus_membership_forced_merge.csv",
    )
    write_texts(dst)
    zip_path = validate_and_package(dst)
    print(json.dumps({
        "status": "PASS",
        "corrected_zip": str(zip_path),
        "loci": 13,
        "membership_rows_each_scenario": 33,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: repair_ld_client_package.py SOURCE_DIR DESTINATION_DIR")
    main(sys.argv[1], sys.argv[2])
