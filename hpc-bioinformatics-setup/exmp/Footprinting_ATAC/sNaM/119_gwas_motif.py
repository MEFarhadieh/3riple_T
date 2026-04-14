#!/usr/bin/env python3
"""
GWAS → RUNX3 Motif Pipeline
============================
گام ۱: parse فایل‌های GWAS
گام ۴: تست motif binding RUNX3

نیازمندی‌ها:
  - base env : pandas, pyfaidx
  - meme_env : fimo (فقط برای گام ۴)

اجرا:
  python3 gwas_runx3_pipeline.py --step 1   # parse فایل‌های GWAS
  python3 gwas_runx3_pipeline.py --step 4a  # استخراج توالی‌ها + موتیف
  python3 gwas_runx3_pipeline.py --step 4b  # مقایسه ref vs alt (بعد از fimo)

  # fimo را دستی در meme_env اجرا کن:
  # conda activate meme_env
  # fimo --thresh 1e-4 --oc runx3_motif/fimo_ref RUNX3.meme runx3_motif/snps_ref.fa
  # fimo --thresh 1e-4 --oc runx3_motif/fimo_alt RUNX3.meme runx3_motif/snps_alt.fa
"""

import argparse
import pandas as pd
import re
import sys
from pathlib import Path

# ─── تنظیمات ─────────────────────────────────────────────────────────────────

GWAS_DIR  = Path("/mnt/archive/farhadie/tn5_bias/skin_Mphage/gwas")
REF_FA    = "/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
MOTIF_DB  = "/mnt/archive/farhadie/ref/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
OUTDIR    = GWAS_DIR / "runx3_motif"
FLANK     = 50  # bp هر طرف SNP

FILES = {
    "atopic_eczema":                 "atopic_eczema.tsv",
    "pain_sensitivity":              "pain_sensitivity.tsv",
    "peripheral_neuropathy":         "peripheral_neuropathy.tsv",
    "pruritus":                      "pruritus.tsv",
    "sensory_peripheral_neuropathy": "sensory_peripheral_neuropathy.tsv",
}

# ═══════════════════════════════════════════════════════════════════════════════
# گام ۱: parse فایل‌های GWAS
# ═══════════════════════════════════════════════════════════════════════════════

def extract_rsid_and_chrpos(val):
    val = str(val).strip()
    m = re.match(r'^(rs\d+)', val)
    if m:
        a = re.search(r'-([A-Za-z?*]+)$', val)
        return m.group(1), None, (a.group(1) if a else None)
    m2 = re.match(r'^chr([XYMxy\d]+):(\d+)-?([A-Za-z?*]*)$', val)
    if m2:
        return None, f"{m2.group(1)}:{m2.group(2)}", (m2.group(3) or None)
    return None, None, None

def parse_location(s):
    if pd.isna(s) or str(s).strip() in ('-', '', 'NA'):
        return None, None
    s = str(s).strip().lstrip('chr')
    if ':' in s:
        p = s.split(':')
        return p[0], p[1]
    return None, None

def parse_gwas_file(filepath, label):
    print(f"  [*] {filepath.name}")
    try:
        df = pd.read_csv(filepath, sep='\t', dtype=str, on_bad_lines='skip')
    except Exception as e:
        print(f"      [!] خطا: {e}")
        return pd.DataFrame()
    print(f"      {len(df):,} ردیف")

    ext = df['riskAllele'].apply(
        lambda x: pd.Series(extract_rsid_and_chrpos(x),
                             index=['rsid', 'chrpos', 'risk_allele']))
    df = pd.concat([df, ext], axis=1)

    loc = df['locations'].apply(
        lambda x: pd.Series(parse_location(x), index=['chr', 'pos']))
    df = pd.concat([df, loc], axis=1)

    mask = df['chr'].isna() & df['chrpos'].notna()
    if mask.any():
        fb = df.loc[mask, 'chrpos'].apply(
            lambda x: pd.Series(parse_location(x), index=['chr', 'pos']))
        df.loc[mask, ['chr', 'pos']] = fb.values

    df['pvalue']     = pd.to_numeric(df['pValue'], errors='coerce')
    df['trait_file'] = label
    df['has_rsid']   = df['rsid'].notna()

    keep = ['rsid', 'chrpos', 'chr', 'pos', 'risk_allele',
            'pvalue', 'beta', 'ci', 'orValue',
            'mappedGenes', 'traitName', 'efoTraits',
            'accessionId', 'pubmedId', 'trait_file', 'has_rsid']
    return df[[c for c in keep if c in df.columns]].copy()

def step1_parse():
    print("\n=== گام ۱: parse فایل‌های GWAS ===\n")
    all_dfs = []
    for label, fname in FILES.items():
        fpath = GWAS_DIR / fname
        if not fpath.exists():
            print(f"  [!] نیست: {fpath}")
            continue
        df = parse_gwas_file(fpath, label)
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        sys.exit("[!] هیچ فایلی parse نشد")

    merged = pd.concat(all_dfs, ignore_index=True)

    print(f"\n{'='*50}")
    print(f"کل ردیف:          {len(merged):,}")
    print(f"rsID دار:         {merged['has_rsid'].sum():,}")
    print(f"chr:pos فقط:      {(~merged['has_rsid']).sum():,}")
    print(f"rsID منحصربه‌فرد: {merged['rsid'].nunique():,}")
    print(f"\nتوزیع trait:")
    print(merged.groupby('trait_file').agg(
        total=('rsid', 'count'),
        unique_rsid=('rsid', 'nunique'),
        best_p=('pvalue', 'min')
    ).to_string())

    merged.to_csv(GWAS_DIR / "snp_list_step1_all.csv", index=False)
    rsid_df = (merged[merged['has_rsid']]
               .drop_duplicates(subset=['rsid', 'trait_file'])
               .sort_values('pvalue'))
    rsid_df.to_csv(GWAS_DIR / "snp_list_step1_rsid.csv", index=False)
    chrpos_df = merged[~merged['has_rsid'] & merged['chrpos'].notna()]
    if not chrpos_df.empty:
        chrpos_df.to_csv(GWAS_DIR / "snp_list_step1_need_rsid.csv", index=False)

    print(f"\n[+] snp_list_step1_all.csv  ({len(merged):,} ردیف)")
    print(f"[+] snp_list_step1_rsid.csv ({len(rsid_df):,} ردیف)")
    if not chrpos_df.empty:
        print(f"[+] snp_list_step1_need_rsid.csv ({len(chrpos_df):,} ردیف)")


# ═══════════════════════════════════════════════════════════════════════════════
# گام ۴-الف: استخراج توالی‌ها + موتیف RUNX3
# ═══════════════════════════════════════════════════════════════════════════════

def extract_runx3_motif():
    content = open(MOTIF_DB).read()
    header  = re.search(r'^(MEME version.+?)(?=MOTIF)', content, re.DOTALL)
    header_str = header.group(1) if header else (
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
    )
    match = re.search(
        r'(MOTIF\s+\S*MA0684\S*\s+RUNX3.+?)(?=MOTIF|\Z)',
        content, re.DOTALL | re.IGNORECASE)
    if not match:
        match = re.search(
            r'(MOTIF\s+\S+\s+RUNX3\b.+?)(?=MOTIF|\Z)',
            content, re.DOTALL | re.IGNORECASE)
    if match:
        out = OUTDIR / "RUNX3.meme"
        out.write_text(header_str + match.group(1))
        print(f"[+] RUNX3 motif ذخیره شد: {out}")
        return True
    else:
        all_runx = re.findall(r'MOTIF\s+(\S+)\s+(RUNX\S*)', content, re.IGNORECASE)
        print("[!] RUNX3 پیدا نشد. موجود:")
        for mid, mname in all_runx:
            print(f"    {mid}  {mname}")
        return False

def extract_sequences():
    import pyfaidx
    df = pd.read_csv(GWAS_DIR / "snp_list_step1_all.csv", dtype=str)
    df = df[df['chr'].notna() & df['pos'].notna()].copy()
    df['pos_int'] = pd.to_numeric(df['pos'], errors='coerce')
    df = df[df['pos_int'].notna()].copy()
    df['pos_int'] = df['pos_int'].astype(int)
    df['risk_allele_clean'] = df['risk_allele'].str.strip().str.upper()
    df = df[df['risk_allele_clean'].str.match(r'^[ACGT]$', na=False)].copy()
    print(f"SNP با آلل معتبر: {len(df):,}")

    genome = pyfaidx.Fasta(REF_FA)
    print(f"chr names نمونه: {list(genome.keys())[:3]}")

    ref_seqs, alt_seqs, snp_info, skipped = [], [], [], 0

    for _, row in df.iterrows():
        chrom_raw  = str(row['chr']).strip()
        chrom_with = f"chr{chrom_raw}" if not chrom_raw.startswith('chr') else chrom_raw
        chrom      = chrom_with if chrom_with in genome else \
                     (chrom_raw if chrom_raw in genome else None)
        if not chrom:
            skipped += 1
            continue

        pos   = row['pos_int']
        start = max(0, pos - FLANK - 1)
        end   = pos + FLANK

        try:
            seq = str(genome[chrom][start:end]).upper()
        except Exception:
            skipped += 1
            continue

        if len(seq) < FLANK:
            skipped += 1
            continue

        snp_idx     = pos - 1 - start
        ref_allele  = seq[snp_idx]
        risk_allele = row['risk_allele_clean']
        rsid        = row['rsid'] if pd.notna(row.get('rsid')) else f"{chrom}_{pos}"
        trait       = str(row['trait_file'])

        ref_seqs.append(f">{rsid}_ref\n{seq}")
        alt_seq = seq[:snp_idx] + risk_allele + seq[snp_idx+1:]
        alt_seqs.append(f">{rsid}_alt\n{alt_seq}")

        snp_info.append({
            'rsid': rsid, 'chrom': chrom, 'pos': pos,
            'ref_allele': ref_allele, 'risk_allele': risk_allele,
            'trait_file': trait,
            'pvalue': row.get('pvalue', '.'),
            'mappedGenes': row.get('mappedGenes', '.'),
            'traitName': row.get('traitName', '.'),
        })

    print(f"استخراج شد: {len(ref_seqs):,}  |  رد شد: {skipped:,}")
    (OUTDIR / "snps_ref.fa").write_text('\n'.join(ref_seqs) + '\n')
    (OUTDIR / "snps_alt.fa").write_text('\n'.join(alt_seqs) + '\n')
    pd.DataFrame(snp_info).to_csv(OUTDIR / "snps_for_fimo.csv", index=False)
    print(f"[+] فایل‌ها آماده: {OUTDIR}")

def step4a_prepare():
    print("\n=== گام ۴-الف: آماده‌سازی توالی‌ها و موتیف ===\n")
    OUTDIR.mkdir(exist_ok=True)
    ok = extract_runx3_motif()
    if not ok:
        sys.exit("[!] موتیف RUNX3 پیدا نشد")
    extract_sequences()

    print("\n" + "="*50)
    print("حالا در meme_env اجرا کن:")
    print(f"  conda activate meme_env")
    print(f"  fimo --thresh 1e-4 --oc {OUTDIR}/fimo_ref {OUTDIR}/RUNX3.meme {OUTDIR}/snps_ref.fa")
    print(f"  fimo --thresh 1e-4 --oc {OUTDIR}/fimo_alt {OUTDIR}/RUNX3.meme {OUTDIR}/snps_alt.fa")
    print(f"  conda deactivate")
    print(f"  python3 gwas_runx3_pipeline.py --step 4b")


# ═══════════════════════════════════════════════════════════════════════════════
# گام ۴-ب: مقایسه ref vs alt
# ═══════════════════════════════════════════════════════════════════════════════

def step4b_compare():
    print("\n=== گام ۴-ب: مقایسه ref vs alt ===\n")

    def load_fimo(path):
        tsv = path / "fimo.tsv"
        if not tsv.exists():
            return pd.DataFrame()
        df = pd.read_csv(tsv, sep='\t', comment='#').dropna(subset=['sequence_name'])
        df['rsid'] = df['sequence_name'].str.replace(r'_(ref|alt)$', '', regex=True)
        return df

    ref_df = load_fimo(OUTDIR / "fimo_ref")
    alt_df = load_fimo(OUTDIR / "fimo_alt")
    print(f"FIMO hits — ref: {len(ref_df):,}  |  alt: {len(alt_df):,}")

    def best_score(fimo_df, allele):
        if fimo_df.empty:
            return pd.DataFrame(columns=['rsid', f'score_{allele}', f'pval_{allele}'])
        g = fimo_df.groupby('rsid')['score'].max().reset_index().rename(
            columns={'score': f'score_{allele}'})
        p = fimo_df.groupby('rsid')['p-value'].min().reset_index().rename(
            columns={'p-value': f'pval_{allele}'})
        return g.merge(p, on='rsid')

    ref_best = best_score(ref_df, 'ref')
    alt_best = best_score(alt_df, 'alt')

    all_rsids = pd.DataFrame(
        {'rsid': list(set(ref_df['rsid'].tolist()) | set(alt_df['rsid'].tolist()))})
    merged = (all_rsids
              .merge(ref_best, on='rsid', how='left')
              .merge(alt_best, on='rsid', how='left'))
    merged['score_ref'] = merged.get('score_ref', 0).fillna(0)
    merged['score_alt'] = merged.get('score_alt', 0).fillna(0)
    merged['delta']     = merged['score_alt'] - merged['score_ref']

    def effect(row):
        r, a = row['score_ref'] > 0, row['score_alt'] > 0
        d = row['delta']
        if r and not a:   return 'disruption'
        if a and not r:   return 'gain'
        if r and a:
            if d < -2:    return 'weakened'
            if d >  2:    return 'strengthened'
            return 'neutral_hit'
        return 'no_hit'

    merged['effect'] = merged.apply(effect, axis=1)

    snp_info = pd.read_csv(OUTDIR / "snps_for_fimo.csv", dtype=str)
    snp_info['pvalue'] = pd.to_numeric(snp_info['pvalue'], errors='coerce')
    result = merged.merge(
        snp_info[['rsid', 'chrom', 'pos', 'ref_allele', 'risk_allele',
                  'pvalue', 'mappedGenes', 'traitName', 'trait_file']],
        on='rsid', how='left')
    result.sort_values('delta', key=abs, ascending=False, inplace=True)

    result.to_csv(OUTDIR / "runx3_annotated.csv", index=False)
    candidates = result[result['effect'].isin(
        ['disruption', 'gain', 'weakened', 'strengthened'])].copy()
    candidates.to_csv(OUTDIR / "runx3_candidates.csv", index=False)

    print(f"\nتوزیع اثرات:")
    print(result['effect'].value_counts().to_string())
    print(f"\nکاندیداهای نهایی: {len(candidates):,}")
    print("\n--- کاندیداهای برتر ---")
    cols = ['rsid', 'chrom', 'pos', 'ref_allele', 'risk_allele',
            'effect', 'delta', 'score_ref', 'score_alt',
            'pvalue', 'mappedGenes', 'traitName', 'trait_file']
    cols = [c for c in cols if c in candidates.columns]
    print(candidates[cols].head(20).to_string(index=False))
    print(f"\n[+] runx3_annotated.csv  — همه SNP ها")
    print(f"[+] runx3_candidates.csv — فقط کاندیداها")


# ═══════════════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GWAS → RUNX3 Motif Pipeline")
    parser.add_argument('--step', required=True,
                        choices=['1', '4a', '4b'],
                        help='گام: 1=parse GWAS  4a=آماده‌سازی  4b=مقایسه')
    args = parser.parse_args()

    if args.step == '1':
        step1_parse()
    elif args.step == '4a':
        step4a_prepare()
    elif args.step == '4b':
        step4b_compare()
