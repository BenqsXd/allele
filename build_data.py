"""Build the extra reference data for Allele: polygenic scores and population allele frequencies.

Usage:
    python build_data.py [--refresh] [--panel <raw DNA file>]
      --panel  a MyHeritage/23andMe raw file whose SNP *positions* (the chip layout, identical for every kit;
               genotypes are ignored) add ancestry-informative SNPs to the ancestry reference

Writes pgs.js:
  PGS_SCORES  published polygenic scores (PGS Catalog, harmonised to GRCh37), each variant with
              its effect allele, weight and effect-allele frequency in 1000 Genomes Europeans
  POP_FREQS   1000 Genomes phase 3 super-population allele frequencies for the same SNPs (ancestry estimate)
"""
import csv, gzip, io, json, os, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "pgs.js")
REFRESH = "--refresh" in sys.argv

# id: (display name, category key, unit or None for disease risk, short description)
SCORES = {
    "PGS000019": ("Coronary artery disease", "cardio", None, "Narrowed heart arteries, heart attack"),
    "PGS000035": ("Atrial fibrillation", "cardio", None, "Irregular heart rhythm"),
    "PGS000301": ("Systolic blood pressure", "cardio", "higher blood pressure", "Upper blood-pressure value"),
    "PGS000310": ("LDL cholesterol", "cardio", "higher LDL cholesterol", "\"Bad\" cholesterol level"),
    "PGS000038": ("Stroke", "cardio", None, "Any stroke"),
    "PGS000043": ("Venous thromboembolism", "clotting", None, "Deep-vein thrombosis and pulmonary embolism"),
    "PGS000805": ("Type 2 diabetes", "metabolic", None, "Adult-onset diabetes"),
    "PGS000298": ("Body mass index", "metabolic", "higher BMI", "Weight relative to height"),
    "PGS000297": ("Height", "trait", "taller", "Adult height"),
    "PGS000876": ("Alzheimer's disease", "neuro", None, "Late-onset Alzheimer's disease"),
    "PGS000903": ("Parkinson's disease", "neuro", None, "Parkinson's disease"),
    "PGS000193": ("Major depression", "neuro", None, "Major depressive disorder"),
    "PGS000004": ("Breast cancer", "cancer", None, "Breast cancer (women)"),
    "PGS000030": ("Prostate cancer", "cancer", None, "Prostate cancer (men)"),
    "PGS000055": ("Colorectal cancer", "cancer", None, "Bowel cancer"),
    "PGS000740": ("Lung cancer", "cancer", None, "Lung cancer"),
    "PGS003745": ("Melanoma", "skin", None, "Melanoma skin cancer"),
    "PGS000869": ("Type 1 diabetes", "immune", None, "Autoimmune diabetes"),
    "PGS000809": ("Multiple sclerosis", "immune", None, "Multiple sclerosis"),
    "PGS002293": ("Psoriasis", "immune", None, "Psoriasis"),
    "PGS002745": ("Rheumatoid arthritis", "immune", None, "Rheumatoid arthritis"),
    "PGS000316": ("Celiac disease", "immune", None, "Coeliac disease"),
    "PGS000799": ("Asthma", "immune", None, "Asthma"),
    "PGS000137": ("Glaucoma", "eye", None, "Primary open-angle glaucoma"),
    "PGS002269": ("Macular degeneration", "eye", None, "Age-related macular degeneration"),
    "PGS019957": ("Myopia", "eye", None, "Short-sightedness"),
    "PGS000859": ("Chronic kidney disease", "kidney", None, "Chronic kidney disease"),
    "PGS003329": ("Gout", "kidney", None, "Gout"),
    "PGS002583": ("Bone mineral density", "bone", "higher bone density", "Heel bone mineral density"),
}
POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]


def fetch(url, path, data=None, headers=None):
    if os.path.exists(path) and not REFRESH:
        return open(path, "rb").read()
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, data=data, headers=headers or {})
            body = urllib.request.urlopen(req, timeout=180).read()
            open(path, "wb").write(body)
            return body
        except Exception as e:
            print("  retry", attempt + 1, url[:80], e)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("download failed: " + url)


def load_score(pgs_id):
    url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh37.txt.gz"
    raw = gzip.decompress(fetch(url, os.path.join(DATA, f"{pgs_id}_GRCh37.txt.gz"))).decode()
    meta = {l[1:].split("=", 1)[0]: l.split("=", 1)[1].strip() for l in raw.splitlines() if l.startswith("#") and "=" in l}
    rows = csv.DictReader(io.StringIO("\n".join(l for l in raw.splitlines() if not l.startswith("#"))), delimiter="\t")
    out = []
    for r in rows:
        rs = (r.get("hm_rsID") or r.get("rsID") or "").strip()
        chrom, pos = (r.get("hm_chr") or "").strip(), (r.get("hm_pos") or "").strip()
        if not (chrom and pos.isdigit()):
            continue
        pk = f"{chrom}:{pos}"  # GRCh37 position: joins to the 1000 Genomes frequencies and to files without rsIDs
        if not rs.startswith("rs"):
            rs = pk
        ea = (r.get("effect_allele") or "").strip().upper()
        oa = (r.get("other_allele") or r.get("hm_inferOtherAllele") or "").strip().upper()
        if not (rs.startswith("rs") or ":" in rs) or len(ea) != 1 or len(oa) != 1 or ea not in "ACGT" or oa not in "ACGT":
            continue
        if (r.get("is_interaction") or "").lower() == "true" or (r.get("is_haplotype") or "").lower() == "true":
            continue
        try:
            w = float(r["effect_weight"])
        except (KeyError, ValueError):
            continue
        out.append([rs, ea, oa, w, pk])
    return meta, out


KG_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"


def frequencies(rsids):
    """1000 Genomes phase 3 super-population allele frequencies for the given rsIDs (biallelic SNVs only)."""
    cache = os.path.join(DATA, "pop_freqs.json")
    if os.path.exists(cache) and not REFRESH:
        freqs = json.load(open(cache))
        if all(r in freqs for r in rsids):
            return freqs
    vcf = os.path.join(DATA, "1kg_phase3_sites.vcf.gz")
    if REFRESH or not os.path.exists(vcf):
        print("downloading 1000 Genomes sites (1.4 GB)")
        urllib.request.urlretrieve(KG_URL, vcf)
    want, freqs, n = set(rsids), {}, 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            n += 1
            head = line.split("	", 8)
            key = f"{head[0]}:{head[1]}"
            if key not in want or len(head[3]) != 1 or len(head[4]) != 1:
                continue
            info = dict(kv.split("=", 1) for kv in head[7].split(";") if "=" in kv)
            ref, alt = head[3], head[4]
            freqs[key] = {p: {alt: round(float(info[p + "_AF"]), 4), ref: round(1 - float(info[p + "_AF"]), 4)} for p in POPS}
    print(f"frequencies: {len(freqs)}/{len(want)} SNPs found in {n:,} 1000 Genomes sites")
    json.dump(freqs, open(cache, "w"))
    return freqs


def main():
    os.makedirs(DATA, exist_ok=True)
    scores = {}
    for pid, (name, cat, unit, desc) in SCORES.items():
        meta, variants = load_score(pid)
        scores[pid] = dict(meta=meta, variants=variants, name=name, cat=cat, unit=unit, desc=desc)
        print(pid, name, len(variants), "variants")
    allpos = sorted({v[4] for s in scores.values() for v in s["variants"]})
    panel = []
    if "--panel" in sys.argv:
        for line in open(sys.argv[sys.argv.index("--panel") + 1], encoding="utf-8", errors="replace"):
            parts = line.replace('"', "").replace("	", ",").split(",")
            if len(parts) >= 3 and parts[0].startswith("rs") and parts[1].isdigit() and parts[2].strip().isdigit():
                panel.append(f"{parts[1]}:{parts[2].strip()}")
        print(f"panel: {len(panel)} chip positions")
    freqs = frequencies(sorted(set(allpos) | set(panel)))

    out_scores = []
    for pid, s in scores.items():
        kept = []
        for rs, ea, oa, w, pk in s["variants"]:
            eur = freqs.get(pk, {}).get("EUR", {})
            if ea not in eur and oa not in eur:
                continue  # no frequency, cannot place on the population scale
            p = eur.get(ea, 1 - eur.get(oa, 0))
            kept.append([rs, ea, oa, w, round(p, 4), pk])
        m = s["meta"]
        out_scores.append(dict(id=pid, name=s["name"], cat=s["cat"], unit=s["unit"], desc=s["desc"],
                               pmid=m.get("citation", ""), weightType=m.get("weight_type", ""),
                               total=len(s["variants"]), variants=kept))
        print(f"{pid}: {len(kept)}/{len(s['variants'])} variants with EUR frequency")

    # Allele frequencies per population for every SNP used, for the ancestry estimate (biallelic only)
    # Ancestry SNPs: every score SNP, plus the 20,000 chip SNPs whose frequencies differ most between populations
    def spread(k):
        f = freqs.get(k, {})
        if not all(x in f for x in POPS):
            return 0
        a = sorted({al for x in POPS for al in f[x]})[0]
        v = [f[x].get(a, 0) for x in POPS]
        return max(v) - min(v)
    informative = sorted((k for k in set(panel) - set(allpos) if spread(k) >= 0.15), key=spread, reverse=True)[:20000]
    pop = {}
    for rs in sorted(set(allpos) | set(informative)):
        f = freqs.get(rs, {})
        if not all(k in f for k in POPS):
            continue
        alleles = sorted({a for k in POPS for a in f[k]})
        if len(alleles) != 2:
            continue
        a = alleles[0]
        pop[rs] = [a, alleles[1]] + [round(f[k].get(a, 0), 3) for k in POPS]
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("// Generated by build_data.py. Polygenic scores: PGS Catalog (https://www.pgscatalog.org), harmonised to GRCh37.\n")
        fh.write("// Allele frequencies: 1000 Genomes phase 3 via Ensembl GRCh37. Do not edit.\n")
        fh.write("const PGS_SCORES=" + json.dumps(out_scores, separators=(",", ":")) + ";\n")
        fh.write("const POP_NAMES=" + json.dumps(POPS) + ";\n")
        fh.write("const POP_FREQS=" + json.dumps(pop, separators=(",", ":")) + ";\n")
    print(f"{len(out_scores)} scores, {len(pop)} ancestry SNPs -> {OUT} ({os.path.getsize(OUT) // 1024} KB)")


if __name__ == "__main__":
    main()
