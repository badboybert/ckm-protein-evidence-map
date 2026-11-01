#!/usr/bin/env python3
# 61_lead_snp_consequences.py -- Paper 3.1 epitope/coding pre-emption supplement (ADDITIVE sidecar; editor-blessed
# 2026-07-09). Same category as the localize/54-59 sensitivity sidecars: reads the FROZEN artifacts, writes ONE
# separate supplementary CSV, and does NOT change tiers/52-48-2 or re-run any phase.
#
# WHAT: systematically annotate the OWN-GENE variant CONSEQUENCE of EVERY hit's lead cis-pQTL SNP, to pre-empt the
# standard Olink-epitope-artifact reviewer criticism. PRIMARY axis = `direct_amino_acid_sequence_altering`
# (PTV / missense / inframe -- the residue changes an antibody could actually see); a `possible_splice_effect` axis
# (donor/acceptor = high, region/polypyrimidine = low) is reported SEPARATELY and is NEVER counted as protein-altering
# (editor-revised 2026-07-09: a splice_region variant that is dominantly synonymous/intronic does not alter the
# sequence). Result: direct sequence-altering lead SNPs are a MINORITY (1 PTV = IDI2-CKD Flagged; ~6 unique
# missense incl. the PCSK9 R46L positive control); both Convergent hits stay interpretable (PAM non-coding; PCSK9 R46L).
# Epitope/coding is a DISCLOSED red-flag, NOT a mechanical downgrade rule (IDI2 precedent).
#
# LOCAL-FIRST: VEP consequences are fetched from Ensembl REST ONCE and cached to a local JSON; every subsequent run
# reads the cache and is fully offline-reproducible. The cache is persisted ONLY after full validation, so a failed
# or partial fetch never leaves a bad cache. (First download may use the network, per the CKM rule.)
#
# OUT: results/repro_corrected/table_s_lead_snp_consequences.csv (one row per gene~outcome hit).
#
# KNOWN LIMITATIONS (3 residual Low findings from Codex cold-review round 16; PI-accepted 2026-07-10 as future/hardening
# only -- ZERO effect on this run's 52-hit / 49-rsID single-chunk output; recorded, not fixed):
#   L1 (display) the bucket order is not EXACTLY standard Ensembl for `inframe_*` (ranked above missense) or for
#      `coding_sequence_variant` (bucketed other_protein_altering though excluded from the direct-AA flag). No own-gene
#      hit here carries those terms, so most_severe/dominant/canonical/driving are unaffected; only future rsIDs with
#      such terms could be mis-displayed.
#   L2 (future chunking) a POST chunk of ONLY unmappable IDs returns HTTP 400 and is (deliberately) treated as a fetch
#      failure -> fail-loud, NOT silent per-SNP vep_unresolved. Kept as-is on purpose: masking every 400 as "unmappable"
#      would hide genuine request errors. The current 49-rsID table is one MIXED chunk (NOTCH2 handled); only a future
#      table >200 rsIDs whose chunk is entirely unmappable would trip this.
#   L3 (cache hardening) canonical/MANE cache validation checks that >=1 transcript carries a flag, not per-rsID; a
#      hand-edited/partially-stripped cache with one flagged transcript could still degrade some SNPs. A machine-written
#      cache is always fully annotated (fetched with canonical+mane), so this is tamper-hardening only.

# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import json, csv, os, sys, urllib.request
from collections import Counter

B      = (_P31R_EXT + "/E/20260331-CKM/Paper 3")
import subprocess
# THE GENERATION IS RESOLVED, NOT HARDCODED (was repro_corrected, the previous generation; the
# rebuild would silently read old tables and write into the old tree). Idiom from 41_2.
_r = subprocess.run([_P31R_RSCRIPT, "-e",
                     ('source("' + _P31R_EXT + '/E/20260331-CKM/Paper 3/R/_phase1_set.R"); cat(PHASE1_SET)')],
                    capture_output=True, text=True, errors="replace")
if _r.returncode != 0 or not _r.stdout.strip():
    sys.exit("61: could not resolve PHASE1_SET via R/_phase1_set.R:\n" + _r.stderr[-900:])
RC     = os.path.join(B, "results", "repro_" + _r.stdout.strip())
TIERED = os.path.join(RC, "table_s_pair_evidence_tiered.csv")
DECODE = os.path.join(RC, "decode_replication", "decode_replication_summary.csv")
CACHE  = os.path.join(RC, "lead_snp_vep_raw.json")               # one-time VEP fetch cache (offline reproducible)
OUT    = os.path.join(RC, "table_s_lead_snp_consequences.csv")   # supplementary table
VEP_URL, VEP_POST_MAX = "https://rest.ensembl.org/vep/human/id", 200   # Ensembl documents a 200-ID POST maximum
VEP_OPTS = {"canonical": 1, "mane": 1}   # request canonical + MANE-Select/Plus-Clinical flags (needed for the canonical column)

def fail(msg):
    sys.stderr.write("LEAD-SNP-ANNOT FAIL: " + msg + "\n"); sys.exit(1)

def atomic_write(path, write_fn):   # write via a .tmp then os.replace; remove the .tmp if interrupted -> no orphan .tmp
    tmp = path + ".tmp"
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as f: write_fn(f)
        os.replace(tmp, path)
    except BaseException:            # incl. KeyboardInterrupt: clean up the partial .tmp, then re-raise (fail-loud)
        if os.path.exists(tmp):
            try: os.remove(tmp)
            except OSError: pass
        raise

# ---- 1. tiered hits + lead SNPs --------------------------------------------------------------------------------
if not os.path.exists(TIERED): fail("missing tiered table: " + TIERED)
hits, rsids = [], set()
with open(TIERED, newline="", encoding="utf-8") as f:
    rd = csv.DictReader(f); _tc = rd.fieldnames or []
    if "gene" not in _tc or not ({"outcome", "outcome_name"} & set(_tc)) or "lead_snps" not in _tc:
        fail("tiered table missing required column(s) (need gene, outcome/outcome_name, lead_snps): %s" % _tc)
    for row in rd:
        oc = row.get("outcome") or row.get("outcome_name")
        snps = [s.strip() for s in (row.get("lead_snps") or "").split(",") if s.strip()]
        if not snps: fail("hit %s~%s has no lead_snps" % (row.get("gene"), oc))
        hits.append((row["gene"], oc, row.get("evidence_tier"), snps))
        rsids.update(snps)
rsids = sorted(rsids)
if not rsids: fail("no lead SNPs parsed from the tiered table")
print("[in ] %d hits, %d unique lead rsIDs" % (len(hits), len(rsids)))

# ---- 2. deCODE status -- read + FULLY validate (schema too) BEFORE any network fetch, so a bad deCODE never
#         leaves a VEP cache side-effect. gene/outcome-level; report-only join. -------------------------------------
if not os.path.exists(DECODE): fail("missing deCODE summary (required input): " + DECODE)
decode = {}   # (gene, outcome_name) -> status
with open(DECODE, newline="", encoding="utf-8") as f:
    r = csv.DictReader(f); cols = r.fieldnames or []
    stat_col = "status" if "status" in cols else next(   # prefer the exact `status` column, then fall back
        (c for c in cols if c not in ("gene", "outcome_id", "outcome_name") and
         any(k in c.lower() for k in ("class", "replicat", "status"))), None)
    if stat_col is None: fail("deCODE summary has no class/status/replication column: " + DECODE)
    ocol = "outcome_name" if "outcome_name" in cols else ("outcome" if "outcome" in cols else None)
    if ocol is None: fail("deCODE summary has no outcome_name/outcome column -- cannot join per (gene, outcome)")
    for row in r:
        dk = (row.get("gene"), row.get(ocol))
        if dk in decode: fail("deCODE summary has a DUPLICATE (gene,outcome) row %s -- non-unique status, refusing" % (dk,))
        decode[dk] = row.get(stat_col, "?")

# every FROZEN hit must have a deCODE row (per-pair summary; "no SomaScan data" = status NO_DECODE_DATA, NOT a missing
# row) -- validate COVERAGE here, BEFORE any network fetch, so a stale/incomplete deCODE never leaves a VEP cache side-effect.
_dmiss = [(g, o) for g, o, _, _ in hits if (g, o) not in decode]
if _dmiss: fail("%d frozen hit(s) have NO deCODE row (stale/incomplete deCODE or outcome-vocab mismatch): %s" % (len(_dmiss), _dmiss[:10]))

# ---- 3. VEP consequences (SELF-DESCRIBING cache-first; one-time CHUNKED fetch; cache persisted ONLY after the
#         per-hit resolvability pre-check, so a semantically-failed run leaves NO cache) ---------------------------
fetched = False
if os.path.exists(CACHE):
    try:
        with open(CACHE, encoding="utf-8") as f: cobj = json.load(f)
    except Exception as e:
        fail("VEP cache unreadable/malformed (%s) -- delete it to re-fetch: %s" % (e, CACHE))
    if not (isinstance(cobj, dict) and isinstance(cobj.get("records"), list) and cobj.get("records")):
        fail("VEP cache is empty/malformed (delete it to re-fetch): " + CACHE)
    if sorted(cobj.get("rsids", [])) != rsids:   # cache fetched for a DIFFERENT rsID set -> stale, do NOT trust it
        fail("VEP cache rsID set != current lead SNPs (stale cache) -- delete %s to re-fetch" % CACHE)
    if cobj.get("vep_opts") != VEP_OPTS:   # cache fetched with DIFFERENT VEP options (e.g. no canonical/MANE) -> stale
        fail("VEP cache options %s != current %s (stale/under-annotated cache) -- delete %s to re-fetch"
             % (cobj.get("vep_opts"), VEP_OPTS, CACHE))
    vep = cobj["records"]
    print("[vep] loaded %d cached VEP records (offline, rsID set + options verified) from %s" % (len(vep), CACHE))
else:
    # NETWORK IS NOT A FALLBACK. This generation is offline by rule; the one-time fetch the
    # local-first rule allows must be asked for, not stumbled into because a cache is missing.
    if os.environ.get("CKM_ALLOW_VEP_FETCH") != "1":
        fail("VEP cache absent and network fetch not authorized -- set CKM_ALLOW_VEP_FETCH=1 to\n"
             "  perform the one-time Ensembl fetch for this generation's lead rsIDs")
    print("[vep] one-time Ensembl VEP REST fetch (GRCh38, opts=%s) for %d rsIDs, chunks of %d ..." % (VEP_OPTS, len(rsids), VEP_POST_MAX))
    vep = []
    for i in range(0, len(rsids), VEP_POST_MAX):
        chunk = rsids[i:i + VEP_POST_MAX]
        req = urllib.request.Request(VEP_URL, data=json.dumps(dict(VEP_OPTS, ids=chunk)).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=180) as resp: part = json.load(resp)
        except Exception as e:
            fail("Ensembl VEP REST fetch failed on chunk %d (%s) -- no cache/output written" % (i // VEP_POST_MAX + 1, e))
        if not isinstance(part, list): fail("Ensembl VEP returned a non-list response chunk")
        vep.extend(part)
    if not vep: fail("Ensembl VEP returned no records")
    fetched = True

# validate (for BOTH the fresh-fetch AND the cache path) that the records actually carry canonical/MANE annotation --
# they MUST, since we fetch with canonical+mane. Catches Ensembl ignoring/changing the options, an under-annotated
# fetch, OR a hand-edited/options-mismatched cache; fail BEFORE building output or persisting a fresh cache, so the
# canonical column can never silently degrade and no under-annotated JSON is ever persisted.
if not any(("canonical" in tc or "mane_select" in tc or "mane_plus_clinical" in tc)
           for r in vep if isinstance(r, dict) for tc in r.get("transcript_consequences", [])):
    fail("VEP records carry NO canonical/MANE flag on any transcript (options not honored / under-annotated / cache "
         "hand-stripped) -- refusing to persist or use (%s)" % ("fresh fetch" if fetched else "delete %s to re-fetch" % CACHE))

# build input -> LIST of records. VEP LEGITIMATELY returns MULTIPLE records for a multiallelic/merged rsID (one per
# allele); collect them ALL and check own-gene allele-consistency below -- do NOT overwrite or fail on the duplicate.
vmap = {}
for rec in vep:
    if isinstance(rec, dict): vmap.setdefault(rec.get("input"), []).append(rec)
resolvable = {rs for rs in rsids if any(r.get("most_severe_consequence") for r in vmap.get(rs, []))}
# an rsID absent/error/unmappable is annotated per-SNP as `vep_unresolved` (editor-approved 2026-07-09); a HIT is
# unannotatable ONLY if NONE of its lead SNPs resolve -- validate that HERE, BEFORE persisting the cache, so a
# semantically-failed run leaves no cache side-effect.
_bad = [(g, o) for g, o, _, sn in hits if not any(rs in resolvable for rs in sn)]
if _bad: fail("%d hit(s) have NO VEP-resolvable lead SNP -- cannot annotate: %s" % (len(_bad), _bad[:10]))
n_unres = len(rsids) - len(resolvable)
if n_unres: print("[vep] %d/%d lead rsID(s) unresolvable by VEP-by-ID -> flagged vep_unresolved per-SNP" % (n_unres, len(rsids)))
if fetched:   # only NOW (every hit has a resolvable lead SNP) persist the SELF-DESCRIBING cache, atomically
    # The endpoint and the moment of the fetch are what identify which external annotation release
    # these consequence labels came from; without them the cache dates itself only by file mtime,
    # which any copy destroys (phase1audit F-07).
    atomic_write(CACHE, lambda f: json.dump({"rsids": rsids, "vep_opts": VEP_OPTS, "records": vep,
        "vep_url": VEP_URL, "fetched_utc": __import__("datetime").datetime.utcnow().isoformat(timespec="seconds") + "Z"}, f))
    print("[vep] cached %d VEP records (+rsID set +opts) -> %s" % (len(vep), CACHE))

# ---- consequence categorization for the HIT'S OWN gene (EPITOPE-oriented; editor-revised 2026-07-09) ------------
# The Olink epitope-artifact concern is whether the lead cis-pQTL SNP changes the assayed protein's AMINO-ACID SEQUENCE
# at the antibody epitope. So the PRIMARY flag is `direct_amino_acid_sequence_altering` (PTV / missense / inframe) --
# NOT a generic "most-severe SO term". A splice-region annotation that is dominantly synonymous or intronic does NOT
# alter the sequence and is reported as a SEPARATE, secondary `possible_splice_effect` (donor/acceptor = high,
# region/polypyrimidine = low), NEVER as "protein-altering". most-severe / dominant-base / canonical-or-MANE own-gene
# consequences are all reported so a minority-isoform effect (e.g. splice_region on 1 of many transcripts) is not
# over-weighted or misread as high-risk.
PTV_T       = {"stop_gained", "frameshift_variant", "start_lost", "stop_lost", "transcript_ablation"}
MISSENSE_T  = {"missense_variant"}
OTHER_AA_T  = {"inframe_insertion", "inframe_deletion", "conservative_inframe_insertion", "conservative_inframe_deletion",
               "protein_altering_variant", "coding_sequence_variant"}
SPLICE_HI_T = {"splice_donor_variant", "splice_acceptor_variant"}                        # canonical splice site: high concern
SPLICE_LO_T = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
               "splice_donor_region_variant", "splice_acceptor_region_variant"}          # near-splice window: secondary/low
SYNON_T     = {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}
DIRECT_AA_T = PTV_T | MISSENSE_T | (OTHER_AA_T - {"coding_sequence_variant"})   # RAW own-gene terms that DIRECTLY alter the
# AA sequence (epitope-relevant). coding_sequence_variant alone does NOT prove an AA change -> excluded from the FLAG
# (still bucketed other_protein_altering for the display column, so it is not mis-called non-coding).
PRIMARY_CHR   = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}   # primary GRCh38 assembly (exclude patch / alt-contig loci)
NON_OVERLAP_T = {"upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"}   # variant does NOT overlap the transcript
def cat_of_terms(terms):    # most-severe bucket for ONE transcript (STANDARD Ensembl VEP severity order; DISPLAY column only)
    ts = set(terms)          # splice_acceptor/donor rank ABOVE stop_gained; splice_region ranks BELOW missense (VEP convention)
    if ts & SPLICE_HI_T:  return "splice_high"
    if ts & PTV_T:        return "PTV"
    if ts & MISSENSE_T:   return "missense"
    if ts & OTHER_AA_T:   return "other_protein_altering"
    if ts & SPLICE_LO_T:  return "splice_low"
    if ts & SYNON_T:      return "synonymous"
    return "non-coding"
def base_cat_of_terms(terms):   # own-gene BASE consequence IGNORING the splice-position annotation (splice reported separately)
    ts = set(terms)              # -> reveals whether a splice_region variant is "really" synonymous / intronic (editor's point)
    if ts & PTV_T:       return "PTV"
    if ts & MISSENSE_T:  return "missense"
    if ts & OTHER_AA_T:  return "other_protein_altering"
    if ts & SYNON_T:     return "synonymous"
    return "non-coding"
# STANDARD Ensembl VEP severity order of the buckets (splice_donor/acceptor rank ABOVE stop_gained; splice_region BELOW
# missense) -- used ONLY for the transparent `most_severe_own_gene_consequence` display + the driving-SNP choice. The
# epitope interpretation does NOT depend on this ordering; it is carried by the INDEPENDENT direct_aa / splice flags,
# which are derived from the RAW own-gene terms (below), not from these collapsed buckets.
STD_SEV = {"splice_high": 8, "PTV": 7, "missense": 6, "other_protein_altering": 5, "splice_low": 4, "synonymous": 3,
           "non-coding": 2, "no_own_gene_transcript": 1, "vep_unresolved": 0}
def snp_annot(recs, gene):   # own-gene annotation for one rsID across ALL its VEP records (multiallelic-safe)
    usable = [r for r in recs if r.get("most_severe_consequence")]
    if not usable:            # rsID not resolvable by VEP-by-ID (e.g. a dbSNP variant that does NOT map to the genome)
        return {"most_severe": "vep_unresolved", "dominant": "vep_unresolved", "canonical": "vep_unresolved",
                "direct_aa": False, "splice": "none", "allele_dep": False}
    prim = [r for r in usable if str(r.get("seq_region_name")) in PRIMARY_CHR]
    usable = prim or usable   # prefer the PRIMARY-assembly mapping(s); drop patch/alt-contig duplicate loci (fall back if none primary)
    tx, base, own_terms = [], [], set()   # per-transcript bucket / splice-stripped base / UNION of raw own-gene terms
    cat_by_allele = {}                     # variant_allele -> set(buckets)  (own-gene, per transcript)
    mane_cats, canon_cats = [], []
    for rec in usable:
        for tc in rec.get("transcript_consequences", []):
            if tc.get("gene_symbol") != gene: continue
            terms = tc.get("consequence_terms", [])
            if terms and set(terms) <= NON_OVERLAP_T: continue   # variant is upstream/downstream of (does NOT OVERLAP) this transcript
            own_terms.update(terms)
            c = cat_of_terms(terms); tx.append(c); base.append(base_cat_of_terms(terms))
            cat_by_allele.setdefault(tc.get("variant_allele"), set()).add(c)
            if tc.get("mane_select") or tc.get("mane_plus_clinical"): mane_cats.append(c)
            if tc.get("canonical") == 1: canon_cats.append(c)
    if not tx:                # own gene has no transcript overlap on ANY allele/record
        return {"most_severe": "no_own_gene_transcript", "dominant": "no_own_gene_transcript",
                "canonical": "no_own_gene_transcript", "direct_aa": False, "splice": "none", "allele_dep": False}
    most_severe = max(tx, key=lambda c: STD_SEV[c])                   # STANDARD VEP most-severe bucket (splice-aware) -> transparency
    dom_base, dom_n = Counter(base).most_common(1)[0]                 # splice-STRIPPED plurality base (reveals synonymous/intronic)
    dominant = "%s(%d/%d)" % (dom_base, dom_n, len(base))
    pool = mane_cats or canon_cats                                    # prefer MANE Select/Plus-Clinical, else canonical transcript
    canonical = (max(pool, key=lambda c: STD_SEV[c]) if pool else "no_canonical_transcript")
    # the two epitope axes come from the RAW own-gene terms -> TRULY INDEPENDENT: a transcript that is BOTH missense AND
    # splice_region sets direct_aa=True AND splice=low; neither term hides the other (unlike a single collapsed bucket).
    direct_aa = bool(own_terms & DIRECT_AA_T)                         # PRIMARY epitope flag: alters the AA sequence?
    splice = "high" if (own_terms & SPLICE_HI_T) else ("low" if (own_terms & SPLICE_LO_T) else "none")
    per_allele = {a: max(cs, key=lambda c: STD_SEV[c]) for a, cs in cat_by_allele.items()}   # most-severe own-gene bucket / allele
    allele_dep = len(set(per_allele.values())) > 1                    # own-gene category differs across the rsID's alleles
    return {"most_severe": most_severe, "dominant": dominant, "canonical": canonical,
            "direct_aa": direct_aa, "splice": splice, "allele_dep": allele_dep}

# ---- 4. per-hit annotation (consequence in the hit's OWN gene) -------------------------------------------------
SPLICE_RANK = {"none": 0, "low": 1, "high": 2}
rows = []
for gene, outcome, tier, snps in hits:   # deCODE per-hit coverage already validated above (pre-fetch)
    ann = {rs: snp_annot(vmap.get(rs, []), gene) for rs in snps}     # per-lead-SNP so a multi-lead hit is auditable
    msc = sorted({rec["most_severe_consequence"] for rs in snps for rec in vmap.get(rs, []) if rec.get("most_severe_consequence")})
    drive = max(snps, key=lambda rs: STD_SEV[ann[rs]["most_severe"]])   # standard-VEP most-severe RESOLVED lead SNP drives the reported category
    d = ann[drive]
    # hit-level flags = OR / max across the hit's lead SNPs (ANY AA-altering lead instrument -> disclose)
    direct_aa = any(ann[rs]["direct_aa"] for rs in snps)
    splice = max((ann[rs]["splice"] for rs in snps), key=lambda s: SPLICE_RANK[s])
    rows.append({
        "gene": gene, "outcome": outcome, "evidence_tier": tier, "lead_snps": ", ".join(snps),   # ", " matches the frozen tiered table's delimiter
        "per_snp_consequence": "; ".join("%s=%s" % (rs, ann[rs]["most_severe"]) for rs in snps),   # per-lead-SNP most-severe
        "most_severe_own_gene_consequence": d["most_severe"],        # standard VEP most-severe, splice-aware (driving SNP)
        "dominant_own_gene_consequence": d["dominant"],              # splice-stripped plurality base + transcript counts (driving SNP)
        "canonical_or_MANE_consequence": d["canonical"],             # MANE Select / canonical transcript category (driving SNP)
        "driving_snp": drive,                                        # which lead SNP drives the reported category
        "direct_amino_acid_sequence_altering": direct_aa,            # PRIMARY epitope flag: PTV/missense/inframe on ANY lead SNP
        "possible_splice_effect": splice,                            # SEPARATE, secondary: high=donor/acceptor, low=region/polypyrimidine
        "multiallelic_allele_dependent": any(ann[rs]["allele_dep"] for rs in snps),   # own-gene cat differs across alleles -> verify lead allele
        "vep_unresolved_lead_snps": sum(1 for rs in snps if ann[rs]["most_severe"] == "vep_unresolved"),
        "most_severe_consequence_variant_global": ";".join(msc) or "NA",   # variant-global, NOT own-gene
        "decode_cross_platform": decode[(gene, outcome)],
    })

# ---- 5. write (atomic tmp->rename, orphan-.tmp-safe) + summary --------------------------------------------------
def _write_csv(f):
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
atomic_write(OUT, _write_csv)
print("[out] wrote %d hit annotations -> %s" % (len(rows), OUT))
print("[sum] most-severe own-gene categories:", dict(Counter(r["most_severe_own_gene_consequence"] for r in rows)))
aa  = [(r["gene"], r["outcome"], r["most_severe_own_gene_consequence"]) for r in rows if r["direct_amino_acid_sequence_altering"]]
spl = [(r["gene"], r["outcome"], r["possible_splice_effect"], r["dominant_own_gene_consequence"])
       for r in rows if r["possible_splice_effect"] != "none" and not r["direct_amino_acid_sequence_altering"]]
print("[sum] DIRECT amino-acid sequence-altering lead SNPs (epitope-relevant) (%d): %s" % (len(aa), aa))
print("[sum] splice-effect-only, NOT AA-altering (dominant base shown) (%d): %s" % (len(spl), spl))
