#!/usr/bin/env python
"""Minimal pure-Python remote tabix slicer (no pysam/htslib).
Fetches a genomic region from a bgzipped+tabix-indexed GWAS sumstats over HTTP
using the .tbi bin index + HTTP Range requests. Built for the Chen-CLSA metaboQTL
deposits on the GWAS-Catalog FTP (hg38, contigs '1'..'22', seq=col1 pos=col2).

Validated against a known sentinel (see __main__): arginine GCST90200372
rs17788484 chr6:131573218 beta=-0.381004 from chen_instruments.tsv.
"""
import gzip, struct, zlib, io, os, sys, urllib.request

def parse_tbi(tbi_path):
    raw = gzip.open(tbi_path, "rb").read()
    assert raw[:4] == b"TBI\x01", "not a tabix index"
    off = 4
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[off:off+32]); off += 32
    names = raw[off:off+l_nm].split(b"\x00"); off += l_nm
    names = [n.decode() for n in names if n]
    idx = {}  # contig -> (bins dict {bin:[(cbeg,cend),...]}, linear[ioff,...])
    for ci in range(n_ref):
        n_bin, = struct.unpack("<i", raw[off:off+4]); off += 4
        bins = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack("<Ii", raw[off:off+8]); off += 8
            chunks = []
            for _ in range(n_chunk):
                cb, ce = struct.unpack("<QQ", raw[off:off+16]); off += 16
                chunks.append((cb, ce))
            bins[b] = chunks
        n_intv, = struct.unpack("<i", raw[off:off+4]); off += 4
        linear = list(struct.unpack("<%dQ" % n_intv, raw[off:off+8*n_intv])); off += 8*n_intv
        idx[names[ci]] = (bins, linear)
    return idx, dict(col_seq=col_seq, col_beg=col_beg, col_end=col_end, skip=skip)

def reg2bins(beg, end):  # 0-based half-open
    end -= 1
    lst = [0]
    for sh, k0 in ((26,1),(23,9),(20,73),(17,585),(14,4681)):
        lst += list(range(k0 + (beg>>sh), k0 + (end>>sh) + 1))
    return lst

def _http_range(url, start, end):
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:
            if attempt == 3: raise
    return b""

def _bgzf_blocks(comp):
    """yield decompressed payload of each bgzf block in `comp` (starting at a block boundary)."""
    p = 0; n = len(comp)
    while p < n - 18:
        if comp[p:p+2] != b"\x1f\x8b": break
        xlen = struct.unpack("<H", comp[p+10:p+12])[0]
        bsize = None; q = p+12; xend = q+xlen
        while q < xend:
            si1, si2, slen = comp[q], comp[q+1], struct.unpack("<H", comp[q+2:q+4])[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack("<H", comp[q+4:q+6])[0]
            q += 4+slen
        if bsize is None: break
        block = comp[p:p+bsize+1]
        try:
            yield zlib.decompress(block[12+xlen:-8], -15)
        except Exception:
            pass
        p += bsize+1

def fetch_region(url, idx, meta, contig, lo, hi, margin=0):
    """return list of row-strings (tab-sep) with col_beg position in [lo-margin, hi+margin].
       lo/hi are 1-based inclusive (genomic). Internally treated 0-based for binning."""
    if contig not in idx:
        contig2 = contig.replace("chr","") if contig.startswith("chr") else "chr"+contig
        contig = contig2 if contig2 in idx else contig
    bins, linear = idx[contig]
    qlo = max(0, lo-1-margin); qhi = hi+margin
    want = set(reg2bins(qlo, qhi))
    # linear-index filter: min offset for the start window
    li = qlo >> 14
    min_lin = linear[li] if li < len(linear) else (linear[-1] if linear else 0)
    chunks = []
    for b in want:
        for cb, ce in bins.get(b, []):
            if ce > min_lin:
                chunks.append((cb, ce))
    if not chunks: return []
    cstart = min(cb>>16 for cb,ce in chunks)
    cend = max(ce>>16 for cb,ce in chunks)
    comp = _http_range(url, cstart, cend + 65536)
    csq, cbg = meta["col_seq"]-1, meta["col_beg"]-1
    out = []
    buf = b""
    for payload in _bgzf_blocks(comp):
        buf += payload
    text = buf.decode("utf-8", "replace")
    for line in text.split("\n"):
        line = line.rstrip("\r")
        if not line or line[0] in ("#",):
            continue
        f = line.split("\t")
        if len(f) <= max(csq, cbg):
            continue
        if f[csq] != contig and f[csq] != contig.replace("chr",""):
            continue
        try:
            pos = int(f[cbg])
        except ValueError:
            continue
        if qlo <= pos-1 < qhi:
            out.append(line)
    return out

def get_header(url):
    """fetch+decompress the first bgzf block, return the header column list.
    Chen deposits do NOT share a uniform column order -> always map by name."""
    comp = _http_range(url, 0, 65535)
    for payload in _bgzf_blocks(comp):
        line = payload.decode("utf-8", "replace").split("\n")[0].rstrip("\r")
        return line.split("\t")
    return []

def open_index(tbi_url, cache_path):
    if not os.path.exists(cache_path):
        urllib.request.urlretrieve(tbi_url, cache_path)
    return parse_tbi(cache_path)

if __name__ == "__main__":
    # VALIDATION: arginine GCST90200372, sentinel rs17788484 chr6:131573218, beta=-0.381004
    base = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90200001-GCST90201000/GCST90200372/harmonised/GCST90200372.h.tsv.gz"
    tbi = base + ".tbi"
    cache = os.path.join(os.path.dirname(__file__), "..", "..", "data", "metabolome", "chen_regional", "GCST90200372.tbi")
    idx, meta = parse_tbi(cache)
    print("meta:", meta)
    # print the harmonized header (first data line region) to learn columns
    rows = fetch_region(base, idx, meta, "6", 131573218-2000, 131573218+2000)
    print(f"rows in ±2kb of sentinel: {len(rows)}")
    for r in rows:
        if "rs17788484" in r:
            print("SENTINEL ROW:", r)
