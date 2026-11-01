"""6.1b class-rule reconciliation (orchestrator, Python reader).

Applies the executed susie_classify() rule of phase4/12_rebuild_pair_evidence.R (deposit v4, as called by
cascade/13 with PE_SUSIE_GRID = phase3_susie/susie_coloc_results.csv = v5 B1) to B1, and compares the result
with susie_class in every local copy of table_s_pair_evidence_tiered.csv. Also applies the grid's
classify_pair() (16_susie_grid.R L587-605) to the local grid output, to show what the OR rule WOULD give.

Rule (12 L518-561), on B1, which has no `setting` and no `classification` column:
  rel = first row for (gene, outcome_name)
  either converged_exp or converged_out not TRUE  -> not_evaluable
  PP.H4_susie NA                                    -> not_evaluable
  PP.H4_susie >= 0.80                               -> supportive
  PP.H3_susie > 0.50                                -> discordant
  otherwise                                         -> inconclusive
"""
import csv, hashlib, os, sys
from collections import Counter

V3 = sys.argv[1]
B1 = os.path.join(V3, 'paper3_deposit_request_v5/files/results/repro_phase1_hg38/phase3_susie/susie_coloc_results.csv')
TIERED = {
    'ledger_canonical(v3.2/_engine_root/phase4_annotation)': os.path.join(V3, 'v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv'),
    'v4_pkg_fig_results': os.path.join(V3, 'paper3_deposit_request_v4/files/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv'),
    'v4_pkg_root_pre_eas': os.path.join(V3, 'paper3_deposit_request_v4/files/results/repro_phase1_hg38/table_s_pair_evidence_tiered.csv.pre_eas_20260923T073343'),
}
GRID = os.path.join(V3, 'V4/source_data/tables_from_rds/phase3_susie__susie_grid_results.csv')


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def truthy(x):
    return str(x).strip().upper() == 'TRUE'


b1 = list(csv.DictReader(open(B1, encoding='utf-8')))
assert len(b1) == 81, len(b1)
print('B1', md5(B1)[:8], 'rows', len(b1))
fallback = {}
for r in b1:
    key = (r['gene'], r['outcome_name'])
    if key in fallback:
        continue  # first row, as rel <- sub[1]
    if not (truthy(r['converged_exp']) and truthy(r['converged_out'])):
        lab = 'not_evaluable'
    else:
        h4, h3 = num(r['PP.H4_susie']), num(r['PP.H3_susie'])
        if h4 is None:
            lab = 'not_evaluable'
        elif h4 >= 0.80:
            lab = 'supportive'
        elif h3 is not None and h3 > 0.50:
            lab = 'discordant'
        else:
            lab = 'inconclusive'
    fallback[key] = lab
print('fallback on B1:', dict(Counter(fallback.values())))

for name, p in TIERED.items():
    if not os.path.exists(p):
        print(name, 'MISSING', p)
        continue
    t = list(csv.DictReader(open(p, encoding='utf-8')))
    oc = 'outcome' if 'outcome' in t[0] else 'outcome_name'
    canon = {(r['gene'], r[oc]): r['susie_class'] for r in t}
    keys_t, keys_b = set(canon), set(fallback)
    diff = [(k, canon[k], fallback[k]) for k in sorted(keys_t & keys_b) if canon[k] != fallback[k]]
    print('\n==', name, md5(p)[:8], 'rows', len(t), 'classes', dict(Counter(canon.values())))
    print('   pairs only in tiered', sorted(keys_t - keys_b)[:5], len(keys_t - keys_b),
          '| only in B1', sorted(keys_b - keys_t)[:5], len(keys_b - keys_t))
    print('   label differences vs fallback-on-B1:', len(diff), diff[:10])

# What the grid's OR rule would give (16_susie_grid.R classify_pair), on the local grid output.
if os.path.exists(GRID):
    g = list(csv.DictReader(open(GRID, encoding='utf-8')))
    print('\nGRID', md5(GRID)[:8], 'rows', len(g), 'settings', dict(Counter(r['setting'] for r in g)))
    by = {}
    for r in g:
        by.setdefault((r['gene'], r.get('outcome_name') or r.get('outcome_id')), {})[r['setting']] = r
    cls_col = 'classification' if 'classification' in g[0] else None
    print('   grid pairs', len(by), '| has classification column:', bool(cls_col))
    if cls_col:
        gc = {k: next(iter(v.values()))[cls_col] for k, v in by.items()}
        print('   grid classification counts', dict(Counter(gc.values())))
        cmp = [(k, gc[k], fallback.get(k)) for k in sorted(gc) if k in fallback and gc[k] != fallback[k]]
        print('   grid label != fallback-on-B1:', len(cmp), cmp[:12])
