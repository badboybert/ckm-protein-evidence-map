#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite (Stage 3b, track D) -- development stage and mode of action vs the genetically implied direction.

ONE implementation, imported by both consumers (never copied):
  * scripts/P31R_build_af1.py      -> Additional file 1, Table S4;
  * R1_stage1/S1R1_06_claims_ledger/ext_af1.py -> ledger rows quoted in the text (L450).

Inputs (all read here, nothing typed):
  * the gene-level annotation master of the pipeline (ChEMBL/Open Targets/DGIdb development stage, tractability, PheWAS);
  * the tiered pair-evidence table (genetic direction per colocalization-supported pair);
  * R2_rewrite/data/P31R_devstage_audit.tsv: the audit of every gene with a direct agent at the approved or clinical stage
    (agents, ChEMBL action types and UniProt location as returned by Open Targets on the audit date; see its columns);
  * R2_rewrite/data/P31R_devstage_ot_snapshot.json: the raw Open Targets responses for those genes (P31R_af1_ot_snapshot.py);
    the audit's agents, stage, action types, effect and withdrawal mark are asserted equal to one rule applied to it
    (snapshot_rule: drug-level maximum stage in any jurisdiction; single-protein, non-nucleic-acid mechanisms only);
  * R2_rewrite/data/P31R_old_Table4_parsed.csv: the old main Table 4 as parsed by the R1 review (cross-check only).

Rules (Methods: "Only agents acting directly on the single protein counted, and each agent's mode of action was compared
with the genetically implied direction"):
  * development stage counts direct single-protein agents only: approved / in clinical development / tractable / no known
    agent or tractability (the pipeline's direct tiers 1-4); a gene whose only agents act on a
    protein group (audit: direct_single_protein == no) is re-assigned from its tractability;
  * the genetically implied direction of a pair is the sign of the MR estimate per SD of circulating protein (risk: higher
    protein, higher risk or lower eGFR; protective: the reverse);
  * an agent that lowers the protein's activity is CONCORDANT with a risk pair and OPPOSITE to a protective pair, and the
    reverse for an agent that raises it;
  * UNKNOWN when the agents act in both directions, or when the protein is a membrane receptor (UniProt: single-pass type I
    membrane protein): the plasma assay then measures a circulating form whose abundance need not track the receptor
    signalling that the agent acts on.
"""
import io, os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V3 = os.path.dirname(ROOT)
MASTER = os.path.join(V3, "v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/phase4_annotation_master.csv")
TIERED = os.path.join(V3, "v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv")
AUDIT = os.path.join(ROOT, "data", "P31R_devstage_audit.tsv")
OLDT4 = os.path.join(ROOT, "data", "P31R_old_Table4_parsed.csv")
OLDT4_MD5 = "d0b61e71523284d2af8d5e6f378572dc"   # md5 of the session-scratchpad original (ev/exec/rev1-bio/Table4_parsed.csv)

SNAPSHOT = os.path.join(ROOT, "data", "P31R_devstage_ot_snapshot.json")   # raw Open Targets responses (P31R_af1_ot_snapshot.py)

STAGE = {1: "approved", 2: "in clinical development", 3: "tractable", 4: "no known agent or tractability"}   # labels as in Additional file 2 ({SM:devstage})
WITHDRAWN = " (agent listed as withdrawn)"
RECEPTOR = re.compile(r"Single-pass type I membrane protein", re.I)
NUCLEIC = {"RNAI INHIBITOR", "ANTISENSE INHIBITOR"}                       # act on the mRNA, not the protein
LOWERS = {"INHIBITOR", "ANTAGONIST", "PROTEOLYTIC ENZYME", "HYDROLYTIC ENZYME"}
RAISES = {"AGONIST", "STABILISER", "EXOGENOUS PROTEIN"}
SALT = re.compile(r"\s+(SODIUM|HYDROCHLORIDE|MALATE|MESYLATE|TOSYLATE|ESYLATE)$")


def snapshot_rule():
    """One rule for every audited gene, applied to the raw Open Targets snapshot: the direct agents are the drugs with a
    mechanism whose target list is exactly this protein (not a nucleic-acid mechanism); the stage is 'approved' when any of
    them has a drug-level maximum clinical stage of APPROVAL (in any jurisdiction, as Open Targets records it), otherwise
    'clinical' when any is in clinical development; the listed agents are those at that stage; an agent whose row on this
    protein is at WITHDRAWAL is marked. Returns {gene: dict(stage, agents, action_types, effect, withdrawn)}."""
    import json
    snap = json.load(io.open(SNAPSHOT, encoding="utf-8"))
    out = {}
    for g, t in snap["targets"].items():
        single = []
        for r in t["drugAndClinicalCandidates"]["rows"]:
            d = r["drug"]
            acts = {m["actionType"] for m in d["mechanismsOfAction"]["rows"] if [x["approvedSymbol"] for x in m["targets"]] == [g]} - NUCLEIC
            if acts:
                single.append(dict(name=SALT.sub("", d["name"]), top=d["maximumClinicalStage"], row=r["maxClinicalStage"], acts=acts))
        if not single:
            out[g] = dict(stage=None, agents=set(), action_types=set(), effect="--", withdrawn=False)
            continue
        appr = [s for s in single if s["top"] == "APPROVAL"]
        at = appr if appr else [s for s in single if s["top"].startswith(("PHASE", "EARLY"))]
        assert at, (g, "no approved or clinical-stage direct agent")
        acts = set().union(*(s["acts"] for s in at))
        assert acts <= LOWERS | RAISES, (g, acts - LOWERS - RAISES, "an action type the effect rule does not cover")
        eff = "mixed" if (acts & LOWERS and acts & RAISES) else ("lowers" if acts & LOWERS else "raises")
        out[g] = dict(stage="approved" if appr else "clinical", agents={s["name"] for s in at}, action_types=acts, effect=eff,
                      withdrawn=any(s["row"] == "WITHDRAWAL" for s in at))
    return out, snap["queried"]


def _md5(p):
    import hashlib
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def _tf(x):
    return str(x).strip().upper() == "TRUE"


def pair_directions():
    pe = pd.read_csv(TIERED)
    # the tier table's mr_beta is eGFR* for eGFR (positive = lower eGFR; ledger L078), so risk <=> mr_beta > 0 for every outcome
    assert ((pe.mr_beta > 0) == (pe.mr_direction == "risk")).all(), "mr_direction disagrees with the sign of mr_beta"
    return pe[["gene", "outcome", "mr_beta", "mr_direction", "evidence_tier"]].copy()


def verdict(effect, direction, receptor):
    if receptor:
        return "unknown (membrane receptor; plasma form measured)"
    if effect == "mixed":
        return "unknown (agents act in both directions)"
    s = (-1 if effect == "lowers" else 1) * (1 if direction == "risk" else -1)
    return "concordant" if s < 0 else "opposite"


def table():
    """Per-gene development-stage table (one row per protein) and a summary dict."""
    m = pd.read_csv(MASTER)
    pe = pair_directions()
    au = pd.read_csv(AUDIT, sep="\t", dtype=str).fillna("")
    assert au.gene.is_unique and set(au.direct_single_protein) <= {"yes", "no"} and set(au.agent_stage) <= {"approved", "clinical"}
    assert set(au.effect_on_protein) <= {"lowers", "raises", "mixed"} and set(au.withdrawal_listed) <= {"yes", "no"}
    # the curated audit agrees with the raw Open Targets snapshot under one rule for every gene (agents, stage, effect, withdrawal)
    rule, queried = snapshot_rule()
    assert set(rule) == set(au.gene), ("snapshot and audit cover different genes", set(rule) ^ set(au.gene))
    for _, r in au.iterrows():
        x = rule[r.gene]
        if r.direct_single_protein == "no":
            assert x["stage"] is None, (r.gene, "audit says no direct agent, snapshot has", x["agents"])
            continue
        assert x["stage"] == r.agent_stage, (r.gene, "stage", r.agent_stage, x["stage"])
        assert set(a.strip() for a in r.agents.split(";")) == x["agents"], (r.gene, "agents", r.agents, sorted(x["agents"]))
        assert set(a.strip() for a in r.action_types.split(";")) == x["action_types"], (r.gene, "action types", r.action_types, x["action_types"])
        assert r.effect_on_protein == x["effect"], (r.gene, "effect", r.effect_on_protein, x["effect"])
        assert (r.withdrawal_listed == "yes") == x["withdrawn"], (r.gene, "withdrawal", r.withdrawal_listed)
        assert queried in r.audit_source, (r.gene, "audit source does not name the snapshot date", queried)
    assert set(m.gene) == set(pe.gene) and m.gene.is_unique, "annotation master and tier table cover different genes"
    # every gene the pipeline put at a direct approved or clinical stage is in the audit, at the same stage
    t12 = m[m.drug_tier_direct <= 2].set_index("gene").drug_tier_direct
    assert set(t12.index) == set(au.gene), ("audit does not cover exactly the direct tier-1/2 genes", set(t12.index) ^ set(au.gene))
    for _, r in au.iterrows():
        assert t12[r.gene] == (1 if r.agent_stage == "approved" else 2), ("audit stage differs from the pipeline", r.gene)
    # cross-check against the old main Table 4 (relocated copy, identity by md5)
    assert _md5(OLDT4) == OLDT4_MD5, "relocated old Table 4 differs from the scratchpad original"
    t4 = pd.read_csv(OLDT4).set_index("Gene")
    assert set(t4.index) == set(m.gene) and (t4.loc[m.gene, "Drug tier"].values == m.drug_tier_direct.values).all(), "old Table 4 drug tier != master"
    rows = []
    aud = au.set_index("gene")
    for _, g in m.sort_values("gene").iterrows():
        pr = pe[pe.gene == g.gene]
        tier_rank = {"Convergent": 0, "Supported": 1, "Flagged": 2}
        best = min(pr.evidence_tier, key=lambda x: tier_rank[x])
        pairs = "; ".join("%s: %s (%s)" % (o, d, t) for o, d, t in sorted(zip(pr.outcome, pr.mr_direction, pr.evidence_tier)))
        sm, ab = _tf(g.tractability_sm), _tf(g.tractability_ab)
        if g.gene in aud.index:
            a = aud.loc[g.gene]
            receptor = bool(RECEPTOR.search(a.uniprot_location))
            if a.direct_single_protein == "yes":
                st = 1 if a.agent_stage == "approved" else 2
                vs = [verdict(a.effect_on_protein, d, receptor) for d in pr.mr_direction]
                vv = vs[0] if len(set(vs)) == 1 else "differs by outcome: " + "; ".join("%s %s" % (o, v) for o, v in zip(pr.outcome, vs))
                agents, action, eff = a.agents, a.action_types, a.effect_on_protein
            else:
                st = 3 if (sm or ab) else 4
                vv, agents, action, eff = "no direct agent", "none acting on this protein alone (%s act on a protein group)" % a.agents.lower(), a.action_types, "--"
            src, note, loc = a.audit_source, a.audit_note, a.uniprot_location
            wd = a.withdrawal_listed == "yes"
            assert not wd or st == 1, (g.gene, "a withdrawal mark on a gene that is not at the approved stage")
            changed = "yes: %s -> %s" % (STAGE[int(g.drug_tier_direct)], STAGE[st]) if st != int(g.drug_tier_direct) else "no"
        else:
            st = int(g.drug_tier_direct)
            assert st in (3, 4)
            vv, agents, action, eff, src, note, loc, changed = "no agent", "--", "--", "--", "pipeline annotation (no direct agent)", "--", "--", "no"
            wd = False
        rows.append(dict(gene=g.gene, pairs=pairs, best_tier=best, stage=STAGE[st] + (WITHDRAWN if wd else ""), stage_n=st, withdrawn=wd, stage_changed=changed, agents=agents,
                         action_types=action, effect=eff, verdict=vv, uniprot_location=loc, tract_sm="yes" if sm else "no",
                         tract_ab="yes" if ab else "no", audit_source=src, audit_note=note,
                         phewas_ot=g.phewas_class_ot, phewas_ieu=g.phewas_class_ieu, phewas_atlas=g.phewas_class_atlas,
                         phewas_consensus=g.phewas_consensus))
    d = pd.DataFrame(rows)
    ap = d[d.stage_n == 1]
    summ = dict(approved=int((d.stage_n == 1).sum()), clinical=int((d.stage_n == 2).sum()), tractable=int((d.stage_n == 3).sum()),
                none=int((d.stage_n == 4).sum()),
                approved_opposite=sorted(ap[ap.verdict == "opposite"].gene), approved_concordant=sorted(ap[ap.verdict == "concordant"].gene),
                approved_unknown=sorted(ap[ap.verdict.str.startswith("unknown")].gene),
                restaged=sorted(d[d.stage_changed != "no"].gene), approved_withdrawn=sorted(ap[ap.withdrawn].gene))
    assert summ["approved"] + summ["clinical"] + summ["tractable"] + summ["none"] == len(d) == len(m)
    assert len(summ["approved_opposite"]) + len(summ["approved_concordant"]) + len(summ["approved_unknown"]) == summ["approved"], \
        "an approved-drug gene has a verdict that differs by outcome; extend the summary"
    return d, summ


def phewas():
    m = pd.read_csv(MASTER)
    flagged = m.phewas_consensus.str.startswith("flag")
    otq = m[m.phewas_class_ot != "not_queryable"]
    ot_flag = otq.phewas_class_ot.isin(["pleiotropy_flagged", "safety_flagged"])
    return dict(consensus_flagged=int(flagged.sum()), genes=len(m), ot_flagged=int(ot_flag.sum()), ot_queryable=len(otq))


if __name__ == "__main__":
    d, s = table()
    print(s)
    print(phewas())
    print(d[d.stage_n <= 2][["gene", "stage", "agents", "effect", "verdict"]].to_string())
