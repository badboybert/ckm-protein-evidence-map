#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Put the reference list in order of first citation, and add the two references
the metabolite Methods deferred.

Why this exists. The manuscript is Vancouver-numbered, and Vancouver numbers
entries in the order they are first cited. This list did not: measured on the
built docx, first appearance ran 1, 2, 3, **33, 34**, 4 ... 27, **20**, 28 ... 32.
Nothing dangled and nothing was uncited -- which is exactly why it survived. Those
are ABSENCE checks, and a list can pass all of them while being in the wrong order.

What it does, in one pass over the built docx:

  1. appends the NEW_REFS entries (the two the metabolite Methods deferred, then the
     L2-L4 update's), applies ENTRY_FIXES, and (PASS A, 2026-09-30) prints a DOI on
     every stage-1 entry that lacked one (DOI_ADDS), asserting the house form on all;
  2. reads every bracketed citation in document order, before the References
     heading, and builds old -> new from first appearance;
  3. asserts the map is a BIJECTION over 1..N before touching anything;
  4. rewrites every in-text citation and re-sorts the entry paragraphs.

Run it AFTER apply_blocks.py -- the inserted blocks carry citations of their own,
and a renumber that ran on the stage-1 file would not see them. `build_all.py`
owns that ordering.

  python renumber_refs.py <built.docx>
  python renumber_refs.py <built.docx> --check    # report, change nothing
"""
from __future__ import print_function
import os
import re
import sys

import docx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from docx_edit import para_text, splice  # noqa: E402

# The two references the metabolite Methods deferred by name. Resolved against
# PubMed, not typed from a draft; the year is the ISSUE year in both cases, which
# is also the entry year here, so nothing had to be adjudicated.
#   [A] PMID 36635386 -- Nat Genet 2023;55(1):44-53   (the CLSA metabolome GWAS)
#   [B] PMID 39837327 -- Cell Genom 2025;5(2):100743  (Han Chinese serum metabolome)
NEW_REFS = [
    (u"Chen Y, Lu T, Pettersson-Kymmer U, et al. Genomic atlas of the plasma metabolome "
     u"prioritizes metabolites implicated in human diseases. Nat Genet. 2023;55(1):44–53. "
     u"doi:10.1038/s41588-022-01270-1."),
    (u"Cheng C, Xu F, Pan XF, et al. Genetic mapping of serum metabolome to chronic diseases "
     u"among Han Chinese. Cell Genom. 2025;5(2):100743. doi:10.1016/j.xgen.2024.100743."),
    # [C] PMID 29412048 -- the urea-cycle step the ARG1 node asserts. PASS A (2026-09-30): PubMed lists FIVE authors
    #     (Caldwell RW, Rodriguez PC, Toque HA, Narayanan SP, Caldwell RB); the list names every author of a <= 6-author
    #     paper (17 entries do), so 'et al.' hid two of them.
    (u"Caldwell RW, Rodriguez PC, Toque HA, Narayanan SP, Caldwell RB. Arginase: a multifaceted enzyme important in "
     u"health and disease. Physiol Rev. 2018;98(2):641\u2013665. doi:10.1152/physrev.00037.2016."),
    # [D] PMID 2303490 -- the renal dipeptidase the kidney-gene anchor asserts (PubMed carries no DOI; the DOI is the
    #     Crossref record of J Biol Chem 265(7):3992-3995, resolved 2026-09-30 with volume/issue/pages/title asserted)
    (u"Adachi H, Tawaragi Y, Inuzuka C, et al. Primary structure of human microsomal dipeptidase "
     u"deduced from molecular cloning. J Biol Chem. 1990;265(7):3992\u20133995. doi:10.1016/s0021-9258(19)39692-9."),
    # ---- the L2-L4 update (2026-09-27): 27 references the native census, the bridge fine-mapping and the metabolome
    # Discussion cite. Every entry typed from its PubMed record (efetch XML: title, author list, ISOAbbreviation, issue
    # year from <JournalIssue><PubDate>, volume, issue, pages, DOI), never from the drafting spec; the one preprint
    # (Kanai 2021, the Biobank Japan CAD GWAS that pheweb.jp cites for 32,512 cases / 146,214 controls) from Crossref.
    # Keys in NEW_REF_KEYS, which l2l4_facts.refs() reads to cite them by stage-1 number.
    (u"Spracklen CN, Horikoshi M, Kim YJ, et al. Identification of type 2 diabetes loci in 433,540 East Asian "
     u"individuals. Nature. 2020;582(7811):240\u2013245. doi:10.1038/s41586-020-2263-3."),                 # PMID 32499647
    (u"Mishra A, Malik R, Hachiya T, et al. Stroke genetics informs drug discovery and risk prediction across "
     u"ancestries. Nature. 2022;611(7934):115\u2013123. doi:10.1038/s41586-022-05165-3."),                 # PMID 36180795
    (u"Chen HL, Chiang HY, Chang DR, et al. Discovery and prioritization of genetic determinants of kidney function in "
     u"297,355 individuals from Taiwan and Japan. Nat Commun. 2024;15:9317. doi:10.1038/s41467-024-53516-7."),  # PMID 39472450
    (u"Kanai M, Ulirsch JC, Karjalainen J, et al. Insights from complex trait fine-mapping across diverse populations. "
     u"medRxiv. 2021; doi:10.1101/2021.09.03.21262975."),                                                   # no PMID
    (u"Miyazawa K, Ito K, Ito M, et al. Cross-ancestry genome-wide analysis of atrial fibrillation unveils disease "
     u"biology and enables cardioembolic risk prediction. Nat Genet. 2023;55(2):187\u2013197. "
     u"doi:10.1038/s41588-022-01284-9."),                                                                  # PMID 36653681
    (u"Enzan N, Miyazawa K, Koyama S, et al. Genome-wide analysis of heart failure yields insights into disease "
     u"heterogeneity and enables prognostic prediction in the Japanese population. Nat Commun. 2025;16:9680. "
     u"doi:10.1038/s41467-025-64659-6."),                                                                  # PMID 41184235
    (u"Sakaue S, Kanai M, Tanigawa Y, et al. A cross-population atlas of genetic associations for 220 human phenotypes. "
     u"Nat Genet. 2021;53(10):1415\u20131424. doi:10.1038/s41588-021-00931-x."),                          # PMID 34594039
    (u"Sin YY, Baron G, Schulze A, Funk CD. Arginase-1 deficiency. J Mol Med (Berl). 2015;93(12):1287\u20131296. "
     u"doi:10.1007/s00109-015-1354-3."),                                                                   # PMID 26467175
    (u"McIntyre T, Curthoys NP. Renal catabolism of glutathione. Characterization of a particulate rat renal dipeptidase "
     u"that catalyzes the hydrolysis of cysteinylglycine. J Biol Chem. 1982;257(20):11915\u201311921. "
     u"doi:10.1016/s0021-9258(18)33653-6."),                                                     # PMID 6126478 (DOI: Crossref)
    (u"Adachi H, Kubota I, Okamura N, et al. Purification and characterization of human microsomal dipeptidase. "
     u"J Biochem. 1989;105(6):957\u2013961. doi:10.1093/oxfordjournals.jbchem.a122787."),                 # PMID 2768222
    (u"Habib GM, Shi ZZ, Cuevas AA, Guo Q, Matzuk MM, Lieberman MW. Leukotriene D4 and cystinyl-bis-glycine metabolism "
     u"in membrane-bound dipeptidase-deficient mice. Proc Natl Acad Sci U S A. 1998;95(9):4859\u20134863. "
     u"doi:10.1073/pnas.95.9.4859."),                                                                      # PMID 9560193
    (u"Hartiala JA, Tang WH, Wang Z, et al. Genome-wide association study and targeted metabolomics identifies "
     u"sex-specific association of CPS1 with coronary artery disease. Nat Commun. 2016;7:10558. "
     u"doi:10.1038/ncomms10558."),                                                                         # PMID 26822151
    (u"Wyss M, Kaddurah-Daouk R. Creatine and creatinine metabolism. Physiol Rev. 2000;80(3):1107\u20131213. "
     u"doi:10.1152/physrev.2000.80.3.1107."),                                                              # PMID 10893433
    (u"Musunuru K, Strong A, Frank-Kamenetsky M, et al. From noncoding variant to phenotype via SORT1 at the 1p13 "
     u"cholesterol locus. Nature. 2010;466(7307):714\u2013719. doi:10.1038/nature09266."),                 # PMID 20686566
    (u"Newton-Cheh C, Johnson T, Gateva V, et al. Genome-wide association study identifies eight loci associated with "
     u"blood pressure. Nat Genet. 2009;41(6):666\u2013676. doi:10.1038/ng.361."),                          # PMID 19430483
    (u"Mallat Z, Gojova A, Marchiol-Fournigault C, et al. Inhibition of transforming growth factor-beta signaling "
     u"accelerates atherosclerosis and induces an unstable plaque phenotype in mice. Circ Res. 2001;89(10):930\u2013934. "
     u"doi:10.1161/hh2201.099415."),                                                                       # PMID 11701621
    (u"Grainger DJ. Transforming growth factor beta and atherosclerosis: so far, so good for the protective cytokine "
     u"hypothesis. Arterioscler Thromb Vasc Biol. 2004;24(3):399\u2013404. doi:10.1161/01.ATV.0000114567.76772.33."),  # PMID 14699019
    (u"IL6R Genetics Consortium Emerging Risk Factors Collaboration; Sarwar N, Butterworth AS, Freitag DF, et al. "
     u"Interleukin-6 receptor pathways in coronary heart disease: a collaborative meta-analysis of 82 studies. Lancet. "
     u"2012;379(9822):1205\u20131213. doi:10.1016/S0140-6736(11)61931-4."),                               # PMID 22421339
    (u"Kurosu H, Choi M, Ogawa Y, et al. Tissue-specific expression of betaKlotho and fibroblast growth factor (FGF) "
     u"receptor isoforms determines metabolic activity of FGF19 and FGF21. J Biol Chem. 2007;282(37):26687\u201326695. "
     u"doi:10.1074/jbc.M704165200."),                                                                      # PMID 17623664
    (u"Kir S, Beddow SA, Samuel VT, et al. FGF19 as a postprandial, insulin-independent activator of hepatic protein "
     u"and glycogen synthesis. Science. 2011;331(6024):1621\u20131624. doi:10.1126/science.1198363."),     # PMID 21436455
    (u"Huang X, Yang C, Luo Y, Jin C, Wang F, McKeehan WL. FGFR4 prevents hyperlipidemia and insulin resistance but "
     u"underlies high-fat diet induced fatty liver. Diabetes. 2007;56(10):2501\u20132510. doi:10.2337/db07-0648."),  # PMID 17664243
    (u"Ge H, Zhang J, Gong Y, et al. Fibroblast growth factor receptor 4 (FGFR4) deficiency improves insulin resistance "
     u"and glucose metabolism under diet-induced obesity conditions. J Biol Chem. 2014;289(44):30470\u201330480. "
     u"doi:10.1074/jbc.M114.592022."),                                                                     # PMID 25204652
    (u"Shemyakin A, K\u00f6vamees O, Rafnsson A, et al. Arginase inhibition improves endothelial function in patients "
     u"with coronary artery disease and type 2 diabetes mellitus. Circulation. 2012;126(25):2943\u20132950. "
     u"doi:10.1161/CIRCULATIONAHA.112.140335."),                                                          # PMID 23183942
    (u"Piatti PM, Monti LD, Valsecchi G, et al. Long-term oral L-arginine administration improves peripheral and "
     u"hepatic insulin sensitivity in type 2 diabetic patients. Diabetes Care. 2001;24(5):875\u2013880. "
     u"doi:10.2337/diacare.24.5.875."),                                                                    # PMID 11347747
    (u"Cervenka I, Agudelo LZ, Ruas JL. Kynurenines: Tryptophan's metabolites in exercise, inflammation, and mental "
     u"health. Science. 2017;357(6349):eaaf9794. doi:10.1126/science.aaf9794."),                           # PMID 28751584
    (u"Br\u00f6er A, Wagner CA, Lang F, Br\u00f6er S. The heterodimeric amino acid transporter 4F2hc/y+LAT2 mediates "
     u"arginine efflux in exchange with glutamine. Biochem J. 2000;349(Pt 3):787\u2013795. doi:10.1042/bj3490787."),  # PMID 10903140
    (u"Tassaneeyakul W, Birkett DJ, McManus ME, et al. Caffeine metabolism by human hepatic cytochromes P450: "
     u"contributions of 1A2, 2E1 and 3A isoforms. Biochem Pharmacol. 1994;47(10):1767\u20131776. "
     u"doi:10.1016/0006-2952(94)90304-2."),                                                                # PMID 8204093
    # FIX round (2026-09-27): two references the corrected prose cites, each taken from its PubMed record (esummary):
    # Kozak & Tate 1982 (reduced Cys-Gly a potent inhibitor of the rat renal microvillus dipeptidase; verify:accuracy) and
    # Ulaganathan 2015 (FGFR4 p.Gly388Arg alters the transmembrane segment and STAT3 signalling; the FGFR4 instrument tag).
    (u"Kozak EM, Tate SS. Glutathione-degrading enzymes of microvillus membranes. J Biol Chem. 1982;257(11):6322\u20136327. "
     u"doi:10.1016/s0021-9258(20)65143-2."),                                                     # PMID 6122685 (DOI: Crossref)
    (u"Ulaganathan VK, Sperl B, Rapp UR, Ullrich A. Germline variant FGFR4 p.G388R exposes a membrane-proximal STAT3 "
     u"binding site. Nature. 2015;528(7583):570\u2013574. doi:10.1038/nature16449."),                         # PMID 26675719
]
# one key per NEW_REFS entry, in the same order (l2l4_facts.refs() cites by key; the four 2026-09-13 entries too)
NEW_REF_KEYS = ["CLSA", "CHENG", "CALDWELL", "ADACHI1990",
                "SPRACKLEN", "MISHRA", "CHEN2024", "KANAI2021", "MIYAZAWA", "ENZAN", "SAKAUE", "SIN", "MCINTYRE",
                "ADACHI1989", "HABIB", "HARTIALA", "WYSS", "MUSUNURU", "NEWTONCHEH", "MALLAT", "GRAINGER", "IL6RGC",
                "KUROSU", "KIR", "HUANG", "GE", "SHEMYAKIN", "PIATTI", "CERVENKA", "BROER", "TASSANEEYAKUL",
                "KOZAK", "ULAGANATHAN"]
# identifiers for every NEW_REFS entry (the refcheck input the build writes; see refcheck_input())
NEW_REF_IDS = {"SPRACKLEN": "PMID 32499647", "MISHRA": "PMID 36180795", "CHEN2024": "PMID 39472450",
               "KANAI2021": "doi:10.1101/2021.09.03.21262975", "MIYAZAWA": "PMID 36653681", "ENZAN": "PMID 41184235",
               "SAKAUE": "PMID 34594039", "SIN": "PMID 26467175", "MCINTYRE": "PMID 6126478", "ADACHI1989": "PMID 2768222",
               "HABIB": "PMID 9560193", "HARTIALA": "PMID 26822151", "WYSS": "PMID 10893433", "MUSUNURU": "PMID 20686566",
               "NEWTONCHEH": "PMID 19430483", "MALLAT": "PMID 11701621", "GRAINGER": "PMID 14699019",
               "IL6RGC": "PMID 22421339", "KUROSU": "PMID 17623664", "KIR": "PMID 21436455", "HUANG": "PMID 17664243",
               "GE": "PMID 25204652", "SHEMYAKIN": "PMID 23183942", "PIATTI": "PMID 11347747",
               "CERVENKA": "PMID 28751584", "BROER": "PMID 10903140", "TASSANEEYAKUL": "PMID 8204093",
               "KOZAK": "PMID 6122685", "ULAGANATHAN": "PMID 26675719"}

# Corrections to individual entries, each with the record that settles it. Nothing
# here is a style preference: every one is a disagreement with the registry, found
# by resolving all 36 entries by title against PubMed. A reference list that is
# contiguous, has nothing dangling and nothing uncited can still point at the wrong
# paper -- those three are ABSENCE checks, and this is the class they cannot see.
ENTRY_FIXES = [
    dict(pmid="28898252",
         why=("the published title reads 'Hemoglobin A1c' (PLoS Med 2017;14(9):e1002383). A "
              "British-spelling pass rewrote a QUOTED title, which is the one place a dialect "
              "sweep must never reach."),
         old=u"Impact of common genetic determinants of haemoglobin A1c on",
         new=u"Impact of common genetic determinants of Hemoglobin A1c on"),
    dict(pmid="39984443",
         why=("PubMed author order is Wang B, Pozarickij A, Mazidi M, Wright N, ... "
              "(Nat Commun 2025;16:1869). The entry promoted the second author to first."),
         old=u"Pozarickij A, Gan W, Lin K, et al.",
         new=u"Wang B, Pozarickij A, Mazidi M, et al."),
    dict(pmid="37220626",
         why=("the Index Medicus form is PubMed's ISOAbbreviation, "
              "'J R Stat Soc Series B Stat Methodol'; 'J R Stat Soc B' abbreviates nothing "
              "an indexer recognises."),
         old=u"J R Stat Soc B. 2020",
         new=u"J R Stat Soc Series B Stat Methodol. 2020"),
    # ---- PACKAGING PASS A (2026-09-30): every record re-fetched from PubMed (efetch) on the day
    dict(pmid="37220626",
         why=("the published title reads 'genetic fine mapping' (no hyphen; PubMed ArticleTitle and Crossref "
              "10.1111/rssb.12388 agree). The main text's own 'fine-mapping' is house spelling and is untouched: "
              "this row reaches the quoted title only."),
         old=u"genetic fine-mapping. J R Stat Soc",
         new=u"genetic fine mapping. J R Stat Soc"),
    dict(pmid="36797296",
         why=("PI-REVIEWABLE: PubMed's Index Medicus initials for the second author are 'EY' (ForeName 'Evan Yi-Wen', "
              "NLM two-initial rule); 'EYW' is the byline form. Aligned to the registry as the list's other entries "
              "are; revert if the PI prefers the byline."),
         old=u"Xu F, Yu EYW, Cai X, et al.",
         new=u"Xu F, Yu EY, Cai X, et al."),
    dict(pmid="29846171",
         why=("PI-REVIEWABLE: the NLM abbreviation (ISOAbbreviation = MedlineTA) is 'Elife'; 'eLife' is the "
              "journal's own styling. Aligned to Index Medicus as the J R Stat Soc row above is; revert if the "
              "PI prefers the journal's styling."),
         old=u"eLife. 2018;7:e34408",
         new=u"Elife. 2018;7:e34408"),
]

# PACKAGING PASS A (2026-09-30): a DOI on every entry, in the house form the 33 DOI-bearing entries already used
# ('... pages. doi:10.xxxx/....' -- a period before 'doi:' and a period after the DOI). The 31 stage-1 entries below
# carried none. Each DOI was re-derived LIVE from the entry's PubMed record (efetch ArticleIdList) after the record was
# asserted to match the PRINTED entry on first author, year, volume and first page (session scratch
# exec/packA1/resolve_dois.py -> dois_resolved.json); the PMID is recorded beside it and is never printed (the journal's
# form is the DOI). The key is a title fragment asserted unique across the whole list at run time.
DOI_ADDS = [
    (u"Cardiovascular-kidney-metabolic health: a", "10.1161/CIR.0000000000001184"),  # PMID 37807924
    (u"Mendelian randomization: genetic anchors", "10.1093/hmg/ddu328"),  # PMID 25064373
    (u"Plasma proteomic associations with genetics", "10.1038/s41586-023-06592-6"),  # PMID 37794186
    (u"The support of human genetic evidence for", "10.1038/ng.3314"),  # PMID 26121088
    (u"Are drug targets with genetic support twice", "10.1371/journal.pgen.1008489"),  # PMID 31830040
    (u"Genetic drug target validation using", "10.1038/s41467-020-16969-0"),  # PMID 32591531
    (u"Strengthening the Reporting of", "10.1001/jama.2021.18236"),  # PMID 34698778
    (u"The MR-Base platform supports systematic", "10.7554/eLife.34408"),  # PMID 29846171
    (u"Bayesian test for colocalisation between", "10.1371/journal.pgen.1004383"),  # PMID 24830394
    (u"A simple new approach to variable selection", "10.1111/rssb.12388"),  # PMID 37220626
    (u"Large-scale cis- and trans-eQTL analyses", "10.1038/s41588-021-00913-z"),  # PMID 34475573
    (u"Large-scale integration of the plasma", "10.1038/s41588-021-00978-w"),  # PMID 34857953
    (u"Plasma proteome analyses in individuals of", "10.1038/s41588-022-01051-w"),  # PMID 35501419
    (u"Comparative studies of 2168 plasma proteins", "10.1038/s41467-025-56935-2"),  # PMID 39984443
    (u"A more accurate method for colocalisation", "10.1371/journal.pgen.1009440"),  # PMID 34587156
    (u"A catalog of genetic loci associated with", "10.1038/s41588-019-0407-x"),  # PMID 31152163
    (u"Discovery and prioritization of variants", "10.1038/s41467-021-24491-0"),  # PMID 34272381
    (u"A fast and efficient colocalization", "10.1038/s41467-020-20885-8"),  # PMID 33536417
    (u"Genetic analysis of over 1 million people", "10.1038/s41588-018-0205-x"),  # PMID 30224653
    (u"Bias due to participant overlap in", "10.1002/gepi.21998"),  # PMID 27625185
    (u"FinnGen provides genetic insights from a", "10.1038/s41586-022-05473-8"),  # PMID 36653562
    (u"Multiancestry genome-wide association study", "10.1038/s41588-018-0058-3"),  # PMID 29531354
    (u"Genetic studies of body mass index yield", "10.1038/nature14177"),  # PMID 25673413
    (u"Discovery and refinement of loci associated", "10.1038/ng.2797"),  # PMID 24097068
    (u"Impact of common genetic determinants of", "10.1371/journal.pmed.1002383"),  # PMID 28898252
    (u"Genome analyses of >200,000 individuals", "10.1016/j.ajhg.2018.09.009"),  # PMID 30388399
    (u"The GTEx Consortium atlas of genetic", "10.1126/science.aaz1776"),  # PMID 32913098
    (u"Tissue-based map of the human proteome.", "10.1126/science.1260419"),  # PMID 25613900
    (u"The ChEMBL Database in 2023: a drug", "10.1093/nar/gkad1004"),  # PMID 37933841
    (u"DGIdb 5.0: rebuilding the drug-gene", "10.1093/nar/gkad1040"),  # PMID 37953380
    (u"The next-generation Open Targets Platform:", "10.1093/nar/gkac1046"),  # PMID 36399499
]
# the house form, asserted on every entry after the adds: '. doi:...' after the pages of an article, '; doi:...' after the
# year of a preprint ('medRxiv. 2021; doi:10.1101/...'), a period after the DOI, nothing else
DOI_FORM = re.compile(r"(?:\. |; )doi:10\.\d{4,9}/\S+\.$")

# ADJUDICATED, deliberately NOT changed. Reference 4 (Ma Q et al., PMID 41387722,
# doi:10.1038/s41467-025-67238-x) prints 2026;17:535. PubMed and Crossref both report
# 12 December 2025 -- the ONLINE date, which is the trap: Crossref also returns
# volume 17, and Nature Communications volume 16 is 2025 (reference 31 in this very
# list is "Nat Commun. 2025;16:1869"), so volume 17 is the 2026 volume. The year is
# kept consistent with the volume, and the DOI resolves it either way.

TOKEN = re.compile(r"\[(\d+(?:\s*[,–‑-]\s*\d+)*)\]")


class Refused(Exception):
    """A precondition does not hold. Never renumber on a guess."""


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def parts_of(tok):
    """The reference numbers inside one bracketed token, in printed order."""
    out = []
    for part in re.split(r"\s*,\s*", tok):
        rng = re.match(r"^(\d+)\s*[–‑-]\s*(\d+)$", part.strip())
        if rng:
            out.extend(range(int(rng.group(1)), int(rng.group(2)) + 1))
        elif part.strip().isdigit():
            out.append(int(part.strip()))
        else:
            raise Refused("unparsable citation part %r" % part)
    return out


def find_list(doc):
    """(heading paragraph, [entry paragraphs]) -- or refuse."""
    paras = doc.paragraphs
    head = [i for i, p in enumerate(paras) if norm(para_text(p)) == "References"]
    if len(head) != 1:
        raise Refused("expected one 'References' heading, found %d" % len(head))
    h = head[0]
    entries, numbers = [], []
    for p in paras[h + 1:]:
        if (p.style.name or "") == "Heading 1":
            break                     # PASS A stage 2: the 'Additional files' section follows the list (build step 4c)
        t = norm(para_text(p))
        if not t:
            continue
        m = re.match(r"^(\d+)\.\s", t)
        if not m:
            raise Refused("paragraph after the reference list is not an entry: %r" % t[:60])
        entries.append(p)
        numbers.append(int(m.group(1)))
    if numbers != list(range(1, len(numbers) + 1)):
        raise Refused("entries are not numbered 1..N: %s" % numbers[:12])
    return paras[h], entries


def body_order(doc, heading):
    """Reference numbers in order of first citation, everything before the list."""
    order, seen = [], set()
    for p in doc.paragraphs:
        if p._p is heading._p:
            break
        for m in TOKEN.finditer(norm(para_text(p))):
            for n in parts_of(m.group(1)):
                if n not in seen:
                    seen.add(n)
                    order.append(n)
    return order


def add_new_refs(doc, entries):
    """Append NEW_REFS as the next numbers, once. Idempotent by title match."""
    have = [norm(para_text(p)) for p in entries]
    added = []
    for text in NEW_REFS:
        key = text.split(".")[1][:40]           # a slice of the title
        if any(key in h for h in have):
            continue
        src = entries[-1]
        new_p = src.insert_paragraph_before(u"%d. %s" % (len(entries) + len(added) + 1, text),
                                            style=src.style)
        src._p.addnext(new_p._p)                # insert_paragraph_before + move = append
        added.append(new_p)
    return entries + added


def apply_entry_fixes(entries):
    """Registry-driven corrections, each asserted unique across the whole list."""
    done = 0
    for fx in ENTRY_FIXES:
        hits = [p for p in entries if fx["old"] in para_text(p)]
        if len(hits) != 1:
            raise Refused("entry fix for PMID %s matches %d entries: %r"
                          % (fx["pmid"], len(hits), fx["old"]))
        p = hits[0]
        t = para_text(p)
        i = t.index(fx["old"])
        splice(p, [(i, i + len(fx["old"]), fx["new"])])
        print("  fixed (PMID %s): %s" % (fx["pmid"], fx["why"][:70]))
        done += 1
    return done


def apply_doi_adds(entries):
    """Append ' doi:<doi>.' to each DOI_ADDS entry (matched by its unique title fragment); then assert that EVERY entry
    prints exactly one DOI in the house form and no PMID string. Idempotent: a fragment whose entry already carries a
    DOI is skipped (the fixed docx is rebuilt from the stage-1 file every time, so this is the normal steady state
    only under --check)."""
    done = 0
    for frag, doi in DOI_ADDS:
        hits = [p for p in entries if frag in para_text(p)]
        if len(hits) != 1:
            raise Refused("DOI add %r matches %d entries" % (frag, len(hits)))
        t = para_text(hits[0])
        if "doi:" in t:
            continue
        if not t.rstrip().endswith("."):
            raise Refused("entry does not end with a period before the DOI is appended: %r" % t[-40:])
        splice(hits[0], [(len(t.rstrip()), len(t), u" doi:%s." % doi)])
        done += 1
    bad = [para_text(p)[:50] for p in entries
           if para_text(p).count("doi:") != 1 or not DOI_FORM.search(para_text(p).strip()) or "PMID" in para_text(p)]
    if bad:
        raise Refused("%d entries without exactly one house-form DOI (or printing a PMID): %s" % (len(bad), bad[:3]))
    return done


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        raise SystemExit(__doc__)
    path = args[0]
    check = "--check" in sys.argv

    doc = docx.Document(path)
    heading, entries = find_list(doc)
    print("reference list: %d entries" % len(entries))

    if not check:
        entries = add_new_refs(doc, entries)
        print("after adding the deferred references: %d entries" % len(entries))

    if not check:
        print("entry fixes: %d" % apply_entry_fixes(entries))
        print("DOIs added: %d (every entry now prints exactly one DOI in the house form)" % apply_doi_adds(entries))

    n = len(entries)
    order = body_order(doc, heading)
    print("first-citation order: %s" % order)

    missing = sorted(set(range(1, n + 1)) - set(order))
    extra = sorted(x for x in order if x < 1 or x > n)
    if missing or extra:
        raise Refused("not a bijection: uncited %s; cited with no entry %s" % (missing, extra))
    if len(order) != n:
        raise Refused("%d distinct numbers cited against %d entries" % (len(order), n))

    old2new = {old: i + 1 for i, old in enumerate(order)}
    moved = sum(1 for k, v in old2new.items() if k != v)
    print("%d of %d entries change number" % (moved, n))
    if check:
        print("ascending already" if moved == 0 else "NOT in first-citation order")
        return 0 if moved == 0 else 1

    # ---- in-text tokens. Rewrite right-to-left so earlier offsets stay valid.
    rewritten = 0
    for p in doc.paragraphs:
        if p._p is heading._p:
            break
        t = para_text(p)
        edits = []
        for m in TOKEN.finditer(t):
            # sorted: a bracket the L2-L4 update writes in stage-1 order prints ascending after the renumber (P18b)
            nums = sorted(old2new[x] for x in parts_of(m.group(1)))
            edits.append((m.start(), m.end(), "[" + ", ".join(str(x) for x in nums) + "]"))
        if edits:
            splice(p, edits)
            rewritten += len(edits)
    print("%d citation tokens rewritten" % rewritten)

    # ---- the entries themselves: new prefix, then re-seated in the new order
    new_order = [None] * n
    for p, old in zip(entries, range(1, n + 1)):
        t = para_text(p)
        m = re.match(r"^\s*(\d+)\.", t)
        if not m or int(m.group(1)) != old:
            raise Refused("entry %d does not start with its own number: %r" % (old, t[:40]))
        splice(p, [(m.start(1), m.end(1), str(old2new[old]))])
        new_order[old2new[old] - 1] = p
    if any(x is None for x in new_order):
        raise Refused("the new order has a hole")

    anchor = heading._p
    for p in new_order:
        anchor.addnext(p._p)
        anchor = p._p
    print("entries re-seated in first-citation order")

    doc.save(path)

    # ---- read it back from disk: the artifact, not the edit
    doc = docx.Document(path)
    heading, entries = find_list(doc)
    got = body_order(doc, heading)
    if got != list(range(1, len(entries) + 1)):
        raise Refused("POST-CHECK FAILED: first-citation order is %s" % got[:12])
    print("post-check: %d entries, first citation ascends 1..%d" % (len(entries), len(entries)))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Refused as e:
        raise SystemExit("REFUSED: %s" % e)
