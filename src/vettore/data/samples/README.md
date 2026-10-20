# Sample constructs

The files in this directory are the sample constructs offered on Vettore's welcome screen
(`GET /api/app/samples`). Each one is a public GenBank record downloaded unchanged from the
NCBI Nucleotide database with E-utilities, for example:

```sh
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\
&id=M77789.2&rettype=gbwithparts&retmode=text&tool=vettore" -o M77789.2.gb
```

Nothing in the records is edited: the sequence, the FEATURES table and the references are the
submitters' own. The display name shown in Vettore (the `label` below) comes from the record's
DEFINITION line, because the GenBank LOCUS name of an older record is an internal identifier
(`SYNPUC19V`) rather than the name of the plasmid.

| File | Accession.version | Length | Shown as | Why it is here | Retrieved |
|---|---|---|---|---|---|
| `M77789.2.gb` | M77789.2 | 2686 bp, circular | pUC19 | The classic high-copy cloning vector: polylinker, lacZα, AmpR, pUC origin | 2026-10-20 |
| `PZ629150.1.gb` | PZ629150.1 | 7079 bp, circular | pAAV-cTNT-SaCas9-U6 | AAV transfer plasmid whose two ITRs are annotated as `repeat_region` with `/rpt_type=inverted` | 2026-10-20 |
| `MT612431.1.gb` | MT612431.1 | 7204 bp, circular | pAAV-EF1A-DIO-Gephyrin.FingR-GFP | Cre-dependent (DIO) AAV transfer plasmid; its ITRs are annotated only by a `/note` ("5' ITR", "3' ITR") | 2026-10-20 |

`PZ629150.1.gb` and `MT612431.1.gb` are the same records already used as test fixtures in
`tests/data/`; the copies here were checked base for base against the live NCBI records on
2026-10-20 (both identical). They are kept as separate files because `tests/` is not part of
the installed package and the welcome screen must work in a packaged build.

## Licence and reuse

GenBank records are produced by their submitters and distributed by the NCBI, which states
that the information in its databases is in the public domain and may be redistributed; the
NCBI asks that the source be cited and makes no warranty about the data
(<https://www.ncbi.nlm.nih.gov/home/about/policies/>). Each record's own citation (authors,
title, journal) is kept inside the file, so a user who opens a sample sees where it came from.

Vettore does not call NCBI at runtime: these files are read from the installed package.
