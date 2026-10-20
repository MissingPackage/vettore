# Test fixtures

The `*.dna` files and `3730.ab1` (ABI Sanger trace) are copied unchanged from the Biopython
test suite (`Tests/`, https://github.com/biopython/biopython),
distributed under the Biopython License Agreement / BSD 3-Clause License.

## File-exchange fixtures (`exchange-*`)

Downloaded on 2026-10-20. INSDC records (GenBank/ENA) are freely available public data.

| File | Source | Changes |
|---|---|---|
| `exchange-pBR322.embl` | ENA record J01749 (pBR322), EMBL flat file from `https://www.ebi.ac.uk/ena/browser/api/embl/J01749` | none |
| `exchange-plasmids.fasta` | NCBI nuccore L09137.2, J01749.1, X06403.1 (pUC19c, pBR322, pACYC184), E-utilities `efetch rettype=fasta` | none |
| `exchange-sample-a.xdna`, `exchange-sample-b.xdna` | Biopython `Tests/Xdna/sample-a.xdna`, `sample-b.xdna` (BSD) | none |
| `exchange-artificial.gck` | Biopython `Tests/Gck/artificial.gck` (BSD) | none |
| `exchange-pUC19-numbered.seq` | NCBI L09137.2 GenBank record | only the numbered ORIGIN block (plain-sequence text with position numbers) |
| `exchange-pUC19-raw.txt` | NCBI L09137.2 FASTA | header line removed (FASTA without header) |
| `exchange-quirks-header.gb` | NCBI J01749.1 GenBank record | hand-edited: LOCUS line replaced by a long name without topology or date; one `/note="non-ASCII test: 37 °C, 5 µl"` added to the tet CDS; final `//` removed; Windows (CRLF) line endings; encoded as cp1252 |
| `exchange-quirks-features.gb` | NCBI J01749.1 GenBank record | hand-edited: eight `misc_feature` entries labelled `test …` added before the last feature: a region across the origin written as `join(4300..4361,1..50)` and as `4300..50`, a site between two bases (`100^101`), regions beyond (`4400..4500`) and past (`4300..4400`) the end, a location in another record (`J01636.1:1..10`), and two features with the same name and unusual qualifiers. The `test …` features mark positions only and carry no biological meaning |
