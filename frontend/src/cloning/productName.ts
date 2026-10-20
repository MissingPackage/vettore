// The name offered for a new product tab. The engine names a combination after its inputs
// ("pUC19 + GFP"), the same whatever the method; the method is added so two products made from
// the same files are told apart, and a number when that name is already open.

/** Short method names, by the engine's operation (Product.operation). */
const METHOD_TAGS: Record<string, string> = {
  "Restriction cloning": "restriction cloning",
  "Gibson assembly": "Gibson",
  "NEBuilder HiFi assembly": "HiFi",
  "In-Fusion cloning": "In-Fusion",
  "Golden Gate assembly": "Golden Gate",
  "Gateway BP reaction": "BP",
  "Gateway LR reaction": "LR",
  "Insert multiple fragments": "ligation",
  "Ligate fragments": "ligation",
  "Anneal oligos": "annealed oligos",
};

/**
 * "pUC19 + GFP" → "pUC19 + GFP (Gibson)", then "pUC19 + GFP (Gibson) 2" if that tab is open.
 * Names that already say what was done ("pUC19 101..900 PCR product", "pUC19 mut") are kept.
 */
export function defaultProductName(name: string, operation: string, openNames: string[]): string {
  const tagged = name.includes(" + ") ? `${name} (${METHOD_TAGS[operation] ?? operation})` : name;
  const taken = new Set(openNames);
  if (!taken.has(tagged)) return tagged;
  for (let i = 2; ; i++) {
    if (!taken.has(`${tagged} ${i}`)) return `${tagged} ${i}`;
  }
}
