// The wrong clones of a cloning step: what a colony screen must tell apart from the product.
// Kept on the product's first History entry, so Find Diagnostic Digest starts from them.
import type { CloningResult, Product, WrongClone } from "./types";

/**
 * The simulation's other products (the other orientation, a by-product), then the wrong
 * clones the engine listed for this method (insert reversed, empty vector, parent vector).
 * A molecule is listed once, under its first name.
 */
export function wrongClonesOf(product: Product, result: CloningResult): WrongClone[] {
  const key = (p: Product) => `${p.doc.circular ? "c" : "l"}:${p.doc.sequence.toUpperCase()}`;
  const seen = new Set([key(product)]);
  const out: WrongClone[] = [];
  const others = result.products.filter((p) => key(p) !== key(product));
  for (const p of [...others, ...(result.screen_against ?? [])]) {
    if (seen.has(key(p))) continue;
    seen.add(key(p));
    out.push({ name: p.label ?? p.doc.name, data: p.doc });
  }
  return out;
}
