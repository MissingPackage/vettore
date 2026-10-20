// Why Save makes a GenBank copy instead of overwriting the file a document came from.
import type { DocSource } from "../documents";

function withArticle(label: string): string {
  return /^[aeiou]/i.test(label) ? `an ${label}` : `a ${label}`;
}

/** One sentence saying why the original file is not overwritten. */
export function copyReason(name: string, source: DocSource): string {
  const file = `“${source.fileName}”`;
  if (source.records > 1) {
    const others = source.records - 1;
    return `“${name}” is one of ${source.records} sequences in ${file}; saving over that file `
      + `would replace the other ${others === 1 ? "one" : others}.`;
  }
  if (source.incomplete) {
    return `Some features of ${file} could not be shown, so saving over it would lose them.`;
  }
  if (source.label === ".dna") {
    return `${file} is a .dna file, which can also hold history, notes and sequencing traces `
      + "that Vettore does not read; saving over it would lose them.";
  }
  if (source.label === "GenBank" || source.label === "FASTA") {
    return `${file} holds ${withArticle(source.label)} record under another file extension.`;
  }
  if (source.label === "plain sequence text") {
    return `${file} is plain sequence text, which Vettore reads but does not write.`;
  }
  return `${file} is ${withArticle(source.label)} file, which Vettore reads but does not write.`;
}
