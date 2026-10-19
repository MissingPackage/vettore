// Export the visible map (circular or linear) as a standalone SVG or PNG.

const STYLE_PROPS = [
  "fill", "fill-opacity", "stroke", "stroke-width", "stroke-opacity", "stroke-dasharray",
  "opacity", "font-family", "font-size", "font-weight", "font-style", "text-anchor",
  "dominant-baseline", "visibility", "display",
];

function findMapSvg(): SVGSVGElement | null {
  const selectors = [".veCircularView svg.circularViewSvg", ".veCircularView svg", ".veLinearView svg"];
  for (const sel of selectors) {
    const candidates = Array.from(document.querySelectorAll<SVGSVGElement>(sel));
    const visible = candidates
      .filter((el) => el.getBoundingClientRect().width > 0)
      .sort((a, b) => b.getBoundingClientRect().width * b.getBoundingClientRect().height -
        a.getBoundingClientRect().width * a.getBoundingClientRect().height);
    if (visible[0]) return visible[0];
  }
  return null;
}

function inlineStyles(source: Element, target: Element) {
  const computed = window.getComputedStyle(source);
  const style = STYLE_PROPS.map((p) => `${p}:${computed.getPropertyValue(p)}`).join(";");
  target.setAttribute("style", style);
  const s = source.children;
  const t = target.children;
  for (let i = 0; i < s.length; i++) inlineStyles(s[i], t[i]);
}

export function mapSvgString(): string {
  const svg = findMapSvg();
  if (!svg) throw new Error("Open the Map view to export the map.");
  const rect = svg.getBoundingClientRect();
  const clone = svg.cloneNode(true) as SVGSVGElement;
  inlineStyles(svg, clone);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("width", String(Math.round(rect.width)));
  clone.setAttribute("height", String(Math.round(rect.height)));
  if (!clone.getAttribute("viewBox")) clone.setAttribute("viewBox", `0 0 ${rect.width} ${rect.height}`);
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("width", "100%");
  bg.setAttribute("height", "100%");
  bg.setAttribute("fill", "white");
  clone.insertBefore(bg, clone.firstChild);
  return new XMLSerializer().serializeToString(clone);
}

export async function exportMap(format: "svg" | "png", scale = 3): Promise<Blob> {
  const svgText = mapSvgString();
  if (format === "svg") return new Blob([svgText], { type: "image/svg+xml" });
  const img = new Image();
  const url = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml" }));
  try {
    await new Promise<void>((resolve, reject) => {
      img.onload = () => resolve();
      img.onerror = () => reject(new Error("Could not render the map image."));
      img.src = url;
    });
    const canvas = document.createElement("canvas");
    canvas.width = img.width * scale;
    canvas.height = img.height * scale;
    const ctx = canvas.getContext("2d")!;
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0);
    return await new Promise<Blob>((resolve, reject) =>
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("PNG encoding failed."))), "image/png"),
    );
  } finally {
    URL.revokeObjectURL(url);
  }
}
