// Files dropped on the window. In the desktop window the engine learns each file's full path
// (pywebview resolves it natively and src/vettore/app.py sends it as DROP_EVENT), so the
// document remembers where it came from and Save writes back in place. In a browser, or for a
// file whose path could not be resolved, the content is uploaded instead.

export const DROP_EVENT = "vettore:dropped-paths";
const WAIT_MS = 1500;

declare global {
  interface Window {
    /** Set by the desktop app once it listens for dropped files. */
    __vettoreNativeDrop?: boolean;
  }
}

export interface DroppedFiles {
  /** Full paths to open from disk. */
  paths: string[];
  /** Files to upload by content. */
  files: File[];
}

export function nativeDropActive(): boolean {
  return !!window.__vettoreNativeDrop && !!window.pywebview;
}

/** Wait for the desktop app to report the dropped files' paths; upload whatever it cannot. */
export function resolveDrop(files: File[]): Promise<DroppedFiles> {
  if (!nativeDropActive()) return Promise.resolve({ paths: [], files });
  return new Promise((resolve) => {
    const finish = (paths: string[], unresolved: string[] | null) => {
      window.removeEventListener(DROP_EVENT, onPaths);
      clearTimeout(timer);
      const left = unresolved === null ? files : files.filter((f) => unresolved.includes(f.name));
      resolve({ paths, files: left });
    };
    const onPaths = (e: Event) => {
      const detail = (e as CustomEvent<{ paths?: string[]; unresolved?: string[] }>).detail ?? {};
      finish(detail.paths ?? [], detail.unresolved ?? []);
    };
    // No answer from the desktop app: upload everything rather than lose the drop.
    const timer = setTimeout(() => finish([], null), WAIT_MS);
    window.addEventListener(DROP_EVENT, onPaths);
  });
}
