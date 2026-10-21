# Vettore design system

How Vettore looks and why: the tokens, the shared classes, and the rules every screen and
dialog follows. Vettore's users are bench scientists, many of them not technical. A screen that
is confusing, or a failure without an explanation, makes them give up, so most of this file is
about being clear rather than decorative.

Files:

| File | What it holds |
|---|---|
| `frontend/src/theme/tokens.css` | Design tokens (CSS variables): colour, type, spacing, radius, elevation, focus, motion |
| `frontend/src/theme/overrides.css` | Blueprint 3 and editor-component chrome brought in line with the tokens |
| `frontend/src/theme/dark-chrome.css` | The optional dark frame (Preferences → Appearance) |
| `frontend/src/app.css` | Vettore's own surfaces (shell, views, dialogs) and the shared classes below |
| `frontend/e2e/design-screens.spec.ts` | Screenshot harness and layout, contrast and focus checks |
| `frontend/public/favicon.svg`, `vettore.png`, `vettore.ico` | App icon |

## Principles

1. **The familiar layout comes first.** Menu bar, document tabs, view tabs, map and sequence
   panes, status bar: their places and order do not change. Styling refines; it never moves
   things around.
2. **Calm.** Neutral greys, one accent blue, and colour only where it carries meaning: red for
   something that failed, amber for something to check, green for something confirmed.
3. **Readable by everyone.** Every text meant to be read reaches WCAG AA contrast (4.5:1) on the
   surface it sits on; control outlines and focus rings reach 3:1. Nothing smaller than 12px.
4. **Say what to do next.** A dialog opens with one sentence on what it does, its defaults give
   a useful result at the first click, problems are reported next to the field they concern,
   and every error ends with the action that fixes it.
5. **One application.** Vettore's own surfaces and the editor component's dialogs, menus and
   tables share the same font, accent, inputs and messages. The map, sequence and alignment
   drawings themselves are the editor's and are not restyled.

## Tokens

Always use a token, never a literal colour. Ratios are WCAG contrast against the surfaces named.

### Colour

| Token | Value | Use |
|---|---|---|
| `--vt-text` | `#374151` | Body text in Vettore's surfaces (10.3:1 on white, 8.0:1 on `surface-sunken`) |
| `--vt-text-strong` | `#111827` | Headings, labels, active tabs, values |
| `--vt-text-muted` | `#5b6472` | Secondary text: hints, counts, captions (4.6:1 or more on every surface) |
| `--vt-text-placeholder` | `#6b7280` | Input placeholders (4.8:1 on white) |
| `--vt-text-faint` | `#9ca3af` | Disabled controls and decoration only, never text to read (2.5:1) |
| `--vt-text-on-accent` | `#ffffff` | Text on accent or danger fills |
| `--vt-text-on-accent-muted` | `rgba(255,255,255,0.9)` | Secondary text on accent fills (about 5:1 on `--vt-accent`) |
| `--vt-surface` | `#ffffff` | Pages, inputs, cards, the active document tab |
| `--vt-surface-subtle` | `#f7f8fa` | Status bar, dialog bodies, sample cards |
| `--vt-surface-sunken` | `#dfe3e8` | Document tab strip |
| `--vt-surface-raised` | `#eceff2` | View-tab track, hovered document tabs |
| `--vt-surface-hover` | `#e9edf2` | Hover on tabs and plain buttons |
| `--vt-overlay-hover` | 8 % ink | Hover tint over any background (icon buttons) |
| `--vt-border` | `#c3c8cf` | Container outlines and dividers (decorative) |
| `--vt-border-subtle` | `#d8dce1` | Lighter dividers inside a surface |
| `--vt-border-strong` | `#868e9a` | Input and select outlines (3.3:1 on white) |
| `--vt-accent` | `#106ba3` | Links, filled primary buttons, active-tab indicators (5.8:1 both ways with white) |
| `--vt-accent-hover` / `--vt-accent-active` | `#0e5a8a` / `#0c4d77` | Pressed states of the accent |
| `--vt-accent-bg` / `--vt-accent-border` | `#eef5fb` / `#b3d4ee` | Info notices, primary tags |
| `--vt-accent-disabled` | accent at 45 % | Disabled primary buttons |
| `--vt-selection` | `#e6f0fb` | Selected table rows |
| `--vt-focus` | `#2d72d2` | Keyboard focus ring (4.7:1 on white) |
| `--vt-focus-halo` | focus at 30 % | Soft ring around a focused text input |
| `--vt-danger` / `-bg` / `-border` | `#b42318` / `#fef3f2` / `#fecdca` | Errors (6.0:1 text on its background) |
| `--vt-danger-halo` | danger at 25 % | Focus ring of an invalid input |
| `--vt-warning` / `-bg` / `-border` | `#b54708` / `#fffaeb` / `#fedf89` | Things to check (5.2:1) |
| `--vt-success` / `-bg` / `-border` | `#027a48` / `#ecfdf3` / `#a6f4c5` | Confirmed results (5.1:1) |

The accent follows the editor component's blue so that its toolbar and Vettore's buttons read
as one application.

### Type

`--vt-font` is `system-ui` first (Segoe UI on Windows, the desktop font on Linux), then common
Linux fonts, so text never falls through to Blueprint's icon font. `--vt-font-mono` starts with
Menlo, which the editor component bundles as its sequence font: sequences typed in Vettore's
fields look like the sequences on the map.

| Token | Size | Use |
|---|---|---|
| `--vt-font-size-xs` | 12px | Hints, captions, badges; the smallest text allowed |
| `--vt-font-size-sm` | 13px | Dense tables, enzyme chips |
| `--vt-font-size` | 14px | Body text and controls (Blueprint's base size) |
| `--vt-font-size-md` | 16px | View headings, empty-state titles |
| `--vt-font-size-lg` | 18px | Dialog titles |
| `--vt-font-size-xl` | 28px | Welcome title |

Line height: `--vt-line-height` (1.43, 20px at 14px) for paragraphs; `--vt-line-height-tight`
(1.25) for labels, tabs and headings. Weights: `--vt-weight-regular` 400, `-medium` 500,
`-strong` 600.

### Space, shape, size, elevation, motion

| Token | Value | Use |
|---|---|---|
| `--vt-space-0` … `--vt-space-7` | 2, 4, 8, 12, 16, 24, 32, 48 px | Gaps and padding; use the scale, not new values |
| `--vt-radius-sm` / `--vt-radius` / `--vt-radius-lg` / `--vt-radius-xl` | 4 / 6 / 8 / 10 px | Controls, toolbar keys / notices, segmented controls / cards, tabs / dialogs |
| `--vt-control-height` | 30px | Every input, select and button in Vettore's forms |
| `--vt-control-height-sm` | 24px | The smallest clickable target (WCAG 2.5.8) |
| `--vt-measure` | 72ch | Longest line for explanatory text |
| `--vt-shadow-1` / `-2` / `-3` | | Raised controls / popovers / dialogs |
| `--vt-duration` | 120ms | Colour transitions on hover |
| `--vt-duration-med` | 180ms | Chrome that moves: the view-tab pill, a new tab, a notice |
| `--vt-duration-data` | 260ms | Data that changes shape: gel bands moving to a new size |
| `--vt-ease-out` | `cubic-bezier(0.2, 0.7, 0.2, 1)` | Easing for everything that moves |

All three durations are 0 when the system asks for reduced motion.

### Motion

Motion is for the frame and for feedback, never for reading. The budget:

- **Chrome, at most 180 ms**: the view-tab pill slides, a new document's tab slides in and its
  outline glows once (so a new product is easy to find), notices slide in from the right,
  welcome actions lift 1 px on hover. Nothing waits for an animation before it can be clicked.
- **Data, at most 260 ms, and only where the movement is the information**: gel bands slide to
  their new size when a lane's enzymes change, so the eye sees which fragment split.
- **Never**: the map, sequence, alignment and tables. Zoom, scroll, selection and view changes
  stay instant. No entrance animations on pages, no parallax, no animated text.

### Dark frame

File → Preferences → Appearance offers Light (the default), Dark menus and tabs, and As the
system is set. `theme/dark-chrome.css` styles the menu bar, toolbar, document tabs, view bar
and status bar under `html[data-chrome="dark"]`, with its own `--vt-chrome-*` tokens (every
text pair at least 6:1, control outlines 3.3:1). The map, sequence, tables, dialogs and menus
stay light: the editor's feature and enzyme colours are drawn for white. A new chrome surface
needs a rule there too.

## Shared classes

Defined in `app.css`. Use them instead of inline styles; they already carry the right colour,
size and spacing. Existing `vt-*` classes keep working.

### Dialog structure

**`.vt-dialog-intro`**: the first line of every dialog, saying what it does in one sentence.

```tsx
<p className="vt-dialog-intro">Cut the vector and the insert with the same enzymes and ligate them.</p>
```

**`.vt-section-title`**: a heading inside a dialog or view ("Primers to order", "Lanes").

**`.vt-field`, `.vt-field-label`, `.vt-field-hint`**: one labelled control. The `Field`
component in `cloning/inputs.tsx` renders this structure; use it in new forms.

**`.vt-field-error`** and **`.vt-field-invalid`**: inline validation directly under the field.
Put `vt-field-invalid` on the `.vt-field` to outline its control in red.

```tsx
<label className={`vt-field${error ? " vt-field-invalid" : ""}`}>
  <span className="vt-field-label">Overlap length (bp)</span>
  <input className="bp3-input" aria-invalid={!!error} aria-describedby="overlap-error" … />
  {error && <span id="overlap-error" className="vt-field-error">Use 15 to 40 bp.</span>}
</label>
```

**`.vt-card`** (same as `.vt-insert-row`): a white group inside a dialog, for one lane, insert
or part.

### Messages

**`.vt-notice`** with **`.vt-notice-danger`**, **`-warning`**, **`-success`** or **`-info`**,
and an optional **`.vt-notice-title`**: a message block in Vettore's own views. In Blueprint
dialogs a `Callout` with the same intent looks the same (see `overrides.css`); either is fine.

```tsx
<div className="vt-notice vt-notice-warning" role="status">
  <span className="vt-notice-title">Check before ordering</span>
  BsaI also cuts inside the insert. Remove the site or choose another enzyme.
</div>
```

**`.vt-error`**: an error paragraph inside a dialog (a danger notice that wraps). In the
document tab strip the same class becomes a one-line pill that is dismissed by a click.

**`.vt-empty-state`** with **`.vt-empty-state-title`**: nothing to show yet. Say why and what to
do.

```tsx
<div className="vt-empty-state">
  <p className="vt-empty-state-title">No product yet</p>
  <p>Choose a vector and at least one enzyme on the left.</p>
</div>
```

**`.vt-busy`**: work in progress, with a spinner drawn before the text (it stands still when
the system asks for reduced motion). Pair it with `role="status"` so screen readers announce it.

```tsx
{busy && <span className="vt-busy" role="status">Simulating…</span>}
```

### Text helpers

| Class | Use |
|---|---|
| `.vt-muted` | Secondary text; also works on `<td>` and `<tr>` |
| `.vt-mono` | Sequences, recognition sites, overhangs |
| `.vt-num` | Numbers in columns: tabular figures, never wraps ("59.7 °C", "3,941 bp") |
| `.vt-ok`, `.vt-warn`, `.vt-bad` | A value that is confirmed, needs checking, or failed (colour plus weight; also on `<td>`) |
| `.vt-kbd` | A key: `<span className="vt-kbd">Ctrl</span>+<span className="vt-kbd">O</span>` |
| `.vt-break` | Long sequences that may break anywhere |
| `.vt-sr-only` | Text for screen readers only, e.g. the name of an icon-only button |
| `.vt-row` | A horizontal group of controls with a 10px gap |

## Blueprint and the editor component

- Use Blueprint components (`Dialog`, `Button`, `Callout`, `HTMLSelect`, `InputGroup`) for
  controls. Their look is adjusted in `overrides.css`: dialog surfaces, input outlines at 3:1,
  filled primary buttons in the accent, tag and callout colours from the tokens, and one
  item per menu row.
- Every rule in `overrides.css` starts with `:root`. Blueprint's stylesheet is bundled twice,
  the second copy after `app.css`, so an override with Blueprint's own specificity would lose.
- Do not restyle the editor's map, sequence rows, alignment tracks or chromatograms.
- Only the page font changes globally; size, line height and colour of the editor's panels stay
  Blueprint's, which those panels are laid out for.

Rules that adapt the editor component and Blueprint to Vettore (in `theme/overrides.css` and
`app.css`):

- Every `.bp3-html-select` draws its own chevron (`::after`), so plain `<select>` wrappers need
  no caret icon; do not add Blueprint's icon on top of it.
- `.vt-shell[data-overlay-view]` hides the editor panels covered by the Enzymes and History
  views, so their sliders cannot catch clicks meant for Vettore's tables.
- `html[data-platform]` selects how menu shortcuts are written: `Ctrl+Shift+S` off a Mac, the
  key symbols on a Mac.
- The editor's toolbar keys have no fill; a tool that is on is a white key with an accent
  outline (`.veToolbar … .bp3-active`).
- The menu search (top of Help) ranks with `menuSearch.tsx`: a name that starts with the typed
  text, then a word that does, then the text anywhere, then the menu it sits in. Ctrl+/ and the
  view bar's "Search commands" button open it with the box focused (the editor's own shortcut
  was Cmd+/ only).
- Messages from the engine appear in the message stack (`files/NoticeStack`, classes
  `vt-toast*`); the `vt-notice*` classes are for messages inside dialogs.

## Dialogs: the pattern

1. Title: what the dialog does, in the menu's words ("Simulate Agarose Gel").
2. `.vt-dialog-intro`: one sentence on what happens and what is needed.
3. Inputs on the left, live result on the right (`.vt-cloning` grid). The body scrolls inside
   the dialog, so the footer buttons stay on screen at 1366×768.
4. Defaults that already give a useful result: the active document preselected, the usual
   overlap length, the usual ladder.
5. Validation next to the field (`.vt-field-error`), not only in the result pane.
6. A busy state while the engine works (`.vt-busy`), and the primary button disabled until
   there is a result.
7. Errors that say what happened and what to do next (see Writing).
8. Footer: `Cancel` on the left of the primary action; the primary action is named for its
   result ("Create Product", "Add 9 features"), not "OK".
9. `data-testid` on the dialog body and on every control an end-to-end test touches.

## Writing

- Plain words a bench scientist uses: "cut", "insert", "primer", "band"; no internal names.
- Errors: what happened, why if known, and what to do. "Forward primer “F1” does not bind
  “pVector”. Check that it matches the template, or choose another template." rather than
  "No product".
- Numbers with units and a space: `59.7 °C`, `3,941 bp`, `1.0 %`. Coordinates are 1-based, as
  biologists count.
- Menu items and dialog titles in Title Case with `…` when they open a dialog; labels, hints
  and messages in sentence case.
- Never name a commercial sequence-editing product; describe the feature in field terms.

## Do and don't

| Do | Don't |
|---|---|
| Use tokens and shared classes | Write literal colours or one-off font sizes |
| Keep one control height (30px) within a form | Mix Blueprint's small 24px inputs with 30px selects |
| Put a hint under a field when the format is not obvious | Hide required information in a placeholder |
| Show why something is disabled (hint or tooltip) | Leave a grey button without explanation |
| Use colour together with text or weight (`.vt-ok` plus "100.00 %") | Rely on colour alone |
| Let long sequences break (`.vt-break`) and keep numbers whole (`.vt-num`) | Let a Tm wrap onto two lines |
| Keep the layout: menu, document tabs, view tabs, panes, status bar | Move or restyle the editor's drawings |

## Accessibility checklist

For every screen or dialog you add or change:

- [ ] Text contrast at least 4.5:1, control outlines and focus rings at least 3:1 (the harness
      checks the common cases with `expectReadable`).
- [ ] Every control reachable with Tab, in reading order, with a visible focus ring.
- [ ] Every input has a visible label (`.vt-field-label` or Blueprint `FormGroup`); icon-only
      buttons have `aria-label` or `.vt-sr-only` text.
- [ ] Targets at least 24×24 px (`--vt-control-height-sm`).
- [ ] Errors are next to the field, linked with `aria-describedby`, and say how to fix them.
- [ ] Busy and result messages use `role="status"`; blocking errors use `role="alert"`.
- [ ] Colour is never the only signal.
- [ ] Nothing overlaps or truncates at 1366×768 and 1920×1080; the dialog footer stays visible.
- [ ] Motion respects `prefers-reduced-motion` (the tokens set every `--vt-duration*` to 0)
      and stays within the budget under Motion.
- [ ] A new chrome surface is readable in the dark frame too (`theme/dark-chrome.css`).

## Checking the design

```sh
cd frontend && npm run build
npx playwright test e2e/design-screens.spec.ts
```

The harness opens the test fixtures and every dialog at 1366×768 and 1920×1080 and saves
screenshots to `frontend/test-results/design/<size>/` (git-ignored; Playwright clears the folder
at the start of every run). It also fails when a dialog or its footer leaves the window, the
view or status bar overflows, secondary text drops below AA, a form mixes control heights, two
menu items share a row, keyboard focus is invisible, or form controls lose the page font.

The app icon's PNG and ICO files are rasterised from `favicon.svg` with Chromium; the ICO holds
PNG-compressed images from 16 to 256 px.
