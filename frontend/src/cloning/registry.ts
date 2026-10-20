// Cloning actions offered in the Actions menu, with the engine endpoint each one calls.
import type { ComponentType } from "react";
import {
  DeleteFragmentForm,
  DesignPrimersForm,
  AnnealOligosForm,
  GoldenGateForm,
  InsertMultipleForm,
  LigateFragmentsForm,
  MutagenesisForm,
  PcrForm,
  RestrictionInsertForm,
  gatewayForm,
  overlapForm,
  taTopoForm,
} from "./forms";
import type { FormProps } from "./types";

export interface CloningAction {
  id: string;
  title: string;
  endpoint: string;
  Form: ComponentType<FormProps>;
  /** One line: what the action does. */
  help: string;
  /** Optional second paragraph: assumptions and how the inputs are read. */
  details?: string;
  /** Engines that are not implemented yet stay out of the menus. */
  ready?: boolean;
}

export const CLONING_ACTIONS: CloningAction[] = [
  {
    id: "restriction",
    ready: true,
    title: "Restriction Cloning: Insert Fragment",
    endpoint: "restriction-insert",
    Form: RestrictionInsertForm,
    help: "Cut the vector and the insert with restriction enzymes and ligate the insert fragment into the vector.",
  },
  {
    id: "delete-fragment",
    ready: true,
    title: "Delete Restriction Fragment",
    endpoint: "delete-fragment",
    Form: DeleteFragmentForm,
    help: "Remove the fragment between two sites and re-ligate compatible ends.",
  },
  {
    id: "pcr",
    ready: true,
    title: "PCR",
    endpoint: "pcr",
    Form: PcrForm,
    help: "Amplify a template with two primers (5′ tails are added to the product).",
  },
  {
    id: "design-primers",
    ready: true,
    title: "Design Primers to Amplify a Region",
    endpoint: "design-pcr-primers",
    Form: DesignPrimersForm,
    help: "Choose primers whose annealing parts reach the target Tm; optional 5′ tails.",
  },
  {
    id: "gibson",
    ready: true,
    title: "Gibson Assembly",
    endpoint: "overlap-assembly",
    Form: overlapForm("gibson"),
    help: "Insert fragments into a linearized vector; primers carry the overlaps.",
  },
  {
    id: "hifi",
    ready: true,
    title: "NEBuilder HiFi DNA Assembly",
    endpoint: "overlap-assembly",
    Form: overlapForm("hifi"),
    help: "Same design as Gibson; NEB recommends 15–30 bp overlaps.",
  },
  {
    id: "infusion",
    ready: true,
    title: "In-Fusion Cloning",
    endpoint: "overlap-assembly",
    Form: overlapForm("infusion"),
    help: "15 bp homology to the vector ends on the insert primers.",
  },
  {
    id: "mutagenesis",
    ready: true,
    title: "Mutagenesis (Q5 SDM style)",
    endpoint: "mutagenesis",
    Form: MutagenesisForm,
    help: "Substitute, insert or delete bases with back-to-back primers, then KLD.",
  },
  {
    id: "anneal-oligos",
    ready: true,
    title: "Anneal Oligos and Insert",
    endpoint: "anneal-oligos",
    Form: AnnealOligosForm,
    help: "Clone a short insert (gRNA, shRNA, linker, tag) from two annealed oligos. Vettore "
      + "reads the overhangs from the cut vector, so the oligos fit any vector and enzyme.",
  },
  {
    id: "golden-gate",
    ready: true,
    title: "Golden Gate Assembly",
    endpoint: "golden-gate",
    Form: GoldenGateForm,
    help: "Type IIS assembly of parts into a vector; overhangs are checked for conflicts.",
  },
  {
    id: "gateway-bp",
    ready: true,
    title: "Gateway BP Reaction",
    endpoint: "gateway",
    Form: gatewayForm("BP"),
    help: "Recombine an attB-flanked insert with a donor vector (attP sites) into an entry clone.",
    details: "Only products where both att sites recombined are shown. The entry clone is listed "
      + "first, the attR by-product after it.",
  },
  {
    id: "gateway-lr",
    ready: true,
    title: "Gateway LR Reaction",
    endpoint: "gateway",
    Form: gatewayForm("LR"),
    help: "Recombine an entry clone (attL sites) with a destination vector (attR sites) into an expression clone.",
    details: "Only products where both att sites recombined are shown. The expression clone is "
      + "listed first, the attP by-product after it; a vector with two Gateway cassettes gives "
      + "one clone per cassette.",
  },
  {
    id: "ta",
    ready: true,
    title: "TA or GC Cloning",
    endpoint: "ta-topo",
    Form: taTopoForm([
      { value: "TA", label: "TA cloning (3′ A insert, 3′ T vector)" },
      { value: "GC", label: "GC cloning (3′ G insert, 3′ C vector)" },
    ]),
    help: "Clone a Taq PCR product into a T-vector (or a G-tailed product into a C-vector); both orientations are shown.",
    details: "The vector file is the vector as double-stranded DNA, without its single-base overhangs: "
      + "a linearised vector file, or a circular map plus the base after which it opens. Vettore adds the "
      + "3′ T (or C) to the vector and the 3′ A (or G) to the PCR product.",
  },
  {
    id: "topo",
    ready: true,
    title: "TOPO Cloning",
    endpoint: "ta-topo",
    Form: taTopoForm([
      { value: "TOPO-TA", label: "TOPO TA (3′ A insert)" },
      { value: "blunt-TOPO", label: "Blunt TOPO (blunt insert)" },
      { value: "directional-TOPO", label: "Directional TOPO (CACC on the forward primer)" },
    ]),
    help: "Clone a PCR product into a topoisomerase-activated vector; directional TOPO needs CACC at the 5′ end of the forward primer.",
    details: "The vector file is the vector as double-stranded DNA, without its single-stranded ends: a "
      + "linearised vector file, or a circular map plus the base after which it opens. For directional "
      + "TOPO the insert is placed in the direction of the vector's numbering, CACC end first.",
  },
  {
    id: "insert-multiple",
    ready: true,
    title: "Insert Multiple Fragments",
    endpoint: "insert-multiple",
    Form: InsertMultipleForm,
    help: "Ligate two or more inserts into a cut vector in one reaction; every product the ends allow is shown.",
    details: "Each insert is a piece cut out with its own enzymes, or a PCR product (blunt). The order "
      + "and orientation of the inserts follow from their ends; the order you list only decides which "
      + "product comes first.",
  },
  {
    id: "ligate",
    ready: true,
    title: "Ligate Fragments",
    endpoint: "ligate",
    Form: LigateFragmentsForm,
    help: "Join two or more fragments end to end into a linear molecule (or close them into a circle).",
    details: "Each fragment is a piece cut out with restriction enzymes, or a whole file or region used "
      + "as a blunt PCR product. Sticky ends join only when their overhangs pair; blunt ends join any "
      + "blunt end.",
  },
];

export const READY_ACTIONS = CLONING_ACTIONS.filter((a) => a.ready);

export const actionById = (id: string) => CLONING_ACTIONS.find((a) => a.id === id) ?? null;
