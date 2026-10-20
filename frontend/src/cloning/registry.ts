// Cloning actions offered in the Actions menu, with the engine endpoint each one calls.
import type { ComponentType } from "react";
import {
  DeleteFragmentForm,
  DesignPrimersForm,
  GoldenGateForm,
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
  help: string;
}

export const CLONING_ACTIONS: CloningAction[] = [
  {
    id: "restriction",
    title: "Restriction Cloning: Insert Fragment",
    endpoint: "restriction-insert",
    Form: RestrictionInsertForm,
    help: "Cut the vector and the insert with restriction enzymes and ligate the insert fragment into the vector.",
  },
  {
    id: "delete-fragment",
    title: "Delete Restriction Fragment",
    endpoint: "delete-fragment",
    Form: DeleteFragmentForm,
    help: "Remove the fragment between two sites and re-ligate compatible ends.",
  },
  {
    id: "pcr",
    title: "PCR",
    endpoint: "pcr",
    Form: PcrForm,
    help: "Amplify a template with two primers (5′ tails are added to the product).",
  },
  {
    id: "design-primers",
    title: "Design Primers to Amplify a Region",
    endpoint: "design-pcr-primers",
    Form: DesignPrimersForm,
    help: "Choose primers whose annealing parts reach the target Tm; optional 5′ tails.",
  },
  {
    id: "gibson",
    title: "Gibson Assembly",
    endpoint: "overlap-assembly",
    Form: overlapForm("gibson"),
    help: "Insert fragments into a linearized vector; primers carry the overlaps.",
  },
  {
    id: "hifi",
    title: "NEBuilder HiFi DNA Assembly",
    endpoint: "overlap-assembly",
    Form: overlapForm("hifi"),
    help: "Same design as Gibson; NEB recommends 15–30 bp overlaps.",
  },
  {
    id: "infusion",
    title: "In-Fusion Cloning",
    endpoint: "overlap-assembly",
    Form: overlapForm("infusion"),
    help: "15 bp homology to the vector ends on the insert primers.",
  },
  {
    id: "mutagenesis",
    title: "Mutagenesis (Q5 SDM style)",
    endpoint: "mutagenesis",
    Form: MutagenesisForm,
    help: "Substitute, insert or delete bases with back-to-back primers, then KLD.",
  },
  {
    id: "golden-gate",
    title: "Golden Gate Assembly",
    endpoint: "golden-gate",
    Form: GoldenGateForm,
    help: "Type IIS assembly of parts into a vector; overhangs are checked for conflicts.",
  },
  {
    id: "gateway-bp",
    title: "Gateway BP Reaction",
    endpoint: "gateway",
    Form: gatewayForm("BP"),
    help: "attB-flanked insert × attP donor vector → entry clone.",
  },
  {
    id: "gateway-lr",
    title: "Gateway LR Reaction",
    endpoint: "gateway",
    Form: gatewayForm("LR"),
    help: "Entry clone × attR destination vector → expression clone.",
  },
  {
    id: "ta",
    title: "TA or GC Cloning",
    endpoint: "ta-topo",
    Form: taTopoForm([
      { value: "TA", label: "TA cloning" },
      { value: "GC", label: "GC cloning" },
    ]),
    help: "A-tailed PCR product into a T-vector; both orientations are possible.",
  },
  {
    id: "topo",
    title: "TOPO Cloning",
    endpoint: "ta-topo",
    Form: taTopoForm([
      { value: "TOPO-TA", label: "TOPO TA" },
      { value: "blunt-TOPO", label: "Blunt TOPO" },
      { value: "directional-TOPO", label: "Directional TOPO (CACC)" },
    ]),
    help: "Topoisomerase-activated vectors; directional TOPO needs CACC on the forward primer.",
  },
];

export const actionById = (id: string) => CLONING_ACTIONS.find((a) => a.id === id) ?? null;
