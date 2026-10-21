// Every menu tool, grouped by module; order within a menu follows this list.
import { tools as alignment } from "./alignment";
import { tools as primerDesign } from "./primerDesign";
import { tools as aav } from "./aav";
import { tools as digestPlanner } from "./digestPlanner";
import { tools as preferences } from "./preferences";
import { tools as help } from "./help";
import type { Tool } from "./types";

export const TOOLS: Tool[] = [...primerDesign, ...alignment, ...digestPlanner, ...aav, ...preferences, ...help];

export const toolById = (id: string) => TOOLS.find((t) => t.id === id) ?? null;
