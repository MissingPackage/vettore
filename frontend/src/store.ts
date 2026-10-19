import { applyMiddleware, combineReducers, createStore } from "redux";
import thunk from "redux-thunk";
import { reducer as form } from "redux-form";
import { tg_modalState } from "@teselagen/ui";
import { vectorEditorMiddleware, vectorEditorReducer } from "@teselagen/ove";

export const store = createStore(
  combineReducers({
    form,
    tg_modalState,
    VectorEditor: vectorEditorReducer(),
  }),
  undefined,
  applyMiddleware(thunk, vectorEditorMiddleware),
);

export type RootState = ReturnType<typeof store.getState>;
