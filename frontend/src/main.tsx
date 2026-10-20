import "./shimGlobal";
import React from "react";
import { createRoot } from "react-dom/client";
import { Provider } from "react-redux";
import "@teselagen/ove/ove.css";
import "./theme/tokens.css";
import "./app.css";
import { store } from "./store";
import App from "./App";

// Menus write shortcuts out as Ctrl+Shift+S except on a Mac, which has its own key symbols
// (theme/overrides.css).
document.documentElement.dataset.platform = /Mac|iPhone|iPad/.test(navigator.platform) ? "mac" : "pc";

createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <Provider store={store}>
      <App />
    </Provider>
  </React.StrictMode>,
);
