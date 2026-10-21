import "./shimGlobal";
import React from "react";
import { createRoot } from "react-dom/client";
import { Provider } from "react-redux";
import "@teselagen/ove/ove.css";
import "./theme/tokens.css";
import "./app.css";
import "./theme/dark-chrome.css";
import { applyChrome, loadPreferences } from "./preferences";
import { store } from "./store";
import App from "./App";

// Menus write shortcuts out as Ctrl+Shift+S except on a Mac, which has its own key symbols
// (theme/overrides.css).
document.documentElement.dataset.platform = /Mac|iPhone|iPad/.test(navigator.platform) ? "mac" : "pc";

// Before the first paint, so a dark frame does not flash light at start (App keeps it in step
// with Preferences afterwards).
applyChrome(loadPreferences().chrome);

createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <Provider store={store}>
      <App />
    </Provider>
  </React.StrictMode>,
);
