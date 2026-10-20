import "./shimGlobal";
import React from "react";
import { createRoot } from "react-dom/client";
import { Provider } from "react-redux";
import "@teselagen/ove/ove.css";
import "./theme/tokens.css";
import "./app.css";
import { store } from "./store";
import App from "./App";

createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <Provider store={store}>
      <App />
    </Provider>
  </React.StrictMode>,
);
