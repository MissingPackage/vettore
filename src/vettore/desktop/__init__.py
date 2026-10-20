"""Desktop integration: closing the window, one running copy, opening files from the system.

Everything here is independent of the window backend (Qt WebEngine on Linux, WebView2 on
Windows): the modules talk to a pywebview window object only through its public API, so they
are unit-tested with a fake window. Nothing here imports the biology engine, so a second start
that only hands its files to the running copy stays quick.
"""
