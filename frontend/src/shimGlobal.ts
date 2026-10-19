// OVE and some of its dependencies expect a Node-style `global`.
(window as unknown as { global: Window }).global ||= window;
export {};
