# Gyroscope Physics

An interactive, rigorous HTML textbook on classical gyroscope physics for readers with graduate-level physics and mathematics background.

The main artifact is [`index.html`](./index.html), a self-contained tutorial with:

- step-by-step derivations and proofs;
- MathJax-rendered equations;
- inline SVG figures and SMIL animations;
- realistic shaded diagrams of gyroscope hardware;
- exercises with solutions.

## Contents

The tutorial covers:

1. rotations, frames, and angular velocity;
2. inertia tensors, kinetic energy, and angular momentum;
3. Euler equations and their variational origin;
4. torque-free rigid bodies and Poinsot geometry;
5. the heavy symmetric top;
6. fast gyroscope asymptotics and slow precession;
7. the geometric mechanics view on `SO(3)`;
8. practical gyroscopes and model limits.

## Reading Locally

Open `index.html` in a browser.

MathJax is loaded from a CDN, so equation rendering needs internet access the first time the page is opened unless MathJax is already cached.

## Validation

The HTML was checked with the `rigorous-explainer` validation scripts:

- TeX delimiter and brace balance;
- raw angle brackets inside math;
- internal anchor links.

The headless browser DOM check was attempted but blocked by the local Chrome sandbox/profile environment during generation.
