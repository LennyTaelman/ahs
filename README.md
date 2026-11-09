# Ahs

Live feedback on the shape of a sung vowel, for practising consistent ahs.

A dot shows where the vowel sits: up and down is open or closed (first vocal tract resonance, F1), left and right is dark or bright (second resonance, F2). A meter shows ring (energy around 2 to 4 kHz against everything below) and a needle shows pitch. Sing an ah you like, tap **Set my ah**, and the pad zooms in on that spot so drift is easy to see.

## Running it

It is a static page with no build step and no server side. The browser only gives microphone access on `https` or `localhost`:

- **GitHub Pages**: repo Settings → Pages → Deploy from a branch → pick the branch, folder `/`. Then open the Pages URL on a phone or laptop.
- **Locally**: `python3 -m http.server` in this folder, then open `http://localhost:8000`.

## Files

- `index.html`: the page and its drawing code.
- `dsp.js`: the analysis. Pitch by YIN, resonances by LPC on the signal decimated to about 11 kHz, ring from an FFT band ratio.
- `test/dsp.test.js`: checks the analysis against synthetic vowels. Run with `node test/dsp.test.js`.

## Limits

- Resonance tracking gets unreliable above roughly 300 Hz, where harmonics are too sparse to outline the resonances.
- The reference is one point. A well-sung ah shifts a little with pitch, so some movement across the range is expected.
- Ring depends on microphone and distance. Compare within a session, with the microphone in one place.
