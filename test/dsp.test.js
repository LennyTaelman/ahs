// Checks the analyzer against synthetic vowels with known pitch and resonances.
// Run: node test/dsp.test.js
const { createAnalyzer } = require('../dsp.js');

function synthVowel(sr, f0, formants, seconds) {
  const n = Math.floor(sr * seconds);
  let x = new Float64Array(n);
  // Glottal-like source: impulse train smoothed by a one-pole lowpass (-6 dB/oct extra tilt).
  let phase = 0, lp = 0;
  for (let i = 0; i < n; i++) {
    phase += f0 / sr;
    let s = 0;
    if (phase >= 1) { phase -= 1; s = 1; }
    lp = 0.9 * lp + s;
    x[i] = lp;
  }
  for (const [f, bw] of formants) {
    const rr = Math.exp((-Math.PI * bw) / sr);
    const c1 = 2 * rr * Math.cos((2 * Math.PI * f) / sr), c2 = -rr * rr;
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) y[i] = x[i] + c1 * (y[i - 1] || 0) + c2 * (y[i - 2] || 0);
    x = y;
  }
  let peak = 0;
  for (const v of x) peak = Math.max(peak, Math.abs(v));
  return Float32Array.from(x, (v) => (0.3 * v) / peak);
}

const vowels = {
  ah: [[700, 80], [1100, 90], [2600, 120], [3300, 150], [4500, 200]],
  ee: [[280, 60], [2250, 100], [2900, 120], [3500, 150], [4500, 200]],
  oo: [[320, 60], [800, 80], [2400, 120], [3300, 150], [4500, 200]],
};

let failures = 0;
function check(label, got, want, tolPct) {
  const err = (100 * Math.abs(got - want)) / want;
  const ok = err <= tolPct;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: got ${got.toFixed(1)}, want ${want} (${err.toFixed(1)}% off, limit ${tolPct}%)`);
}

for (const sr of [44100, 48000]) {
  const an = createAnalyzer(sr);
  for (const [name, fm] of Object.entries(vowels)) {
    for (const f0 of [98, 147, 220, 294]) {
      const r = an.analyze(synthVowel(sr, f0, fm, 0.5));
      const tag = `${sr} ${name} @${f0}Hz`;
      check(`${tag} pitch`, r.f0, f0, 0.5);
      // Ah is the vowel this tool is for and must be tight. Closed vowels have a low F1
      // that sits between sparse harmonics, so LPC is pulled toward the nearest harmonic;
      // that is a known limit of the method, and the bounds below only guard against regressions.
      const tol = name === 'ah' ? (f0 > 250 ? 20 : 10) : f0 >= 200 ? 45 : 15;
      check(`${tag} F1`, r.f1, fm[0][0], tol);
      check(`${tag} F2`, r.f2, fm[1][0], tol);
    }
  }
}

// Ring should rise when the upper resonances are stronger (narrower, clustered).
{
  const sr = 48000, an = createAnalyzer(sr);
  const dull = an.analyze(synthVowel(sr, 147, [[700, 80], [1100, 90], [2600, 400], [3300, 500], [4500, 500]], 0.5));
  const ringy = an.analyze(synthVowel(sr, 147, [[700, 80], [1100, 90], [2700, 60], [3000, 60], [3300, 80]], 0.5));
  const ok = ringy.ring > dull.ring + 6;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ring: dull ${dull.ring.toFixed(1)} dB, ringy ${ringy.ring.toFixed(1)} dB`);
}

// Silence and noise must not report a pitch.
{
  const sr = 48000, an = createAnalyzer(sr);
  const silent = an.analyze(new Float32Array(8192));
  const noise = an.analyze(Float32Array.from({ length: 8192 }, () => 0.2 * (Math.random() * 2 - 1)));
  const ok = silent.f0 === 0 && (noise.f0 === 0 || noise.clarity < 0.8);
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} unvoiced: silence f0=${silent.f0}, noise clarity=${noise.clarity.toFixed(2)}`);
}

console.log(failures ? `\n${failures} check(s) failed` : '\nall checks passed');
process.exit(failures ? 1 : 0);
