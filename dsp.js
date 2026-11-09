// Signal analysis for Ahs: pitch (YIN), vowel resonances F1/F2 (LPC), and ring.
// Pure functions, no browser APIs, so the same file runs in the page and in node tests.
(function (root) {
  'use strict';

  // Windowed-sinc lowpass used before decimating by M.
  function designLowpass(M) {
    const L = 16 * M + 1;
    const fc = 0.45 / M;
    const h = new Float64Array(L);
    let sum = 0;
    for (let n = 0; n < L; n++) {
      const k = n - (L - 1) / 2;
      const sinc = k === 0 ? 2 * fc : Math.sin(2 * Math.PI * fc * k) / (Math.PI * k);
      h[n] = sinc * (0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (L - 1)));
      sum += h[n];
    }
    for (let n = 0; n < L; n++) h[n] /= sum;
    return h;
  }

  // YIN pitch detector. Returns {f0, clarity}; f0 is 0 when nothing periodic is found.
  function yin(x, start, sr, d) {
    const tauMin = Math.floor(sr / 700);
    const tauMax = d.length - 1;
    const W = tauMax;
    d[0] = 1;
    let running = 0;
    for (let tau = 1; tau <= tauMax; tau++) {
      let s = 0;
      for (let j = 0; j < W; j++) {
        const diff = x[start + j] - x[start + j + tau];
        s += diff * diff;
      }
      running += s;
      d[tau] = running > 0 ? (s * tau) / running : 1;
    }
    let best = -1;
    for (let tau = tauMin; tau < tauMax; tau++) {
      if (d[tau] < 0.15) {
        while (tau + 1 < tauMax && d[tau + 1] < d[tau]) tau++;
        best = tau;
        break;
      }
    }
    if (best < 0) {
      let m = Infinity;
      for (let tau = tauMin; tau < tauMax; tau++) if (d[tau] < m) { m = d[tau]; best = tau; }
      if (m > 0.4) return { f0: 0, clarity: 1 - Math.min(1, m) };
    }
    let t = best;
    if (best > 1 && best < tauMax - 1) {
      const a = d[best - 1], b = d[best], c = d[best + 1];
      const den = a - 2 * b + c;
      if (den !== 0) t = best + (0.5 * (a - c)) / den;
    }
    return { f0: sr / t, clarity: 1 - Math.min(1, d[best]) };
  }

  // Levinson-Durbin: autocorrelation r[0..p] -> prediction polynomial a[0..p], a[0] = 1.
  function levinson(r, p) {
    const a = new Float64Array(p + 1);
    const tmp = new Float64Array(p + 1);
    a[0] = 1;
    let e = r[0];
    for (let i = 1; i <= p; i++) {
      let acc = r[i];
      for (let j = 1; j < i; j++) acc += a[j] * r[i - j];
      const k = -acc / e;
      tmp.set(a);
      for (let j = 1; j < i; j++) a[j] = tmp[j] + k * tmp[i - j];
      a[i] = k;
      e *= 1 - k * k;
      if (e <= 0) break;
    }
    return a;
  }

  // Durand-Kerner roots of the monic polynomial z^p + a[1] z^(p-1) + ... + a[p].
  function polyRoots(a) {
    const p = a.length - 1;
    const re = new Float64Array(p), im = new Float64Array(p);
    let pr = 1, pi = 0;
    for (let k = 0; k < p; k++) {
      re[k] = 0.9 * pr; im[k] = 0.9 * pi;
      const nr = pr * 0.4 - pi * 0.9, ni = pr * 0.9 + pi * 0.4;
      const mag = Math.hypot(nr, ni);
      pr = nr / mag; pi = ni / mag;
    }
    for (let iter = 0; iter < 80; iter++) {
      let move = 0;
      for (let k = 0; k < p; k++) {
        let nr = 1, ni = 0;
        for (let j = 1; j <= p; j++) {
          const t = nr * re[k] - ni * im[k] + a[j];
          ni = nr * im[k] + ni * re[k];
          nr = t;
        }
        let dr = 1, di = 0;
        for (let j = 0; j < p; j++) {
          if (j === k) continue;
          const xr = re[k] - re[j], xi = im[k] - im[j];
          const t = dr * xr - di * xi;
          di = dr * xi + di * xr;
          dr = t;
        }
        const den = dr * dr + di * di || 1e-300;
        const qr = (nr * dr + ni * di) / den, qi = (ni * dr - nr * di) / den;
        re[k] -= qr; im[k] -= qi;
        move += Math.abs(qr) + Math.abs(qi);
      }
      if (move < 1e-10) break;
    }
    return { re, im };
  }

  // In-place radix-2 FFT.
  function fft(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) {
        let t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = (-2 * Math.PI) / len;
      const wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2;
          const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti;
          re[a] += tr; im[a] += ti;
          const t = cr * wr - ci * wi;
          ci = cr * wi + ci * wr;
          cr = t;
        }
      }
    }
  }

  function createAnalyzer(sr) {
    const M = Math.max(1, Math.round(sr / 11025));
    const fs = sr / M;
    const h = designLowpass(M);
    const N = 512;
    const order = 2 * Math.round(fs / 2000);
    const tauMax = Math.floor(sr / 65);
    const need = Math.max((N - 1) * M + h.length, 2 * tauMax + 2);

    const yinBuf = new Float64Array(tauMax + 1);
    const dec = new Float64Array(N);
    const frame = new Float64Array(N);
    const r = new Float64Array(order + 1);
    const fre = new Float64Array(N), fim = new Float64Array(N);

    // Analyzes the last `need` samples of buf.
    function analyze(buf) {
      const end = buf.length;
      const out = { rms: 0, f0: 0, clarity: 0, f1: 0, f2: 0, ring: -60, formants: [] };

      const pStart = end - (2 * tauMax + 1);
      let sq = 0;
      for (let i = pStart; i < end; i++) sq += buf[i] * buf[i];
      out.rms = Math.sqrt(sq / (end - pStart));
      if (out.rms < 1e-6) return out;

      const y = yin(buf, pStart, sr, yinBuf);
      out.f0 = y.f0; out.clarity = y.clarity;

      // Lowpass and decimate the tail to ~11 kHz.
      const dStart = end - ((N - 1) * M + h.length);
      for (let n = 0; n < N; n++) {
        let acc = 0;
        const base = dStart + n * M;
        for (let k = 0; k < h.length; k++) acc += h[k] * buf[base + k];
        dec[n] = acc;
      }

      // Ring: energy in 2-4 kHz relative to energy below 2 kHz.
      for (let n = 0; n < N; n++) {
        fre[n] = dec[n] * (0.5 - 0.5 * Math.cos((2 * Math.PI * n) / (N - 1)));
        fim[n] = 0;
      }
      fft(fre, fim);
      const binHz = fs / N;
      let lo = 0, hi = 0;
      const top = Math.min(4000, fs / 2 - 300);
      for (let k = Math.ceil(60 / binHz); k * binHz < top; k++) {
        const p = fre[k] * fre[k] + fim[k] * fim[k];
        if (k * binHz < 2000) lo += p; else hi += p;
      }
      out.ring = 10 * Math.log10((hi + 1e-20) / (lo + 1e-20));

      // LPC resonances.
      let prev = dec[0];
      frame[0] = 0;
      for (let n = 1; n < N; n++) {
        const w = 0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (N - 1));
        frame[n] = (dec[n] - 0.97 * prev) * w;
        prev = dec[n];
      }
      for (let lag = 0; lag <= order; lag++) {
        let acc = 0;
        for (let n = lag; n < N; n++) acc += frame[n] * frame[n - lag];
        r[lag] = acc;
      }
      if (r[0] <= 0) return out;
      r[0] *= 1.0001;
      const roots = polyRoots(levinson(r, order));
      const fm = [];
      for (let k = 0; k < order; k++) {
        if (roots.im[k] <= 0) continue;
        const f = (Math.atan2(roots.im[k], roots.re[k]) * fs) / (2 * Math.PI);
        const bw = (-Math.log(Math.hypot(roots.re[k], roots.im[k])) * fs) / Math.PI;
        if (f > 200 && f < fs / 2 - 200 && bw < 700) fm.push(f);
      }
      fm.sort((a, b) => a - b);
      out.formants = fm;
      // F1 is the lowest resonance in its plausible range, F2 the next one clearly above it.
      const f1 = fm.find((f) => f >= 230 && f <= 1100);
      const f2 = f1 && fm.find((f) => f > f1 + 180 && f < 2800);
      if (f1 && f2) { out.f1 = f1; out.f2 = f2; }
      return out;
    }

    return { need, analyze, decimatedRate: fs, order };
  }

  const api = { createAnalyzer };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.AhsDSP = api;
})(typeof self !== 'undefined' ? self : this);
