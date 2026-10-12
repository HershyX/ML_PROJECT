// Small numerical toolkit used by the analysis worker. Everything runs in the browser.

export function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Parse a CSV expression matrix. Returns gene-major raw values (NaN for missing). */
export function parseMatrix(text) {
  if (text.charCodeAt(0) === 0xfeff) text = text.slice(1);
  const lines = text.split(/\r?\n/).filter((l) => l.length);
  if (lines.length < 3) throw new Error('The file needs a header row and at least two data rows.');
  const split = (l) => l.split(',').map((s) => s.replace(/^"|"$/g, '').trim());
  const header = split(lines[0]);
  const rows = lines.slice(1).map(split);
  const colNames = header.slice(1);
  const rowNames = rows.map((r) => r[0]);
  const genesInRows = rows.length >= colNames.length;
  const genes = genesInRows ? rowNames : colNames;
  const samples = genesInRows ? colNames : rowNames;
  const n = samples.length, p = genes.length;
  if (n < 10) throw new Error('Need at least 10 samples to cluster.');
  const raw = new Float64Array(n * p);
  const val = (s) => { const v = parseFloat(s); return Number.isFinite(v) ? v : NaN; };
  if (genesInRows) {
    for (let g = 0; g < p; g++) for (let s = 0; s < n; s++) raw[g * n + s] = val(rows[g][s + 1]);
  } else {
    for (let s = 0; s < n; s++) for (let g = 0; g < p; g++) raw[g * n + s] = val(rows[s][g + 1]);
  }
  return { genes, samples, n, p, raw, genesInRows };
}

/** Impute with gene mean, then standardise per gene. */
export function preprocess(raw, n, p) {
  const imp = new Float64Array(raw);
  const Zg = new Float32Array(n * p);
  const mean = new Float64Array(p), variance = new Float64Array(p), missingPerGene = new Int32Array(p);
  let nMissing = 0;
  for (let g = 0; g < p; g++) {
    let s = 0, c = 0;
    for (let i = 0; i < n; i++) { const v = raw[g * n + i]; if (v === v) { s += v; c++; } else missingPerGene[g]++; }
    const mu = c ? s / c : 0;
    let ss = 0;
    for (let i = 0; i < n; i++) {
      let v = raw[g * n + i];
      if (v !== v) { v = mu; imp[g * n + i] = mu; }
      ss += (v - mu) * (v - mu);
    }
    const varg = ss / n, sd = Math.sqrt(varg) || 1;
    mean[g] = mu; variance[g] = varg; nMissing += missingPerGene[g];
    for (let i = 0; i < n; i++) Zg[g * n + i] = (imp[g * n + i] - mu) / sd;
  }
  return { imp, Zg, mean, variance, missingPerGene, nMissing };
}

/** PCA through the sample Gram matrix with power iteration and deflation. */
export function pca(Zg, n, p, m = 10) {
  const Zs = new Float32Array(n * p); // sample-major copy for contiguous dot products
  for (let g = 0; g < p; g++) for (let i = 0; i < n; i++) Zs[i * p + g] = Zg[g * n + i];
  const G = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let s = 0; const a = i * p, b = j * p;
      for (let g = 0; g < p; g++) s += Zs[a + g] * Zs[b + g];
      G[i * n + j] = G[j * n + i] = s;
    }
  }
  let trace = 0; for (let i = 0; i < n; i++) trace += G[i * n + i];
  const rnd = mulberry32(7), vecs = [], vals = [];
  for (let c = 0; c < m; c++) {
    let v = new Float64Array(n); for (let i = 0; i < n; i++) v[i] = rnd() - 0.5;
    let lam = 0;
    for (let it = 0; it < 400; it++) {
      for (const u of vecs) { let d = 0; for (let i = 0; i < n; i++) d += v[i] * u[i]; for (let i = 0; i < n; i++) v[i] -= d * u[i]; }
      const w = new Float64Array(n);
      for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += G[i * n + j] * v[j]; w[i] = s; }
      let nr = 0; for (let i = 0; i < n; i++) nr += w[i] * w[i]; nr = Math.sqrt(nr) || 1;
      let diff = 0; for (let i = 0; i < n; i++) { w[i] /= nr; diff += Math.abs(w[i] - v[i]); }
      v = w; lam = nr;
      if (diff < 1e-9) break;
    }
    vecs.push(v); vals.push(lam);
  }
  const scores = new Float64Array(n * m);
  for (let c = 0; c < m; c++) { const sc = Math.sqrt(Math.max(vals[c], 0)); for (let i = 0; i < n; i++) scores[i * m + c] = vecs[c][i] * sc; }
  return { scores, m, explained: vals.map((v) => (v / trace) * 100) };
}

const sqd = (P, m, i, j) => { let s = 0; for (let c = 0; c < m; c++) { const d = P[i * m + c] - P[j * m + c]; s += d * d; } return s; };

export function kmeans(P, n, m, k, seed = 1, nInit = 10, maxIter = 100) {
  const rnd = mulberry32(seed); let best = null;
  for (let r = 0; r < nInit; r++) {
    const C = new Float64Array(k * m);
    const first = Math.floor(rnd() * n);
    for (let c = 0; c < m; c++) C[c] = P[first * m + c];
    const d2 = new Float64Array(n).fill(Infinity);
    for (let q = 1; q < k; q++) {
      let tot = 0;
      for (let i = 0; i < n; i++) {
        let d = 0; for (let c = 0; c < m; c++) { const x = P[i * m + c] - C[(q - 1) * m + c]; d += x * x; }
        if (d < d2[i]) d2[i] = d; tot += d2[i];
      }
      let t = rnd() * tot, pick = n - 1;
      for (let i = 0; i < n; i++) { t -= d2[i]; if (t <= 0) { pick = i; break; } }
      for (let c = 0; c < m; c++) C[q * m + c] = P[pick * m + c];
    }
    const lab = new Int32Array(n).fill(-1);
    let inertia = 0;
    for (let it = 0; it < maxIter; it++) {
      let changed = false; inertia = 0;
      for (let i = 0; i < n; i++) {
        let bd = Infinity, bq = 0;
        for (let q = 0; q < k; q++) {
          let d = 0; for (let c = 0; c < m; c++) { const x = P[i * m + c] - C[q * m + c]; d += x * x; }
          if (d < bd) { bd = d; bq = q; }
        }
        if (lab[i] !== bq) { lab[i] = bq; changed = true; }
        inertia += bd;
      }
      if (!changed) break;
      const cnt = new Int32Array(k); C.fill(0);
      for (let i = 0; i < n; i++) { cnt[lab[i]]++; for (let c = 0; c < m; c++) C[lab[i] * m + c] += P[i * m + c]; }
      for (let q = 0; q < k; q++) if (cnt[q]) for (let c = 0; c < m; c++) C[q * m + c] /= cnt[q];
    }
    if (!best || inertia < best.inertia) best = { labels: Int32Array.from(lab), inertia };
  }
  return best;
}

/** Ward linkage tree via Lance-Williams updates. Returns merge list. */
export function wardTree(P, n, m) {
  const D = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) D[i * n + j] = D[j * n + i] = sqd(P, m, i, j);
  const size = new Float64Array(n).fill(1), active = new Uint8Array(n).fill(1), merges = [];
  for (let step = 0; step < n - 1; step++) {
    let bd = Infinity, a = -1, b = -1;
    for (let i = 0; i < n; i++) {
      if (!active[i]) continue;
      for (let j = i + 1; j < n; j++) if (active[j] && D[i * n + j] < bd) { bd = D[i * n + j]; a = i; b = j; }
    }
    merges.push([a, b]);
    const si = size[a], sj = size[b], dab = D[a * n + b];
    for (let q = 0; q < n; q++) {
      if (!active[q] || q === a || q === b) continue;
      const sk = size[q];
      const nd = ((si + sk) * D[a * n + q] + (sj + sk) * D[b * n + q] - sk * dab) / (si + sj + sk);
      D[a * n + q] = D[q * n + a] = nd;
    }
    size[a] = si + sj; active[b] = 0;
  }
  return merges;
}

export function cutTree(merges, n, k) {
  const par = Int32Array.from({ length: n }, (_, i) => i);
  const find = (x) => { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
  for (let s = 0; s < n - k; s++) par[find(merges[s][1])] = find(merges[s][0]);
  const map = new Map(), lab = new Int32Array(n);
  for (let i = 0; i < n; i++) { const r = find(i); if (!map.has(r)) map.set(r, map.size); lab[i] = map.get(r); }
  return lab;
}

/** Diagonal covariance Gaussian mixture fitted by EM, started from given labels. */
export function gmm(P, n, m, k, init, iters = 80) {
  const mu = new Float64Array(k * m), vr = new Float64Array(k * m), w = new Float64Array(k);
  const gv = new Float64Array(m);
  for (let c = 0; c < m; c++) { let s = 0, ss = 0; for (let i = 0; i < n; i++) { const v = P[i * m + c]; s += v; ss += v * v; } gv[c] = Math.max(ss / n - (s / n) ** 2, 1e-9); }
  const R = new Float64Array(n * k);
  for (let i = 0; i < n; i++) R[i * k + init[i]] = 1;
  const mstep = () => {
    for (let q = 0; q < k; q++) {
      let nk = 0; for (let i = 0; i < n; i++) nk += R[i * k + q];
      nk = Math.max(nk, 1e-9); w[q] = nk / n;
      for (let c = 0; c < m; c++) {
        let s = 0; for (let i = 0; i < n; i++) s += R[i * k + q] * P[i * m + c];
        const mean = s / nk; let v = 0;
        for (let i = 0; i < n; i++) { const d = P[i * m + c] - mean; v += R[i * k + q] * d * d; }
        mu[q * m + c] = mean; vr[q * m + c] = Math.max(v / nk, gv[c] * 1e-3);
      }
    }
  };
  mstep();
  for (let it = 0; it < iters; it++) {
    for (let i = 0; i < n; i++) {
      let mx = -Infinity; const lp = new Float64Array(k);
      for (let q = 0; q < k; q++) {
        let s = Math.log(w[q] + 1e-300);
        for (let c = 0; c < m; c++) { const v = vr[q * m + c], d = P[i * m + c] - mu[q * m + c]; s += -0.5 * (Math.log(2 * Math.PI * v) + (d * d) / v); }
        lp[q] = s; if (s > mx) mx = s;
      }
      let tot = 0; for (let q = 0; q < k; q++) { lp[q] = Math.exp(lp[q] - mx); tot += lp[q]; }
      for (let q = 0; q < k; q++) R[i * k + q] = lp[q] / tot;
    }
    mstep();
  }
  const lab = new Int32Array(n);
  for (let i = 0; i < n; i++) { let b = 0; for (let q = 1; q < k; q++) if (R[i * k + q] > R[i * k + b]) b = q; lab[i] = b; }
  return lab;
}

export function silhouette(P, n, m, lab, k) {
  const cnt = new Int32Array(k); for (let i = 0; i < n; i++) cnt[lab[i]]++;
  let total = 0;
  for (let i = 0; i < n; i++) {
    const sums = new Float64Array(k);
    for (let j = 0; j < n; j++) if (j !== i) sums[lab[j]] += Math.sqrt(sqd(P, m, i, j));
    const own = lab[i];
    if (cnt[own] <= 1) continue;
    const a = sums[own] / (cnt[own] - 1);
    let b = Infinity;
    for (let q = 0; q < k; q++) if (q !== own && cnt[q]) b = Math.min(b, sums[q] / cnt[q]);
    total += (b - a) / Math.max(a, b);
  }
  return total / n;
}

export function ari(a, b, n) {
  const ka = Math.max(...a) + 1, kb = Math.max(...b) + 1;
  const t = new Float64Array(ka * kb), ra = new Float64Array(ka), rb = new Float64Array(kb);
  for (let i = 0; i < n; i++) { t[a[i] * kb + b[i]]++; ra[a[i]]++; rb[b[i]]++; }
  const c2 = (x) => (x * (x - 1)) / 2;
  let sij = 0, sa = 0, sb = 0;
  for (const x of t) sij += c2(x); for (const x of ra) sa += c2(x); for (const x of rb) sb += c2(x);
  const exp = (sa * sb) / c2(n), max = (sa + sb) / 2;
  return max === exp ? 1 : (sij - exp) / (max - exp);
}

export function relabelBySize(lab, k) {
  const cnt = new Array(k).fill(0); for (const l of lab) cnt[l]++;
  const order = cnt.map((c, i) => [c, i]).sort((x, y) => y[0] - x[0]);
  const map = new Array(k); order.forEach(([, i], r) => { map[i] = r; });
  return Int32Array.from(lab, (l) => map[l]);
}

/** Per sample fraction of its reference cluster mates that stay together across reseeded runs. */
export function stability(P, n, m, k, ref, runs = 12) {
  const out = new Float64Array(n);
  const members = Array.from({ length: k }, () => []);
  for (let i = 0; i < n; i++) members[ref[i]].push(i);
  for (let r = 0; r < runs; r++) {
    const L = kmeans(P, n, m, k, 100 + r, 1).labels;
    for (let i = 0; i < n; i++) {
      const mates = members[ref[i]]; if (mates.length <= 1) { out[i] += 1; continue; }
      let same = 0; for (const j of mates) if (j !== i && L[j] === L[i]) same++;
      out[i] += same / (mates.length - 1);
    }
  }
  for (let i = 0; i < n; i++) out[i] /= runs;
  return out;
}

/** Exact t-SNE (n is small). Distances are normalised so the perplexity search is stable. */
export function tsne(P, n, m, { perplexity = 30, iters = 500, seed = 1 } = {}) {
  const rnd = mulberry32(seed);
  const D = new Float64Array(n * n); let mean = 0;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = sqd(P, m, i, j); D[i * n + j] = D[j * n + i] = d; mean += 2 * d; }
  mean /= n * (n - 1); for (let i = 0; i < D.length; i++) D[i] /= mean;
  const Pm = new Float64Array(n * n), target = Math.log(perplexity);
  for (let i = 0; i < n; i++) {
    let beta = 1, lo = -Infinity, hi = Infinity;
    const row = new Float64Array(n);
    for (let t = 0; t < 60; t++) {
      let sum = 0, dp = 0;
      for (let j = 0; j < n; j++) { if (j === i) { row[j] = 0; continue; } row[j] = Math.exp(-D[i * n + j] * beta); sum += row[j]; dp += D[i * n + j] * row[j]; }
      sum = Math.max(sum, 1e-300);
      const H = Math.log(sum) + (beta * dp) / sum, diff = H - target;
      if (Math.abs(diff) < 1e-5) break;
      if (diff > 0) { lo = beta; beta = hi === Infinity ? beta * 2 : (beta + hi) / 2; } else { hi = beta; beta = lo === -Infinity ? beta / 2 : (beta + lo) / 2; }
    }
    let sum = 0; for (let j = 0; j < n; j++) sum += row[j]; sum = Math.max(sum, 1e-300);
    for (let j = 0; j < n; j++) Pm[i * n + j] = row[j] / sum;
  }
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const v = Math.max((Pm[i * n + j] + Pm[j * n + i]) / (2 * n), 1e-12); Pm[i * n + j] = Pm[j * n + i] = v; }
  const Y = new Float64Array(n * 2), iY = new Float64Array(n * 2), gains = new Float64Array(n * 2).fill(1), grad = new Float64Array(n * 2);
  const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
  for (let i = 0; i < Y.length; i++) Y[i] = gauss() * 1e-2;
  const num = new Float64Array(n * n);
  for (let it = 0; it < iters; it++) {
    const exag = it < 100 ? 12 : 1, mom = it < 250 ? 0.5 : 0.8;
    let qs = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const dx = Y[i * 2] - Y[j * 2], dy = Y[i * 2 + 1] - Y[j * 2 + 1];
      const v = 1 / (1 + dx * dx + dy * dy); num[i * n + j] = num[j * n + i] = v; qs += 2 * v;
    }
    grad.fill(0);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const mult = (exag * Pm[i * n + j] - num[i * n + j] / qs) * num[i * n + j];
      const dx = Y[i * 2] - Y[j * 2], dy = Y[i * 2 + 1] - Y[j * 2 + 1];
      grad[i * 2] += 4 * mult * dx; grad[i * 2 + 1] += 4 * mult * dy;
      grad[j * 2] -= 4 * mult * dx; grad[j * 2 + 1] -= 4 * mult * dy;
    }
    for (let i = 0; i < Y.length; i++) {
      gains[i] = Math.sign(grad[i]) !== Math.sign(iY[i]) ? gains[i] + 0.2 : gains[i] * 0.8;
      if (gains[i] < 0.01) gains[i] = 0.01;
      iY[i] = mom * iY[i] - 200 * gains[i] * grad[i]; Y[i] += iY[i];
    }
    let cx = 0, cy = 0; for (let i = 0; i < n; i++) { cx += Y[i * 2]; cy += Y[i * 2 + 1]; }
    cx /= n; cy /= n; for (let i = 0; i < n; i++) { Y[i * 2] -= cx; Y[i * 2 + 1] -= cy; }
  }
  return Y;
}
