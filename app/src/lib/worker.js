import { parseMatrix, preprocess, pca, kmeans, wardTree, cutTree, gmm, silhouette, ari, relabelBySize, stability, tsne } from './ml.js';

const METHODS = ['kmeans', 'ward', 'gmm'];
const KS = [2, 3, 4, 5, 6];
const NPC = 10;
const progress = (step, pct) => postMessage({ type: 'progress', step, pct });

function pearson(a, b) {
  const n = a.length; let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) { const x = a[i] - ma, y = b[i] - mb; sab += x * y; saa += x * x; sbb += y * y; }
  return sab / (Math.sqrt(saa * sbb) || 1);
}

function clusterGeneZ(Zg, n, p, lab, k) {
  const out = new Float32Array(p * k), cnt = new Int32Array(k);
  for (let i = 0; i < n; i++) cnt[lab[i]]++;
  for (let g = 0; g < p; g++) for (let i = 0; i < n; i++) out[g * k + lab[i]] += Zg[g * n + i];
  for (let g = 0; g < p; g++) for (let c = 0; c < k; c++) out[g * k + c] /= cnt[c] || 1;
  return out;
}

function markers(geneZ, genes, p, k) {
  const res = [];
  for (let c = 0; c < k; c++) {
    const diffs = new Float64Array(p);
    for (let g = 0; g < p; g++) {
      let other = 0; for (let q = 0; q < k; q++) if (q !== c) other += geneZ[g * k + q];
      diffs[g] = geneZ[g * k + c] - other / (k - 1);
    }
    const idx = Array.from(diffs.keys()).sort((a, b) => diffs[b] - diffs[a]);
    res.push({ top: idx.slice(0, 8).map((g) => ({ gene: genes[g], z: +geneZ[g * k + c].toFixed(2) })), top50: idx.slice(0, 50).map((g) => genes[g]) });
  }
  return res;
}

function analyze(text) {
  progress('Reading the file', 3);
  const { genes, samples, n, p, raw, genesInRows } = parseMatrix(text);
  progress('Standardising genes', 12);
  const { imp, Zg, mean, variance, missingPerGene, nMissing } = preprocess(raw, n, p);

  progress('Running PCA', 20);
  const { scores, m, explained } = pca(Zg, n, p, NPC);

  const methods = {}, transfers = [];
  const geneZ = {};
  METHODS.forEach((meth) => { methods[meth] = {}; });
  const pctPer = 50 / KS.length;
  KS.forEach((k, ki) => {
    progress(`Clustering with k = ${k}`, 35 + ki * pctPer);
    const km = relabelBySize(kmeans(scores, n, m, k, 1, 20).labels, k);
    const labs = { kmeans: km };
    if (!methods._tree) methods._tree = wardTree(scores, n, m);
    labs.ward = relabelBySize(cutTree(methods._tree, n, k), k);
    labs.gmm = relabelBySize(gmm(scores, n, m, k, km), k);
    for (const meth of METHODS) {
      const lab = labs[meth], sizes = new Array(k).fill(0);
      for (const l of lab) sizes[l]++;
      const gz = clusterGeneZ(Zg, n, p, lab, k);
      geneZ[`${meth}-${k}`] = gz; transfers.push(gz.buffer);
      methods[meth][k] = {
        labels: lab, sizes, silhouette: +silhouette(scores, n, m, lab, k).toFixed(4),
        markers: markers(gz, genes, p, k),
        stability: meth === 'kmeans' ? stability(scores, n, m, k, lab, 12) : null,
      };
    }
    const ar = METHODS.map((a) => METHODS.map((b) => +ari(labs[a], labs[b], n).toFixed(3)));
    methods.kmeans[k].ari = ar;
  });
  delete methods._tree;

  progress('Summarising the data', 90);
  const order = Array.from(variance.keys()).sort((a, b) => variance[b] - variance[a]);
  const heatIdx = order.slice(0, 40);
  const sampleSd = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let mu = 0; for (let g = 0; g < p; g++) mu += imp[g * n + i]; mu /= p;
    let ss = 0; for (let g = 0; g < p; g++) { const d = imp[g * n + i] - mu; ss += d * d; }
    sampleSd[i] = Math.sqrt(ss / (p - 1));
  }
  const col = (g) => Array.from(imp.subarray(g * n, (g + 1) * n));
  const nb = Math.min(12, p), nc = Math.min(15, p);
  const corr = Array.from({ length: nc }, (_, a) => Array.from({ length: nc }, (_, b) => +pearson(imp.subarray(a * n, (a + 1) * n), imp.subarray(b * n, (b + 1) * n)).toFixed(3)));
  let dup = 0; const seen = new Set();
  for (let i = 0; i < n; i++) {
    let h = 2166136261;
    for (let g = 0; g < p; g++) { h = Math.imul(h ^ Math.round(imp[g * n + i] * 1e4), 16777619); }
    const key = h + ':' + imp[i] + ':' + imp[(p - 1) * n + i];
    if (seen.has(key)) dup++; else seen.add(key);
  }
  const missRank = Array.from(missingPerGene.keys()).filter((g) => missingPerGene[g] > 0).sort((a, b) => missingPerGene[b] - missingPerGene[a]);
  let rmin = Infinity, rmax = -Infinity, rsum = 0;
  for (let i = 0; i < imp.length; i++) { const v = imp[i]; if (v < rmin) rmin = v; if (v > rmax) rmax = v; rsum += v; }
  let lowVar = 0; for (let g = 0; g < p; g++) if (variance[g] <= 1e-5) lowVar++;

  const pc = Array.from({ length: Math.min(5, m) }, (_, c) => Float64Array.from({ length: n }, (_, i) => scores[i * m + c]));
  transfers.push(Zg.buffer);
  const result = {
    meta: { n, p, nMissing, genesInRows, duplicates: dup, lowVar, valueMin: rmin, valueMax: rmax, valueMean: rsum / imp.length, missingGenes: missRank.length, missingRows: 0 },
    genes, samples, pc, explained, sampleSd, geneMean: Array.from(mean), geneVar: Array.from(variance),
    heatIdx, Zg, geneZ, methods,
    eda: {
      topVar: order.slice(0, 20).map((g) => ({ gene: genes[g], v: variance[g] })),
      missing: missRank.slice(0, 15).map((g) => ({ gene: genes[g], c: missingPerGene[g] })),
      box: Array.from({ length: nb }, (_, g) => ({ gene: genes[g], v: col(g) })),
      corr, corrGenes: genes.slice(0, nc),
    },
  };
  return { result, transfers, scores, m, n };
}

self.onmessage = (e) => {
  try {
    const { result, transfers, scores, m, n } = analyze(e.data.text);
    progress('Done', 100);
    postMessage({ type: 'result', result }, transfers);
    progress('Computing t-SNE map', 100);
    const Y = tsne(scores, n, m);
    postMessage({ type: 'embedding', tsne: Y }, [Y.buffer]);
  } catch (err) {
    postMessage({ type: 'error', message: err.message || String(err) });
  }
};
