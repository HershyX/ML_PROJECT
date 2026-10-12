export const PALETTE = {
  light: ['#0072b2', '#e69f00', '#009e73', '#cc79a7', '#56b4e9', '#d55e00'],
  dark: ['#56a8e8', '#f0b429', '#2fc79a', '#e08fc0', '#8fd0f5', '#ff8a4c'],
};

export const groupName = (i) => 'Group ' + 'ABCDEF'[i];

export const METHOD_LABEL = { kmeans: 'Standard grouping', ward: 'Family-tree grouping', gmm: 'Probability-based grouping' };
export const METHOD_TECH = { kmeans: 'k-means clustering', ward: 'hierarchical clustering (Ward)', gmm: 'Gaussian mixture model' };
export const METHOD_DESC = {
  kmeans: 'Places a centre point in each group and puts every patient with the nearest one. Fast and the most widely used.',
  ward: 'Starts with every patient alone and keeps joining the two most similar groups, like building a family tree.',
  gmm: 'Treats each group as a soft cloud of patients, so a borderline patient can belong partly to two groups.',
};

/** Plain-language reading of a silhouette score. */
export function sepWord(s) {
  if (s >= 0.5) return { label: 'Clearly separated', tone: 'good', text: 'The groups are distinct from each other.' };
  if (s >= 0.25) return { label: 'Moderately separated', tone: 'ok', text: 'The groups are visible but touch at the edges.' };
  return { label: 'Overlapping', tone: 'warn', text: 'The groups blend into each other at the edges, so treat them as a starting point rather than firm categories.' };
}

export function agreeWord(a) {
  const pct = Math.round(Math.max(0, a) * 100);
  if (a >= 0.8) return { pct, label: 'Strong agreement', tone: 'good' };
  if (a >= 0.5) return { pct, label: 'Partial agreement', tone: 'ok' };
  return { pct, label: 'Little agreement', tone: 'warn' };
}

export function zWord(z) {
  const a = Math.abs(z);
  if (a < 0.3) return 'About average';
  const dir = z > 0 ? 'higher' : 'lower';
  return a > 1 ? `Much ${dir}` : `Somewhat ${dir}`;
}

export function bestK(data, method = 'kmeans') {
  let best = 2, bs = -Infinity;
  for (const k of [2, 3, 4, 5, 6]) { const s = data.methods[method][k].silhouette; if (s > bs) { bs = s; best = k; } }
  return best;
}

export function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

export function readHash() {
  const out = {};
  new URLSearchParams(location.hash.slice(1)).forEach((v, k) => { out[k] = v; });
  return out;
}

export function writeHash(obj) {
  const p = new URLSearchParams();
  Object.entries(obj).forEach(([k, v]) => { if (v !== '' && v != null) p.set(k, v); });
  history.replaceState(null, '', '#' + p.toString());
}

export function geneRow(data, g) { return data.Zg.subarray(g * data.meta.n, (g + 1) * data.meta.n); }

/** Genes whose average differs most between a set of patients and everyone else. */
export function topGenesFor(data, members, top = 8) {
  const { n, p } = data.meta, inSel = new Uint8Array(n);
  members.forEach((i) => { inSel[i] = 1; });
  const ns = members.length, nr = n - ns;
  if (!ns || !nr) return [];
  const diff = new Float64Array(p);
  for (let g = 0; g < p; g++) {
    let a = 0, b = 0; const row = geneRow(data, g);
    for (let i = 0; i < n; i++) { if (inSel[i]) a += row[i]; else b += row[i]; }
    diff[g] = a / ns - b / nr;
  }
  return Array.from(diff.keys()).sort((x, y) => Math.abs(diff[y]) - Math.abs(diff[x])).slice(0, top).map((g) => ({ g, gene: data.genes[g], d: diff[g] }));
}
