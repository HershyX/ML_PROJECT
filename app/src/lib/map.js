import { useMemo } from 'react';
import { geneRow, groupName } from './util.js';

const mix = (a, b, t) => a.map((v, i) => Math.round(v + (b[i] - v) * t));
const COLD = [70, 120, 235], MID = [200, 204, 214], HOT = [232, 100, 40];
export const zColor = (z) => { const t = Math.max(-3, Math.min(3, z)) / 3; const c = t < 0 ? mix(MID, COLD, -t) : mix(MID, HOT, t); return `rgb(${c[0]},${c[1]},${c[2]})`; };

/** Everything a patient map needs for the current layout and colouring. */
export function useMap(s, over = {}) {
  const { data, tsne, colors, labels } = s;
  const k = over.k ?? s.k, method = over.method ?? s.method, view = over.view ?? s.view;
  const colorBy = over.colorBy ?? s.colorBy, px = s.px, py = s.py, gene = over.gene ?? s.gene;
  return useMemo(() => {
    const n = data.meta.n;
    let xs, ys, xl, yl;
    if (view === 'tsne' && tsne) {
      xs = Float64Array.from({ length: n }, (_, i) => tsne[i * 2]); ys = Float64Array.from({ length: n }, (_, i) => tsne[i * 2 + 1]);
      xl = 'Closeness map, direction 1'; yl = 'Closeness map, direction 2';
    } else {
      xs = data.pc[px]; ys = data.pc[py];
      xl = `Main pattern ${px + 1} (${data.explained[px].toFixed(0)}% of differences)`; yl = `Main pattern ${py + 1} (${data.explained[py].toFixed(0)}%)`;
    }
    const M = data.methods[method][k], lab = M.labels;
    let cols, hollow = null, legend, scale = null;
    const clusterCols = Array.from(lab, (l) => colors[l % colors.length]);
    const clusterLegend = M.sizes.map((c, i) => ({ color: colors[i % colors.length], label: groupName(i), n: c }));
    if (colorBy === 'gene' && gene && data.genes.includes(gene)) {
      const row = geneRow(data, data.genes.indexOf(gene));
      cols = Array.from(row, zColor); legend = null; scale = gene;
    } else if (colorBy === 'labels' && labels) {
      cols = labels.values.map((v) => colors[labels.cats.indexOf(v) % colors.length]);
      legend = labels.cats.map((c) => ({ color: colors[labels.cats.indexOf(c) % colors.length], label: c, n: labels.values.filter((v) => v === c).length }));
    } else if (colorBy === 'stability' && method === 'kmeans' && data.methods.kmeans[k].stability) {
      cols = clusterCols; hollow = Array.from(M.stability, (v) => v < 0.75); legend = clusterLegend;
    } else { cols = clusterCols; legend = clusterLegend; }
    return { xs, ys, xl, yl, cols, hollow, legend, scale, labels: lab, M };
  }, [data, tsne, colors, labels, k, method, view, colorBy, px, py, gene]);
}
