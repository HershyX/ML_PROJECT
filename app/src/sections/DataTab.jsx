import React, { useState, useMemo } from 'react';
import Plot from '../lib/Plot.jsx';
import { Badge, Section } from '../lib/ui.jsx';

export default function DataTab({ s }) {
  const { data, colors, theme } = s;
  const [open, setOpen] = useState(null);
  const { n, p, nMissing, duplicates, lowVar, missingGenes, valueMin, valueMax } = data.meta;
  const total = n * p, pct = (nMissing / total) * 100;

  const checks = [
    { id: 'miss', title: 'Blank readings', value: nMissing.toLocaleString(), tone: pct < 1 ? 'good' : pct < 5 ? 'ok' : 'warn', badge: pct < 1 ? 'Minor' : pct < 5 ? 'Worth noting' : 'Needs attention',
      text: nMissing === 0 ? 'No readings are missing.' : `${nMissing.toLocaleString()} of ${total.toLocaleString()} readings (${pct.toFixed(3)}%) were blank, touching ${missingGenes} genes. Each blank was filled with that gene's average so no patient had to be dropped.` },
    { id: 'dup', title: 'Duplicate patients', value: duplicates, tone: duplicates === 0 ? 'good' : 'warn', badge: duplicates === 0 ? 'None found' : 'Check the file',
      text: duplicates === 0 ? 'Every patient has a different set of readings, so no one is counted twice.' : `${duplicates} patients have readings identical to another patient and may be counted twice.` },
    { id: 'flat', title: 'Genes that barely change', value: lowVar, tone: lowVar === 0 ? 'good' : 'ok', badge: lowVar === 0 ? 'None found' : 'Worth noting',
      text: lowVar === 0 ? 'Every gene differs between patients, so every gene has something to tell us.' : `${lowVar} genes show almost no difference between patients and add little to the grouping.` },
    { id: 'size', title: 'Size of the study', value: `${n} patients`, tone: n >= 100 ? 'good' : 'ok', badge: n >= 100 ? 'Good size' : 'Small', text: `${p.toLocaleString()} genes were measured per patient. Readings range from ${valueMin.toFixed(1)} to ${valueMax.toFixed(1)}.` },
  ];

  const geneMeans = useMemo(() => data.geneMean, [data]);
  const cards = [
    { id: 'means', title: 'How active are the genes overall?', what: 'Each bar counts genes with a similar average activity.', mean: 'Most genes sit in the middle. A long tail to one side means a few genes are much more active than the rest.',
      el: (h) => <Plot theme={theme} height={h} data={[{ type: 'histogram', x: geneMeans, nbinsx: 40, marker: { color: colors[0] }, hovertemplate: '%{y} genes<extra></extra>' }]} layout={{ xaxis: { title: { text: 'Average gene activity' } }, yaxis: { title: { text: 'Number of genes' } }, bargap: 0.04 }} /> },
    { id: 'var', title: 'Which genes differ most between patients?', what: 'The 20 genes whose activity varies the most from one patient to the next.', mean: 'Highly variable genes are the ones most likely to separate patients into groups.',
      el: (h) => <Plot theme={theme} height={h} data={[{ type: 'bar', orientation: 'h', y: data.eda.topVar.map((x) => x.gene).reverse(), x: data.eda.topVar.map((x) => x.v).reverse(), marker: { color: colors[4] }, hovertemplate: '%{y}<extra></extra>' }]} layout={{ margin: { l: 70 }, xaxis: { title: { text: 'How much it varies' } } }} /> },
    { id: 'sd', title: 'How much does each patient stand out?', what: 'Each bar counts patients whose readings are spread out by a similar amount.', mean: 'A few patients far from the crowd may be unusual cases or measurement problems and are worth a closer look.',
      el: (h) => <Plot theme={theme} height={h} data={[{ type: 'histogram', x: data.sampleSd, nbinsx: 30, marker: { color: colors[2] }, hovertemplate: '%{y} patients<extra></extra>' }]} layout={{ xaxis: { title: { text: 'Spread of a patient\'s readings' } }, yaxis: { title: { text: 'Number of patients' } }, bargap: 0.04 }} /> },
    { id: 'box', title: 'Range of activity for 12 example genes', what: 'Each box shows where the middle half of patients fall. Dots are unusually high or low patients.', mean: 'Dots far from the box are not errors by default. In biology, extreme values can be the most interesting patients.',
      el: (h) => <Plot theme={theme} height={h} data={data.eda.box.map((b, i) => ({ type: 'box', y: b.v, name: b.gene, marker: { color: colors[i % 6] }, boxpoints: 'outliers', hoverinfo: 'y+name' }))} layout={{ showlegend: false, yaxis: { title: { text: 'Gene activity' } } }} /> },
    { id: 'corr', title: 'Which genes move together?', what: 'Each square compares two genes. Orange means they rise and fall together, blue means one rises as the other falls.', mean: 'Genes that move together often work in the same biological process, which is why a few patterns can summarise thousands of genes.',
      el: (h) => <Plot theme={theme} height={h} data={[{ type: 'heatmap', z: data.eda.corr, x: data.eda.corrGenes, y: data.eda.corrGenes, zmin: -1, zmax: 1, colorscale: [[0, '#4678eb'], [0.5, '#e9ebf0'], [1, '#e86428']], colorbar: { thickness: 10, len: 0.8 } }]} layout={{ margin: { l: 64, b: 64 }, yaxis: { autorange: 'reversed' } }} /> },
  ];
  if (data.eda.missing.length) cards.push({ id: 'gaps', title: 'Where are the blank readings?', what: 'The genes with the most blank readings.', mean: 'A handful of genes account for the gaps, so the overall picture is not affected.',
    el: (h) => <Plot theme={theme} height={h} data={[{ type: 'bar', y: data.eda.missing.map((x) => x.gene).reverse(), x: data.eda.missing.map((x) => x.c).reverse(), orientation: 'h', marker: { color: colors[1] } }]} layout={{ margin: { l: 70 }, xaxis: { title: { text: 'Blank readings' } } }} /> });

  const cur = cards.find((c) => c.id === open);
  return (
    <Section title="Data check" lead="Before trusting any grouping, we check the data is sound. Here is what we found, in plain terms.">
      <div className="tiles">{checks.map((c) => <div className="card tile" key={c.id}><span className="muted small">{c.title}</span><b>{c.value}</b><Badge tone={c.tone}>{c.badge}</Badge><p>{c.text}</p></div>)}</div>
      <h2 className="sub">Charts. Click any chart to enlarge it.</h2>
      <div className="cards">{cards.map((c) => (
        <div className="card chart" key={c.id}><h3>{c.title}</h3><p className="muted small">{c.what}</p><button className="plain" onClick={() => setOpen(c.id)} aria-label={`Enlarge: ${c.title}`}>{c.el(240)}</button><p className="meaning"><b>What it means.</b> {c.mean}</p></div>
      ))}</div>
      {cur && <div className="modal-bg" onClick={() => setOpen(null)}><div className="modal big" onClick={(e) => e.stopPropagation()}><div className="modal-head"><h2>{cur.title}</h2><button className="btn" onClick={() => setOpen(null)}>Close</button></div><p className="muted">{cur.what}</p>{cur.el(460)}<p className="meaning"><b>What it means.</b> {cur.mean}</p></div></div>}
    </Section>
  );
}
