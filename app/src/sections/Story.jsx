import React, { useEffect, useRef, useState, useMemo } from 'react';
import Plot from '../lib/Plot.jsx';
import Scatter from '../lib/Scatter.jsx';
import { useMap } from '../lib/map.js';
import { Legend, Badge } from '../lib/ui.jsx';
import { bestK, geneRow, groupName, sepWord, zWord } from '../lib/util.js';

export default function Story({ s }) {
  const { data, colors, setSection, setK, setColorBy } = s;
  const n = data.meta.n, p = data.meta.p;
  const kb = bestK(data), best = data.methods.kmeans[kb], sep = sepWord(best.silhouette);
  const [step, setStep] = useState(0);
  const refs = useRef([]);

  const top10 = data.explained.reduce((a, b) => a + b, 0);
  const steps = useMemo(() => [
    { t: 'The question', h: `${n} patients, ${p.toLocaleString()} gene readings each`, b: 'Researchers measured how active thousands of genes are in tumour samples from many patients. The hope is that patients fall into a few natural groups, such as disease subtypes, that could one day guide care. Nobody can read a table this wide, so we let a computer look for patterns.' },
    { t: 'Simplifying', h: 'A few main patterns explain most differences', b: `Thousands of readings per patient is too much to picture. A statistical method boils them down into a handful of main patterns. The first pattern alone captures ${data.explained[0].toFixed(0)}% of how patients differ from each other, and the first ten together capture ${top10.toFixed(0)}%.` },
    { t: 'The patient map', h: 'One dot for every patient', b: 'Using the two strongest patterns, each patient becomes a single dot. Patients with similar gene activity sit close together, patients who differ sit far apart. No groups have been drawn yet. This is only the layout.' },
    { t: 'Finding groups', h: `The computer finds ${kb} groups`, b: `We asked the computer to sort patients into 2, 3, 4, 5 and 6 groups without telling it what to look for. ${kb} groups gave the clearest result, so that is what the colours show.` },
    { t: 'How confident can we be?', h: sep.label, b: `${sep.text} The separation score is ${best.silhouette.toFixed(2)} on a scale where 1 means perfectly distinct groups and 0 means complete overlap. Rings mark borderline patients who change group when the analysis is repeated. This is a research lead, not a diagnosis.` },
    { t: 'What makes each group different', h: 'Genes that stand out in each group', b: 'For every group we list the genes that are most active compared with all the other groups. A biologist would look at these first to understand what the group might represent.' },
    { t: 'Your turn', h: 'Open the map and explore', b: 'Click any patient, draw a box around a region to see what those patients share, try a different grouping method, or look up a gene you know.' },
  ], [data, n, p, kb, best, sep, top10]);

  useEffect(() => {
    const io = new IntersectionObserver((es) => es.forEach((e) => { if (e.isIntersecting) setStep(+e.target.dataset.i); }), { rootMargin: '-45% 0px -45% 0px' });
    refs.current.forEach((el) => el && io.observe(el)); return () => io.disconnect();
  }, [steps]);

  const mapPlain = useMap(s, { k: kb, method: 'kmeans', view: 'pca', colorBy: 'cluster' });
  const mapStab = useMap(s, { k: kb, method: 'kmeans', view: 'pca', colorBy: 'stability' });
  const grey = useMemo(() => mapPlain.cols.map(() => colors[0]), [mapPlain, colors]);
  const heatZ = useMemo(() => data.heatIdx.map((g) => Array.from(geneRow(data, g))), [data]);
  const order = useMemo(() => Array.from({ length: n }, (_, i) => i), [n]);

  return (
    <div className="story">
      <div className="story-text">
        <div className="intro"><span className="kicker">Genomic Subtype Explorer</span><h1>Do patients fall into natural groups?</h1><p className="lead">A guided tour of a gene activity study, written for clinicians. Scroll to follow the story, or jump straight to the map.</p><button className="btn primary" onClick={() => setSection('explore')}>Skip to the patient map</button></div>
        {steps.map((st, i) => (
          <section key={i} ref={(el) => (refs.current[i] = el)} data-i={i} className={'step' + (step === i ? ' on' : '')}>
            <span className="kicker">Step {i + 1} of {steps.length}. {st.t}</span>
            <h2>{st.h}</h2><p>{st.b}</p>
            {i === steps.length - 1 && <div className="row"><button className="btn primary" onClick={() => { setK(kb); setColorBy('cluster'); setSection('explore'); }}>Open the patient map</button><button className="btn" onClick={() => setSection('data')}>Check the data quality</button></div>}
          </section>
        ))}
      </div>
      <div className="story-vis"><div className="card vis">
        {step === 0 && <><div className="vtitle">Gene activity for the 40 most variable genes. Each row is a gene, each column a patient. Orange is more active than average, blue is less.</div><Plot theme={s.theme} height={360} data={[{ z: heatZ, type: 'heatmap', zmin: -3, zmax: 3, colorscale: [[0, '#4678eb'], [0.5, '#c8ccd6'], [1, '#e86428']], showscale: false, x: order, y: data.heatIdx.map((g) => data.genes[g]), hovertemplate: '%{y}<extra></extra>' }]} layout={{ margin: { l: 56, r: 4, t: 4, b: 28 }, xaxis: { title: { text: 'Patients' }, showticklabels: false }, yaxis: { autorange: 'reversed', tickfont: { size: 9 } } }} /></>}
        {step === 1 && <><div className="vtitle">How much of the differences between patients each main pattern captures</div><Plot theme={s.theme} height={340} data={[{ type: 'bar', x: data.explained.map((_, i) => 'Pattern ' + (i + 1)), y: data.explained, marker: { color: data.explained.map((_, i) => (i < 2 ? colors[0] : colors[4])) }, hovertemplate: '%{y:.1f}% of differences<extra></extra>' }]} layout={{ yaxis: { ticksuffix: '%' } }} /></>}
        {step === 2 && <><div className="vtitle">The patient map before any grouping</div><Scatter xs={mapPlain.xs} ys={mapPlain.ys} colors={grey} dim height={380} xLabel={mapPlain.xl} yLabel={mapPlain.yl} theme={s.theme} /></>}
        {step === 3 && <><div className="vtitle">The patient map coloured by group</div><Scatter xs={mapPlain.xs} ys={mapPlain.ys} colors={mapPlain.cols} height={380} xLabel={mapPlain.xl} yLabel={mapPlain.yl} theme={s.theme} /><Legend items={mapPlain.legend} /></>}
        {step === 4 && <><div className="vtitle">Rings are borderline patients</div><Scatter xs={mapStab.xs} ys={mapStab.ys} colors={mapStab.cols} hollow={mapStab.hollow} height={340} xLabel={mapStab.xl} yLabel={mapStab.yl} theme={s.theme} /><Legend items={mapStab.legend} /><div className="row" style={{ marginTop: 10 }}><Badge tone={sep.tone}>{sep.label}</Badge><span className="muted">{mapStab.hollow.filter(Boolean).length} of {n} patients are borderline</span></div></>}
        {step === 5 && <div className="mk-grid">{best.markers.map((m, i) => <div className="mk" key={i}><h3 style={{ color: colors[i] }}>{groupName(i)} <span>{best.sizes[i]} patients</span></h3>{m.top.slice(0, 5).map((x) => <div className="gline" key={x.gene}><b>{x.gene}</b><span>{zWord(x.z)}</span></div>)}</div>)}</div>}
        {step === 6 && <><div className="vtitle">The patient map, ready to explore</div><Scatter xs={mapPlain.xs} ys={mapPlain.ys} colors={mapPlain.cols} height={380} xLabel={mapPlain.xl} yLabel={mapPlain.yl} theme={s.theme} /><Legend items={mapPlain.legend} /></>}
      </div></div>
    </div>
  );
}
