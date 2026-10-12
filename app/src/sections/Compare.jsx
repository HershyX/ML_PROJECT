import React from 'react';
import Plot from '../lib/Plot.jsx';
import Scatter from '../lib/Scatter.jsx';
import { useMap } from '../lib/map.js';
import { Badge, Legend, Section } from '../lib/ui.jsx';
import { METHOD_LABEL, METHOD_DESC, METHOD_TECH, sepWord, agreeWord } from '../lib/util.js';

const IDS = ['kmeans', 'ward', 'gmm'];

function Mini({ s, id }) {
  const m = useMap(s, { method: id, colorBy: 'cluster', view: 'pca' });
  const sep = sepWord(m.M.silhouette);
  return (
    <div className="card mini">
      <h3>{METHOD_LABEL[id]}</h3>
      <p className="muted small">{METHOD_DESC[id]}</p>
      <Scatter xs={m.xs} ys={m.ys} colors={m.cols} height={230} mini theme={s.theme} />
      <Legend items={m.legend} />
      <div className="row"><Badge tone={sep.tone}>{sep.label}</Badge><span className="muted small">score {m.M.silhouette.toFixed(2)}</span></div>
      <p className="tech">Technical name: {METHOD_TECH[id]}</p>
    </div>
  );
}

export default function Compare({ s }) {
  const { data, k, setK, colors, theme } = s;
  const ari = data.methods.kmeans[k].ari;
  const pairs = [[0, 1], [0, 2], [1, 2]].map(([a, b]) => ({ a, b, ...agreeWord(ari[a][b]) }));
  const avg = pairs.reduce((t, x) => t + x.pct, 0) / 3;
  const verdict = avg >= 80 ? 'The three methods mostly put the same patients together, which makes these groups more trustworthy.' : avg >= 50 ? 'The methods agree on the broad picture but disagree about some patients. The overall groups are reasonable, the exact membership is less certain.' : 'The methods disagree about who belongs where. Take these groups as one possible view of the data, not a firm answer.';
  const ks = [2, 3, 4, 5, 6];
  const consistency = ks.map((v) => { const st = data.methods.kmeans[v].stability; return Math.round((st.reduce((a, b) => a + b, 0) / st.length) * 100); });

  return (
    <Section title="Do different methods agree?" lead="There is more than one way to sort patients into groups. If very different methods find the same groups, we can trust the result more.">
      <div className="card toolbar"><div className="field"><label>Number of groups</label><div className="seg">{ks.map((v) => <button key={v} className={v === k ? 'on' : ''} onClick={() => setK(v)}>{v}</button>)}</div></div></div>
      <div className="three">{IDS.map((id) => <Mini key={id} s={s} id={id} />)}</div>
      <div className="two">
        <div className="card">
          <h3>How much do the methods agree, with {k} groups?</h3>
          <p className="verdict">{verdict}</p>
          <div className="agree">{pairs.map((x) => <div key={x.a + '-' + x.b} className="arow"><span>{METHOD_LABEL[IDS[x.a]]} and {METHOD_LABEL[IDS[x.b]]}</span><div className="atrack"><b className={x.tone} style={{ width: x.pct + '%' }} /></div><em>{x.pct}%</em><Badge tone={x.tone}>{x.label}</Badge></div>)}</div>
          <p className="tech">Technical name: adjusted Rand index.</p>
        </div>
        <div className="card">
          <h3>How clearly separated are the groups at each size?</h3>
          <p className="muted small">Higher is better. Anything below 0.25 means the groups overlap noticeably.</p>
          <Plot theme={theme} height={240} data={IDS.map((id, i) => ({ type: 'bar', name: METHOD_LABEL[id], x: ks.map((v) => v + ' groups'), y: ks.map((v) => data.methods[id][v].silhouette), marker: { color: colors[i] }, hovertemplate: '%{y:.2f}<extra>' + METHOD_LABEL[id] + '</extra>' }))} layout={{ barmode: 'group', legend: { orientation: 'h', y: -0.25 }, yaxis: { title: { text: 'Separation score' } }, margin: { b: 70 } }} />
        </div>
      </div>
      <div className="card">
        <h3>How consistent is the standard grouping when repeated?</h3>
        <p className="muted small">The analysis was restarted 12 times with different random starts. This shows what share of each group's patients stay together on average.</p>
        <Plot theme={theme} height={220} data={[{ type: 'bar', x: ks.map((v) => v + ' groups'), y: consistency, marker: { color: ks.map((v) => (v === k ? colors[0] : colors[4])) }, hovertemplate: '%{y}% stay together<extra></extra>' }]} layout={{ yaxis: { range: [0, 100], ticksuffix: '%' } }} />
      </div>
    </Section>
  );
}
