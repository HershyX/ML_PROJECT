import React, { useMemo, useRef } from 'react';
import Scatter from '../lib/Scatter.jsx';
import { useMap } from '../lib/map.js';
import { Legend, Badge, Hint } from '../lib/ui.jsx';
import { METHOD_LABEL, METHOD_DESC, METHOD_TECH, groupName, sepWord, zWord, download, geneRow, topGenesFor } from '../lib/util.js';

export default function Explore({ s }) {
  const { data, k, setK, method, setMethod, view, setView, colorBy, setColorBy, gene, setGene, pinned, setPinned, selection, setSelection, colors, px, setPx, py, setPy, tsne, labels, labelRef } = s;
  const m = useMap(s);
  const cv = useRef(null);
  const n = data.meta.n;
  const selSet = useMemo(() => (selection.length ? new Set(selection) : null), [selection]);
  const sep = sepWord(m.M.silhouette);

  const colorOptions = [['cluster', 'Group'], ['stability', 'Borderline patients'], ['gene', 'Activity of one gene'], ['labels', 'My own labels']];
  const savePng = () => cv.current.toBlob((b) => download('patient-map.png', b));
  const saveCsv = () => {
    const rows = ['patient,group'].concat(data.samples.map((id, i) => `${id},${groupName(m.labels[i])}`));
    download(`patient-groups-${k}.csv`, new Blob([rows.join('\n')], { type: 'text/csv' }));
  };

  return (
    <div className="page wide">
      <div className="page-head"><div><h1>Patient map</h1><p className="lead">Each dot is a patient. Patients with similar gene activity sit close together. Click a dot to see that patient, or draw a box around several dots to see what they share.</p></div></div>
      <div className="explore">
        <div className="card map-card">
          <div className="toolbar">
            <div className="field"><label>Number of groups <Hint text="How many groups the computer should sort patients into. Try a few and watch how the picture changes. Keys 2 to 6 also work." /></label>
              <div className="seg">{[2, 3, 4, 5, 6].map((v) => <button key={v} className={v === k ? 'on' : ''} onClick={() => setK(v)}>{v}</button>)}</div></div>
            <div className="field"><label>Grouping method <Hint text={`${METHOD_DESC[method]} Technical name: ${METHOD_TECH[method]}.`} /></label>
              <select value={method} onChange={(e) => setMethod(e.target.value)}>{Object.entries(METHOD_LABEL).map(([id, l]) => <option key={id} value={id}>{l}</option>)}</select></div>
            <div className="field"><label>Layout <Hint text="Overall pattern shows the biggest differences between patients. Close neighbours puts patients with the most similar readings next to each other and is good for small groups." /></label>
              <div className="seg"><button className={view === 'pca' ? 'on' : ''} onClick={() => setView('pca')}>Overall pattern</button><button className={view === 'tsne' ? 'on' : ''} disabled={!tsne} onClick={() => setView('tsne')}>{tsne ? 'Close neighbours' : 'Preparing'}</button></div></div>
            <div className="field"><label>Colour dots by</label>
              <select value={colorBy} onChange={(e) => { setColorBy(e.target.value); if (e.target.value === 'labels' && !labels) labelRef.current.click(); }}>{colorOptions.map(([id, l]) => <option key={id} value={id}>{l}</option>)}</select></div>
            {colorBy === 'gene' && <div className="field"><label>Gene</label><input list="genes" value={gene} onChange={(e) => setGene(e.target.value.toUpperCase())} placeholder="For example EGFR" /><datalist id="genes">{data.genes.slice(0, 5000).map((g) => <option key={g} value={g} />)}</datalist></div>}
          </div>
          <Scatter xs={m.xs} ys={m.ys} colors={m.cols} hollow={m.hollow} selected={selSet} pinned={pinned} names={data.samples} onPin={(i) => { setPinned(i); setSelection([]); }} onSelect={(idx) => { setSelection(idx); setPinned(-1); }} xLabel={m.xl} yLabel={m.yl} height={520} canvasRef={cv} theme={s.theme} />
          {m.scale ? <div className="scale"><span>Less active than average</span><i /><span>More active than average</span><b>{m.scale}</b></div> : <Legend items={m.legend} />}
          {colorBy === 'stability' && method !== 'kmeans' && <p className="muted small">Borderline marking is only available for the standard grouping method.</p>}
          <div className="map-foot">
            <div className="row"><Badge tone={sep.tone}>{sep.label}</Badge><span className="muted small">{sep.text}</span></div>
            <div className="row"><button className="btn" onClick={savePng}>Save picture</button><button className="btn" onClick={saveCsv}>Download groups as spreadsheet</button></div>
          </div>
          <details className="adv"><summary>Advanced: choose which patterns to plot</summary>
            <div className="row"><label>Across <select value={px} onChange={(e) => setPx(+e.target.value)}>{data.pc.map((_, i) => <option key={i} value={i}>Main pattern {i + 1}</option>)}</select></label><label>Up <select value={py} onChange={(e) => setPy(+e.target.value)}>{data.pc.map((_, i) => <option key={i} value={i}>Main pattern {i + 1}</option>)}</select></label></div></details>
        </div>
        <Drawer s={s} m={m} />
      </div>
    </div>
  );
}

function Drawer({ s, m }) {
  const { data, k, method, pinned, setPinned, selection, setSelection, colors } = s;
  const n = data.meta.n;
  const lab = m.labels;

  const profile = useMemo(() => {
    if (pinned < 0) return null;
    const ups = [], downs = [];
    const all = [];
    for (let g = 0; g < data.meta.p; g++) all.push([g, geneRow(data, g)[pinned]]);
    all.sort((a, b) => b[1] - a[1]);
    all.slice(0, 4).forEach(([g, z]) => ups.push({ gene: data.genes[g], z })); all.slice(-4).reverse().forEach(([g, z]) => downs.push({ gene: data.genes[g], z }));
    const d = (j) => { let t = 0; for (let c = 0; c < data.pc.length; c++) { const x = data.pc[c][pinned] - data.pc[c][j]; t += x * x; } return t; };
    const near = Array.from({ length: n }, (_, j) => j).filter((j) => j !== pinned).sort((a, b) => d(a) - d(b)).slice(0, 5);
    return { ups, downs, near };
  }, [data, pinned, n]);

  const sel = useMemo(() => {
    if (!selection.length) return null;
    const comp = new Array(k).fill(0); selection.forEach((i) => { comp[lab[i]]++; });
    return { comp, genes: topGenesFor(data, selection, 8) };
  }, [selection, data, lab, k]);

  if (profile) {
    const stab = data.methods.kmeans[k].stability?.[pinned];
    return (
      <aside className="card drawer">
        <div className="dhead"><h3>Patient {data.samples[pinned]}</h3><button className="btn sm" onClick={() => setPinned(-1)}>Close</button></div>
        <div className="row"><i className="sw" style={{ background: colors[lab[pinned] % 6] }} /><b>{groupName(lab[pinned])}</b>{method === 'kmeans' && stab != null && <Badge tone={stab < 0.75 ? 'warn' : 'good'}>{stab < 0.75 ? 'Borderline' : 'Firmly in this group'}</Badge>}</div>
        <h4>Genes unusually active in this patient</h4>
        {profile.ups.map((x) => <div className="gline" key={x.gene}><b>{x.gene}</b><span>{zWord(x.z)}</span></div>)}
        <h4>Genes unusually quiet in this patient</h4>
        {profile.downs.map((x) => <div className="gline" key={x.gene}><b>{x.gene}</b><span>{zWord(x.z)}</span></div>)}
        <h4>Most similar patients</h4>
        <div className="chips">{profile.near.map((j) => <button key={j} className="gchip" onClick={() => setPinned(j)}><i className="sw sm" style={{ background: colors[lab[j] % 6] }} />{data.samples[j]}</button>)}</div>
      </aside>
    );
  }
  if (sel) {
    return (
      <aside className="card drawer">
        <div className="dhead"><h3>{selection.length} patients selected</h3><button className="btn sm" onClick={() => setSelection([])}>Clear</button></div>
        <h4>Which groups they belong to</h4>
        {sel.comp.map((c, i) => c > 0 && <div className="cbar" key={i}><span>{groupName(i)}</span><div><b style={{ width: (c / selection.length) * 100 + '%', background: colors[i % 6] }} /></div><em>{c}</em></div>)}
        <h4>What these patients share</h4>
        <p className="muted small">Genes that differ most between the selected patients and everyone else.</p>
        {sel.genes.map((x) => <div className="gline" key={x.gene}><b>{x.gene}</b><span>{zWord(x.d)}</span></div>)}
      </aside>
    );
  }
  const M = data.methods[method][k];
  return (
    <aside className="card drawer">
      <h3>How to read this map</h3>
      <ul className="tips"><li>Each dot is one patient.</li><li>Dots close together have similar gene activity.</li><li>Colours show the group the computer placed each patient in.</li><li>Click a dot for one patient, or drag a box over several.</li></ul>
      <h3>The {k} groups</h3>
      {M.markers.map((mk, i) => (
        <div className="grp" key={i}><div className="row"><i className="sw" style={{ background: colors[i % 6] }} /><b>{groupName(i)}</b><span className="muted">{M.sizes[i]} patients, {Math.round((M.sizes[i] / n) * 100)}%</span></div>
          <p className="muted small">Stands out for {mk.top.slice(0, 3).map((x) => x.gene).join(', ')}</p></div>
      ))}
    </aside>
  );
}
