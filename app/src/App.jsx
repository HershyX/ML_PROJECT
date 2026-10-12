import React, { useEffect, useMemo, useRef, useState, useCallback } from 'react';
import Worker from './lib/worker.js?worker';
import { PALETTE, bestK, readHash, writeHash } from './lib/util.js';
import Story from './sections/Story.jsx';
import Explore from './sections/Explore.jsx';
import DataTab from './sections/DataTab.jsx';
import Compare from './sections/Compare.jsx';
import Genes from './sections/Genes.jsx';
import { Glossary } from './lib/ui.jsx';

const SECTIONS = [
  ['story', 'Guide', 'M4 5h16M4 12h10M4 19h16'],
  ['explore', 'Patient map', 'M12 3a9 9 0 100 18 9 9 0 000-18zm0 5l3 7-7-3z'],
  ['data', 'Data check', 'M5 20V10m7 10V4m7 16v-7'],
  ['compare', 'Methods', 'M4 6h7v12H4zm9 0h7v12h-7z'],
  ['genes', 'Genes', 'M7 3c0 6 10 6 10 12M17 3c0 6-10 6-10 12M8 7h8M8 15h8'],
];

const Icon = ({ d }) => <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round"><path d={d} /></svg>;

export default function App() {
  const h0 = useMemo(readHash, []);
  const [section, setSection] = useState(h0.sec || 'story');
  const [k, setK] = useState(+h0.k || 0);
  const [method, setMethod] = useState(h0.m || 'kmeans');
  const [view, setView] = useState(h0.v || 'pca');
  const [px, setPx] = useState(+h0.px || 0);
  const [py, setPy] = useState(h0.py !== undefined ? +h0.py : 1);
  const [colorBy, setColorBy] = useState(h0.c || 'cluster');
  const [gene, setGene] = useState(h0.gene || '');
  const [pinned, setPinned] = useState(h0.pin !== undefined ? +h0.pin : -1);
  const [selection, setSelection] = useState([]);
  const [theme, setTheme] = useState('auto');
  const [data, setData] = useState(null);
  const [tsne, setTsne] = useState(null);
  const [status, setStatus] = useState({ step: 'Loading data', pct: 0 });
  const [error, setError] = useState('');
  const [dataset, setDataset] = useState('TCGA RNA-seq (built-in)');
  const [labels, setLabels] = useState(null);
  const [toast, setToast] = useState('');
  const [help, setHelp] = useState(false);
  const worker = useRef(null);
  const fileRef = useRef(null), labelRef = useRef(null);

  const flash = useCallback((t) => { setToast(t); setTimeout(() => setToast(''), 2200); }, []);

  const run = useCallback((text, name) => {
    worker.current?.terminate();
    const w = new Worker(); worker.current = w;
    setData(null); setTsne(null); setError(''); setPinned(-1); setSelection([]); setLabels(null); setStatus({ step: 'Starting', pct: 1 }); setDataset(name);
    w.onmessage = (e) => {
      const m = e.data;
      if (m.type === 'progress') setStatus({ step: m.step, pct: m.pct });
      else if (m.type === 'result') { setData(m.result); setK((cur) => cur || bestK(m.result)); setGene((g) => g || m.result.methods.kmeans[3].markers[0].top[0].gene); }
      else if (m.type === 'embedding') setTsne(m.tsne);
      else if (m.type === 'error') setError(m.message);
    };
    w.postMessage({ text });
  }, []);

  useEffect(() => {
    fetch('RNA_seq_1.csv').then((r) => r.text()).then((t) => run(t, 'TCGA RNA-seq (built-in)')).catch(() => setError('Could not load the built-in dataset.'));
    return () => worker.current?.terminate();
  }, [run]);

  useEffect(() => {
    writeHash({ sec: section, k: k || '', m: method, v: view, px, py, c: colorBy, gene, pin: pinned >= 0 ? pinned : '' });
  }, [section, k, method, view, px, py, colorBy, gene, pinned]);

  useEffect(() => {
    const r = document.documentElement;
    if (theme === 'auto') r.removeAttribute('data-theme'); else r.setAttribute('data-theme', theme);
  }, [theme]);

  const dark = theme === 'dark' || (theme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches);
  const colors = dark ? PALETTE.dark : PALETTE.light;

  useEffect(() => {
    const onKey = (e) => {
      const t = e.target.tagName; if (t === 'INPUT' || t === 'SELECT' || t === 'TEXTAREA') return;
      if (e.key >= '2' && e.key <= '6') setK(+e.key);
      else if (e.key === '/') { e.preventDefault(); setSection('genes'); setTimeout(() => document.getElementById('gene-search')?.focus(), 50); }
      else if (e.key === 'Escape') { setSelection([]); setPinned(-1); }
    };
    window.addEventListener('keydown', onKey); return () => window.removeEventListener('keydown', onKey);
  }, []);

  const onUpload = (e) => { const f = e.target.files[0]; if (!f) return; f.text().then((t) => run(t, f.name)); e.target.value = ''; };
  const onLabels = (e) => {
    const f = e.target.files[0]; if (!f || !data) return;
    f.text().then((t) => {
      const map = new Map();
      t.split(/\r?\n/).forEach((l) => { const [a, b] = l.split(',').map((x) => x.replace(/^"|"$/g, '').trim()); if (a && b) map.set(a, b); });
      const values = data.samples.map((s) => map.get(s) || 'unknown');
      const matched = values.filter((v) => v !== 'unknown').length;
      if (!matched) { flash('No sample IDs in that file matched this dataset.'); return; }
      setLabels({ name: f.name, values, cats: [...new Set(values)] }); setColorBy('labels'); flash(`Matched ${matched} of ${data.samples.length} samples`);
    });
    e.target.value = '';
  };
  const share = () => navigator.clipboard.writeText(location.href).then(() => flash('Link copied'), () => flash('Copy failed'));

  const s = { data, tsne, k: k || 3, setK, method, setMethod, view, setView, px, setPx, py, setPy, colorBy, setColorBy, gene, setGene, pinned, setPinned, selection, setSelection, colors, dark, theme, labels, setSection, flash, labelRef };

  return (
    <div className="app">
      <nav className="rail" aria-label="Sections">
        <div className="logo" title="Genomic Subtype Explorer">GS</div>
        {SECTIONS.map(([id, label, d]) => (
          <button key={id} className={'railbtn' + (section === id ? ' on' : '')} onClick={() => setSection(id)} aria-label={label}><Icon d={d} /><span>{label}</span></button>
        ))}
      </nav>
      <div className="main">
        <header className="top">
          <div className="dsname"><span className="dot" />{dataset}{data && <em>{data.meta.n} patients, {data.meta.p.toLocaleString()} genes</em>}</div>
          <div className="topbtns">
            <button className="btn" onClick={() => setHelp(true)}>Plain words guide</button>
            <button className="btn" onClick={() => fileRef.current.click()}>Use your own data</button>
            <button className="btn" disabled={!data} onClick={() => labelRef.current.click()}>Add patient labels</button>
            <button className="btn" onClick={share}>Copy link</button>
            <button className="btn" onClick={() => setTheme(dark ? 'light' : 'dark')} aria-label="Toggle theme">{dark ? 'Light' : 'Dark'}</button>
            <input ref={fileRef} type="file" accept=".csv,text/csv" hidden onChange={onUpload} />
            <input ref={labelRef} type="file" accept=".csv,text/csv" hidden onChange={onLabels} />
          </div>
        </header>
        <main className="content">
          {error && <div className="state err"><h2>Something went wrong</h2><p>{error}</p><p className="muted">Your file should be a spreadsheet saved as CSV, with gene names down the first column and one column per patient.</p></div>}
          {!error && !data && (
            <div className="state"><div className="spinner" /><h2>{status.step}</h2><div className="progress"><div style={{ width: status.pct + '%' }} /></div><p className="muted">Everything is calculated on your own computer, so patient data never leaves this page. This can take about half a minute.</p></div>
          )}
          {!error && data && (
            <>
              {section === 'story' && <Story s={s} />}
              {section === 'explore' && <Explore s={s} />}
              {section === 'data' && <DataTab s={s} />}
              {section === 'compare' && <Compare s={s} />}
              {section === 'genes' && <Genes s={s} />}
            </>
          )}
          <p className="notice">For research and teaching only. These groups are statistical patterns, not diagnoses, and should not guide the care of any individual patient.</p>
        </main>
      </div>
      {help && <Glossary onClose={() => setHelp(false)} />}
      {data && !tsne && section !== 'story' && <div className="chip-status"><div className="spinner sm" />Preparing the close neighbours layout</div>}
      {toast && <div className="toast">{toast}</div>}
    </div>
  );
}
