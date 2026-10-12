import React, { useMemo, useState } from 'react';
import Plot from '../lib/Plot.jsx';
import { Badge, Section } from '../lib/ui.jsx';
import { METHOD_LABEL, groupName, zWord, geneRow } from '../lib/util.js';

export default function Genes({ s }) {
  const { data, k, method, gene, setGene, colors, theme } = s;
  const M = data.methods[method][k];
  const gi = data.genes.indexOf(gene.trim().toUpperCase());
  const [grp, setGrp] = useState(0);
  const [enrich, setEnrich] = useState({ state: 'idle' });
  const g = Math.min(grp, k - 1);

  const z = gi >= 0 ? Array.from({ length: k }, (_, c) => data.geneZ[`${method}-${k}`][gi * k + c]) : null;
  const byGroup = useMemo(() => {
    if (gi < 0) return null;
    const row = geneRow(data, gi), out = Array.from({ length: k }, () => []);
    M.labels.forEach((l, i) => out[l].push(row[i]));
    return out;
  }, [gi, data, M, k]);

  const lookup = async () => {
    setEnrich({ state: 'loading' });
    try {
      const body = new FormData();
      body.append('list', M.markers[g].top50.join('\n')); body.append('description', 'Genomic Subtype Explorer group');
      const add = await (await fetch('https://maayanlab.cloud/Enrichr/addList', { method: 'POST', body })).json();
      const lib = 'GO_Biological_Process_2023';
      const res = await (await fetch(`https://maayanlab.cloud/Enrichr/enrich?userListId=${add.userListId}&backgroundType=${lib}`)).json();
      setEnrich({ state: 'done', rows: res[lib].slice(0, 8).map((r) => ({ term: r[1].replace(/\s*\(GO:\d+\)/, ''), p: r[6], n: r[5].length })) });
    } catch (e) { setEnrich({ state: 'error' }); }
  };

  return (
    <Section title="Genes" lead="Look up any gene to see how active it is in each group, or see which genes define a group and what they do in the body.">
      <div className="two">
        <div className="card">
          <h3>Look up a gene</h3>
          <input id="gene-search" className="big-input" list="genelist" value={gene} onChange={(e) => setGene(e.target.value.toUpperCase())} placeholder="Type a gene name, for example EGFR" autoComplete="off" />
          <datalist id="genelist">{data.genes.map((x) => <option key={x} value={x} />)}</datalist>
          {gi < 0 && gene && <p className="muted">No gene called "{gene.replace(/</g, '')}" in this dataset.</p>}
          {z && (
            <>
              <p className="muted small">Average activity {data.geneMean[gi].toFixed(2)}. Compared with the average patient, each group looks like this.</p>
              <div className="zbars">{z.map((v, c) => { const half = Math.min(Math.abs(v), 2.5) / 2.5 * 50; return (
                <div className="zrow" key={c}><span>{groupName(c)}</span><div className="ztrack"><u /><b style={{ background: colors[c % 6], left: v >= 0 ? '50%' : 50 - half + '%', width: half + '%' }} /></div><em>{zWord(v)}</em></div>); })}
              </div>
              <Plot theme={theme} height={230} data={byGroup.map((v, c) => ({ type: 'box', y: v, name: groupName(c), marker: { color: colors[c % 6] }, boxpoints: false, hoverinfo: 'y' }))} layout={{ showlegend: false, yaxis: { title: { text: 'Activity vs average' } } }} />
            </>
          )}
        </div>
        <div className="card">
          <h3>Genes that define each group</h3>
          <p className="muted small">Using {METHOD_LABEL[method].toLowerCase()} with {k} groups. Click a gene to look it up.</p>
          <div className="seg wrap">{M.markers.map((_, c) => <button key={c} className={c === g ? 'on' : ''} onClick={() => { setGrp(c); setEnrich({ state: 'idle' }); }}><i className="sw sm" style={{ background: colors[c % 6] }} />{groupName(c)}</button>)}</div>
          <div className="chips top">{M.markers[g].top.map((x) => <button key={x.gene} className="gchip" onClick={() => setGene(x.gene)}>{x.gene}<em>{zWord(x.z)}</em></button>)}</div>
          <h4>What do these genes do?</h4>
          <p className="muted small">This asks Enrichr, a free public research database, which body processes the group's top 50 genes are linked to. Only the 50 gene names are sent, never patient IDs or readings.</p>
          {enrich.state === 'idle' && <button className="btn primary" onClick={lookup}>Look up processes for {groupName(g)}</button>}
          {enrich.state === 'loading' && <p className="muted">Asking Enrichr</p>}
          {enrich.state === 'error' && <><p className="muted">Could not reach Enrichr. Check your internet connection and try again.</p><button className="btn" onClick={lookup}>Try again</button></>}
          {enrich.state === 'done' && <ol className="proc">{enrich.rows.map((r) => <li key={r.term}><b>{r.term}</b><span className="muted small">{r.n} of the genes match</span><Badge tone={r.p < 0.01 ? 'good' : r.p < 0.05 ? 'ok' : 'warn'}>{r.p < 0.01 ? 'Strong match' : r.p < 0.05 ? 'Likely match' : 'Weak match'}</Badge></li>)}</ol>}
        </div>
      </div>
    </Section>
  );
}
