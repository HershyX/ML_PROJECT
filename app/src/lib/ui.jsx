import React from 'react';

export function Legend({ items }) {
  if (!items) return null;
  return <div className="legend">{items.map((it) => <span key={it.label}><i style={{ background: it.color }} />{it.label}{it.n != null && <em>{it.n} patients</em>}</span>)}</div>;
}

export function Badge({ tone = 'ok', children }) { return <span className={'badge ' + tone}>{children}</span>; }

export function Hint({ text }) { return <span className="hint-q" tabIndex={0} data-tip={text} aria-label={text}>?</span>; }

export function Section({ title, lead, children, actions }) {
  return (
    <div className="page">
      <div className="page-head"><div><h1>{title}</h1>{lead && <p className="lead">{lead}</p>}</div>{actions}</div>
      {children}
    </div>
  );
}

const TERMS = [
  ['Gene activity', 'How strongly a gene is switched on in a patient sample. Labs measure it from tissue. Technical name: gene expression.'],
  ['Group', 'A set of patients whose gene activity looks alike. A group may point to a disease subtype, but it needs medical checking. Technical name: cluster.'],
  ['Patient map', 'A picture where each dot is a patient. Dots that sit close together have similar gene activity. Technical name: PCA plot.'],
  ['Main pattern', 'A combination of genes that moves together and explains a large share of the differences between patients. Technical name: principal component.'],
  ['Closeness map', 'A second layout that keeps close neighbours together, good for spotting small groups. Technical name: t-SNE.'],
  ['Separation score', 'A number from below 0 to 1 that says how clearly the groups are apart. Near 1 means distinct groups, near 0 means they overlap. Technical name: silhouette score.'],
  ['Borderline patient', 'A patient who lands in a different group when the analysis is repeated with a different random start. Technical name: low stability.'],
  ['Genes that stand out', 'Genes that are much more active in one group than in the others, useful leads for a biologist. Technical name: marker genes.'],
  ['Agreement between methods', 'Whether different grouping methods put the same patients together. High agreement makes a result more trustworthy. Technical name: adjusted Rand index.'],
  ['Blank readings', 'Measurements that are missing in the file. They are filled in with that gene\'s average so no patient is dropped. Technical name: imputation.'],
];

export function Glossary({ onClose }) {
  return (
    <div className="modal-bg" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()} role="dialog" aria-label="Plain words guide">
        <div className="modal-head"><h2>Plain words guide</h2><button className="btn" onClick={onClose}>Close</button></div>
        <p className="muted">What each term on this page means, with the technical name in case you need it for a report.</p>
        <dl>{TERMS.map(([t, d]) => <div key={t}><dt>{t}</dt><dd>{d}</dd></div>)}</dl>
      </div>
    </div>
  );
}
