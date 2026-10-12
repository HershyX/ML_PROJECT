import React, { useEffect, useRef } from 'react';
import Plotly from 'plotly.js-cartesian-dist-min';

const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

export default function Plot({ data, layout = {}, height = 260, theme }) {
  const ref = useRef(null);
  useEffect(() => {
    const ink = cssVar('--muted'), grid = cssVar('--soft');
    const base = {
      height, margin: { l: 44, r: 12, t: 8, b: 40 }, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: ink, family: 'Inter, system-ui, sans-serif', size: 12 },
      xaxis: { gridcolor: grid, zerolinecolor: grid, ...layout.xaxis }, yaxis: { gridcolor: grid, zerolinecolor: grid, ...layout.yaxis },
    };
    Plotly.react(ref.current, data, { ...base, ...layout, xaxis: base.xaxis, yaxis: base.yaxis, margin: { ...base.margin, ...layout.margin } }, { displayModeBar: false, responsive: true });
  }, [data, layout, height, theme]);
  useEffect(() => { const el = ref.current; return () => Plotly.purge(el); }, []);
  return <div ref={ref} style={{ width: '100%' }} />;
}
