import React, { useEffect, useRef, useState, useCallback } from 'react';

const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

/** Canvas scatter with tweened positions, hover tooltip, click to pin and drag to select. */
export default function Scatter({ xs, ys, colors, hollow, selected, pinned = -1, names, onPin, onSelect, xLabel, yLabel, height = 520, mini = false, dim = false, canvasRef, theme }) {
  const local = useRef(null), cv = canvasRef || local;
  const cur = useRef(null), tgt = useRef(null), raf = useRef(0), box = useRef(null);
  const props = useRef({}); props.current = { colors, hollow, selected, pinned, mini, dim, xLabel, yLabel, height };
  const [tip, setTip] = useState(null);
  const n = xs.length;

  const pad = mini ? 8 : 30;
  const toPx = useCallback((u, v, w, h) => [pad + u * (w - 2 * pad), h - pad - v * (h - 2 * pad)], [pad]);

  const draw = useCallback(() => {
    const c = cv.current; if (!c || !cur.current) return;
    const dpr = window.devicePixelRatio || 1, w = c.clientWidth, h = props.current.height;
    if (c.width !== Math.round(w * dpr) || c.height !== Math.round(h * dpr)) { c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); }
    const g = c.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
    const { colors: col, hollow: hol, selected: sel, pinned: pin, mini: mn, dim: dm } = props.current;
    const r = mn ? 2.6 : 4.2, hasSel = sel && sel.size > 0;
    if (!mn) {
      g.fillStyle = css('--muted'); g.font = '12px Inter, system-ui, sans-serif'; g.globalAlpha = 0.9;
      g.textAlign = 'right'; g.fillText(props.current.xLabel || '', w - 10, h - 8);
      g.save(); g.translate(14, 10); g.rotate(-Math.PI / 2); g.textAlign = 'right'; g.fillText(props.current.yLabel || '', 0, 0); g.restore();
    }
    for (let i = 0; i < n; i++) {
      const [x, y] = toPx(cur.current[i * 2], cur.current[i * 2 + 1], w, h);
      const on = !hasSel || sel.has(i);
      g.globalAlpha = dm ? 0.55 : on ? 0.85 : 0.18;
      g.beginPath(); g.arc(x, y, r, 0, 7);
      const colr = dm ? css('--muted') : col[i];
      if (hol && hol[i]) { g.strokeStyle = colr; g.lineWidth = 1.6; g.stroke(); g.globalAlpha *= 0.25; g.fillStyle = colr; g.fill(); }
      else { g.fillStyle = colr; g.fill(); }
    }
    g.globalAlpha = 1;
    if (pin >= 0) {
      const [x, y] = toPx(cur.current[pin * 2], cur.current[pin * 2 + 1], w, h);
      g.strokeStyle = css('--ink'); g.lineWidth = 2; g.beginPath(); g.arc(x, y, r + 5, 0, 7); g.stroke();
    }
    if (box.current) {
      const { x0, y0, x1, y1 } = box.current;
      g.fillStyle = css('--accent'); g.globalAlpha = 0.12; g.fillRect(Math.min(x0, x1), Math.min(y0, y1), Math.abs(x1 - x0), Math.abs(y1 - y0));
      g.globalAlpha = 0.8; g.strokeStyle = css('--accent'); g.lineWidth = 1.5; g.strokeRect(Math.min(x0, x1), Math.min(y0, y1), Math.abs(x1 - x0), Math.abs(y1 - y0)); g.globalAlpha = 1;
    }
  }, [n, toPx, cv]);

  // Retarget positions and animate towards them.
  useEffect(() => {
    const norm = (a) => { let mn = Infinity, mx = -Infinity; for (const v of a) { if (v < mn) mn = v; if (v > mx) mx = v; } const s = mx - mn || 1; return (v) => (v - mn) / s; };
    const nx = norm(xs), ny = norm(ys), t = new Float64Array(n * 2);
    for (let i = 0; i < n; i++) { t[i * 2] = nx(xs[i]); t[i * 2 + 1] = ny(ys[i]); }
    tgt.current = t;
    if (!cur.current || cur.current.length !== t.length) { cur.current = Float64Array.from(t); draw(); return; }
    const from = Float64Array.from(cur.current), t0 = performance.now(); cancelAnimationFrame(raf.current);
    const step = (now) => {
      const k = Math.min(1, (now - t0) / 450), e = 1 - Math.pow(1 - k, 3);
      for (let i = 0; i < t.length; i++) cur.current[i] = from[i] + (t[i] - from[i]) * e;
      draw(); if (k < 1) raf.current = requestAnimationFrame(step);
    };
    raf.current = requestAnimationFrame(step);
    return () => cancelAnimationFrame(raf.current);
  }, [xs, ys, n, draw]);

  useEffect(() => { draw(); }, [colors, hollow, selected, pinned, dim, theme, height, xLabel, yLabel, draw]);
  useEffect(() => { const f = () => draw(); window.addEventListener('resize', f); return () => window.removeEventListener('resize', f); }, [draw]);

  const pos = (e) => { const r = cv.current.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
  const nearest = (mx, my) => {
    const c = cv.current, w = c.clientWidth, h = props.current.height; let best = -1, bd = 144;
    for (let i = 0; i < n; i++) { const [x, y] = toPx(cur.current[i * 2], cur.current[i * 2 + 1], w, h); const d = (x - mx) ** 2 + (y - my) ** 2; if (d < bd) { bd = d; best = i; } }
    return best;
  };
  const down = (e) => { if (mini) return; const [x, y] = pos(e); box.current = { x0: x, y0: y, x1: x, y1: y }; };
  const move = (e) => {
    const [x, y] = pos(e);
    if (box.current) { box.current.x1 = x; box.current.y1 = y; draw(); setTip(null); return; }
    const i = nearest(x, y); setTip(i >= 0 ? { x, y, i } : null);
  };
  const up = (e) => {
    const b = box.current; box.current = null; if (mini) return;
    const [x, y] = pos(e);
    if (b && (Math.abs(b.x1 - b.x0) > 5 || Math.abs(b.y1 - b.y0) > 5)) {
      const c = cv.current, w = c.clientWidth, h = props.current.height, out = [];
      const lx = Math.min(b.x0, b.x1), hx = Math.max(b.x0, b.x1), ly = Math.min(b.y0, b.y1), hy = Math.max(b.y0, b.y1);
      for (let i = 0; i < n; i++) { const [px, py] = toPx(cur.current[i * 2], cur.current[i * 2 + 1], w, h); if (px >= lx && px <= hx && py >= ly && py <= hy) out.push(i); }
      onSelect && onSelect(out);
    } else { const i = nearest(x, y); onPin && onPin(i); }
    draw();
  };

  return (
    <div className="plotbox">
      <canvas ref={cv} style={{ height, cursor: mini ? 'default' : 'crosshair' }} onMouseDown={down} onMouseMove={mini ? undefined : move} onMouseUp={up} onMouseLeave={() => { box.current = null; setTip(null); draw(); }} />
      {tip && !mini && <div className="tip" style={{ left: tip.x, top: tip.y }}>{names ? names[tip.i] : tip.i}</div>}
    </div>
  );
}
