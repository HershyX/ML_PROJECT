import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { EffectComposer } from "three/addons/postprocessing/EffectComposer.js";
import { RenderPass } from "three/addons/postprocessing/RenderPass.js";
import { UnrealBloomPass } from "three/addons/postprocessing/UnrealBloomPass.js";
import { OutputPass } from "three/addons/postprocessing/OutputPass.js";

const PALETTE = ["#4cc9f0", "#f72585", "#ffd166", "#06d6a0", "#b388ff", "#ff8c42"];
const GROUP_COLORS = PALETTE.map((c) => new THREE.Color(c));
const groupName = (g) => `Group ${"ABCDEF"[g]}`;
const $ = (id) => document.getElementById(id);
const ease = (u) => (u < 0.5 ? 4 * u * u * u : 1 - Math.pow(-2 * u + 2, 3) / 2);
const BASE_SIZE = 2.6;
const FOV = 55;

const state = { layout: "tsne", method: "kmeans", k: 4, gene: null, selected: null, hover: null, isolate: null };

// Data

const data = await (await fetch("data.json")).json();
const N = data.train.ids.length;
const H = data.unseen.ids.length;
const GENES = data.genes.names;
const geneCol = new Map(GENES.map((g, i) => [g, i]));
const Z = new Int8Array(Uint8Array.from(atob(data.genes.z), (c) => c.charCodeAt(0)).buffer);
// rows are the grouped patients first, then the unseen ones
const zOf = (row, col) => Z[row * GENES.length + col] / data.genes.scale;
const MAX_IMPORTS = 300;
const CAP = H + MAX_IMPORTS;
const withRoom = (values) => {
  const a = new Float32Array(CAP * 3);
  a.set(values);
  return a;
};
const LAYOUTS = {
  train: { tsne: Float32Array.from(data.train.tsne), pca: Float32Array.from(data.train.pca) },
  unseen: { tsne: withRoom(data.unseen.tsne), pca: withRoom(data.unseen.pca) },
};
// Unseen patients: the held-out ones from data.json, then any imported from a file
const U = {
  ids: [...data.unseen.ids],
  voteIdx: [...data.unseen.vote_idx],
  voteDist: [...data.unseen.vote_dist],
  imported: [], // { z, missing } for each imported patient, z over model.genes
  count: H,
};
let model = null; // loaded from model.json on the first import

function unseenZ(h, gene) {
  if (h < H) return zOf(N + h, geneCol.get(gene));
  return U.imported[h - H].z[model.index.get(gene)] ?? 0;
}
let cfg = data.configs[`${state.method}-${state.k}`];
let blend = 0; // 0 is the t-SNE layout, 1 is PCA

// Scene

const canvas = $("scene");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
const scene = new THREE.Scene();
scene.background = new THREE.Color("#04050a");
const camera = new THREE.PerspectiveCamera(FOV, 1, 0.1, 4000);
const HOME = new THREE.Vector3(0, 35, 175);
camera.position.copy(HOME);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.35;
controls.minDistance = 15;
controls.maxDistance = 500;

const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));
const bloom = new UnrealBloomPass(new THREE.Vector2(1, 1), 0.8, 0.55, 0.05);
composer.addPass(bloom);
composer.addPass(new OutputPass());

const material = new THREE.ShaderMaterial({
  uniforms: { uScale: { value: 1 } },
  vertexShader: `
    attribute vec3 aColor;
    attribute float aSize;
    attribute float aAlpha;
    uniform float uScale;
    varying vec3 vColor;
    varying float vAlpha;
    void main() {
      vColor = aColor;
      vAlpha = aAlpha;
      vec4 mv = modelViewMatrix * vec4(position, 1.0);
      gl_PointSize = aSize * uScale / -mv.z;
      gl_Position = projectionMatrix * mv;
    }`,
  fragmentShader: `
    varying vec3 vColor;
    varying float vAlpha;
    void main() {
      float d = length(gl_PointCoord - 0.5) * 2.0;
      if (d > 1.0) discard;
      float core = smoothstep(0.4, 0.0, d);
      float halo = pow(1.0 - d, 2.4);
      gl_FragColor = vec4(vColor * (halo * 0.6 + core * 1.2) * vAlpha, 1.0);
    }`,
  transparent: true,
  depthWrite: false,
  blending: THREE.AdditiveBlending,
});

function makePoints(count) {
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(count * 3), 3));
  g.setAttribute("aColor", new THREE.BufferAttribute(new Float32Array(count * 3), 3));
  g.setAttribute("aSize", new THREE.BufferAttribute(new Float32Array(count).fill(BASE_SIZE), 1));
  g.setAttribute("aAlpha", new THREE.BufferAttribute(new Float32Array(count).fill(1), 1));
  const points = new THREE.Points(g, material);
  points.frustumCulled = false;
  scene.add(points);
  const a = g.attributes;
  return { pos: a.position.array, col: a.aColor.array, size: a.aSize.array, alpha: a.aAlpha.array, attr: a };
}

function touch(p, ...names) {
  for (const n of names) p.attr[n].needsUpdate = true;
}

function starDust() {
  const count = 2600;
  const pos = new Float32Array(count * 3);
  const v = new THREE.Vector3();
  for (let i = 0; i < count; i++) {
    v.randomDirection().multiplyScalar(600 + Math.random() * 900).toArray(pos, i * 3);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  const m = new THREE.PointsMaterial({ color: 0x6b7394, size: 1.3, sizeAttenuation: false, transparent: true, opacity: 0.5 });
  scene.add(new THREE.Points(g, m));
}
starDust();

const stars = makePoints(N);
const unseen = makePoints(CAP);
unseen.alpha.fill(0);
const arrived = new Set();
const flying = new Set();

const lines = (() => {
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(30 * 3 * 2), 3));
  g.setAttribute("color", new THREE.BufferAttribute(new Float32Array(30 * 3 * 2), 3));
  const m = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false });
  const seg = new THREE.LineSegments(g, m);
  seg.frustumCulled = false;
  scene.add(seg);
  return seg;
})();

// Animation helpers

const tweens = new Set();
function tween(ms, fn, done) {
  tweens.add({ start: performance.now(), ms, fn, done });
}

function layoutPos(kind, i, out = new THREE.Vector3()) {
  const L = LAYOUTS[kind];
  const j = i * 3;
  return out.set(
    L.tsne[j] + (L.pca[j] - L.tsne[j]) * blend,
    L.tsne[j + 1] + (L.pca[j + 1] - L.tsne[j + 1]) * blend,
    L.tsne[j + 2] + (L.pca[j + 2] - L.tsne[j + 2]) * blend,
  );
}

function currentPos(kind, i, out = new THREE.Vector3()) {
  return out.fromArray((kind === "train" ? stars : unseen).pos, i * 3);
}

function writePositions() {
  const v = new THREE.Vector3();
  for (let i = 0; i < N; i++) layoutPos("train", i, v).toArray(stars.pos, i * 3);
  for (let h = 0; h < U.count; h++) if (!flying.has(h)) layoutPos("unseen", h, v).toArray(unseen.pos, h * 3);
  touch(stars, "position");
  touch(unseen, "position");
}

function setLayout(layout) {
  if (layout === state.layout) return;
  state.layout = layout;
  const from = blend;
  const to = layout === "pca" ? 1 : 0;
  tween(1500, (u) => {
    blend = from + (to - from) * ease(u);
    writePositions();
  });
  renderNote();
}

// Model results in the browser

function vote(h) {
  const idx = U.voteIdx[h];
  const dist = U.voteDist[h];
  const w = new Array(state.k).fill(0);
  idx.forEach((j, m) => (w[cfg.labels[j]] += 1 / (dist[m] + 1e-9)));
  const total = w.reduce((a, b) => a + b, 0);
  const p = w.map((x) => x / total);
  const g = p.indexOf(Math.max(...p));
  const label = p[g] >= 0.8 ? "Clear match" : p[g] >= 0.5 ? "Moderate match" : "Borderline";
  return { p, g, conf: p[g], label };
}

function geneColor(z, out) {
  const t = Math.max(-1, Math.min(1, z / 2.5));
  const low = new THREE.Color("#3b5bdb");
  const mid = new THREE.Color("#4a4e63");
  const high = new THREE.Color("#ffb347");
  return t < 0 ? out.copy(mid).lerp(low, -t) : out.copy(mid).lerp(high, t);
}

const WHITE = new THREE.Color("#ffffff");
function targetColor(kind, i, out = new THREE.Color()) {
  if (state.gene != null) return geneColor(kind === "train" ? zOf(i, state.gene) : unseenZ(i, GENES[state.gene]), out);
  if (kind === "train") return out.copy(GROUP_COLORS[cfg.labels[i]]);
  return out.copy(WHITE).lerp(GROUP_COLORS[vote(i).g], 0.35);
}

function recolor(ms = 700) {
  const c = new THREE.Color();
  const jobs = [["train", stars, N], ["unseen", unseen, U.count]].map(([kind, p, count]) => {
    const from = p.col.slice();
    const to = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) targetColor(kind, i, c).toArray(to, i * 3);
    return { p, from, to };
  });
  tween(ms, (u) => {
    for (const { p, from, to } of jobs) {
      for (let j = 0; j < to.length; j++) p.col[j] = from[j] + (to[j] - from[j]) * u;
      touch(p, "aColor");
    }
  });
}

function focusSet() {
  const s = state.selected;
  if (!s) return null;
  return new Set(s.kind === "train" ? [s.i, ...data.train.neighbors[s.i]] : U.voteIdx[s.i]);
}

function applyEmphasis() {
  const s = state.selected;
  const hv = state.hover;
  const focus = focusSet();
  for (let i = 0; i < N; i++) {
    let a = 1;
    if (state.isolate != null && cfg.labels[i] !== state.isolate) a = 0.1;
    if (focus && !focus.has(i)) a = Math.min(a, 0.3);
    stars.alpha[i] = a;
    let size = focus?.has(i) ? 3.6 : BASE_SIZE;
    if (s?.kind === "train" && s.i === i) size = 7;
    if (hv?.kind === "train" && hv.i === i) size = Math.max(size, 5.5);
    stars.size[i] = size;
  }
  for (let h = 0; h < U.count; h++) {
    unseen.alpha[h] = arrived.has(h) || flying.has(h) ? 1 : 0;
    let size = 4.6;
    if (s?.kind === "unseen" && s.i === h) size = 8.5;
    if (hv?.kind === "unseen" && hv.i === h) size = Math.max(size, 6.5);
    unseen.size[h] = size;
  }
  touch(stars, "aAlpha", "aSize");
  touch(unseen, "aAlpha", "aSize");
}

function updateLines() {
  const s = state.selected;
  const g = lines.geometry;
  if (!s) return g.setDrawRange(0, 0);
  const others = s.kind === "train" ? data.train.neighbors[s.i] : U.voteIdx[s.i];
  const pos = g.attributes.position.array;
  const col = g.attributes.color.array;
  const a = currentPos(s.kind, s.i);
  const b = new THREE.Vector3();
  const c = new THREE.Color();
  others.forEach((j, m) => {
    currentPos("train", j, b);
    a.toArray(pos, m * 6);
    b.toArray(pos, m * 6 + 3);
    c.copy(GROUP_COLORS[cfg.labels[j]]).multiplyScalar(0.55);
    c.toArray(col, m * 6);
    c.toArray(col, m * 6 + 3);
  });
  g.attributes.position.needsUpdate = true;
  g.attributes.color.needsUpdate = true;
  g.setDrawRange(0, others.length * 2);
}

// Camera

function flyTo(target, dist = 55) {
  const fromT = controls.target.clone();
  const fromC = camera.position.clone();
  const dir = fromC.clone().sub(fromT).normalize();
  const toC = target.clone().add(dir.multiplyScalar(dist));
  controls.autoRotate = false;
  tween(1200, (u) => {
    const e = ease(u);
    controls.target.lerpVectors(fromT, target, e);
    camera.position.lerpVectors(fromC, toC, e);
  });
}

let idleTimer;
controls.addEventListener("start", () => {
  controls.autoRotate = false;
  clearTimeout(idleTimer);
});
controls.addEventListener("end", () => {
  clearTimeout(idleTimer);
  idleTimer = setTimeout(() => (controls.autoRotate = !state.selected), 7000);
});

// Unseen patients flying in

const trail = (() => {
  const count = 60;
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(count * 3), 3));
  const colors = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) new THREE.Color("#ffffff").multiplyScalar((i / count) ** 2).toArray(colors, i * 3);
  g.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  const m = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false });
  const line = new THREE.Line(g, m);
  line.frustumCulled = false;
  line.visible = false;
  scene.add(line);
  return { line, count, points: [] };
})();

function pushTrail(v) {
  trail.points.push(v.clone());
  if (trail.points.length > trail.count) trail.points.shift();
  const pos = trail.line.geometry.attributes.position.array;
  const pts = trail.points;
  for (let i = 0; i < trail.count; i++) pts[Math.max(0, i - (trail.count - pts.length))].toArray(pos, i * 3);
  trail.line.geometry.attributes.position.needsUpdate = true;
}

function pulse(at, color) {
  const m = new THREE.MeshBasicMaterial({ color, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false, side: THREE.DoubleSide });
  const ring = new THREE.Mesh(new THREE.RingGeometry(0.85, 1, 64), m);
  ring.position.copy(at);
  scene.add(ring);
  tween(1100, (u) => {
    ring.quaternion.copy(camera.quaternion);
    ring.scale.setScalar(1 + 16 * ease(u));
    m.opacity = 0.9 * (1 - u);
  }, () => {
    scene.remove(ring);
    ring.geometry.dispose();
    m.dispose();
  });
}

function sendIn(h, { follow = true } = {}) {
  if (arrived.has(h) || flying.has(h)) {
    if (follow && arrived.has(h)) select({ kind: "unseen", i: h });
    return;
  }
  flying.add(h);
  const start = new THREE.Vector3().randomDirection().multiplyScalar(320);
  const side = new THREE.Vector3().randomDirection().multiplyScalar(110);
  targetColor("unseen", h).toArray(unseen.col, h * 3);
  touch(unseen, "aColor");
  if (follow) {
    trail.points = [];
    trail.line.visible = true;
    select(null);
    flyTo(layoutPos("unseen", h), 120);
  }
  applyEmphasis();
  const end = new THREE.Vector3();
  const ctrl = new THREE.Vector3();
  const p = new THREE.Vector3();
  tween(follow ? 2600 : 1800 + Math.random() * 900, (u) => {
    const e = 1 - Math.pow(1 - u, 3);
    layoutPos("unseen", h, end);
    ctrl.copy(start).lerp(end, 0.5).add(side);
    p.copy(start).multiplyScalar((1 - e) ** 2).addScaledVector(ctrl, 2 * (1 - e) * e).addScaledVector(end, e * e);
    p.toArray(unseen.pos, h * 3);
    touch(unseen, "position");
    if (follow) pushTrail(p);
  }, () => {
    flying.delete(h);
    arrived.add(h);
    pulse(layoutPos("unseen", h), GROUP_COLORS[vote(h).g]);
    if (follow) {
      trail.line.visible = false;
      select({ kind: "unseen", i: h });
    } else {
      applyEmphasis();
    }
    renderUnseen();
  });
}

// Picking

function pick(x, y) {
  let best = null;
  let bestD = 14 * 14;
  const v = new THREE.Vector3();
  const test = (kind, p, i) => {
    v.fromArray(p.pos, i * 3).project(camera);
    if (v.z > 1) return;
    const d = ((v.x + 1) / 2 * innerWidth - x) ** 2 + ((1 - v.y) / 2 * innerHeight - y) ** 2;
    if (d < bestD) {
      bestD = d;
      best = { kind, i };
    }
  };
  for (let i = 0; i < N; i++) if (stars.alpha[i] > 0.2) test("train", stars, i);
  for (const h of arrived) test("unseen", unseen, h);
  return best;
}

const tooltip = $("tooltip");
let down = null;
canvas.addEventListener("pointerdown", (e) => {
  down = { x: e.clientX, y: e.clientY };
  $("hint").style.opacity = 0;
});
canvas.addEventListener("pointermove", (e) => {
  const hit = pick(e.clientX, e.clientY);
  if (hit?.kind !== state.hover?.kind || hit?.i !== state.hover?.i) {
    state.hover = hit;
    applyEmphasis();
  }
  if (!hit) return (tooltip.hidden = true);
  const id = hit.kind === "train" ? data.train.ids[hit.i] : U.ids[hit.i];
  const g = hit.kind === "train" ? cfg.labels[hit.i] : vote(hit.i).g;
  tooltip.innerHTML = `<b>${id}</b> <span>${hit.kind === "unseen" ? "unseen, lands in " : ""}${groupName(g)}</span>`;
  tooltip.style.left = `${e.clientX + 14}px`;
  tooltip.style.top = `${e.clientY + 14}px`;
  tooltip.hidden = false;
});
canvas.addEventListener("pointerup", (e) => {
  if (!down || Math.hypot(e.clientX - down.x, e.clientY - down.y) > 5) return;
  const hit = pick(e.clientX, e.clientY);
  select(hit, { fly: !!hit });
});
canvas.addEventListener("pointerleave", () => (tooltip.hidden = true));

// Selection and the patient card

const card = $("card");

function select(sel, { fly = true } = {}) {
  state.selected = sel;
  applyEmphasis();
  renderCard();
  if (sel && fly) flyTo(currentPos(sel.kind, sel.i), 50);
  if (!sel) controls.autoRotate = false;
}

function geneBars(zFor, g) {
  const rows = cfg.markers[g].map(([gene]) => {
    const z = zFor(gene);
    const w = Math.min(Math.abs(z) / 4, 1) * 50;
    const left = z >= 0 ? 50 : 50 - w;
    const color = z >= 0 ? "#ffb347" : "#3b5bdb";
    return `<div class="bar"><span class="name">${gene}</span>
      <span class="track centered"><span class="fill" style="left:${left}%;width:${w}%;background:${color}"></span></span>
      <span class="val">${z > 0 ? "+" : ""}${z.toFixed(1)}</span></div>`;
  });
  return `<h3>Genes that mark ${groupName(g)}</h3><div class="bars">${rows.join("")}</div>
    <p class="caption">This patient's level, in standard deviations from the average patient.</p>`;
}

function patientButton(j) {
  const g = cfg.labels[j];
  return `<li><button data-train="${j}"><span class="dot" style="--c:${PALETTE[g]}"></span>${data.train.ids[j]}<small>${groupName(g)}</small></button></li>`;
}

function renderCard() {
  const s = state.selected;
  card.hidden = !s;
  if (!s) return;
  let html = `<button class="close" aria-label="Close">&times;</button>`;
  if (s.kind === "train") {
    const g = cfg.labels[s.i];
    const size = cfg.labels.filter((x) => x === g).length;
    html += `<p class="eyebrow">Patient</p><h2>${data.train.ids[s.i]}</h2>
      <p class="group-line"><span class="dot" style="--c:${PALETTE[g]}"></span>${groupName(g)}, one of ${size} patients</p>
      <h3>Most similar patients</h3><ul class="neighbors">${data.train.neighbors[s.i].map(patientButton).join("")}</ul>
      ${geneBars((gene) => zOf(s.i, geneCol.get(gene)), g)}`;
  } else {
    const v = vote(s.i);
    const bars = v.p.map((p, g) => `<div class="bar"><span class="name">${groupName(g)}</span>
      <span class="track"><span class="fill" style="left:0;width:${p * 100}%;background:${PALETTE[g]}"></span></span>
      <span class="val">${Math.round(p * 100)}%</span></div>`).join("");
    const extra = s.i >= H ? U.imported[s.i - H] : null;
    const missing = extra?.missing
      ? `<p class="caption">${extra.missing.toLocaleString()} genes were missing from the file and were filled with the dataset median.</p>`
      : "";
    html += `<p class="eyebrow">${extra ? "Imported patient" : "Unseen patient"}</p><h2>${U.ids[s.i]}</h2>
      <p class="verdict">Lands in <span style="color:${PALETTE[v.g]}">${groupName(v.g)}</span></p>
      <p>${v.label}. ${Math.round(v.conf * 100)}% of the vote from the 15 most similar patients, weighted by how close they are. The lines show who voted.</p>
      ${missing}
      <h3>Vote</h3><div class="bars">${bars}</div>
      ${geneBars((gene) => unseenZ(s.i, gene), v.g)}`;
  }
  card.innerHTML = html;
}

card.addEventListener("click", (e) => {
  if (e.target.closest(".close")) return select(null);
  const b = e.target.closest("[data-train]");
  if (b) select({ kind: "train", i: +b.dataset.train });
});

// Panels

function renderMetrics() {
  const m = cfg.metrics;
  $("metrics").innerHTML = [
    [m.silhouette.toFixed(3), "Silhouette", "Higher is better, from -1 to 1"],
    [m.davies_bouldin.toFixed(2), "Davies-Bouldin", "Lower is better"],
    [m.calinski_harabasz.toFixed(1), "Calinski-H.", "Higher is better"],
  ].map(([v, name, tip]) => `<div class="metric" title="${tip}"><b>${v}</b><span>${name}</span></div>`).join("");
  renderNote();
}

function renderNote() {
  const s = cfg.metrics.silhouette;
  const sep = s >= 0.5 ? "Groups are clearly separated." : s >= 0.25 ? "Groups are moderately separated." : "Groups overlap a lot, so their edges are soft.";
  const layout = state.layout === "pca"
    ? `PCA layout: the first 3 components, ${Math.round(data.meta.pca3_variance * 100)}% of the variation.`
    : "t-SNE layout: nearby points have similar gene activity. Distances between groups mean little.";
  $("note").textContent = `${sep} ${layout}`;
}

function renderLegend() {
  const legend = $("legend");
  if (state.gene != null) {
    legend.innerHTML = `<div class="gradient"><span>Low</span><i></i><span>High</span><b>${GENES[state.gene]}</b></div>`;
    return;
  }
  legend.innerHTML = Array.from({ length: state.k }, (_, g) => {
    const count = cfg.labels.filter((x) => x === g).length;
    const dim = state.isolate != null && state.isolate !== g ? "dim" : "";
    return `<button class="${dim}" data-g="${g}"><span class="dot" style="--c:${PALETTE[g]}"></span>${groupName(g)}<small>${count}</small></button>`;
  }).join("");
}

$("legend").addEventListener("click", (e) => {
  const b = e.target.closest("[data-g]");
  if (!b) return;
  const g = +b.dataset.g;
  state.isolate = state.isolate === g ? null : g;
  renderLegend();
  applyEmphasis();
});

function renderUnseen() {
  const n = arrived.size;
  $("tally").innerHTML = n
    ? (() => {
        const counts = { "Clear match": 0, "Moderate match": 0, Borderline: 0 };
        for (const h of arrived) counts[vote(h).label]++;
        return `${n} of ${U.count} sent in: ${counts["Clear match"]} clear, ${counts["Moderate match"]} moderate, ${counts.Borderline} borderline.`;
      })()
    : "";
}

// Group labels floating in the scene

const labelLayer = $("labels");
let labelEls = [];
function buildLabels() {
  labelLayer.innerHTML = "";
  labelEls = Array.from({ length: state.k }, (_, g) => {
    const el = document.createElement("div");
    el.className = "glabel";
    el.style.color = PALETTE[g];
    el.innerHTML = `${groupName(g)}<small>${cfg.labels.filter((x) => x === g).length} patients</small>`;
    labelLayer.appendChild(el);
    return el;
  });
}

function updateLabels() {
  const sums = Array.from({ length: state.k }, () => [0, 0, 0, 0]);
  for (let i = 0; i < N; i++) {
    const s = sums[cfg.labels[i]];
    s[0] += stars.pos[i * 3];
    s[1] += stars.pos[i * 3 + 1];
    s[2] += stars.pos[i * 3 + 2];
    s[3]++;
  }
  const v = new THREE.Vector3();
  labelEls.forEach((el, g) => {
    const s = sums[g];
    v.set(s[0] / s[3], s[1] / s[3] + 6, s[2] / s[3]).project(camera);
    const hide = v.z > 1 || state.gene != null || (state.isolate != null && state.isolate !== g);
    el.style.opacity = hide ? 0 : state.selected ? 0.35 : 1;
    el.style.transform = `translate(-50%, -100%) translate(${(v.x + 1) / 2 * innerWidth}px, ${(1 - v.y) / 2 * innerHeight}px)`;
  });
}

// Controls

function setConfig() {
  cfg = data.configs[`${state.method}-${state.k}`];
  state.isolate = null;
  buildLabels();
  renderLegend();
  renderMetrics();
  recolor();
  applyEmphasis();
  renderCard();
  renderUnseen();
}

function bindSeg(id, onPick) {
  $(id).addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b) return;
    for (const x of $(id).children) x.classList.toggle("on", x === b);
    onPick(b.dataset.v);
  });
}

$("k").innerHTML = [2, 3, 4, 5, 6].map((k) => `<button data-v="${k}" class="${k === state.k ? "on" : ""}">${k}</button>`).join("");
bindSeg("layout", setLayout);
bindSeg("method", (v) => {
  state.method = v;
  setConfig();
});
bindSeg("k", (v) => {
  state.k = +v;
  setConfig();
});

$("gene-list").innerHTML = GENES.map((g) => `<option value="${g}">`).join("");
$("gene").addEventListener("input", (e) => {
  const typed = e.target.value.trim();
  const col = geneCol.get(typed) ?? geneCol.get(typed.toUpperCase());
  if (col === undefined) return;
  state.gene = col;
  $("gene-clear").hidden = false;
  renderLegend();
  recolor();
});
$("gene-clear").addEventListener("click", () => {
  state.gene = null;
  $("gene").value = "";
  $("gene-clear").hidden = true;
  renderLegend();
  recolor();
});

$("unseen-pick").innerHTML = data.unseen.ids.map((id, h) => `<option value="${h}">${id}</option>`).join("");
$("send").addEventListener("click", () => sendIn(+$("unseen-pick").value));
$("random").addEventListener("click", () => {
  const left = [...Array(U.count).keys()].filter((h) => !arrived.has(h) && !flying.has(h));
  if (!left.length) return;
  const h = left[Math.floor(Math.random() * left.length)];
  $("unseen-pick").value = h;
  sendIn(h);
});
$("send-all").addEventListener("click", () => {
  select(null);
  flyTo(new THREE.Vector3(), 190);
  [...Array(U.count).keys()].filter((h) => !arrived.has(h) && !flying.has(h))
    .forEach((h, m) => setTimeout(() => sendIn(h, { follow: false }), m * 45));
});

addEventListener("keydown", (e) => {
  if (e.target.closest("input, select")) return;
  if (e.key === "Escape") {
    state.isolate = null;
    renderLegend();
    select(null);
    flyTo(new THREE.Vector3(), HOME.length());
  }
  if (e.key >= "2" && e.key <= "6") $("k").querySelector(`[data-v="${e.key}"]`).click();
});

// Importing patients from a file

async function loadModel() {
  if (model) return model;
  const m = await (await fetch("model.json")).json();
  const f32 = (s) => new Float32Array(Uint8Array.from(atob(s), (c) => c.charCodeAt(0)).buffer);
  model = {
    genes: m.genes,
    upper: m.genes.map((g) => g.toUpperCase()),
    index: new Map(m.genes.map((g, i) => [g, i])),
    k: m.k,
    votes: m.votes,
    mapNeighbours: m.map_neighbours,
    median: f32(m.median),
    mean: f32(m.mean),
    scale: f32(m.scale),
    pcaMean: f32(m.pca_mean),
    components: f32(m.components),
    trainPca: f32(m.train_pca),
    center: m.pca3_center,
    spread: m.pca3_scale,
  };
  return model;
}

// Same layout as the dataset: gene names in the first column, one column per patient
function parsePatients(text, fileName) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/).filter((l) => l.trim());
  if (lines.length < 2) throw new Error("The file has no gene rows.");
  const delim = lines[0].includes("\t") ? "\t" : ",";
  const split = (l) => l.split(delim).map((c) => c.trim().replace(/^"|"$/g, ""));
  const header = split(lines[0]);
  if (header.length < 2) throw new Error("Expected a gene column followed by at least one patient column.");
  const base = fileName.replace(/\.[^.]+$/, "");
  const names = header.slice(1).map((h, j) =>
    !h || /^(value|values|expression)$/i.test(h) ? (header.length === 2 ? base : `${base} ${j + 1}`) : h);
  const values = names.map(() => new Map());
  for (const line of lines.slice(1)) {
    const cells = split(line);
    const gene = cells[0].toUpperCase();
    for (let j = 1; j < header.length; j++) {
      const v = parseFloat(cells[j]);
      if (Number.isFinite(v)) values[j - 1].set(gene, v);
    }
  }
  return names.map((name, j) => ({ name, values: values[j] }));
}

// The browser version of PatientClassifier.classify in ui/pipeline.py
function addImported({ name, values }) {
  const m = model;
  const G = m.genes.length;
  const z = new Float32Array(G);
  let found = 0;
  for (let g = 0; g < G; g++) {
    const v = values.get(m.upper[g]);
    if (v !== undefined) found++;
    z[g] = ((v ?? m.median[g]) - m.mean[g]) / m.scale[g];
  }
  if (found < 10) throw new Error(`Only ${found} genes in ${name} match the dataset. Check the gene names in the first column.`);

  const pc = new Float64Array(m.k);
  for (let c = 0; c < m.k; c++) {
    let sum = 0;
    for (let g = 0, row = c * G; g < G; g++) sum += (z[g] - m.pcaMean[g]) * m.components[row + g];
    pc[c] = sum;
  }
  const dist = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    let sum = 0;
    for (let c = 0; c < m.k; c++) sum += (pc[c] - m.trainPca[i * m.k + c]) ** 2;
    dist[i] = Math.sqrt(sum);
  }
  const nearest = [...dist.keys()].sort((a, b) => dist[a] - dist[b]).slice(0, m.votes);

  const h = U.count++;
  let id = name;
  for (let n = 2; U.ids.includes(id); n++) id = `${name} (${n})`;
  U.ids.push(id);
  U.voteIdx.push(nearest);
  U.voteDist.push(nearest.map((i) => dist[i]));
  U.imported.push({ z, missing: G - found });

  const tsne = [0, 0, 0];
  for (const i of nearest.slice(0, m.mapNeighbours)) {
    for (let a = 0; a < 3; a++) tsne[a] += LAYOUTS.train.tsne[i * 3 + a] / m.mapNeighbours;
  }
  LAYOUTS.unseen.tsne.set(tsne, h * 3);
  LAYOUTS.unseen.pca.set([0, 1, 2].map((a) => (pc[a] - m.center[a]) / m.spread), h * 3);
  layoutPos("unseen", h).toArray(unseen.pos, h * 3);
  $("unseen-pick").insertAdjacentHTML("beforeend", `<option value="${h}">${id} (imported)</option>`);
  return h;
}

function importStatus(text, isError = false) {
  const el = $("import-status");
  el.textContent = text;
  el.classList.toggle("error", isError);
}

$("import").addEventListener("click", () => $("import-file").click());
$("import-file").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  e.target.value = "";
  if (!file) return;
  try {
    importStatus(`Reading ${file.name}...`);
    const [patients] = await Promise.all([file.text().then((t) => parsePatients(t, file.name)), loadModel()]);
    if (U.count + patients.length > CAP) throw new Error(`Up to ${MAX_IMPORTS} imported patients at a time. Reload the page to clear them.`);
    const added = patients.map(addImported);
    importStatus(`Imported ${added.length} patient${added.length > 1 ? "s" : ""} from ${file.name}.`);
    if (added.length === 1) {
      $("unseen-pick").value = added[0];
      sendIn(added[0]);
    } else {
      select(null);
      flyTo(new THREE.Vector3(), 190);
      added.forEach((h, m) => setTimeout(() => sendIn(h, { follow: false }), m * 45));
    }
    renderUnseen();
  } catch (err) {
    importStatus(err.message, true);
  }
});

// Start

$("sub").textContent = `${N} TCGA patients grouped from ${data.meta.genes.toLocaleString()} genes, reduced with PCA (${data.meta.n_pca} components) then t-SNE. ${H} more were kept out to test the model.`;
$("unseen-intro").textContent = `These ${H} patients were left out when the groups were made. Send one in and watch where it lands, or import your own patient file.`;

function resize() {
  const w = innerWidth;
  const h = innerHeight;
  renderer.setSize(w, h);
  composer.setSize(w, h);
  bloom.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  material.uniforms.uScale.value = (h * renderer.getPixelRatio()) / (2 * Math.tan(THREE.MathUtils.degToRad(FOV / 2)));
}
addEventListener("resize", resize);
resize();

writePositions();
setConfig();

function frame(now) {
  for (const t of [...tweens]) {
    const u = Math.min(1, Math.max(0, (now - t.start) / t.ms));
    t.fn(u);
    if (u >= 1) {
      tweens.delete(t);
      t.done?.();
    }
  }
  controls.update();
  updateLines();
  updateLabels();
  composer.render();
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
