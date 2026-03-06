/* ── app.js — Interactive Mantle Flow Web App (Phase 1) ───────────────────
 *
 * Displays a D3 Natural-Earth-projection globe, sends a POST to /compute
 * with the selected plate-model parameters, and renders the returned
 * pressure field as a colour overlay.
 * ─────────────────────────────────────────────────────────────────────── */

// ── Utilities ──────────────────────────────────────────────────────────────

/** Convert 0-360 longitude to -180 to 180 (required by D3 / GeoJSON). */
function normLon(lon) {
  return lon > 180 ? lon - 360 : lon;
}

// ── Euler-pole velocity (replicates Python euler_pole.py) ──────────────────

const _EP_R = 6371.0;  // Earth radius km (= mm/yr when rate in deg/Myr)

function _sph2cart(lat, lon, r) {
  const la = lat * Math.PI / 180, lo = lon * Math.PI / 180;
  return [r * Math.cos(la) * Math.cos(lo),
          r * Math.cos(la) * Math.sin(lo),
          r * Math.sin(la)];
}
function _cross3(a, b) {
  return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}
function _localEN(lat, lon, vx, vy, vz) {
  const la = lat * Math.PI / 180, lo = lon * Math.PI / 180;
  const east  = -vx * Math.sin(lo) + vy * Math.cos(lo);
  const north = -vx * Math.sin(la) * Math.cos(lo)
                -vy * Math.sin(la) * Math.sin(lo)
                +vz * Math.cos(la);
  return [east, north];
}
/** Returns [vew, vns] in mm/yr at (ptLat, ptLon) for given Euler pole. */
function eulerVelocity(poleLat, poleLon, rateDegMyr, ptLat, ptLon) {
  const omega = _sph2cart(poleLat, poleLon, -rateDegMyr * Math.PI / 180);
  const r     = _sph2cart(ptLat,   ptLon,   _EP_R);
  const v     = _cross3(omega, r);
  return _localEN(ptLat, ptLon, v[0], v[1], v[2]);
}

/**
 * Build a GeoJSON LineString from two lon/lat pairs (0-360 convention).
 */
function lineFeature(lona, lata, lonb, latb, meta) {
  return {
    type: 'Feature',
    geometry: {
      type: 'LineString',
      coordinates: [[normLon(lona), lata], [normLon(lonb), latb]],
    },
    properties: meta,
  };
}

// ── Colour scales ───────────────────────────────────────────────────────────

/**
 * Map a pressure value (MPa) to a CSS colour string using RdBu diverging scale.
 *   positive (compression) → red,  zero → white,  negative (suction) → blue
 * @param {number} p     pressure in MPa
 * @param {number} vmax  scale limit (MPa), typically 95th-percentile of |p|
 */
function pressureToColor(p, vmax) {
  if (!isFinite(p)) return 'rgba(0,0,0,0)';
  // t=0 → red, t=0.5 → white, t=1 → blue  (matches d3.interpolateRdBu)
  const t = 0.5 - 0.5 * Math.max(-1, Math.min(1, p / vmax));
  return d3.interpolateRdBu(t);
}

// ── Map setup ───────────────────────────────────────────────────────────────

const mapWrap = document.getElementById('map-wrap');

/** Compute SVG size from the container. */
function svgSize() {
  const rect = mapWrap.getBoundingClientRect();
  return { w: rect.width || 900, h: rect.height || 520 };
}

const svg = d3.select('#map-svg');

// Defs: arrowhead markers
const defs0 = svg.append('defs');

// Asthenospheric flow arrows (dark grey)
defs0.append('marker')
  .attr('id', 'arrowhead')
  .attr('viewBox', '0 -4 8 8')
  .attr('refX', 7).attr('refY', 0)
  .attr('markerWidth', 3).attr('markerHeight', 3)
  .attr('orient', 'auto')
  .append('path').attr('d', 'M0,-4L8,0L0,4').attr('fill', '#555');

// Plate velocity arrows (arrowhead inherits parent line stroke via context-stroke)
defs0.append('marker')
  .attr('id', 'plate-arrowhead')
  .attr('viewBox', '0 -4 8 8')
  .attr('refX', 7).attr('refY', 0)
  .attr('markerWidth', 3).attr('markerHeight', 3)
  .attr('orient', 'auto')
  .append('path').attr('d', 'M0,-4L8,0L0,4').attr('fill', 'context-stroke');

// Trench velocity arrows (orange)
defs0.append('marker')
  .attr('id', 'trench-arrowhead')
  .attr('viewBox', '0 -4 8 8')
  .attr('refX', 7).attr('refY', 0)
  .attr('markerWidth', 3).attr('markerHeight', 3)
  .attr('orient', 'auto')
  .append('path').attr('d', 'M0,-4L8,0L0,4').attr('fill', '#00e8ff');

let proj, path;

function buildProjection() {
  const { w, h } = svgSize();
  proj = d3.geoNaturalEarth1()
    .rotate([-180, 0])   // Pacific-centred: 180° at centre
    .scale(w / 6.4)
    .translate([w / 2, h / 2]);
  path = d3.geoPath().projection(proj);
  // Update clip-path if it already exists (on resize).
  svg.select('#sphere-clip path').attr('d', path);
}

buildProjection();

// Offscreen canvas used as a pixel buffer for pressure rendering.
const pressureCanvas = document.getElementById('pressure-canvas');

// SVG defs: clip-path to keep pressure image inside the sphere outline.
// path is already set by buildProjection() above, so we can set d immediately.
const defs = svg.append('defs');
defs.append('clipPath').attr('id', 'sphere-clip')
  .append('path').datum({ type: 'Sphere' }).attr('d', path);

// Layer groups — order determines z-index (first = bottom).
const gSphere    = svg.append('g').attr('id', 'g-sphere');
const gGraticule = svg.append('g').attr('id', 'g-graticule');
const gLand      = svg.append('g').attr('id', 'g-land');
// Pressure image sits between land and boundaries, clipped to the sphere.
const gPressure  = svg.append('g').attr('id', 'g-pressure')
                      .attr('clip-path', 'url(#sphere-clip)');
const gBounds    = svg.append('g').attr('id', 'g-bounds');
const gPlateVel  = svg.append('g').attr('id', 'g-plate-vel');
const gTrenchVel = svg.append('g').attr('id', 'g-trench-vel');
const gVelocity  = svg.append('g').attr('id', 'g-velocity');
const gEulerPole = svg.append('g').attr('id', 'g-euler-pole');

// Sphere outline
gSphere.append('path')
  .datum({ type: 'Sphere' })
  .attr('class', 'sphere')
  .attr('d', path);

// Graticule
gGraticule.append('path')
  .datum(d3.geoGraticule()())
  .attr('class', 'graticule')
  .attr('d', path);

// ── Load Natural Earth coastlines ──────────────────────────────────────────

const WORLD_URL = 'https://cdn.jsdelivr.net/npm/world-atlas@2/countries-110m.json';

d3.json(WORLD_URL).then(world => {
  gLand.selectAll('path')
    .data(topojson.feature(world, world.objects.land).features)
    .join('path')
      .attr('class', 'land')
      .attr('d', path);
}).catch(() => {
  // Silently ignore if offline — the map still works without coastlines.
});

// ── Resize handler ─────────────────────────────────────────────────────────

window.addEventListener('resize', () => {
  buildProjection();
  redrawStatic();
});

function redrawStatic() {
  path = d3.geoPath().projection(proj);
  gSphere.select('path').attr('d', path);
  gGraticule.select('path').attr('d', path);
  gLand.selectAll('path').attr('d', path);
  if (lastResult) renderResult(lastResult);
  else if (_lastBounds) renderBoundaries(_lastBounds);
  // Re-project geometry velocity arrows on resize.
  _reprojectArrows(gPlateVel,  '.plate-vel-arrow');
  _reprojectArrows(gTrenchVel, '.trench-vel-arrow');
  if (_lastEulerPole) updateEulerPoleMarker(_lastEulerPole.lat, _lastEulerPole.lon);
}

/** Re-project stored arrow data after a projection change (resize). */
function _reprojectArrows(layer, cls) {
  layer.selectAll(cls).each(function(d) {
    const xy0 = proj([normLon(d.lon), d.lat]);
    if (!xy0) return;
    const [x0, y0] = xy0;
    d3.select(this).attr('x1', x0).attr('y1', y0)
      .attr('x2', x0 + d._dx).attr('y2', y0 + d._dy);
  });
}

// ── Geometry (boundaries shown before computation) ─────────────────────────

async function loadGeometry(model) {
  try {
    let url = `/geometry?model=${encodeURIComponent(model)}`;
    if (model === 'Simple') {
      url += `&width=${document.getElementById('inp-width').value}`;
      url += `&length=${document.getElementById('inp-length').value}`;
    }
    const resp = await fetch(url);
    const data = await resp.json();

    _isEarth      = data.is_earth      || false;
    _spDomainIdx  = data.sp_domain_idx || null;
    _plateVelGrid = data.plate_velocities  || [];
    _trenchVelGrid = data.trench_velocities || [];   // keep vew/vns for re-render on scale change

    // Show / hide velocity controls and plate-dimension controls
    document.getElementById('vel-idealized').style.display  = _isEarth ? 'none'  : 'block';
    document.getElementById('vel-earth-note').style.display = _isEarth ? 'block' : 'none';
    document.getElementById('earth-flux').style.display     = _isEarth ? 'block' : 'none';
    document.getElementById('simple-dims').style.display    = (model === 'Simple') ? 'block' : 'none';

    if (!_isEarth) {
      if (data.sp_euler_pole) {
        document.getElementById('inp-euler-lat') .value = data.sp_euler_pole.lat .toFixed(1);
        document.getElementById('inp-euler-lon') .value = data.sp_euler_pole.lon .toFixed(1);
        document.getElementById('inp-euler-rate').value = data.sp_euler_pole.rate.toFixed(2);
        _lastEulerPole = { lat: data.sp_euler_pole.lat, lon: data.sp_euler_pole.lon };
        updateEulerPoleMarker(data.sp_euler_pole.lat, data.sp_euler_pole.lon);
      }
      if (data.trench_vel_default) {
        document.getElementById('inp-trench-vew').value = data.trench_vel_default.vew.toFixed(1);
        document.getElementById('inp-trench-vns').value = data.trench_vel_default.vns.toFixed(1);
      }
    } else {
      _lastEulerPole = null;
      updateEulerPoleMarker(null, null);
    }

    gPlateVel.style('display', null);
    gTrenchVel.style('display', null);
    if (data.boundaries)       renderBoundaries(data.boundaries);
    if (data.plate_velocities) renderPlateVelocities(data.plate_velocities);
    if (data.trench_velocities) renderTrenchVelocities(data.trench_velocities);
  } catch (_) { /* ignore */ }
}

// ── Controls ───────────────────────────────────────────────────────────────

const btnRun  = document.getElementById('btn-run');
const btnStop = document.getElementById('btn-stop');
const status  = document.getElementById('status');
const overlay = document.getElementById('loading-overlay');
const progBar = document.getElementById('progress-bar');
const progSub = document.getElementById('loading-sub');

btnStop.addEventListener('click', async () => {
  btnStop.disabled = true;
  btnStop.textContent = 'Stopping…';
  try {
    await fetch('/compute/cancel', { method: 'POST' });
  } catch (_) { /* ignore */ }
});

// Live slider labels
const inpDepth = document.getElementById('inp-depth');
inpDepth.addEventListener('input', () => {
  document.getElementById('depth-val').textContent = inpDepth.value;
});

function setStatus(msg, cls) {
  status.innerHTML = msg;
  status.className = cls || '';
}

// Animate the progress bar while waiting.
// Eases toward 90 % over ~estimatedMs, never quite reaching 100 until done.
let _progTimer = null;
function startProgress(estimatedMs) {
  progBar.style.width = '0%';
  progBar.style.transition = 'none';
  const start = Date.now();
  clearInterval(_progTimer);
  _progTimer = setInterval(() => {
    const elapsed = Date.now() - start;
    // Asymptotic curve: reaches ~90 % at estimatedMs
    const pct = 90 * (1 - Math.exp(-3 * elapsed / estimatedMs));
    progBar.style.transition = 'width 0.4s ease';
    progBar.style.width = pct.toFixed(1) + '%';
    const remaining = Math.max(0, (estimatedMs - elapsed) / 1000);
    progSub.textContent = remaining > 1
      ? `~${remaining.toFixed(0)} s remaining`
      : 'Finalising…';
  }, 400);
}
function finishProgress() {
  clearInterval(_progTimer);
  progBar.style.transition = 'width 0.2s ease';
  progBar.style.width = '100%';
}

// Draw the plate-velocity colour bar once on load.
(function() {
  const cvs = document.getElementById('plate-vel-cbar');
  if (!cvs) return;
  const ctx = cvs.getContext('2d');
  for (let px = 0; px < 180; px++) {
    ctx.fillStyle = d3.interpolateViridis(px / 179);
    ctx.fillRect(px, 0, 1, 8);
  }
})();

// ── Live velocity state ─────────────────────────────────────────────────────

let _plateVelGrid  = [];   // [{lon, lat, domain, vew, vns}] from last geometry load
let _trenchVelGrid = [];   // [{lon, lat}] positions only
let _lastPlateVels = [];   // last rendered plate velocities (for CSV export)
let _spDomainIdx   = null; // 1-indexed domain of SP
let _isEarth       = false;
let _lastEulerPole = null; // {lat, lon} for resize redraw

/** User-set velocity scale: arrows representing this speed (mm/yr) are drawn at max length. */
function getVelScale()   { return parseFloat(document.getElementById('inp-vel-scale').value)   ||  40; }
/** User-set pressure scale: colour range is ±this value (MPa). */
function getPressScale() { return parseFloat(document.getElementById('inp-press-scale').value) ||  25; }

// ── Live velocity arrow updates ─────────────────────────────────────────────

function recomputePlateVelArrows() {
  if (_isEarth || !_plateVelGrid.length || _spDomainIdx === null) return;
  const poleLat  = parseFloat(document.getElementById('inp-euler-lat') .value);
  const poleLon  = parseFloat(document.getElementById('inp-euler-lon') .value);
  const poleRate = parseFloat(document.getElementById('inp-euler-rate').value);
  if (!isFinite(poleLat) || !isFinite(poleLon) || !isFinite(poleRate)) return;
  _lastEulerPole = { lat: poleLat, lon: poleLon };
  updateEulerPoleMarker(poleLat, poleLon);
  const updated = _plateVelGrid.map(pt => {
    if (pt.domain === _spDomainIdx) {
      const [vew, vns] = eulerVelocity(poleLat, poleLon, poleRate, pt.lat, pt.lon);
      return { ...pt, vew, vns };
    }
    return pt;
  });
  renderPlateVelocities(updated);
}

function recomputeTrenchVelArrows() {
  if (_isEarth || !_trenchVelGrid.length) return;
  const vew = parseFloat(document.getElementById('inp-trench-vew').value);
  const vns = parseFloat(document.getElementById('inp-trench-vns').value);
  if (!isFinite(vew) || !isFinite(vns)) return;
  renderTrenchVelocities(_trenchVelGrid.map(pt => ({ ...pt, vew, vns })));
}

['inp-euler-lat', 'inp-euler-lon', 'inp-euler-rate'].forEach(id => {
  document.getElementById(id).addEventListener('input', recomputePlateVelArrows);
});
['inp-trench-vew', 'inp-trench-vns'].forEach(id => {
  document.getElementById(id).addEventListener('input', recomputeTrenchVelArrows);
});

// ── Load geometry immediately on page load, and whenever the model changes. ─

const selModel = document.getElementById('sel-model');
loadGeometry(selModel.value);
selModel.addEventListener('change', () => {
  gPressure.selectAll('image').remove();
  gVelocity.selectAll('.vel-arrow').remove();
  gPlateVel.selectAll('.plate-vel-arrow').remove();
  gTrenchVel.selectAll('.trench-vel-arrow').remove();
  gEulerPole.selectAll('.euler-pole-marker').remove();
  _plateVelGrid = []; _trenchVelGrid = []; _spDomainIdx = null; _lastEulerPole = null;
  _lastPlateVels = [];
  document.getElementById('colorbar').style.display = 'none';
  document.getElementById('simple-dims').style.display = (selModel.value === 'Simple') ? 'block' : 'none';
  document.getElementById('plate-vel-bar').style.display = 'none';
  document.getElementById('plate-vel-arrow-key').style.display = 'none';
  const _lkl = document.getElementById('lk-plate-vel-line');
  if (_lkl) _lkl.setAttribute('stroke', '#c0c8e0');
  document.getElementById('btn-pv-png').disabled = true;
  document.getElementById('btn-pv-csv').disabled = true;
  document.getElementById('depth-control').style.display = 'none';
  inpDepth.value = 330; document.getElementById('depth-val').textContent = '330';
  lastResult = null;
  loadGeometry(selModel.value);
});

// ── Run computation ────────────────────────────────────────────────────────

let lastResult  = null;
let _lastBounds = null;

// Per-model timing estimates (seconds) split into solve vs grid phases.
const MODEL_EST = {
  'Simple':                         { solve:  2, grid:  2 },
  'Slab2.0Final_NoJapTail_nnr_FS':  { solve:  4, grid:  5 },
};

const loadingLabel = document.getElementById('loading-label');

btnRun.addEventListener('click', async () => {
  const model = selModel.value;
  const fullPayload = {
    model,
    viscosity:        parseFloat(document.getElementById('inp-visc').value),
    basal_bc:         document.getElementById('sel-bc').value,
    press_depth:      parseFloat(inpDepth.value) * 1000,
    grid_spacing_deg: parseFloat(document.getElementById('sel-res').value),
  };
  if (_isEarth) {
    const fluxMode = document.getElementById('sel-slab-flux').value;
    if (fluxMode === 'both' || fluxMode === 'penetrating') {
      fullPayload.flux_slab  = 2;
      fullPayload.flux_alpha = 0.5;
      fullPayload.flux_width = parseFloat(document.getElementById('inp-flux-width').value) * 1000;
      if (fluxMode === 'penetrating') fullPayload.flux_penetrating = true;
    }
    // else defaults: flux_slab=0 (no flux)
  } else {
    fullPayload.euler_lat  = parseFloat(document.getElementById('inp-euler-lat') .value);
    fullPayload.euler_lon  = parseFloat(document.getElementById('inp-euler-lon') .value);
    fullPayload.euler_rate = parseFloat(document.getElementById('inp-euler-rate').value);
    fullPayload.trench_vew = parseFloat(document.getElementById('inp-trench-vew').value);
    fullPayload.trench_vns = parseFloat(document.getElementById('inp-trench-vns').value);
    if (model === 'Simple') {
      fullPayload.width_deg  = parseFloat(document.getElementById('inp-width') .value);
      fullPayload.length_deg = parseFloat(document.getElementById('inp-length').value);
    }
  }
  // Grid-only payload: the fields that grid_only() reads.
  const gridPayload = {
    viscosity:        fullPayload.viscosity,
    press_depth:      fullPayload.press_depth,
    grid_spacing_deg: fullPayload.grid_spacing_deg,
  };

  const est = MODEL_EST[model] || { solve: 10, grid: 10 };

  btnRun.disabled = true;
  btnRun.textContent = '⏳  Computing…';
  btnStop.disabled = false;
  btnStop.textContent = '■ Stop';
  overlay.style.display = 'flex';
  setStatus('');

  const post = (url, body) => fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });

  try {
    // ── Phase 1: BEM matrix solve ──────────────────────────────────────────
    loadingLabel.textContent = 'Running BEM solver…';
    startProgress(est.solve * 1000);

    const solveResp = await post('/compute/solve', fullPayload);
    const solveData = await solveResp.json();
    if (solveData.cancelled) { gPlateVel.style('display', null); gTrenchVel.style('display', null); setStatus('Cancelled.'); return; }
    if (!solveResp.ok || solveData.error)
      throw new Error(solveData.error || `HTTP ${solveResp.status}`);

    // ── Phase 2: grid output ───────────────────────────────────────────────
    loadingLabel.textContent = 'Outputting grid…';
    startProgress(est.grid * 1000);

    const gridResp = await post('/compute/grid', gridPayload);
    const data = await gridResp.json();
    if (data.cancelled) { gPlateVel.style('display', null); gTrenchVel.style('display', null); setStatus('Cancelled.'); return; }
    if (!gridResp.ok || data.error)
      throw new Error(data.error || `HTTP ${gridResp.status}`);

    finishProgress();
    await new Promise(r => setTimeout(r, 200));

    lastResult = data;
    renderResult(data);
    const cacheNote = solveData.cached ? ' (cached)' : '';
    setStatus(`Done in ${(solveData.timing_s + data.timing_s).toFixed(1)} s${cacheNote}`, 'ok');
    document.getElementById('btn-png').disabled = false;
    document.getElementById('btn-csv').disabled = false;
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'error');
  } finally {
    overlay.style.display = 'none';
    progBar.style.width = '0%';
    loadingLabel.textContent = 'Running BEM solver…';
    btnRun.disabled = false;
    btnRun.textContent = '▶  Run';
    btnStop.disabled = false;
    btnStop.textContent = '■ Stop';
  }
});

// ── Geometry velocity renderers ─────────────────────────────────────────────

/**
 * Plate velocity arrows at boundary segment midpoints, coloured by speed.
 * Uses d3.interpolatePlasma: dark=slow, bright=fast.
 * @param {Array} vels  [{lon, lat, vew, vns, iwall}, …]  (mm/yr)
 */
function renderPlateVelocities(vels) {
  if (!vels || !vels.length) return;
  gPlateVel.selectAll('.plate-vel-arrow').remove();

  const vmax = getVelScale();   // fixed user-set scale (mm/yr)
  const { w, h } = svgSize();
  // Scale: arrow representing vmax mm/yr ≈ 3 % of the shorter map dimension
  const arrowScale = Math.min(w, h) * 0.03 / vmax;
  const colorScale = s => d3.interpolateViridis(Math.min(1, s / vmax));

  const arrows = vels.map(v => {
    const spd = Math.hypot(v.vew, v.vns);
    if (!isFinite(spd) || spd === 0) return null;
    const xy = proj([normLon(v.lon), v.lat]);
    if (!xy) return null;
    const dx =  v.vew * arrowScale;
    const dy = -v.vns * arrowScale;
    return { lon: v.lon, lat: v.lat, _dx: dx, _dy: dy, spd };
  }).filter(Boolean);

  // Store for CSV export (use the input vels, which have vew/vns)
  _lastPlateVels = vels.filter(v => isFinite(v.vew) && isFinite(v.vns));
  document.getElementById('btn-pv-png').disabled = false;
  document.getElementById('btn-pv-csv').disabled = false;

  // Update plate-vel bar and legend arrow key
  const pvBar = document.getElementById('plate-vel-bar');
  if (pvBar) {
    pvBar.style.display = 'block';
    document.getElementById('plate-vel-vmax').textContent = `${vmax}`;
  }
  const pvKey = document.getElementById('plate-vel-arrow-key');
  if (pvKey) {
    pvKey.style.display = 'flex';
    document.getElementById('plate-vel-arrow-label').textContent = `${vmax} mm/yr`;
  }

  gPlateVel.selectAll('.plate-vel-arrow')
    .data(arrows)
    .join('line')
      .attr('class', 'plate-vel-arrow')
      .attr('x1', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[0] : 0; })
      .attr('y1', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[1] : 0; })
      .attr('x2', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[0] + d._dx : 0; })
      .attr('y2', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[1] + d._dy : 0; })
      .attr('stroke', d => colorScale(d.spd))
      .attr('marker-end', 'url(#plate-arrowhead)');
}

/**
 * Trench / slab prescribed velocity arrows (orange).
 * @param {Array} vels  [{lon, lat, vew, vns}, …]  (mm/yr)
 */
function renderTrenchVelocities(vels) {
  if (!vels || !vels.length) return;
  gTrenchVel.selectAll('.trench-vel-arrow').remove();

  const refVmax = getVelScale();   // fixed user-set scale (mm/yr)
  const { w, h } = svgSize();
  const arrowScale = Math.min(w, h) * 0.03 / refVmax;

  const arrows = vels.map(v => {
    const spd = Math.hypot(v.vew, v.vns);
    if (!isFinite(spd) || spd === 0) return null;
    const xy = proj([normLon(v.lon), v.lat]);
    if (!xy) return null;
    const dx =  v.vew * arrowScale;
    const dy = -v.vns * arrowScale;
    return { lon: v.lon, lat: v.lat, _dx: dx, _dy: dy };
  }).filter(Boolean);

  gTrenchVel.selectAll('.trench-vel-arrow')
    .data(arrows)
    .join('line')
      .attr('class', 'trench-vel-arrow')
      .attr('x1', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[0] : 0; })
      .attr('y1', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[1] : 0; })
      .attr('x2', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[0] + d._dx : 0; })
      .attr('y2', d => { const xy = proj([normLon(d.lon), d.lat]); return xy ? xy[1] + d._dy : 0; })
      .attr('marker-end', 'url(#trench-arrowhead)');
}

// ── Render ─────────────────────────────────────────────────────────────────

function renderResult(data) {
  // Hide geometry velocity arrows — too cluttered with the pressure field.
  gPlateVel.style('display', 'none');
  gTrenchVel.style('display', 'none');
  renderPressure(data);
  renderBoundaries(data.boundaries);
  renderVelocity(data.velocity);
  renderColorbar();
  updateVelScaleLegend();
  // Switch legend arrow to mantle-flow colour.
  const line = document.getElementById('lk-plate-vel-line');
  if (line) line.setAttribute('stroke', '#555');
  document.getElementById('depth-control').style.display = 'block';
}

/**
 * Revert the map to the geometry/velocity view after a pressure result has been shown.
 * No-ops if no result is currently displayed.
 */
function revertToGeometry() {
  if (!lastResult) return;
  lastResult = null;
  gPressure.selectAll('image').remove();
  gVelocity.selectAll('.vel-arrow').remove();
  gPlateVel.style('display', null);
  gTrenchVel.style('display', null);
  document.getElementById('colorbar').style.display = 'none';
  document.getElementById('btn-png').disabled = true;
  document.getElementById('btn-csv').disabled = true;
  // Restore legend arrow to plate-velocity colour.
  const line = document.getElementById('lk-plate-vel-line');
  if (line) line.setAttribute('stroke', '#c0c8e0');
  const dc = document.getElementById('depth-control');
  dc.style.display = 'none';
  inpDepth.value = 330;
  document.getElementById('depth-val').textContent = '330';
  setStatus('');
}

/**
 * Render the pressure field onto the canvas using bilinear interpolation.
 * Samples the model grid at SAMPLE_DEG intervals for smooth rendering.
 */
function renderPressure(data) {
  const { lons, lats, pressure } = data;
  const { w, h } = svgSize();

  const vmax = getPressScale();   // user-set symmetric scale (MPa)

  pressureCanvas.width  = w;
  pressureCanvas.height = h;
  const ctx = pressureCanvas.getContext('2d');
  ctx.clearRect(0, 0, w, h);
  ctx.globalAlpha = 0.85;

  const nlon = lons.length, nlat = lats.length;
  const dlon = lons[1] - lons[0];   // grid spacing in lon (positive)
  const dlat = lats[1] - lats[0];   // grid spacing in lat (may be negative)

  const latMin = Math.min(lats[0], lats[nlat - 1]);
  const latMax = Math.max(lats[0], lats[nlat - 1]);

  // Sample at SAMPLE_DEG steps with bilinear interpolation between model grid nodes.
  const SAMPLE_DEG = 1.5;
  // Size of each drawn square in pixels (slightly oversized to avoid gaps).
  const halfEdge = Math.ceil((SAMPLE_DEG / 360) * w) + 2;

  for (let slon = lons[0]; slon <= lons[nlon - 1]; slon += SAMPLE_DEG) {
    const glon = normLon(slon);
    if (Math.abs(glon) < 5) continue;  // skip near the map seam

    // Fractional column index
    const jf = (slon - lons[0]) / dlon;
    const j0 = Math.max(0, Math.min(nlon - 2, Math.floor(jf)));
    const tj = jf - j0;

    for (let slat = latMin; slat <= latMax; slat += SAMPLE_DEG) {
      // Fractional row index (handles both +/- dlat)
      const if_ = (slat - lats[0]) / dlat;
      const i0 = Math.max(0, Math.min(nlat - 2, Math.floor(if_)));
      const ti = if_ - i0;
      const i1 = i0 + 1, j1 = j0 + 1;

      const p00 = pressure[i0][j0], p01 = pressure[i0][j1];
      const p10 = pressure[i1][j0], p11 = pressure[i1][j1];
      if (!isFinite(p00) || !isFinite(p01) || !isFinite(p10) || !isFinite(p11)) continue;

      // Bilinear interpolation
      const p = (1 - ti) * ((1 - tj) * p00 + tj * p01)
              +      ti  * ((1 - tj) * p10 + tj * p11);

      const xy = proj([glon, slat]);
      if (!xy) continue;

      ctx.fillStyle = pressureToColor(p, vmax);
      ctx.fillRect(xy[0] - halfEdge, xy[1] - halfEdge, halfEdge * 2, halfEdge * 2);
    }
  }

  // Inject canvas pixels as an SVG <image> (above land, below plate boundaries).
  const dataUrl = pressureCanvas.toDataURL();
  let img = gPressure.select('image');
  if (img.empty()) img = gPressure.append('image');
  img.attr('href', dataUrl).attr('x', 0).attr('y', 0).attr('width', w).attr('height', h);
}

/**
 * Draw original plate boundaries and subduction zone teeth.
 *   type 0 → plate edge (blue)
 *   type 1 → slab wall / trench (red) + teeth triangles
 *   type 2 → strike-slip (green)
 */
function trenchColor(b) {
  const penetratingMode = document.getElementById('sel-slab-flux').value === 'penetrating';
  if (penetratingMode && b.penetrating) return '#ffffff';
  return '#ff6040';
}

function renderBoundaries(bounds) {
  if (!bounds || !bounds.length) return;
  _lastBounds = bounds;

  // Sort so penetrating type-1 segments render last (on top in SVG).
  const sorted = [...bounds].sort((a, b) => {
    const aTop = a.type === 1 && a.penetrating &&
                 document.getElementById('sel-slab-flux').value === 'penetrating';
    const bTop = b.type === 1 && b.penetrating &&
                 document.getElementById('sel-slab-flux').value === 'penetrating';
    return aTop - bTop;
  });
  const features = sorted.map(b =>
    lineFeature(b.lona, b.lata, b.lonb, b.latb, { btype: b.type, penetrating: b.penetrating })
  );

  gBounds.selectAll('.boundary').remove();
  gBounds.selectAll('.boundary')
    .data(features)
    .join('path')
      .attr('class', d => `boundary boundary-${d.properties.btype}`)
      .style('stroke', d => d.properties.btype === 1 ? trenchColor(d.properties) : null)
      .attr('d', path);

  // ── Subduction teeth on type-1 segments ───────────────────────────────────
  gBounds.selectAll('.trench-tooth').remove();

  const TOOTH_SPACING = 22;  // target px between tooth centres
  const TOOTH_H = 7;         // height perpendicular to segment
  const TOOTH_W = 5;         // half-base width along segment
  const MIN_LEN = 8;         // skip segments shorter than this (px)

  const penetratingMode = document.getElementById('sel-slab-flux').value === 'penetrating';
  const type1 = sorted.filter(b => b.type === 1);
  type1.forEach(b => {
    const pa = proj([normLon(b.lona), b.lata]);
    const pb = proj([normLon(b.lonb), b.latb]);
    if (!pa || !pb) return;

    const dx = pb[0] - pa[0], dy = pb[1] - pa[1];
    const len = Math.hypot(dx, dy);
    if (len < MIN_LEN) return;

    const ux = dx / len, uy = dy / len;
    // polarity > 0 → teeth to the left of a→b direction: (-uy, ux)
    // polarity < 0 → teeth to the right: (uy, -ux)
    const pol = (b.polarity != null && b.polarity < 0) ? 1 : -1;
    const nx = uy * pol, ny = -ux * pol;

    const nTeeth = Math.max(1, Math.round(len / TOOTH_SPACING));
    const fill = trenchColor(b);
    for (let k = 0; k < nTeeth; k++) {
      const t  = (k + 0.5) / nTeeth;
      const mx = pa[0] + dx * t, my = pa[1] + dy * t;
      const apex  = [mx + nx * TOOTH_H,  my + ny * TOOTH_H];
      const base1 = [mx - ux * TOOTH_W,  my - uy * TOOTH_W];
      const base2 = [mx + ux * TOOTH_W,  my + uy * TOOTH_W];
      gBounds.append('polygon')
        .attr('class', 'trench-tooth')
        .style('fill', fill)
        .attr('points',
          `${base1[0]},${base1[1]} ${base2[0]},${base2[1]} ${apex[0]},${apex[1]}`);
    }
  });
}

/**
 * Draw or update the Euler pole marker (×) on the map.
 * Pass null to remove the marker.
 */
function updateEulerPoleMarker(lat, lon) {
  gEulerPole.selectAll('.euler-pole-marker').remove();
  if (lat == null || !isFinite(lat) || !isFinite(lon)) return;
  const xy = proj([normLon(lon), lat]);
  if (!xy) return;
  const s = 7;
  gEulerPole.append('line').attr('class', 'euler-pole-marker')
    .attr('x1', xy[0]-s).attr('y1', xy[1]).attr('x2', xy[0]+s).attr('y2', xy[1]);
  gEulerPole.append('line').attr('class', 'euler-pole-marker')
    .attr('x1', xy[0]).attr('y1', xy[1]-s).attr('x2', xy[0]).attr('y2', xy[1]+s);
}


/**
 * Approximate point-to-segment distance in degrees (cos-lat corrected).
 * Used to exclude velocity arrows that are too close to slab walls.
 */
function _minDistToSlabWalls(lat, lon, slabWalls) {
  const RAD = Math.PI / 180;
  let minDist = Infinity;
  for (const b of slabWalls) {
    const midLat = (lat + b.lata + b.latb) / 3;
    const coslat = Math.cos(midLat * RAD);
    // Normalise longitudes relative to b.lona to handle antimeridian wrap.
    const lonb = b.lona + (((b.lonb - b.lona) + 540) % 360 - 180);
    const lonp = b.lona + (((lon   - b.lona) + 540) % 360 - 180);
    const px = (lonp - b.lona) * coslat, py = lat - b.lata;
    const dx = (lonb - b.lona) * coslat, dy = b.latb - b.lata;
    const lenSq = dx * dx + dy * dy;
    let d;
    if (lenSq < 1e-10) {
      d = Math.hypot(px, py);
    } else {
      const t = Math.max(0, Math.min(1, (px * dx + py * dy) / lenSq));
      d = Math.hypot(px - t * dx, py - t * dy);
    }
    if (d < minDist) minDist = d;
  }
  return minDist;
}

/**
 * Draw velocity arrows at subsampled grid points.
 * Arrows scaled to map units so they are readable.
 * Points within 5° of any slab wall are excluded (vectors are unrealistically large there).
 */
function renderVelocity(velocity) {
  if (!velocity || !velocity.length) return;

  const slabWalls = (_lastBounds || []).filter(b => b.type === 1);

  // velocity data is in cm/yr; scale input is in mm/yr (1 cm/yr = 10 mm/yr).
  const velScaleCmyr = getVelScale() / 10;
  const { w, h } = svgSize();
  const arrowScale = Math.min(w, h) * 0.03 / velScaleCmyr;  // pixels per (cm/yr)

  const arrows = velocity
    .filter(v => isFinite(v.vew) && isFinite(v.vns))
    .filter(v => _minDistToSlabWalls(v.lat, v.lon, slabWalls) >= 5)
    .map(v => {
      const [x0, y0] = proj([normLon(v.lon), v.lat]) || [null, null];
      if (x0 === null) return null;
      const dx =  v.vew * arrowScale;
      const dy = -v.vns * arrowScale;   // screen y is inverted
      return { x0, y0, x1: x0 + dx, y1: y0 + dy };
    })
    .filter(Boolean);

  gVelocity.selectAll('.vel-arrow')
    .data(arrows)
    .join('line')
      .attr('class', 'vel-arrow')
      .attr('x1', d => d.x0).attr('y1', d => d.y0)
      .attr('x2', d => d.x1).attr('y2', d => d.y1);
}

/**
 * Draw the pressure colour-bar using the current user-set scale.
 */
function renderColorbar() {
  const vmax = getPressScale();
  const cbDiv = document.getElementById('colorbar');
  cbDiv.style.display = 'block';

  let canvas = cbDiv.querySelector('canvas');
  if (!canvas) {
    canvas = document.createElement('canvas');
    cbDiv.insertBefore(canvas, cbDiv.firstChild);
  }
  canvas.width = 180;
  canvas.height = 14;
  const ctx = canvas.getContext('2d');
  for (let px = 0; px < 180; px++) {
    const p = ((px / 179) * 2 - 1) * vmax;  // -vmax → +vmax
    ctx.fillStyle = pressureToColor(p, vmax);
    ctx.fillRect(px, 0, 1, 14);
  }

  document.getElementById('cb-min').textContent = `${(-vmax).toFixed(1)}`;
  document.getElementById('cb-max').textContent = `${vmax.toFixed(1)}`;
}

/**
 * Update all velocity legend labels to reflect the current user-set scale.
 */
function updateVelScaleLegend() {
  const s = getVelScale();
  const el1 = document.getElementById('plate-vel-vmax');
  if (el1) el1.textContent = `${s}`;
  const el2 = document.getElementById('plate-vel-arrow-label');
  if (el2) el2.textContent = `${s} mm/yr`;
}

// ── Export ──────────────────────────────────────────────────────────────────

const EXPORT_SCALE = 2;  // output at 2× screen resolution

/**
 * Inject an inline <style> + explicit dimensions into the SVG so that class-
 * based styles survive when the SVG is loaded as a standalone <img> (external
 * stylesheets are not available in that context).  Returns the injected element
 * so the caller can remove it after serialisation.
 */
function _svgExportSetup(svgEl, pw, ph) {
  svgEl.setAttribute('width',  pw);
  svgEl.setAttribute('height', ph);
  const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
  style.textContent =
    '.sphere{fill:#0d1b2a}' +
    '.graticule{fill:none;stroke:rgba(255,255,255,0.08);stroke-width:0.5}' +
    '.land{fill:#2a3a4a;stroke:#4a6070;stroke-width:0.4}' +
    '.boundary{fill:none;stroke-width:1.4}' +
    '.boundary-0{stroke:#60a0ff}' +
    '.boundary-1{stroke:#ff6040;stroke-width:2.2}' +
    '.boundary-2{stroke:#a0c860}' +
    '.vel-arrow{stroke:#555;stroke-width:2.5;fill:none;marker-end:url(#arrowhead)}' +
    '.plate-vel-arrow{stroke-width:2.5;fill:none}' +
    '.trench-vel-arrow{stroke:#00e8ff;stroke-width:2.5;fill:none;marker-end:url(#trench-arrowhead)}' +
    '.trench-tooth{fill:#ff6040;stroke:none}' +
    '.euler-pole-marker{stroke:#ffaa00;stroke-width:2;fill:none}';
  svgEl.insertBefore(style, svgEl.firstChild);
  return style;
}

function _svgExportTeardown(svgEl, style) {
  svgEl.removeChild(style);
  svgEl.removeAttribute('width');
  svgEl.removeAttribute('height');
}

/** Serialise svgEl to a blob URL and resolve with the loaded Image. */
function _svgToImage(svgEl) {
  const data = new XMLSerializer().serializeToString(svgEl);
  const url  = URL.createObjectURL(new Blob([data], { type: 'image/svg+xml;charset=utf-8' }));
  return new Promise(resolve => {
    const img = new Image();
    img.onload = () => { URL.revokeObjectURL(url); resolve(img); };
    img.src = url;
  });
}

function _downloadCanvas(canvas, filename) {
  const a = document.createElement('a');
  a.download = filename;
  a.href = canvas.toDataURL('image/png');
  a.click();
}

/** Download pressure map + mantle flow vectors as PNG at 2× resolution. */
async function exportPNG() {
  if (!lastResult) return;
  const { w, h } = svgSize();
  const pw = w * EXPORT_SCALE, ph = h * EXPORT_SCALE;
  const svgEl = document.getElementById('map-svg');

  // Hide the SVG-embedded pressure <image>; we draw pressureCanvas directly
  // so the browser doesn't need to load a data-URL from within an <img> SVG.
  const pressImg = gPressure.select('image');
  pressImg.style('display', 'none');
  const style = _svgExportSetup(svgEl, pw, ph);
  const img = await _svgToImage(svgEl);
  _svgExportTeardown(svgEl, style);
  pressImg.style('display', null);

  const canvas = document.createElement('canvas');
  canvas.width = pw; canvas.height = ph;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#0d1b2a';
  ctx.fillRect(0, 0, pw, ph);
  ctx.drawImage(pressureCanvas, 0, 0, pw, ph);  // pressure field (upscaled)
  ctx.drawImage(img, 0, 0);                      // land, bounds, flow arrows at 2×
  _downloadCanvas(canvas, 'mantle-flow.png');
}

/** Download plate velocity map as PNG at 2× resolution (geometry view). */
async function exportPlateVelPNG() {
  if (!_lastPlateVels.length) return;
  const { w, h } = svgSize();
  const pw = w * EXPORT_SCALE, ph = h * EXPORT_SCALE;
  const svgEl = document.getElementById('map-svg');
  const inPressureView = lastResult !== null;

  // Switch to geometry view momentarily for the snapshot
  if (inPressureView) {
    gPressure.style('display', 'none');
    gVelocity.style('display', 'none');
    gPlateVel.style('display', null);
    gTrenchVel.style('display', null);
  }
  const style = _svgExportSetup(svgEl, pw, ph);
  const img = await _svgToImage(svgEl);
  _svgExportTeardown(svgEl, style);
  if (inPressureView) {
    gPressure.style('display', null);
    gVelocity.style('display', null);
    gPlateVel.style('display', 'none');
    gTrenchVel.style('display', 'none');
  }

  const canvas = document.createElement('canvas');
  canvas.width = pw; canvas.height = ph;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#0d1b2a';
  ctx.fillRect(0, 0, pw, ph);
  ctx.drawImage(img, 0, 0);
  _downloadCanvas(canvas, 'plate-velocities.png');
}

/** Download pressure grid as CSV (lon, lat, pressure_MPa). */
function exportCSV() {
  if (!lastResult) return;
  const { lons, lats, pressure } = lastResult;
  const rows = ['lon,lat,pressure_MPa'];
  for (let i = 0; i < lats.length; i++) {
    for (let j = 0; j < lons.length; j++) {
      const p = pressure[i][j];
      if (isFinite(p)) rows.push(`${lons[j].toFixed(2)},${lats[i].toFixed(2)},${p.toFixed(4)}`);
    }
  }
  const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
  const a = document.createElement('a');
  a.download = 'mantle-flow.csv';
  a.href = URL.createObjectURL(blob);
  a.click();
  URL.revokeObjectURL(a.href);
}

/** Download plate velocity grid as CSV (lon, lat, vew, vns, speed — mm/yr). */
function exportPlateVelCSV() {
  if (!_lastPlateVels.length) return;
  const rows = ['lon,lat,vew_mmyr,vns_mmyr,speed_mmyr'];
  _lastPlateVels.forEach(v => {
    const spd = Math.hypot(v.vew, v.vns);
    rows.push(`${normLon(v.lon).toFixed(2)},${v.lat.toFixed(2)},${v.vew.toFixed(2)},${v.vns.toFixed(2)},${spd.toFixed(2)}`);
  });
  const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
  const a = document.createElement('a');
  a.download = 'plate-velocities.csv';
  a.href = URL.createObjectURL(blob);
  a.click();
  URL.revokeObjectURL(a.href);
}

document.getElementById('btn-png').addEventListener('click', exportPNG);
document.getElementById('btn-csv').addEventListener('click', exportCSV);
document.getElementById('btn-pv-png').addEventListener('click', exportPlateVelPNG);
document.getElementById('btn-pv-csv').addEventListener('click', exportPlateVelCSV);

// ── Scale input listeners ────────────────────────────────────────────────────

document.getElementById('inp-press-scale').addEventListener('input', () => {
  renderColorbar();
  if (lastResult) renderPressure(lastResult);
});

document.getElementById('inp-vel-scale').addEventListener('input', () => {
  updateVelScaleLegend();
  if (_isEarth) {
    // Earth models: re-render stored geometry data directly
    if (_plateVelGrid.length) renderPlateVelocities(_plateVelGrid);
    if (_trenchVelGrid.length) renderTrenchVelocities(_trenchVelGrid);
  } else {
    recomputePlateVelArrows();
    recomputeTrenchVelArrows();
  }
  if (lastResult) renderVelocity(lastResult.velocity);
});

// ── Viscosity input: order-of-magnitude stepping ────────────────────────────

const viscInput = document.getElementById('inp-visc');

/** Format a viscosity value as clean scientific notation, e.g. "4e+20". */
function _fmtVisc(v) {
  return v.toExponential().replace(/\.?0+(e)/, '$1');
}

/** Step size = largest power of 10 that fits in the current value. */
function _viscStep() {
  const v = parseFloat(viscInput.value);
  return (isFinite(v) && v > 0) ? Math.pow(10, Math.floor(Math.log10(v))) : 1e20;
}

/** Keep the step attribute in sync so the spinner buttons use the right increment. */
function _syncViscStep() { viscInput.step = _viscStep(); }

viscInput.addEventListener('focus', _syncViscStep);
viscInput.addEventListener('input', _syncViscStep);

viscInput.addEventListener('keydown', e => {
  if (e.key !== 'ArrowUp' && e.key !== 'ArrowDown') return;
  e.preventDefault();
  const v = parseFloat(viscInput.value);
  if (!isFinite(v) || v <= 0) return;
  const newV = v + (e.key === 'ArrowUp' ? 1 : -1) * _viscStep();
  if (newV > 0) {
    viscInput.value = _fmtVisc(newV);
    viscInput.dispatchEvent(new Event('input'));
    viscInput.dispatchEvent(new Event('change'));
  }
});

_syncViscStep();  // initialise step on page load

// ── Re-run grid phase using cached BEM solve ────────────────────────────────
// Used when only viscosity or plot depth changes — no need to redo the solve.
async function rerunGrid() {
  if (!lastResult) return;
  const model = selModel.value;
  const est = MODEL_EST[model] || { solve: 10, grid: 10 };
  const gridPayload = {
    viscosity:        parseFloat(document.getElementById('inp-visc').value),
    press_depth:      parseFloat(inpDepth.value) * 1000,
    grid_spacing_deg: parseFloat(document.getElementById('sel-res').value),
  };

  btnRun.disabled = true;
  btnStop.disabled = false;
  btnStop.textContent = '■ Stop';
  overlay.style.display = 'flex';
  loadingLabel.textContent = 'Updating grid…';
  startProgress(est.grid * 1000);
  setStatus('');

  try {
    const resp = await fetch('/compute/grid', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(gridPayload),
    });
    const data = await resp.json();
    if (data.cancelled) { setStatus('Cancelled.'); return; }
    if (!resp.ok || data.error) throw new Error(data.error || `HTTP ${resp.status}`);

    finishProgress();
    await new Promise(r => setTimeout(r, 200));
    lastResult = data;
    renderResult(data);
    setStatus(`Done in ${data.timing_s.toFixed(1)} s`, 'ok');
    document.getElementById('btn-png').disabled = false;
    document.getElementById('btn-csv').disabled = false;
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'error');
  } finally {
    overlay.style.display = 'none';
    progBar.style.width = '0%';
    loadingLabel.textContent = 'Running BEM solver…';
    btnRun.disabled = false;
    btnStop.disabled = false;
  }
}

// ── Plate dimensions (Simple model): changing shape → full geometry reload ───
function _reloadSimpleGeometry() {
  gBounds.selectAll('*').remove();
  _lastBounds = null;
  gPressure.selectAll('image').remove();
  gVelocity.selectAll('.vel-arrow').remove();
  gPlateVel.selectAll('.plate-vel-arrow').remove();
  gTrenchVel.selectAll('.trench-vel-arrow').remove();
  gEulerPole.selectAll('.euler-pole-marker').remove();
  _plateVelGrid = []; _trenchVelGrid = []; _spDomainIdx = null; _lastEulerPole = null;
  _lastPlateVels = [];
  document.getElementById('colorbar').style.display = 'none';
  document.getElementById('plate-vel-bar').style.display = 'none';
  document.getElementById('plate-vel-arrow-key').style.display = 'none';
  const _lkl = document.getElementById('lk-plate-vel-line');
  if (_lkl) _lkl.setAttribute('stroke', '#c0c8e0');
  document.getElementById('btn-pv-png').disabled = true;
  document.getElementById('btn-pv-csv').disabled = true;
  document.getElementById('depth-control').style.display = 'none';
  inpDepth.value = 330; document.getElementById('depth-val').textContent = '330';
  lastResult = null;
  setStatus('');
  loadGeometry(selModel.value);
}
['inp-width', 'inp-length'].forEach(id => {
  document.getElementById(id).addEventListener('change', _reloadSimpleGeometry);
});

// ── Revert to geometry view when model parameters change ────────────────────
// Basal BC and slab flux affect the solve cache key → full revert required.
document.getElementById('sel-bc').addEventListener('change', revertToGeometry);
document.getElementById('sel-slab-flux').addEventListener('change', () => {
  const active = document.getElementById('sel-slab-flux').value !== 'none';
  document.getElementById('inp-flux-width').disabled = !active;
  if (_lastBounds) renderBoundaries(_lastBounds);
  revertToGeometry();
});
document.getElementById('inp-flux-width').addEventListener('change', revertToGeometry);
// Viscosity / depth only affect the grid phase → re-run grid if result exists.
document.getElementById('inp-visc').addEventListener('change', () => {
  if (lastResult) rerunGrid(); else revertToGeometry();
});
inpDepth.addEventListener('change', () => {
  if (lastResult) rerunGrid();
});
// Velocity inputs (live — use 'input' so dragging/typing reverts immediately)
['inp-euler-lat', 'inp-euler-lon', 'inp-euler-rate',
 'inp-trench-vew', 'inp-trench-vns'].forEach(id => {
  document.getElementById(id).addEventListener('input', revertToGeometry);
});

// Initialise legend labels immediately on page load.
updateVelScaleLegend();
