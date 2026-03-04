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

// Defs: arrowhead marker for velocity arrows
svg.append('defs').append('marker')
  .attr('id', 'arrowhead')
  .attr('viewBox', '0 -4 8 8')
  .attr('refX', 7).attr('refY', 0)
  .attr('markerWidth', 4).attr('markerHeight', 4)
  .attr('orient', 'auto')
  .append('path')
    .attr('d', 'M0,-4L8,0L0,4')
    .attr('fill', '#ffee80');

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
const gVelocity  = svg.append('g').attr('id', 'g-velocity');

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
  else {
    // Re-project boundaries even without a full result.
    gBounds.selectAll('.boundary').attr('d', path);
  }
}

// ── Geometry (boundaries shown before computation) ─────────────────────────

async function loadGeometry(model) {
  try {
    const resp = await fetch(`/geometry?model=${encodeURIComponent(model)}`);
    const data = await resp.json();
    if (data.boundaries) renderBoundaries(data.boundaries);
  } catch (_) { /* ignore — map still usable */ }
}

// ── Controls ───────────────────────────────────────────────────────────────

const btnRun = document.getElementById('btn-run');
const status  = document.getElementById('status');
const overlay = document.getElementById('loading-overlay');
const progBar = document.getElementById('progress-bar');
const progSub = document.getElementById('loading-sub');

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

// Load geometry immediately on page load, and whenever the model changes.
const selModel = document.getElementById('sel-model');
loadGeometry(selModel.value);
selModel.addEventListener('change', () => {
  // Clear any existing pressure overlay when switching models.
  gPressure.selectAll('image').remove();
  gVelocity.selectAll('.vel-arrow').remove();
  document.getElementById('colorbar').style.display = 'none';
  lastResult = null;
  loadGeometry(selModel.value);
});

// ── Run computation ────────────────────────────────────────────────────────

let lastResult = null;

// Rough estimate for each model (seconds), used to size the progress animation.
const MODEL_EST = {
  'LargeSP_RetreatingTrench':                    15,
  'SmallSP_RetreatingTrench':                    15,
  'LargeSP_RetreatingTrenchSlabGap':             15,
  'Slab2.0Final_NoJapTail_nnr_FS':               60,
  'Slab2.0Final_NoJapTailNoPhil_nnr_FS':         60,
};

btnRun.addEventListener('click', async () => {
  const model = selModel.value;
  const payload = {
    model,
    viscosity: parseFloat(document.getElementById('inp-visc').value),
    basal_bc:  document.getElementById('sel-bc').value,
  };

  const estMs = (MODEL_EST[model] || 30) * 1000;

  btnRun.disabled = true;
  btnRun.textContent = '⏳  Computing…';
  overlay.style.display = 'flex';
  startProgress(estMs);
  setStatus('');

  try {
    const resp = await fetch('/compute', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });

    const data = await resp.json();

    if (!resp.ok || data.error) {
      throw new Error(data.error || `HTTP ${resp.status}`);
    }

    finishProgress();
    // Brief pause so the bar visibly hits 100 %
    await new Promise(r => setTimeout(r, 300));

    lastResult = data;
    renderResult(data);
    setStatus(`Done in ${data.timing_s.toFixed(1)} s`, 'ok');
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'error');
  } finally {
    overlay.style.display = 'none';
    progBar.style.width = '0%';
    btnRun.disabled = false;
    btnRun.textContent = '▶  Run';
  }
});

// ── Render ─────────────────────────────────────────────────────────────────

function renderResult(data) {
  renderPressure(data);
  renderBoundaries(data.boundaries);
  renderVelocity(data.velocity);
  renderColorbar(data);
}

/**
 * Render the pressure field onto the canvas using bilinear interpolation.
 * Samples the model grid at SAMPLE_DEG intervals for smooth rendering.
 */
function renderPressure(data) {
  const { lons, lats, pressure } = data;
  const { w, h } = svgSize();

  // vmax = 95th percentile of |p| to avoid extreme outliers dominating the scale.
  const absVals = [];
  for (let i = 0; i < lats.length; i++)
    for (let j = 0; j < lons.length; j++)
      if (isFinite(pressure[i][j])) absVals.push(Math.abs(pressure[i][j]));
  absVals.sort((a, b) => a - b);
  const vmax = absVals[Math.floor(absVals.length * 0.95)] || 1;
  lastResult._vmax = vmax;

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
 * Draw original plate boundaries (edges and slab walls).
 *   type 0 → plate edge (blue)
 *   type 1 → slab wall / trench (red)
 *   type 2 → strike-slip (green)
 */
function renderBoundaries(bounds) {
  if (!bounds || !bounds.length) return;

  const features = bounds.map(b =>
    lineFeature(b.lona, b.lata, b.lonb, b.latb, { btype: b.type })
  );

  gBounds.selectAll('.boundary')
    .data(features)
    .join('path')
      .attr('class', d => `boundary boundary-${d.properties.btype}`)
      .attr('d', path);
}

/**
 * Draw velocity arrows at subsampled grid points.
 * Arrows scaled to map units so they are readable.
 */
function renderVelocity(velocity) {
  if (!velocity || !velocity.length) return;

  // Find max speed for normalisation.
  let vmax = 0;
  velocity.forEach(v => {
    const spd = Math.hypot(v.vew, v.vns);
    if (isFinite(spd)) vmax = Math.max(vmax, spd);
  });
  if (vmax === 0) return;

  const { w, h } = svgSize();
  const arrowScale = Math.min(w, h) * 0.04 / vmax;  // pixels per (m/s)

  const arrows = velocity
    .filter(v => isFinite(v.vew) && isFinite(v.vns))
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
 * Draw a horizontal colour-bar legend.
 */
function renderColorbar(data) {
  const vmax = data._vmax || 1;
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

  document.getElementById('cb-min').textContent = `${(-vmax).toFixed(1)} MPa`;
  document.getElementById('cb-max').textContent = `${vmax.toFixed(1)} MPa`;
}
