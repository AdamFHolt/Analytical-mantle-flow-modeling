#!/usr/bin/env python3
"""Flask server for the interactive mantle flow web app."""
import sys
import pathlib
import threading
import multiprocessing as mp

# Add flow_computations/ to the module search path so functions.py and euler_pole.py are importable.
FC_DIR = pathlib.Path(__file__).parent.parent / 'flow_computations'
sys.path.insert(0, str(FC_DIR))

from flask import Flask, request, jsonify, send_from_directory
import compute as _compute_mod

app = Flask(__name__, static_folder='static', static_url_path='')

# ── Process management (for cancellation) ─────────────────────────────────────
_mp_ctx      = mp.get_context('spawn')
_proc_lock   = threading.Lock()
_active_proc = None   # type: mp.Process | None


def _solve_worker(fc_dir, payload, q):
    sys.path.insert(0, fc_dir)
    import compute
    try:
        result = compute.solve_only(payload)
        q.put({'ok': True, 'result': result, 'cache': compute._cache})
    except Exception as exc:
        import traceback
        q.put({'ok': False, 'error': str(exc), 'tb': traceback.format_exc()})


def _grid_worker(fc_dir, payload, cache, q):
    sys.path.insert(0, fc_dir)
    import compute
    compute._cache = cache
    try:
        result = compute.grid_only(payload)
        q.put({'ok': True, 'result': result, 'cache': compute._cache})
    except Exception as exc:
        import traceback
        q.put({'ok': False, 'error': str(exc), 'tb': traceback.format_exc()})


def _run_in_proc(target, args):
    """Spawn a child process, block until it finishes or is cancelled."""
    global _active_proc
    q = _mp_ctx.Queue()
    proc = _mp_ctx.Process(target=target, args=args + (q,))
    with _proc_lock:
        _active_proc = proc
    try:
        proc.start()
        proc.join()
    finally:
        with _proc_lock:
            if _active_proc is proc:
                _active_proc = None
    return proc, q


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    return send_from_directory('static', 'index.html')


@app.route('/compute', methods=['POST'])
def compute_route():
    try:
        result = _compute_mod.run_computation(request.json or {})
        return jsonify(result)
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


@app.route('/compute/solve', methods=['POST'])
def compute_solve():
    try:
        payload = request.json or {}
        proc, q = _run_in_proc(_solve_worker, (str(FC_DIR), payload))
        if proc.exitcode != 0:
            return jsonify({'error': 'Computation cancelled.', 'cancelled': True}), 499
        item = q.get_nowait()
        if not item['ok']:
            return jsonify({'error': item['error'], 'traceback': item.get('tb', '')}), 500
        _compute_mod._cache = item['cache']
        return jsonify(item['result'])
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


@app.route('/compute/grid', methods=['POST'])
def compute_grid():
    try:
        payload = request.json or {}
        if _compute_mod._cache is None:
            return jsonify({'error': 'No cached solve — call /compute/solve first'}), 400
        proc, q = _run_in_proc(_grid_worker, (str(FC_DIR), payload, _compute_mod._cache))
        if proc.exitcode != 0:
            return jsonify({'error': 'Computation cancelled.', 'cancelled': True}), 499
        item = q.get_nowait()
        if not item['ok']:
            return jsonify({'error': item['error'], 'traceback': item.get('tb', '')}), 500
        # Persist polygon cache and any other state updated in the subprocess.
        if 'cache' in item and _compute_mod._cache is not None:
            _compute_mod._cache['poly'] = item['cache'].get('poly', {})
        return jsonify(item['result'])
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


@app.route('/compute/cancel', methods=['POST'])
def compute_cancel():
    with _proc_lock:
        proc = _active_proc
    if proc and proc.is_alive():
        proc.terminate()
        return jsonify({'cancelled': True})
    return jsonify({'cancelled': False})


@app.route('/geometry')
def geometry():
    try:
        model = request.args.get('model', 'Simple')
        w = request.args.get('width',  None)
        l = request.args.get('length', None)
        return jsonify(_compute_mod.get_geometry(
            model,
            width_deg  = float(w) if w else None,
            length_deg = float(l) if l else None,
        ))
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


if __name__ == '__main__':
    app.run(debug=True, port=5001, use_reloader=False)
