#!/usr/bin/env python3
"""Flask server for the interactive mantle flow web app."""
import sys
import pathlib

# Add flow_computations/ to the module search path so functions.py and euler_pole.py are importable.
FC_DIR = pathlib.Path(__file__).parent.parent / 'flow_computations'
sys.path.insert(0, str(FC_DIR))

from flask import Flask, request, jsonify, send_from_directory
from compute import run_computation, get_geometry

app = Flask(__name__, static_folder='static', static_url_path='')


@app.route('/')
def index():
    return send_from_directory('static', 'index.html')


@app.route('/compute', methods=['POST'])
def compute():
    try:
        result = run_computation(request.json or {})
        return jsonify(result)
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


@app.route('/geometry')
def geometry():
    try:
        model = request.args.get('model', 'LargeSP_RetreatingTrench')
        return jsonify(get_geometry(model))
    except Exception as exc:
        import traceback
        return jsonify({'error': str(exc), 'traceback': traceback.format_exc()}), 500


if __name__ == '__main__':
    app.run(debug=True, port=5001)
