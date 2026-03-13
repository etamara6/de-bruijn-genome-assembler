"""
api.py — Flask REST API for the De Bruijn genome assembler.

Exposes:
    POST /api/assemble  — assemble contigs from a list of reads
    GET  /api/health    — liveness check
"""

from __future__ import annotations

import logging
import traceback
from dataclasses import asdict

from flask import Flask, jsonify, request, Response
from flask_cors import CORS

from assembler import GenomeAssembler, AssemblyResult



app = Flask(__name__)
CORS(app)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



def _error_response(message: str, status_code: int) -> tuple[Response, int]:
    """Build a uniform JSON error response.

    Args:
        message: Human-readable error description.
        status_code: HTTP status code to return.

    Returns:
        A ``(Response, int)`` tuple suitable for returning from a Flask view.
    """
    return jsonify({"error": message}), status_code




@app.route("/api/health", methods=["GET"])
def health() -> tuple[Response, int]:
    """Liveness check endpoint.

    Returns:
        JSON ``{"status": "ok"}`` with HTTP 200.
    """
    return jsonify({"status": "ok"}), 200


@app.route("/api/assemble", methods=["POST"])
def assemble() -> tuple[Response, int]:
    """Assemble contigs from a list of DNA reads.

    Request body (JSON):
        reads (list[str]): One or more DNA sequences.
        k     (int):       K-mer length (default 4).

    Returns:
        JSON with keys:
            contigs (list[str]):      Assembled contigs, longest-first.
            stats   (dict):           Assembly statistics.
                total_kmers (int)
                graph_nodes (int)
                graph_edges (int)
                elapsed_ms  (float)
                longest_contig (str)
                n50 (int)

    HTTP status codes:
        200 — success
        400 — bad request (missing / malformed input)
        422 — unprocessable entity (algorithm error such as no Eulerian path)
        500 — unexpected server error
    """
    payload = request.get_json(silent=True)
    if payload is None:
        return _error_response("Request body must be valid JSON.", 400)

    reads = payload.get("reads")
    if reads is None:
        return _error_response("Missing required field: 'reads'.", 400)
    if not isinstance(reads, list):
        return _error_response("'reads' must be a JSON array of strings.", 400)
    if not all(isinstance(r, str) for r in reads):
        return _error_response("Every element of 'reads' must be a string.", 400)

    k = payload.get("k", 4)
    if not isinstance(k, int) or isinstance(k, bool):
        return _error_response("'k' must be an integer.", 400)

    logger.info("Received assembly request: %d reads, k=%d", len(reads), k)

    try:
        assembler = GenomeAssembler(reads=reads, k=k)
        result: AssemblyResult = assembler.assemble()
    except ValueError as exc:
        logger.warning("Validation error during assembly: %s", exc)
        return _error_response(str(exc), 400)
    except RuntimeError as exc:
        logger.warning("Algorithm error during assembly: %s", exc)
        return _error_response(str(exc), 422)
    except Exception:
        logger.error("Unexpected error:\n%s", traceback.format_exc())
        return _error_response("An unexpected server error occurred.", 500)

    response_body = {
        "contigs": result.contigs,
        "stats": {
            "total_kmers": result.total_kmers,
            "graph_nodes": result.graph_nodes,
            "graph_edges": result.graph_edges,
            "elapsed_ms": round(result.elapsed_ms, 4),
            "longest_contig": result.longest_contig,
            "n50": result.n50,
        },
    }

    logger.info(
        "Assembly complete: %d contig(s), longest=%d bp, %.2f ms",
        len(result.contigs),
        len(result.longest_contig),
        result.elapsed_ms,
    )

    return jsonify(response_body), 200



if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)