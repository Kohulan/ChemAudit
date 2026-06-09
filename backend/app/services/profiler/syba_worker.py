"""Long-lived SYBA scoring worker (GPL-3.0 isolation).

Runs as a **separate process** so the GPL-3.0 ``syba`` package is never imported
into the Apache-2.0 application process. It loads the SYBA model once, then serves
one request per line of stdin and writes one response per line of stdout:

    request : {"smiles": "<SMILES>"}\\n
    response: {"score": <float>}\\n   or   {"error": "<message>"}\\n

On startup it emits a handshake line: {"ready": true} once the model is loaded,
or {"ready": false, "fatal": true, "error": "..."} if SYBA cannot be imported or
the model cannot be fitted (a permanent condition — the manager will not retry).

Keeping the model resident amortises the multi-second model-load cost across all
subsequent requests, instead of paying it on every call.
"""

import json
import sys


def _emit(obj: dict) -> None:
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


def main() -> None:
    try:
        # GPL-3.0 import confined to this isolated worker process.
        from syba.syba import SybaClassifier

        classifier = SybaClassifier()
        classifier.fitDefaultScore()
    except Exception as exc:  # pragma: no cover - exercised only without syba
        _emit({"ready": False, "fatal": True, "error": str(exc)[:200]})
        return

    _emit({"ready": True})

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
            score = classifier.predict(request.get("smiles", ""))
            _emit({"score": float(score)})
        except Exception as exc:  # pragma: no cover - per-request failure path
            _emit({"error": str(exc)[:200]})


if __name__ == "__main__":
    main()
