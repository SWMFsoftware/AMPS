#!/usr/bin/env python3
"""Local HTTP front end for the Stage-9 verification interface.

The server binds to loopback by default and exposes only named bundles below a
configured output root.  It never accepts arbitrary output paths or arbitrary
commands from the browser.  Scientific request validation is delegated to
``stage9_interface.py`` so CLI and browser submissions have identical guards.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, Optional
from urllib.parse import parse_qs, urlparse

# Keep UI operation from leaving temporary Python bytecode in the source tree.
sys.dont_write_bytecode = True
import stage9_interface as core


UI_ROOT = core.ROOT / "ui" / "stage9"
_BUNDLE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,95}$")


class Stage9Server(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], output_root: pathlib.Path,
                 executable: Optional[pathlib.Path]):
        super().__init__(address, Handler)
        self.output_root = output_root.resolve()
        self.output_root.mkdir(parents=True, exist_ok=True)
        self.executable = executable.resolve() if executable else None
        self.runners: Dict[str, subprocess.Popen] = {}
        self.runner_logs: Dict[str, Any] = {}

    def bundle(self, name: Any) -> pathlib.Path:
        if not isinstance(name, str) or not _BUNDLE_NAME.fullmatch(name):
            raise core.Stage9Error("bundle_name must contain only letters, digits, '.', '_' or '-'")
        return self.output_root / name


class Handler(BaseHTTPRequestHandler):
    server: Stage9Server

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        sys.stderr.write("stage9-ui: " + format % args + "\n")

    def _json(self, status: HTTPStatus, value: Any) -> None:
        data = json.dumps(value, indent=2, sort_keys=True).encode("utf-8")
        self.send_response(status.value)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _body(self) -> Any:
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 10 * 1024 * 1024:
                raise core.Stage9Error("request body must contain at most 10 MiB of JSON")
            return json.loads(self.rfile.read(length))
        except (ValueError, json.JSONDecodeError) as exc:
            raise core.Stage9Error(f"invalid JSON body: {exc}") from exc

    def _static(self, name: str, content_type: str) -> None:
        path = UI_ROOT / name
        if not path.is_file():
            self.send_error(HTTPStatus.NOT_FOUND.value)
            return
        data = path.read_bytes()
        self.send_response(HTTPStatus.OK.value)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        try:
            if parsed.path in ("/", "/index.html"):
                return self._static("index.html", "text/html; charset=utf-8")
            if parsed.path == "/app.js":
                return self._static("app.js", "text/javascript; charset=utf-8")
            if parsed.path == "/styles.css":
                return self._static("styles.css", "text/css; charset=utf-8")
            if parsed.path == "/api/capabilities":
                return self._json(HTTPStatus.OK, core._read_json(core.DEFAULT_REGISTRY))
            if parsed.path == "/api/server":
                return self._json(HTTPStatus.OK, {
                    "submission_enabled": self.server.executable is not None,
                    "submission_disabled_reason": (None if self.server.executable is not None
                        else "server was started without --executable; preview and compile-only are available"),
                    "output_root": str(self.server.output_root),
                })
            if parsed.path == "/api/materials":
                return self._json(HTTPStatus.OK, core._read_json(core.DEFAULT_MATERIALS))
            if parsed.path == "/api/example":
                return self._json(HTTPStatus.OK,
                                  core._read_json(core.ROOT / "examples/stage9_ui_request.json"))
            if parsed.path == "/api/job":
                query = parse_qs(parsed.query)
                bundle = self.server.bundle(query.get("bundle_name", [None])[0])
                state = bundle / "job_state.json"
                if not state.exists():
                    raise core.Stage9Error("bundle has no job state")
                return self._json(HTTPStatus.OK, core._read_json(state))
            self.send_error(HTTPStatus.NOT_FOUND.value)
        except core.Stage9Error as exc:
            self._json(HTTPStatus.BAD_REQUEST, {"status": "error", "error": str(exc)})

    def do_POST(self) -> None:  # noqa: N802
        try:
            body = core._object(self._body(), "HTTP request")
            if self.path == "/api/preview":
                compiled = core.compile_request(core._object(body.get("request"), "request"))
                return self._json(HTTPStatus.OK, compiled["manifest"])
            if self.path == "/api/compile":
                bundle_name = body.get("bundle_name")
                bundle = self.server.bundle(bundle_name)
                compiled = core.compile_request(core._object(body.get("request"), "request"))
                core.write_bundle(compiled, bundle)
                validation_records = []
                if self.server.executable is not None:
                    for entry in compiled["manifest"]["generated_cases"]:
                        validation_dir = bundle / "preflight_validation" / entry["case_id"]
                        validation_dir.mkdir(parents=True, exist_ok=True)
                        log_path = validation_dir / "validate.log"
                        with log_path.open("w", encoding="utf-8") as log:
                            result = subprocess.run(
                                [str(self.server.executable), "--validate-case=" +
                                 str(bundle / entry["relative_path"])],
                                cwd=validation_dir, stdout=log,
                                stderr=subprocess.STDOUT, text=True, check=False)
                        record = {"case_id": entry["case_id"],
                                  "return_code": result.returncode,
                                  "log": str(log_path.relative_to(bundle))}
                        validation_records.append(record)
                        if result.returncode != 0:
                            core._write_json(bundle / "compile_validation.json", {
                                "status": "failed_cpp_case_validation",
                                "records": validation_records})
                            raise core.Stage9Error(
                                f"generated case '{entry['case_id']}' failed production "
                                f"C++ validation; see {log_path}")
                    validation_status = "passed_all_generated_cases"
                else:
                    validation_status = "not_performed_server_has_no_executable"
                core._write_json(bundle / "compile_validation.json", {
                    "status": validation_status, "records": validation_records})
                return self._json(HTTPStatus.CREATED, {
                    "status": "compiled", "bundle_name": bundle_name,
                    "bundle": str(bundle),
                    "manifest_content_hash": compiled["manifest"]["manifest_content_hash"],
                    "cpp_case_validation": validation_status,
                })
            if self.path == "/api/job/start":
                bundle_name = body.get("bundle_name")
                bundle = self.server.bundle(bundle_name)
                if not (bundle / "manifest.json").is_file():
                    raise core.Stage9Error("compile the named bundle before submission")
                if self.server.executable is None:
                    raise core.Stage9Error("server started without --executable; submission is disabled")
                old = self.server.runners.get(bundle_name)
                if old is not None and old.poll() is None:
                    raise core.Stage9Error("the named bundle is already running")
                resume = bool(body.get("resume", False))
                command = [sys.executable, str(core.ROOT / "scripts/stage9_interface.py"),
                           "run", "--bundle", str(bundle), "--executable",
                           str(self.server.executable)]
                if resume:
                    command.append("--resume")
                log = (bundle / "runner.log").open("a", encoding="utf-8")
                process = subprocess.Popen(command, cwd=core.ROOT, stdout=log,
                                           stderr=subprocess.STDOUT, text=True)
                self.server.runners[bundle_name] = process
                self.server.runner_logs[bundle_name] = log
                return self._json(HTTPStatus.ACCEPTED, {
                    "status": "started", "bundle_name": bundle_name,
                    "pid": process.pid, "resume": resume,
                })
            if self.path == "/api/job/cancel":
                bundle_name = body.get("bundle_name")
                marker = core.request_cancel(self.server.bundle(bundle_name))
                return self._json(HTTPStatus.ACCEPTED, {
                    "status": "cancellation_requested", "marker": str(marker)})
            self.send_error(HTTPStatus.NOT_FOUND.value)
        except core.Stage9Error as exc:
            self._json(HTTPStatus.BAD_REQUEST, {"status": "error", "error": str(exc)})
        except Exception as exc:  # Return a stable API failure while retaining server logs.
            self.log_error("unhandled request error: %s", exc)
            self._json(HTTPStatus.INTERNAL_SERVER_ERROR,
                       {"status": "error", "error": "internal server error; see server log"})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1",
                        help="listen address; loopback is the safe default")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--output-root", type=pathlib.Path,
                        default=core.ROOT / "test_output" / "stage9_ui_jobs")
    parser.add_argument("--executable", type=pathlib.Path,
                        help="shieldSim executable; omit to keep Run disabled")
    args = parser.parse_args()
    server = Stage9Server((args.host, args.port), args.output_root, args.executable)
    print(f"Stage-9 verification UI: http://{args.host}:{args.port}/")
    print(f"Bundles: {server.output_root}")
    print("Submission:", str(server.executable) if server.executable else "disabled")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
