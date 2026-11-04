"""The process coordinator should not retain a second numerical runtime."""

import json
from pathlib import Path
import subprocess
import sys


def test_process_coordinator_stays_lightweight_and_file_errors_continue(sample_file, tmp_path):
    sources = []
    for i in range(3):
        source = tmp_path / f"good-{i}.mzML"
        source.write_bytes(sample_file.read_bytes())
        sources.append(source)
    damaged = tmp_path / "damaged.mzML"
    damaged.write_text("<mzML><broken")
    sources.insert(1, damaged)
    script = """
import json, sys
from masscubex import Config
from masscubex.cli import run_batch
heavy = ('numpy', 'numba', 'scipy', 'pyarrow')
before = [name for name in heavy if name in sys.modules]
rows = list(run_batch(sys.argv[2:], sys.argv[1], Config(), workers=2))
after = [name for name in heavy if name in sys.modules]
print(json.dumps({'before': before, 'after': after, 'rows': rows}))
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path / "outputs"), *map(str, sources)],
        check=True, capture_output=True, text=True, timeout=120,
    )
    report = json.loads(result.stdout)
    assert report["before"] == report["after"] == []
    rows = {Path(row["input"]).name: row for row in report["rows"]}
    assert set(rows) == {p.name for p in sources}
    assert rows[damaged.name]["status"] == "error"
    assert all(rows[p.name]["status"] == "ok" for p in sources if p != damaged)
