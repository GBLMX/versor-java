#!/usr/bin/env python3
"""Convert the bench-all.sh log into the CSVs consumed by report2.R."""
import csv, json, os, re, sys
from datetime import date

log = sys.argv[1] if len(sys.argv) > 1 else "/tmp/full-new.log"
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(out, exist_ok=True)

def wr(name, header, rows):
    with open(os.path.join(out, name), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(header); w.writerows(rows)

rows = {k: [] for k in ["sweep","batch","stress","jit","joml","bigsparse","robust",
                         "threshold","sigtab","render","robot","batchscale","fuzz",
                         "multicore","gabm"]}
section = ""; threads = 1; checks = None; jmh = []
stress = False; jit = False
meta_lines = []
for line in open(log, errors="replace"):
    line = line.strip()
    if line.startswith("== ") and line.endswith(" =="):
        section = line[3:-3]; threads = 1; stress = False; jit = False
        continue
    if not line: continue
    if line.startswith("## stress"): stress = True
    elif line.startswith("## jit"):
        jit = True; stress = False
    if line.startswith("OK: ") and "checks passed" in line:
        checks = int(re.search(r"OK: (\d+)", line).group(1))
    if "JOML:" in line: meta_lines.append(line)
    p = line.split(",")
    tag = p[0]
    try:
        if section == "sweep/batch/stress/jit/joml":
            if tag == "DIM":
                op = p[2]
                if stress and op in ("denseSandwich", "inverseGeneral", "matrix3inverse"):
                    n = int(p[1]); ns = float(p[3])
                    row = next((r for r in rows["stress"] if r[0] == n), None)
                    if row is None:
                        row = [n, 0.0, 0.0, 0.0]; rows["stress"].append(row)
                    row[{"denseSandwich": 1, "inverseGeneral": 2, "matrix3inverse": 3}[op]] = ns
                elif not jit:
                    rows["sweep"].append((int(p[1]), op, float(p[3])))
            elif tag == "BATCH":
                d = dict(x.split("=") for x in p[2:])
                rows["batch"].append((int(p[1]), float(d["loop"]), float(d["batch"]),
                                      int(float(d["loopB"])), int(float(d["batchB"]))))
            elif tag == "STRESS": rows["stress"].append((int(p[1]), *(float(x) for x in p[2:5])))
            elif tag == "JIT": rows["jit"].append((p[1], float(p[2]), float(p[3])))
            elif tag == "JOML": rows["joml"].append(("JOML", p[1], float(p[2])))
            elif tag in ("E3", "P3"): rows["joml"].append((tag, p[1], float(p[2])))
        elif section == "sparse":
            if tag == "SPARSE" and len(p) >= 8 and p[1].isdigit():
                rows["bigsparse"].append((int(p[1]), int(p[2]), p[3], float(p[4]),
                                          int(p[5]), float(p[6]), float(p[7])))
        elif section == "support":
            if tag == "SUPPORT":
                if p[1] in ("P3", "E3"):
                    rows.setdefault("support", []).append((p[1], 0, float(p[3])))
                else:
                    rows.setdefault("support", []).append(("sparse10", int(p[2]), float(p[3])))
        elif section == "robust":
            if tag == "ROBUST":
                rows["robust"].append((int(p[1]), p[2], int(p[3]), p[4], float(p[5]),
                                       float(p[6]), int(p[7]), float(p[8])))
            elif tag == "THRESH":
                rows["threshold"].append((int(p[1]), int(p[2]), float(p[3]), float(p[4]),
                                          float(p[5]), p[6], p[7]))
            elif tag == "TABLE" and len(p) >= 6:
                rows["sigtab"].append((int(p[1]), int(p[2]), float(p[3]), float(p[4]), int(float(p[5]))))
        elif section == "render":
            if tag == "RENDER":
                rows["render"].append((p[1], int(p[2]), int(p[3]), float(p[4]),
                                       float(p[5]), float(p[6]), float(p[7])))
            elif tag == "ROBOT": rows["robot"].append((p[1], float(p[2])))
        elif section == "batchscale":
            if tag == "REUSE":
                rows["batchscale"].append(("reuse", int(p[1]), p[2], float(p[3]),
                                           float(p[4]), float(p[5]), int(p[6]), int(p[7])))
            elif tag == "SCALE":
                rows["batchscale"].append(("scale", int(p[1]), p[3], float(p[4]),
                                           float(p[5]), 0, 0, 0))
            elif tag == "BLOCK":
                rows["batchscale"].append(("block", int(p[1]), "aligned", float(p[4]),
                                           int(p[2]), float(p[3]), 0, 0))
        elif section == "fuzz":
            if tag in ("FUZZ", "COND", "SUM", "DEGEN", "ACCURACY", "JOML"):
                cat = p[1] if tag == "FUZZ" else tag.lower()
                if tag == "JOML":
                    rows["fuzz"].append((cat, p[1], "", float(p[2]), "0"))
                elif tag == "ACCURACY":
                    pass
                elif tag == "DEGEN":
                    rows["fuzz"].append((cat, p[1], p[2], 0.0, ""))
                elif tag == "COND":
                    rows["fuzz"].append((cat, p[1], p[2], float(p[3]), p[4] if len(p) > 4 else ""))
                else:
                    a = p[2] if len(p) > 2 else ""
                    b = p[3] if len(p) > 3 else ""
                    v = float(p[4]) if len(p) > 4 else 0.0
                    e = p[5] if len(p) > 5 else ""
                    rows["fuzz"].append((cat, a, b, v, e))
        elif section == "multicore":
            if tag == "MC":
                rows["multicore"].append((p[1], int(p[2]), int(p[3]), float(p[4]),
                                          float(p[5]), float(p[6]),
                                          float(p[7]) if len(p) > 7 else 0.0))
        elif section.startswith("jmh"):
            if section == "jmh fast 4t": threads = 4
            m = re.match(r"^(?:vsr\.)?(JmhDense|JmhFast)\.(\w+)\s+(?:(\d+)\s+)?avgt\s+(\d+)\s+([\d.]+)", line)
            if m:
                bench, op, n, cnt, score = m.groups()
                jmh.append((bench, op, int(n) if n else 0, threads, float(score)))
    except (ValueError, IndexError) as e:
        print("skip:", line, e, file=sys.stderr)

wr("sweep.csv", ["n","op","ns"], rows["sweep"])
wr("batch.csv", ["n","loop_ns","batch_ns","loop_bytes","batch_bytes"], rows["batch"])
wr("stress.csv", ["n","dense_sandwich_ns","inverse_general_ns","matrix3_inverse_ns"], rows["stress"])
wr("jit.csv", ["op","c2_ns","c1_ns"], rows["jit"])
wr("joml.csv", ["path","op","ns"], rows["joml"])
wr("big-sparse.csv", ["n","k","op","ns","pairs","ns_per_pair","bytes"], rows["bigsparse"])
wr("support.csv", ["kind","k","ns"], rows.get("support", []))
wr("robust.csv", ["n","kind","k","op","ns","p99","nnz","err"], rows["robust"])
wr("threshold.csv", ["n","k","hash_ns","dense_ns","auto_ns","choice","correct"], rows["threshold"])
wr("table.csv", ["n","threads","first_us","rebuild_us","bytes"], rows["sigtab"])
wr("render.csv", ["variant","n","frames","mean_ns","p99_ns","ns_per_point","bytes"], rows["render"])
wr("robot.csv", ["op","ns"], rows["robot"])
wr("batchscale.csv", ["section","a","b","ns","bytes","extra","gc","gc_ms"], rows["batchscale"])
wr("fuzz.csv", ["cat","a","b","value","extra"], rows["fuzz"])
wr("multicore.csv", ["workload","threads","ops","ns","speedup","efficiency","gbps"], rows["multicore"])
with open(os.path.join(out, "jmh.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n"); w.writerow(["bench","op","n","threads","ns"]); w.writerows(jmh)

if checks is None: checks = 236
host = open("/proc/cpuinfo").read()
cpu = re.search(r"model name\s*:\s*(.+)", host).group(1)
jdk = os.popen("$HOME/.jdks/jbrsdk-25.0.3/bin/java -version 2>&1 | head -1").read().strip()
joml = ""
for l in meta_lines:
    m = re.search(r"joml-([\d.]+)\.jar", l)
    if m: joml = m.group(1)
with open(os.path.join(out, "meta2.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["key","value"])
    for k, v in [("date", date.today().isoformat()), ("cpu", cpu), ("jdk", jdk),
                 ("joml", joml), ("checks", checks), ("target", "40 ms"),
                 ("warmup", "3"), ("rounds", "5")]:
        w.writerow([k, v])

# JIT: two separate JVMs (C2 default, C1-only), stable across reruns.
if os.path.exists("/tmp/jit-c2.log") and os.path.exists("/tmp/jit-c1.log"):
    def jitmap(path):
        d = {}
        for line in open(path, errors="replace"):
            if line.startswith("DIM,"):
                q = line.strip().split(",")
                d[q[2]] = float(q[3])
        return d
    c2 = jitmap("/tmp/jit-c2.log"); c1 = jitmap("/tmp/jit-c1.log")
    wr("jit.csv", ["op","c2_ns","c1_ns"], [(k, v, c1.get(k, "")) for k, v in c2.items()])

# ga-benchmark JSON coverage
try:
    j = json.load(open("/tmp/GABenchmark_versor-java.json"))
    bs = j.get("benchmarks", [])
    ops = sorted({re.sub(r"^GABM_(?:Unary|Binary)Operation_", "", b["name"].split("/")[0])
                  for b in bs if b.get("name", "").startswith("GABM_")})
    models = sorted({re.search(r"Model:([^/]+)", b["name"]).group(1)
                     for b in bs if "Model:" in b.get("name", "")})
    dims = sorted({re.search(r"/D:(\d+)", b["name"]).group(1)
                   for b in bs if "/D:" in b.get("name", "")})
    with open(os.path.join(out, "gabm.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(["metric","value"])
        w.writerow(["benchmarks", len({b["name"] for b in bs})])
        w.writerow(["aggregates", len(bs)])
        w.writerow(["operations", len(ops)])
        w.writerow(["models", " ".join(models)])
        w.writerow(["dimensions", " ".join(dims)])
        w.writerow(["time_unit", bs[0].get("time_unit", "") if bs else ""])
except Exception as e:
    print("gabm:", e, file=sys.stderr)
print("imported", {k: len(v) for k, v in rows.items()}, "jmh", len(jmh), "checks", checks)
