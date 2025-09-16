#!/usr/bin/env python3
import argparse, csv, json, math

def read_truth_origin(path):
    y = []
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            y.append(0 if row["true_parent"]=="A" else 1)
    return y

def read_truth_breaks(path):
    return json.load(open(path))["breakpoints_1based"]

def read_red_origins(path):
    pred = []
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            p = row["parent"]
            pred.append(0 if p=="A" else (1 if p=="B" else -1))
    return pred

def deriv_breakpoints(labels):
    """1-based indices where label changes (ignore INS=-1 by forward-fill)."""
    lab = []
    cur = None
    for v in labels:
        if v == -1:
            lab.append(cur if cur is not None else 0)  # forward-fill for stability
        else:
            lab.append(v); cur = v
    bps = []
    for i in range(1, len(lab)):
        if lab[i] != lab[i-1]:
            bps.append(i+1)  # 1-based start of new segment
    return bps

def per_base_accuracy(truth, pred):
    ok = 0; tot = 0
    cur = None
    for t, p in zip(truth, pred):
        if p == -1:
            # treat insertion as current parent (neutral); do not penalize twice
            p2 = cur if cur is not None else t
        else:
            p2 = p; cur = p
        ok += (t == p2)
        tot += 1
    return ok / max(1, tot)

def switch_count(x): return max(0, len(deriv_breakpoints(x)))

def match_breakpoints(true_bps, pred_bps, tol):
    """Greedy match within ±tol. Return TP, FP, FN and per-true error in bp."""
    used = [False]*len(pred_bps)
    tp=0; errs=[]
    for tb in true_bps:
        best=None; best_err=None; best_j=None
        for j, pb in enumerate(pred_bps):
            if used[j]: continue
            err = abs(pb - tb)
            if err <= tol and (best is None or err < best_err):
                best, best_err, best_j = tb, err, j
        if best is not None:
            used[best_j] = True
            tp += 1
            errs.append(best_err)
    fp = used.count(True)
    fn = len(true_bps) - tp
    prec = tp / max(1, (tp + (len(pred_bps)-tp)))  # tp / (tp+fp)
    rec  = tp / max(1, (tp + fn))
    f1 = 0.0 if (prec+rec)==0 else 2*prec*rec/(prec+rec)
    mean_err = (sum(errs)/len(errs)) if errs else None
    return dict(tp=tp, fp=(len(pred_bps)-tp), fn=fn, precision=prec, recall=rec, f1=f1, mean_abs_error_bp=mean_err)

def jaccard_segments(truth, pred, label=0):
    """Jaccard index over positions assigned to 'label' (0=A or 1=B). INS=-1 handled as current label."""
    T = []
    for t in truth: T.append(1 if t==label else 0)
    P = []
    cur=None
    for p in pred:
        if p==-1:
            P.append(1 if (cur==label) else 0)
        else:
            P.append(1 if p==label else 0); cur=p
    inter = sum(1 for a,b in zip(T,P) if a==1 and b==1)
    union = sum(1 for a,b in zip(T,P) if a==1 or b==1)
    return inter / max(1, union)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-origin", required=True, help="truth.origin.tsv from generator")
    ap.add_argument("--truth-breaks", required=True, help="truth.breakpoints.json from generator")
    ap.add_argument("--red-origins", required=True, help="*.origins.tsv produced by RED")
    ap.add_argument("--tol", type=int, default=10, help="breakpoint tolerance (bp)")
    args = ap.parse_args()

    truth = read_truth_origin(args.truth_origin)
    true_bps = read_truth_breaks(args.truth_breaks)
    pred = read_red_origins(args.red_origins)
    assert len(truth)==len(pred), "length mismatch: truth vs RED origins"

    acc = per_base_accuracy(truth, pred)
    pred_bps = deriv_breakpoints(pred)
    sw_true = max(0, len(true_bps))
    sw_pred = max(0, len(pred_bps))
    sw_err  = abs(sw_true - sw_pred)

    bp_metrics = match_breakpoints(true_bps, pred_bps, tol=args.tol)
    jA = jaccard_segments(truth, pred, label=0)
    jB = jaccard_segments(truth, pred, label=1)

    print("=== RED evaluation ===")
    print(f"Per-base parent accuracy: {acc:.4f}")
    print(f"Switch count (true={sw_true}, pred={sw_pred}) | abs error: {sw_err}")
    print(f"Breakpoint tolerance: ±{args.tol} bp")
    print(f"Breakpoints: TP={bp_metrics['tp']} FP={bp_metrics['fp']} FN={bp_metrics['fn']}")
    print(f"Precision={bp_metrics['precision']:.3f}  Recall={bp_metrics['recall']:.3f}  F1={bp_metrics['f1']:.3f}")
    if bp_metrics['mean_abs_error_bp'] is not None:
        print(f"Mean |Δbp| among matched breakpoints: {bp_metrics['mean_abs_error_bp']:.1f} bp")
    print(f"Jaccard(A)={jA:.3f}  Jaccard(B)={jB:.3f}")

if __name__ == "__main__":
    main()

