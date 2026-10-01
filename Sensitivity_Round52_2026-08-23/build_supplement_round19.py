"""Assemble supplement_revision_round19_clean.md from the Round 18 supplement plus Round 52 additions.

Sources (read only):
- draft/supplement_revision_round18_clean.md (base)
- C:/Code/UG_ERP_Project/results/Behavior/{E1,E2}_TwoStage_Bates/GLMM_Rejection/posthoc_Main_Effect_emotion.md (2026-08-22)
- project-details/.../hDDM_E{1,2}/results_hddm_final/audit/focal_parameter_audit_{va,vaz}.csv
- draft/sensitivity/actor-random-intercept_2026-08-23/*.csv
- draft/sensitivity/erp-registered-E2_2026-08-23/*.csv (only the +/-100 uV rows; registered windows are not reported by author decision)
- draft/phase6_revision/round52_model_specification_table.md
"""
import re
import pandas as pd

D = "E:/MyPaper_withARS/draft/"
BASE = D + "supplement_revision_round18_clean.md"
OUT = D + "supplement_revision_round19_clean.md"
CODE = "C:/Code/UG_ERP_Project/results/"
HDDM = "E:/MyPaper_withARS/project-details/UG_ERP_Project/hDDM__details/"
ACT = D + "sensitivity/actor-random-intercept_2026-08-23/"
ERP = D + "sensitivity/erp-registered-E2_2026-08-23/"

text = open(BASE, encoding="utf-8").read()
assert text.startswith("# Supplementary Materials")
# Table S3: restore the half-point medians of the canonical CSV (Round 51 Table 1 precision)
for old, new in [("| 32 [19, 57] |", "| 32.5 [19, 57] |"), ("| 76 [53, 99] |", "| 76.5 [53, 99] |")]:
    assert text.count(old) == 1, old
    text = text.replace(old, new)

# ---------- S11 replacement: va and vaz side by side ----------
LAB = {"aff": "Affiliative", "dis": "Disgust", "dom": "Dominance", "rew": "Reward"}
def audit(exp, model):
    df = pd.read_csv(f"{HDDM}hDDM_{exp}/results_hddm_final/audit/focal_parameter_audit_{model}.csv")
    df = df.rename(columns={"Unnamed: 0": "param"})
    out = {}
    for _, r in df.iterrows():
        m = re.match(r"([vaz])_C\(emotion, Treatment\('neu'\)\)\[T\.(\w+)\]", r["param"])
        if m:
            out[(m.group(1), m.group(2))] = (r["mean"], r["hdi_2.5%"], r["hdi_97.5%"])
    return out
def fmt(v):
    return f"{v[0]:+.3f} [{v[1]:.3f}, {v[2]:.3f}]"
s11 = ["**Supplementary Table S11. Reported `va` model and exploratory `vaz` model side by side.** Group-level proposer-expression contrasts relative to Neutral (posterior means with 95% HDIs from the canonical posterior summaries) for drift rate (*v*), decision threshold (*a*), and, in `vaz`, starting point (*z*; positive = toward the acceptance boundary). `va` is the preregistered parameterization (Experiment 2 registered *v* and *a* as the HDDM dependent variables and no starting-point hypothesis); `vaz` is the DIC-minimal model in Experiment 2 (ΔDIC = 12.85) and is reported as exploratory. Conclusions that differ between the two columns are model-dependent.", "",
       "| Experiment | Contrast vs Neutral | `va` Δ*v* | `vaz` Δ*v* | `va` Δ*a* | `vaz` Δ*a* | `vaz` Δ*z* |", "|---|---|---:|---:|---:|---:|---:|"]
for exp in ["E1", "E2"]:
    va, vaz = audit(exp, "va"), audit(exp, "vaz")
    for k in ["dis", "dom", "aff", "rew"]:
        s11.append(f"| {exp[1]} | {LAB[k]} | {fmt(va[('v',k)])} | {fmt(vaz[('v',k)])} | {fmt(va[('a',k)])} | {fmt(vaz[('a',k)])} | {fmt(vaz[('z',k)])} |")
s11.append("")
s11.append("*Note.* Drift direction and rank order (Disgust < Dominance < 0 < Reward ≈ Affiliative) are identical in both models and every drift HDI excludes zero. Threshold contrasts are model-dependent: in Experiment 2 the Disgust and Dominance threshold HDIs exclude zero under `va` but include zero under `vaz`, and the Affiliative threshold contrast becomes credibly negative under `vaz`. All focal parameters converged (maximum R̂ ≤ 1.004 in both models and experiments). Main-text figures report posterior medians of the same posteriors. Source: `results_hddm_final/audit/focal_parameter_audit_{va,vaz}.csv`.")
i = text.index("**Supplementary Table S11.")
j = text.index("## Supplementary figures")
text = text[:i] + "\n".join(s11) + "\n\n" + text[j:]

# ---------- S12 evidence-status map ----------
s12 = """**Supplementary Table S12. Preregistration-to-result map.** Each registered hypothesis or analysis component of both experiments, the estimand actually analyzed, the outcome, the deviation, and the final classification (`confirmatory` = registered hypothesis analyzed as registered or with a disclosed harmonized specification preserving the registered estimand; `harmonized/post hoc` = Experiment 2 outcome analyzed with a window fixed after Experiment 1; `exploratory` = not registered; `sensitivity` = robustness check; `unevaluated` = registered but not analyzed, reason given).

| ID | Experiment | Registered prediction | Analyzed estimand | Outcome | Deviation | Classification |
|---|---|---|---|---|---|---|
| E1-H1 | 1 | Unfair-offer rejection ordering Dominance > Affiliative > Reward | Expression × offer-type GLMM; expression main-effect pairs (FDR over 10 pairs) | Interaction χ²(4) = 6.34, *p* = .175 → offer-specific ordering not supported; across offer types Dominance > Affiliative and Dominance > Reward; Affiliative vs Reward n.s. | Inference standard harmonized to Experiment 2 (Type III + FDR) | confirmatory (offer-specific prediction unsupported); across-offer contrasts use the registered factor at an unregistered level |
| E1-H2a | 1 | FRN: unfair > fair | FRN 252–304 ms, collapsed localizer, baseline covariate | *F*(1, 14142.02) = 0.16, *p* = .686 | Baseline covariate primary; registered subtraction as sensitivity | confirmatory (unsupported) |
| E1-H2b | 1 | FRN: expression × offer interaction (unfair Reward > Affiliative > Dominance) | Same model | Interaction *p* = .677; expression main effect *p* = .109 (.014 under subtraction) | As above | confirmatory (unsupported); baseline-dependent expression effect disclosed |
| E1-H3a | 1 | P3–LPP: fair > unfair | 400–600 ms fixed window | *F*(1, 29.12) = 45.19, *p* < .001; +0.89 µV | None for window | confirmatory (supported) |
| E1-H3b | 1 | P3–LPP: expression × offer interaction | Same model | *p* = .322 | As above | confirmatory (unsupported) |
| E1-H4 | 1 | Face-locked P1/N170/EPN/LPP: all smiles > Neutral; Reward > Affiliative | Face-locked LMMs | Not supported; Disgust-led N170/EPN effects; face-LPP Dominance > Reward/Affiliative | None | confirmatory (unsupported); face-LPP ordering exploratory |
| E1-DDM | 1 | Not registered | HDDM `va`, unfair offers | Drift contrasts all exclude zero; threshold credible for Disgust only | — | exploratory |
| E1-RT | 1 | RT registered as DV without directional hypothesis | Log-RT LMM | Expression × offer interaction reliable; fair-offer contrasts | — | exploratory |
| E2-H1a | 2 | Choice: expression × fairness interaction; Reward reduces unfair rejection relative to Disgust | GLMM as Experiment 1 | Interaction χ²(4) = 8.71, *p* = .069 → not supported; across-offer contrasts recurred | None | confirmatory (unsupported) |
| E2-H1b | 2 | RT: Reward-unfair slower than Disgust-unfair | Log-RT LMM, unfair simple contrasts | Disgust − Reward = −1.81%, *p* = .965 | RT window 150–3000 → 300–3000 ms | confirmatory (unsupported) |
| E2-H2 | 2 | HDDM: *v* and *a* as DVs; Reward decreases and Disgust increases drift (registered in an ambiguity framing) | `va` model; Reward and Disgust Δ*v* vs Neutral | Reward +0.221, Disgust −0.616; both HDIs exclude zero | Directional restatement of the registered framing | confirmatory (supported); five-expression pattern, other threshold contrasts, *z*, and model comparison exploratory |
| E2-H3a | 2 | N400 (350–450 ms): fairness and expression main effects, no interaction | Not analyzed | — | Component set harmonized to Experiment 1 (FRN and 400–600 ms P3–LPP) | unevaluated (see Experiment 2 Deviations) |
| E2-H3b | 2 | LPP 500–800 ms: expression × fairness interaction (Reward-unfair > Disgust-unfair) | 400–600 ms harmonized window | Fair > unfair replicated; interaction *p* = .051 | Window harmonized to Experiment 1 | harmonized/post hoc conceptual replication; registered 500–800 ms window unevaluated |
| E2-H4a | 2 | N2 (~200–300 ms; registered as exploratory): unfair > fair, modulated by expression | Early fronto-central 192–244 ms window (Experiment 1 localizer rule) | Offer *p* = .152; expression *p* = .001 (compound face + offer window) | Window not registered for Experiment 2 | harmonized/post hoc; registered N2 unevaluated |
| E2-H4b | 2 | CPP (response-locked; registered as exploratory) | Not analyzed | — | No response-locked analysis developed | unevaluated; future work |
| E2-prep | 2 | ±100 µV peak-to-peak rejection | ±200 µV pipeline default; ±100 µV subset as sensitivity (Table S16) | No conclusion changed | Disclosed | sensitivity |
"""

# ---------- S13 model specification ----------
spec = open(D + "phase6_revision/round52_model_specification_table.md", encoding="utf-8").read()
rows = [l for l in spec.splitlines() if l.startswith("| ") and not l.startswith("| Analysis")]
s13 = ["**Supplementary Table S13. Final specification of every reported mixed model.** Formulas are the final models retained by the preregistered `buildmer` backward likelihood-ratio protocol (elimination α = .20) applied to the design-justified maximal model; `RE_<factor>_k` denotes the *k*-th sum-to-zero contrast column of a within-participant factor, each with its own one-dimensional random-slope term (zero-correlation parameterization), with correlations tested in a final step. `dims` = retained / nominal random-effect dimensions. All final models met the acceptance gate (optimizer code 0, no singular variance component, positive-definite Hessian, finite standard errors and likelihood). `Alday` = baseline-as-covariate primary analysis; `Traditional` = subtraction sensitivity analysis; `FacePhase_E1` = face-locked models. Source: `REPRODUCIBILITY.txt` of each output folder (2026-08-22 run).", "",
       "| Analysis | Engine | Optimizer | Singular | PD Hessian | dims | Inference | *n* | Final formula |", "|---|---|---|---|---|---|---|---|---|"]
for l in rows:
    c = [x.strip() for x in l.strip().strip("|").split("|")]
    label = c[0].replace("Behavior/", "").replace("EEG/", "").replace("/Stage1_TrialLevel", "").replace("_TwoStage_Bates", "")
    s13.append(f"| {label} | {c[1]} | {c[3]} | {c[4]} | {c[5]} | {c[6]} | {c[7].replace('_',' ')} | {c[8]} | {c[9]} |")
s13 = "\n".join(s13) + "\n"

# ---------- S14 actor sensitivity ----------
diag = pd.read_csv(ACT + "diagnostics.csv"); diag = diag[~diag.model.str.contains("attempt")]
omni = pd.read_csv(ACT + "omnibus_comparison.csv")
vc = pd.read_csv(ACT + "variance_components.csv")
con = pd.read_csv(ACT + "emotion_contrasts_comparison.csv")
ed = pd.read_csv(ACT + "erp_diagnostics.csv"); eo = pd.read_csv(ACT + "erp_omnibus_comparison.csv"); ef = pd.read_csv(ACT + "erp_fair_minus_unfair_comparison.csv")
def p(x):
    return "< .001" if x < .001 else f"= {x:.3f}".replace("0.", ".")
s14 = ["**Supplementary Table S14. Actor-random-intercept sensitivity analysis.** Every focal model of Table S13 for choice, response time, and the two confirmatory offer-locked components was refitted with its exact final formula (canonical refit, reproducing the reported statistics) and again with a crossed random intercept for the 60 proposer actors added. Actor-level expression slopes were additionally attempted and were singular (boundary) except for the Experiment 2 response-time model, where they were negligible and left all conclusions unchanged. The Experiment 1 early-window (FRN) actor model was also singular (actor variance estimated at zero) and is therefore not interpreted. Variances are on the logit scale (choice), log-RT scale (RT), and µV² (ERP). Canonical refits reproduce the reported statistics to the second decimal (for example, Experiment 2 offer type χ²(1) = 196.02 here versus 196.04 in the canonical output), the small differences reflecting optimizer tolerance. Source: `draft/sensitivity/actor-random-intercept_2026-08-23/`.", "",
       "| Outcome | Experiment | Model | Actor intercept variance (SD) | Participant intercept variance (SD) | Expression | Offer type | Expression × offer |", "|---|---|---|---|---|---|---|---|"]
for _, r in diag[diag.model.isin(["canonical_refit", "actor_intercept"]) & diag.valid].iterrows():
    o = omni[(omni.experiment == r.experiment) & (omni.outcome == r.outcome) & (omni.model == r.model)].set_index("term")
    pv = vc[(vc.experiment == r.experiment) & (vc.outcome == r.outcome) & (vc.model == r.model) & vc.group.str.startswith("participant") & vc.term.str.contains("Intercept")]
    pvar = pv.variance.iloc[0] if len(pv) else float("nan")
    av = "—" if pd.isna(r.actor_intercept_var) else f"{r.actor_intercept_var:.4f} ({r.actor_intercept_var**0.5:.3f})"
    def cell(term):
        x = o.loc[term]
        return (f"χ²({int(x.df1)}) = {x.statistic:.2f}, *p* {p(x.p)}" if r.outcome == "choice" else f"*F*({int(x.df1)}, {x.df2:.2f}) = {x.statistic:.2f}, *p* {p(x.p)}")
    s14.append(f"| {'Choice' if r.outcome=='choice' else 'log RT'} | {r.experiment[1]} | {r.model.replace('_',' ')} | {av} | {pvar:.4f} ({pvar**0.5:.3f}) | {cell('emotion')} | {cell('offer_type')} | {cell('emotion:offer_type')} |")
for _, r in ed[ed.valid].iterrows():
    o = eo[(eo.experiment == r.experiment) & (eo.component == r.component) & (eo.model == r.model)].set_index("term")
    f = ef[(ef.experiment == r.experiment) & (ef.component == r.component) & (ef.model == r.model)].iloc[0]
    av = "—" if pd.isna(r.actor_intercept_var) else f"{r.actor_intercept_var:.4f} ({r.actor_intercept_var**0.5:.3f})"
    comp = "P3–LPP 400–600 ms" if r.component == "LPP_pre" else "Early fronto-central window"
    def cell(term):
        x = o.loc[term]; return f"*F*({int(x.df1)}, {x.df2:.2f}) = {x.F:.2f}, *p* {p(x.p)}"
    s14.append(f"| {comp} (fair − unfair = {f.estimate:+.3f} µV [{f['lower.CL']:.3f}, {f['upper.CL']:.3f}]) | {r.experiment[1]} | {r.model.replace('_',' ')} | {av} | residual {r.residual_var:.2f} | {cell('emotion')} | {cell('offer_type')} | {cell('emotion:offer_type')} |")
s14 += ["", "Focal across-offer choice contrasts (link scale, FDR over the ten-pair family):", "", "| Experiment | Contrast | Canonical *b* (SE), *p*~FDR~ | With actor intercept *b* (SE), *p*~FDR~ |", "|---|---|---|---|"]
cm = con[(con.family == "emotion_marginal") & (con.outcome == "choice")]
NAME = {"dis - dom": "Disgust − Dominance", "aff - dom": "Affiliative − Dominance", "dom - enj": "Dominance − Reward", "aff - enj": "Affiliative − Reward"}
for exp in ["E1", "E2"]:
    for k, nm in NAME.items():
        a = cm[(cm.experiment == exp) & (cm.contrast == k) & (cm.model == "canonical_refit")].iloc[0]
        b = cm[(cm.experiment == exp) & (cm.contrast == k) & (cm.model == "actor_intercept")].iloc[0]
        s14.append(f"| {exp[1]} | {nm} | {a.estimate:.3f} ({a.SE:.3f}), {p(a.p_fdr)} | {b.estimate:.3f} ({b.SE:.3f}), {p(b.p_fdr)} |")
s14 = "\n".join(s14) + "\n"

# ---------- S15 complete ten-pair family + interaction contrasts ----------
def pairs_md(exp):
    t = open(f"{CODE}Behavior/{exp}_TwoStage_Bates/GLMM_Rejection/posthoc_Main_Effect_emotion.md", encoding="utf-8").read()
    rows = [l for l in t.splitlines() if l.startswith("|") and not l.startswith("|:") and not l.startswith("|contrast")]
    out = []
    for l in rows:
        c = [x.strip() for x in l.strip().strip("|").split("|")]
        out.append((c[0], c[1], c[2], c[7], c[13], c[14], c[15]))
    return out
L2 = {"neu": "Neutral", "aff": "Affiliative", "dis": "Disgust", "dom": "Dominance", "enj": "Reward"}
s15 = ["**Supplementary Table S15. Complete proposer-expression main-effect family (all ten pairs) and expression-by-offer interaction contrasts for the choice models.** Pairwise contrasts are estimated marginal means on the logit scale from the final rejection GLMMs (Table S13), with *p* values FDR-adjusted over the ten-pair family (the values reported in the main text) and odds ratios with 95% CIs. Interaction contrasts give the fair-minus-unfair difference of each expression contrast (unadjusted 95% CIs); their width indicates the smallest offer-specific modulation the design could have detected. Source: `GLMM_Rejection/posthoc_Main_Effect_emotion.md` (2026-08-22) and the canonical refits in `draft/sensitivity/actor-random-intercept_2026-08-23/choice_interaction_contrasts_canonical.csv`.", "",
       "| Experiment | Contrast | *b* | SE | *p*~FDR~ | OR | 95% CI |", "|---|---|---:|---:|---:|---:|---|"]
for exp in ["E1", "E2"]:
    for c in pairs_md(exp):
        a, b = c[0].split(" - ")
        s15.append(f"| {exp[1]} | {L2[a]} − {L2[b]} | {c[1]} | {c[2]} | {c[3]} | {c[4]} | [{c[5]}, {c[6]}] |")
ic = pd.read_csv(ACT + "choice_interaction_contrasts_canonical.csv")
s15 += ["", "| Experiment | Expression contrast | Fair − unfair difference | SE | 95% CI | *p* (unadjusted) |", "|---|---|---:|---:|---|---:|"]
for exp in ["E1", "E2"]:
    for _, r in ic[ic.experiment == exp].iterrows():
        a, b = r.emotion_pairwise.split(" - ")
        s15.append(f"| {exp[1]} | {L2[a]} − {L2[b]} | {r.estimate:.3f} | {r.SE:.3f} | [{r['asymp.LCL']:.3f}, {r['asymp.UCL']:.3f}] | {r['p.value']:.3f} |")
s15 = "\n".join(s15) + "\n"

# ---------- S16 +/-100 uV sensitivity ----------
om = pd.read_csv(ERP + "erp_sensitivity_omnibus.csv"); fu = pd.read_csv(ERP + "erp_sensitivity_fair_minus_unfair.csv"); dg = pd.read_csv(ERP + "erp_sensitivity_diagnostics.csv")
keep = ["LPP_pre canonical set (validation)", "LPP_pre +/-100 uV subset (R1.6)", "FRN_pre canonical set (validation)", "FRN_pre +/-100 uV subset (R1.6)"]
s16 = ["**Supplementary Table S16. Experiment 2 offer-locked components under the preregistered ±100 µV artifact threshold.** The primary analysis used the pipeline's ±200 µV peak-to-peak rejection. For this sensitivity analysis the maximum peak-to-peak amplitude over all EEG channels and the full −0.5 to 1.5 s epoch was recomputed from the canonical single-trial epochs (the recomputation reproduced the canonical ±200 µV rejection set exactly and the canonical 400–600 ms amplitudes to < 0.0001 µV), and the final models of Table S13 were refitted on the epochs that also satisfy the registered ±100 µV criterion. Fair − unfair contrasts are evaluated at the mean-centered baseline. Source: `draft/sensitivity/erp-registered-E2_2026-08-23/`.", "",
       "| Component | Epoch set | *n* trials | Offer type | fair − unfair (µV) [95% CI] | Expression | Expression × offer |", "|---|---|---:|---|---|---|---|"]
for a in keep:
    o = om[om.analysis == a].set_index("term"); f = fu[fu.analysis == a].iloc[0]; n = int(dg[dg.analysis == a].n_trials.iloc[0])
    comp = "P3–LPP 400–600 ms" if a.startswith("LPP") else "Early fronto-central 192–244 ms"
    es = "±200 µV (primary)" if "canonical" in a else "±100 µV (registered)"
    def cell(term):
        x = o.loc[term]; return f"*F*({int(x.df1)}, {x.df2:.2f}) = {x.F:.2f}, *p* {p(x.p)}"
    s16.append(f"| {comp} | {es} | {n} | {cell('offer_type')} | {f.estimate:+.3f} [{f['lower.CL']:.3f}, {f['upper.CL']:.3f}] | {cell('emotion')} | {cell('emotion:offer_type')} |")
s16 = "\n".join(s16) + "\n"

j = text.index("## Supplementary figures")
text = text[:j] + s12 + "\n" + s13 + "\n" + s14 + "\n" + s15 + "\n" + s16 + "\n" + text[j:]

# ---------- Table of contents ----------
toc = ["## Contents", ""]
for m in re.finditer(r"^(\*\*Supplementary (Table|Figure) (S\d+[ab]?)\. ([^*]+?)\*\*)", text, re.M):
    toc.append(f"- Supplementary {m.group(2)} {m.group(3)}. {m.group(4).strip().rstrip('.')}")
toc_block = "\n".join(["- Supplementary methods: Experiment 1 face-locked ERP analysis"] + toc[2:])
text = text.replace("# Supplementary Materials\n", "# Supplementary Materials\n\n## Contents\n\n" + toc_block + "\n", 1)
open(OUT, "w", encoding="utf-8").write(text)
print("written", OUT, len(text.split()), "words")
