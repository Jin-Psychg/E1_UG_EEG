"""Create a minimal, auditable actor-preserving copy of the manuscript engine."""
from pathlib import Path
import difflib
import hashlib
import json

code = Path(__file__).resolve().parent
source = code.parent / "EqualOffer_ManuscriptAligned" / "manuscript_engine.R"
original = source.read_text(encoding="utf-8")
text = original
changes = [
    ('fit_stage1_bates <- function(', 'fit_stage1_actor <- function('),
    ('                             timeout_sec = 1800) {',
     '                             timeout_sec = 1800, actor_audit_dir = NULL) {'),
    ('  fixed_str <- sprintf("%s ~ %s", outcome_var, fixed_terms)',
     '  stopifnot(model_family == "binomial", "actor_id" %in% names(data))\n'
     '  fixed_str <- sprintf("%s ~ %s + (1 | actor_id)", outcome_var, fixed_terms)'),
    ('  fixed_include_str <- paste(c(include_terms, sprintf("(1 | %s)", group_name)),',
     '  fixed_include_str <- paste(c(include_terms, sprintf("(1 | %s)", group_name), "(1 | actor_id)"),'),
    ('      ddf = "Wald"\n    )',
     '      ddf = "Wald",\n      args = list(control = custom_control_glmmTMB)\n    )'),
    ('  reduced_fixed_str <- paste(deparse(reduced_fixed), collapse = " ")',
     '  reduced_fixed_str <- paste(paste(deparse(reduced_fixed), collapse = " "), "+ (1 | actor_id)")'),
    ('  diag_log$step3_final_formula <- paste(deparse(formula(m_reduced)), collapse = " ")',
     '  if (!is.null(actor_audit_dir)) saveRDS(m_reduced, file.path(actor_audit_dir, "step3_selected.rds"))\n'
     '  diag_log$step3_final_formula <- paste(deparse(formula(m_reduced)), collapse = " ")'),
    ('  diag_log$step2_status <- if (is.null(step2_diag)) {',
     '  if (!is.null(actor_audit_dir) && !is.null(m_zcp)) saveRDS(m_zcp, file.path(actor_audit_dir, "step2_zcp.rds"))\n'
     '  diag_log$step2_status <- if (is.null(step2_diag)) {'),
    ('  diag_log$step4_status <- attr(m_final, "step4_selected") %||% "ZCP_kept"',
     '  if (!is.null(actor_audit_dir)) saveRDS(m_final, file.path(actor_audit_dir, "step4_selected.rds"))\n'
     '  diag_log$step4_status <- attr(m_final, "step4_selected") %||% "ZCP_kept"'),
]
for old, new in changes:
    if text.count(old) != 1:
        raise ValueError(f"Expected exactly one patch location: {old}")
    text = text.replace(old, new, 1)
dest = code / "actor_engine.R"
if dest.exists():
    raise FileExistsError("Preserve the existing engine")
dest.write_text(text, encoding="utf-8")
(code / "actor_engine.diff").write_text("".join(difflib.unified_diff(
    original.splitlines(True), text.splitlines(True), fromfile=str(source), tofile=str(dest))), encoding="utf-8")
(code / "engine_provenance.json").write_text(json.dumps({
    "source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "adapted_sha256": hashlib.sha256(dest.read_bytes()).hexdigest(), "patch_count": len(changes)}, indent=2))
print(f"Created actor engine with {len(changes)} explicitly checked adaptations")
