# Actor-preserving random-structure sensitivity

Recorded 2026-09-29 before executing the new structure selection. Jin approved
the proposed second stage after unchanged-formula numerical retries. This is
retrospective exploratory sensitivity analysis, not a new primary analysis.

## Fixed scientific and reporting contract

Use verified fair-only trial exports from EqualOffer_Aligned_20260928_v2, matching
canonical trial inputs. Preserve reaction coding, RT limits, IDs, factor levels
and emotion * allocation fixed effects. allocation is 5:5 versus 6:4. Preserve
the original primary fits, figures, BH families and primary decomposition gate.
Actor adjustment means baseline actor variation only, not actor-specific emotion
slopes or a guarantee of stimulus-population generalization.

## Selection rule, set before new output

Re-use the original manuscript five-step engine in a separately versioned copy.
Adapt only actor preservation and diagnostic logging:

1. Include (1 | actor_id) in the maximal and zero-correlation models and force it
   in buildmer along with all design fixed effects and the participant intercept.
2. Keep emotion, allocation and emotion-by-allocation participant slope columns
   grouped by factor effect, using the original sum coding. Apply the original
   buildmer order/backward grouped LRT criterion, alpha=.20. No hand deletion of
   individual contrast-coded columns or fixed effects.
3. Retain actor intercept in every correlation candidate; use the original LRT
   alpha=.20 and AIC rule against the valid zero-correlation reference. Do not
   choose a correlated model when the designated reference is invalid.
4. Use the previously recorded extended optimizer budgets (BFGS maxit=2000;
   nlminb iter.max=2000, eval.max=4000) for direct fits/final validation, and pass
   the BFGS control explicitly into buildmer. Record this optimizer-control
   adaptation; it is not an identical rerun of the original default controls.
5. The original numerical gate is unchanged. Report covariance eigen-directions
   separately because its individual-variance rule does not establish full rank.
   If the selection/reference/final gate fails, retain all logs and stop that
   experiment's dependent comparison/deletion stages. No unrecorded fallback.

## Matched comparisons and influence

For an accepted actor model, remove only the actor intercept to create the matched
no-actor sensitivity model; do not reselect its participant structure. Confirm
both models have identical fixed terms and participant random terms and trial rows.
Interpret the paired comparison only when both fits pass the original gate.
Also retain the original primary fit as a separate reference, making clear that
its participant structure may differ from the newly selected sensitivity model.

Then delete each participant in turn from the accepted ACTOR model only, keeping
that full-sample-selected structure fixed. Use both extended-budget optimizers
and the original acceptance gate, preserving failures. The matched no-actor
model is a full-sample comparison; its deletion sequence is not added here.

Export all marginal emotion pairs (BH within experiment/model, ten tests), the
single marginal allocation comparison, and Type III omnibus tests. Marginal means
average equally on the log-odds scale; Wald intervals are unadjusted. Sensitivity
p-values cannot unlock the original primary decomposition gate. Compare estimates,
ORs and uncertainty rather than selecting by significance. No actor variance
significance test or new simple-effect family is introduced.

## Verification and bounds

Verify source hashes, actor/fixed-term preservation and whole slope blocks.
Cross-check contrast arithmetic and BH correction; independently reconstruct
canonical data and refit accepted full-sample models. Archive engine adaptation
diff, stage models, selection diagnostics, optimizer attempts, code and sessionInfo.
No primary/result overwrite, raw edit, installation or manuscript/theory change.
One selection per experiment; stage timeout 1800 seconds and original maximal-fit
timeout; 300-second guards on subsequent fits, finite iterations, one thread.
Native routines may return control late; monitor logs and process activity.
