import FentonProduction

/-!
Acceptance contracts over mathematical real numbers. The Python branch
correspondence is tested by audit.py; it is not a formal Python semantics.
The maximum below is over evaluated samples/candidates, not continuous phase.
-/
namespace ProductionAcceptance
noncomputable section

def budget (reference : ℝ) : ℝ := 1/10^7 + |reference|/10^8
def withinBudget (value reference : ℝ) : Prop := |value-reference| ≤ budget reference

theorem budget_positive (reference : ℝ) : 0 < budget reference := by
  unfold budget
  positivity

theorem budget_rejects_excess (value reference : ℝ)
    (h : budget reference < |value-reference|) : ¬ withinBudget value reference := by
  exact not_le_of_gt h

def verified (residual condition surface : ℝ) (ratios : List ℝ) : Prop :=
  residual ≤ 1/10^10 ∧ condition ≤ 10^11 ∧ surface ≤ 1/10^6 ∧
  ∀ ratio ∈ ratios, ratio ≤ 1

theorem verification_requires_residual (r c s : ℝ) (ratios : List ℝ)
    (h : verified r c s ratios) : r ≤ 1/10^10 := h.1

theorem verification_requires_condition (r c s : ℝ) (ratios : List ℝ)
    (h : verified r c s ratios) : c ≤ 10^11 := h.2.1

theorem verification_requires_surface (r c s : ℝ) (ratios : List ℝ)
    (h : verified r c s ratios) : s ≤ 1/10^6 := h.2.2.1

theorem verification_requires_every_comparison (r c s ratio : ℝ) (ratios : List ℝ)
    (h : verified r c s ratios) (hm : ratio ∈ ratios) : ratio ≤ 1 := h.2.2.2 ratio hm

theorem failed_comparison_rejects (r c s ratio : ℝ) (ratios : List ℝ)
    (hm : ratio ∈ ratios) (hf : 1 < ratio) : ¬ verified r c s ratios := by
  intro h
  exact (not_le_of_gt hf) (verification_requires_every_comparison r c s ratio ratios h hm)

theorem small_collocation_residual_does_not_override_resolution :
    ¬ verified 0 100 (1/10^5) [] := by
  norm_num [verified]

def evaluatedPeak : List ℝ → ℝ
  | [] => 0
  | x::xs => max |x| (evaluatedPeak xs)

theorem evaluated_peak_nonnegative (xs : List ℝ) : 0 ≤ evaluatedPeak xs := by
  cases xs with
  | nil => simp [evaluatedPeak]
  | cons x xs => exact le_trans (abs_nonneg x) (le_max_left _ _)

theorem evaluated_peak_dominates (xs : List ℝ) (x : ℝ) (h : x ∈ xs) :
    |x| ≤ evaluatedPeak xs := by
  induction xs with
  | nil => simp at h
  | cons y ys ih =>
    simp only [List.mem_cons] at h
    rcases h with h | h
    · subst x
      exact le_max_left _ _
    · exact le_trans (ih h) (le_max_right _ _)

theorem refining_peak_cannot_decrease (best candidate : ℝ) :
    best ≤ max best |candidate| := le_max_left _ _

theorem refined_peak_dominates_candidate (best candidate : ℝ) :
    |candidate| ≤ max best |candidate| := le_max_right _ _

theorem evaluated_peak_is_not_signed_crest :
    evaluatedPeak [(-1/10), (-9/10)] = (9/10 : ℝ) := by
  norm_num [evaluatedPeak]

end
end ProductionAcceptance

#print axioms ProductionAcceptance.budget_positive
#print axioms ProductionAcceptance.budget_rejects_excess
#print axioms ProductionAcceptance.verification_requires_residual
#print axioms ProductionAcceptance.verification_requires_condition
#print axioms ProductionAcceptance.verification_requires_surface
#print axioms ProductionAcceptance.verification_requires_every_comparison
#print axioms ProductionAcceptance.failed_comparison_rejects
#print axioms ProductionAcceptance.small_collocation_residual_does_not_override_resolution
#print axioms ProductionAcceptance.evaluated_peak_nonnegative
#print axioms ProductionAcceptance.evaluated_peak_dominates
#print axioms ProductionAcceptance.refining_peak_cannot_decrease
#print axioms ProductionAcceptance.refined_peak_dominates_candidate
#print axioms ProductionAcceptance.evaluated_peak_is_not_signed_crest
