import FentonProduction
import Mathlib.Analysis.Real.Pi.Bounds

/- Real-arithmetic enclosure rules and the finite Fourier output contract.
   Runtime Decimal operations and the Taylor evaluator are separately audited. -/
namespace FiniteSystem
open scoped BigOperators
noncomputable section

def enclosed (lo hi x : ℝ) : Prop := lo ≤ x ∧ x ≤ hi

theorem pi_enclosed : enclosed 3.14159265358979323846 3.14159265358979323847 Real.pi :=
  ⟨le_of_lt Real.pi_gt_d20, le_of_lt Real.pi_lt_d20⟩

theorem exact_enclosed (x : ℝ) : enclosed x x x := ⟨le_rfl, le_rfl⟩

theorem enclosure_order (lo hi x : ℝ) (h : enclosed lo hi x) : lo ≤ hi :=
  le_trans h.1 h.2

theorem widen (a b lo hi x : ℝ) (h : enclosed a b x) (hl : lo ≤ a) (hu : b ≤ hi) :
    enclosed lo hi x := ⟨le_trans hl h.1, le_trans h.2 hu⟩

theorem add_enclosed (a b c d x y : ℝ) (hx : enclosed a b x) (hy : enclosed c d y) :
    enclosed (a+c) (b+d) (x+y) := ⟨add_le_add hx.1 hy.1, add_le_add hx.2 hy.2⟩

theorem neg_enclosed (a b x : ℝ) (h : enclosed a b x) : enclosed (-b) (-a) (-x) :=
  ⟨neg_le_neg h.2, neg_le_neg h.1⟩

theorem sub_enclosed (a b c d x y : ℝ) (hx : enclosed a b x) (hy : enclosed c d y) :
    enclosed (a-d) (b-c) (x-y) := ⟨sub_le_sub hx.1 hy.2, sub_le_sub hx.2 hy.1⟩

theorem mul_nonnegative_enclosed (a b c d x y : ℝ)
    (ha : 0 ≤ a) (hc : 0 ≤ c) (hx : enclosed a b x) (hy : enclosed c d y) :
    enclosed (a*c) (b*d) (x*y) := by
  constructor
  · exact mul_le_mul hx.1 hy.1 hc (le_trans ha hx.1)
  · exact mul_le_mul hx.2 hy.2 (le_trans hc hy.1) (le_trans ha (le_trans hx.1 hx.2))

theorem mul_fixed_enclosed (a b x c lo hi : ℝ) (hx : enclosed a b x)
    (ha : enclosed lo hi (a*c)) (hb : enclosed lo hi (b*c)) : enclosed lo hi (x*c) := by
  rcases le_total 0 c with hc | hc
  · exact ⟨le_trans ha.1 (mul_le_mul_of_nonneg_right hx.1 hc),
           le_trans (mul_le_mul_of_nonneg_right hx.2 hc) hb.2⟩
  · exact ⟨le_trans hb.1 (mul_le_mul_of_nonpos_right hx.2 hc),
           le_trans (mul_le_mul_of_nonpos_right hx.1 hc) ha.2⟩

theorem mul_corners_enclosed (a b c d x y lo hi : ℝ)
    (hx : enclosed a b x) (hy : enclosed c d y)
    (hac : enclosed lo hi (a*c)) (had : enclosed lo hi (a*d))
    (hbc : enclosed lo hi (b*c)) (hbd : enclosed lo hi (b*d)) :
    enclosed lo hi (x*y) := by
  have ha : enclosed lo hi (a*y) := by
    simpa [mul_comm] using mul_fixed_enclosed c d y a lo hi hy
      (by simpa [mul_comm] using hac) (by simpa [mul_comm] using had)
  have hb : enclosed lo hi (b*y) := by
    simpa [mul_comm] using mul_fixed_enclosed c d y b lo hi hy
      (by simpa [mul_comm] using hbc) (by simpa [mul_comm] using hbd)
  exact mul_fixed_enclosed a b x y lo hi hx ha hb

theorem exp_enclosed (a b x : ℝ) (h : enclosed a b x) :
    enclosed (Real.exp a) (Real.exp b) (Real.exp x) :=
  ⟨Real.exp_le_exp.mpr h.1, Real.exp_le_exp.mpr h.2⟩

theorem sqrt_enclosed (a b x : ℝ) (h : enclosed a b x) :
    enclosed (Real.sqrt a) (Real.sqrt b) (Real.sqrt x) :=
  ⟨Real.sqrt_le_sqrt h.1, Real.sqrt_le_sqrt h.2⟩

theorem square_nonnegative (x : ℝ) : 0 ≤ x^2 := sq_nonneg x

theorem square_enclosed (x radius : ℝ) (h : |x| ≤ radius) :
    enclosed 0 (radius^2) (x^2) := by
  constructor
  · exact sq_nonneg x
  · have hr : 0 ≤ radius := le_trans (abs_nonneg x) h
    have hx := (abs_le.mp h)
    nlinarith [sq_nonneg (radius-x), sq_nonneg (radius+x)]

theorem residual_from_enclosure (lo hi r tolerance : ℝ)
    (h : enclosed lo hi r) (hl : -tolerance ≤ lo) (hu : hi ≤ tolerance) :
    |r| ≤ tolerance := abs_le.mpr ⟨le_trans hl h.1, le_trans h.2 hu⟩

theorem output_from_enclosure (lo hi x value error : ℝ)
    (h : enclosed lo hi x) (hl : value-error ≤ lo) (hu : hi ≤ value+error) :
    |value-x| ≤ error := by
  apply abs_le.mpr
  constructor <;> linarith [h.1, h.2]

theorem output_budget (lo hi x value error budget : ℝ)
    (h : enclosed lo hi x) (hl : value-error ≤ lo) (hu : hi ≤ value+error)
    (he : error ≤ budget) : |value-x| ≤ budget :=
  le_trans (output_from_enclosure lo hi x value error h hl hu) he

theorem scaled_output_budget (lo hi x value error budget denominator : ℝ)
    (hd : 0 < denominator) (h : enclosed (lo/denominator) (hi/denominator) x)
    (hl : value-error ≤ lo) (hu : hi ≤ value+error) (he : error ≤ budget) :
    |value/denominator-x| ≤ budget/denominator := by
  apply output_budget (lo/denominator) (hi/denominator) x (value/denominator)
    (error/denominator) (budget/denominator) h
  · simpa [sub_div] using div_le_div_of_nonneg_right hl (le_of_lt hd)
  · simpa [add_div] using div_le_div_of_nonneg_right hu (le_of_lt hd)
  · exact div_le_div_of_nonneg_right he (le_of_lt hd)

theorem scaled_residual_budget (lo hi residual limit denominator : ℝ)
    (hd : 0 < denominator) (h : enclosed (lo/denominator) (hi/denominator) residual)
    (hl : -limit ≤ lo) (hu : hi ≤ limit) : |residual| ≤ limit/denominator := by
  apply residual_from_enclosure (lo/denominator) (hi/denominator) residual (limit/denominator) h
  · simpa only [neg_div] using div_le_div_of_nonneg_right hl (le_of_lt hd)
  · exact div_le_div_of_nonneg_right hu (le_of_lt hd)

theorem finite_residual_contract {n : ℕ} (lo hi residual : Fin n → ℝ) (tolerance : ℝ)
    (h : ∀ i, enclosed (lo i) (hi i) (residual i))
    (hl : ∀ i, -tolerance ≤ lo i) (hu : ∀ i, hi i ≤ tolerance) :
    ∀ i, |residual i| ≤ tolerance := fun i =>
  residual_from_enclosure (lo i) (hi i) (residual i) tolerance (h i) (hl i) (hu i)

theorem finite_equation_count (n : ℕ) : 8+(n+1)+(n+1) = 2*n+10 := by omega

def potential {n : ℕ} (coefficients : Fin n → ℝ) (terminal : ℝ) : ℝ :=
  (∑ j, coefficients j ^ 2)/4 + terminal^2/16

theorem potential_nonnegative {n : ℕ} (coefficients : Fin n → ℝ) (terminal : ℝ) :
    0 ≤ potential coefficients terminal := by
  unfold potential
  positivity

theorem terminal_half_amplitude (terminal : ℝ) : (terminal/2)^2/4 = terminal^2/16 := by ring

theorem potential_add_terminal {n : ℕ} (coefficients : Fin n → ℝ) (terminal : ℝ) :
    potential coefficients terminal = (∑ j, coefficients j^2)/4 + (terminal/2)^2/4 := by
  unfold potential
  ring

theorem dimensional_potential_nonnegative {n : ℕ} (coefficients : Fin n → ℝ)
    (terminal scale : ℝ) (hs : 0 ≤ scale) : 0 ≤ scale*potential coefficients terminal :=
  mul_nonneg hs (potential_nonnegative coefficients terminal)

theorem reconstructed_impulse (z1 z5 z8 density speed k : ℝ) :
    density*speed/k*(z8+z1*z5) = density*speed/k*z8 + density*speed/k*z1*z5 := by ring

theorem kinetic_flux_identity (c current flux pulse : ℝ) :
    (c*pulse-current*flux)/2 = c*pulse/2-current*flux/2 := by ring

theorem energy_enclosure (kl kh pl ph kinetic potential : ℝ)
    (hk : enclosed kl kh kinetic) (hp : enclosed pl ph potential) :
    enclosed (kl+pl) (kh+ph) (kinetic+potential) := add_enclosed _ _ _ _ _ _ hk hp

theorem surface_value_bound {n : ℕ} (a t : Fin n → ℝ) (h : ∀ j, |t j| ≤ 1) :
    |∑ j, a j*t j| ≤ ∑ j, |a j| := by
  calc
    |∑ j, a j*t j| ≤ ∑ j, |a j*t j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ j, |a j| := Finset.sum_le_sum (fun j _ => by
      rw [abs_mul]
      simpa using mul_le_mul_of_nonneg_left (h j) (abs_nonneg (a j)))

end
end FiniteSystem

#print axioms FiniteSystem.pi_enclosed
#print axioms FiniteSystem.exact_enclosed
#print axioms FiniteSystem.enclosure_order
#print axioms FiniteSystem.widen
#print axioms FiniteSystem.add_enclosed
#print axioms FiniteSystem.neg_enclosed
#print axioms FiniteSystem.sub_enclosed
#print axioms FiniteSystem.mul_nonnegative_enclosed
#print axioms FiniteSystem.mul_fixed_enclosed
#print axioms FiniteSystem.mul_corners_enclosed
#print axioms FiniteSystem.exp_enclosed
#print axioms FiniteSystem.sqrt_enclosed
#print axioms FiniteSystem.square_nonnegative
#print axioms FiniteSystem.square_enclosed
#print axioms FiniteSystem.residual_from_enclosure
#print axioms FiniteSystem.output_from_enclosure
#print axioms FiniteSystem.output_budget
#print axioms FiniteSystem.scaled_output_budget
#print axioms FiniteSystem.scaled_residual_budget
#print axioms FiniteSystem.finite_residual_contract
#print axioms FiniteSystem.finite_equation_count
#print axioms FiniteSystem.potential_nonnegative
#print axioms FiniteSystem.terminal_half_amplitude
#print axioms FiniteSystem.potential_add_terminal
#print axioms FiniteSystem.dimensional_potential_nonnegative
#print axioms FiniteSystem.reconstructed_impulse
#print axioms FiniteSystem.kinetic_flux_identity
#print axioms FiniteSystem.energy_enclosure
#print axioms FiniteSystem.surface_value_bound
