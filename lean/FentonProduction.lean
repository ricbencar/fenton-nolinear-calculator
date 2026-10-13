/-
# Production formula audit: Fenton Fourier solver

Source identities are pinned in sources.json. README sections 3-8 and
Fenton's published mathematical definitions are the specifications.

This file proves identities in real arithmetic for a manually reviewed
translation of the source equations. It does not formalise Python, NumPy,
Numba, IEEE-754, quadrature error, Newton convergence or uniqueness.
The attached executable is tested independently by audit.py.
Every reported check refers to the pinned source files in this package.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.Analysis.Complex.Trigonometric
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith



namespace FentonProduction
open scoped BigOperators
noncomputable section


def basisS (a b : ℝ) := Real.sinh (a+b) / Real.cosh b
def basisC (a b : ℝ) := Real.cosh (a+b) / Real.cosh b
def stableS (a b : ℝ) := Real.sinh a + Real.cosh a * Real.tanh b
def stableC (a b : ℝ) := Real.cosh a + Real.sinh a * Real.tanh b

theorem stableS_eq_basisS (a b : ℝ) : stableS a b = basisS a b := by
  unfold stableS basisS
  rw [Real.sinh_add, Real.tanh_eq_sinh_div_cosh]
  have hb : Real.cosh b ≠ 0 := ne_of_gt (Real.cosh_pos b)
  field_simp

theorem stableC_eq_basisC (a b : ℝ) : stableC a b = basisC a b := by
  unfold stableC basisC
  rw [Real.cosh_add, Real.tanh_eq_sinh_div_cosh]
  have hb : Real.cosh b ≠ 0 := ne_of_gt (Real.cosh_pos b)
  field_simp


def pySinh (a : ℝ) := (Real.exp a - 1 / Real.exp a) / 2
def pyCosh (a : ℝ) := (Real.exp a + 1 / Real.exp a) / 2

theorem pySinh_eq (a : ℝ) : pySinh a = Real.sinh a := by
  simp [pySinh, Real.sinh_eq, Real.exp_neg, one_div]

theorem pyCosh_eq (a : ℝ) : pyCosh a = Real.cosh a := by
  simp [pyCosh, Real.cosh_eq, Real.exp_neg, one_div]

theorem bed_streamline_basis (b : ℝ) : stableS (-b) b = 0 := by
  rw [stableS_eq_basisS]
  simp [basisS]



theorem exponential_branch (a : ℝ) :
    Real.sinh a + Real.cosh a * 1 = Real.exp a := by
  simp [Real.sinh_add_cosh]


def scaledS (a b : ℝ) := (Real.exp a-Real.exp (-a-2*b))/(1+Real.exp (-2*b))
def scaledC (a b : ℝ) := (Real.exp a+Real.exp (-a-2*b))/(1+Real.exp (-2*b))

theorem scaledS_matches_basis (a b : ℝ) : scaledS a b = basisS a b := by
  have e1 : Real.exp (-a-2*b) = (Real.exp a*(Real.exp b)^2)⁻¹ := by
    rw [show -a-2*b = -(a+(b+b)) by ring, Real.exp_neg, Real.exp_add, Real.exp_add]
    simp [pow_two]
  have e2 : Real.exp (-2*b) = ((Real.exp b)^2)⁻¹ := by
    rw [show -2*b = -(b+b) by ring, Real.exp_neg, Real.exp_add]
    simp [pow_two]
  unfold scaledS basisS
  rw [e1, e2, Real.sinh_eq, Real.cosh_eq]
  simp only [Real.exp_add, Real.exp_neg]
  have ha : Real.exp a ≠ 0 := ne_of_gt (Real.exp_pos a)
  have hb : Real.exp b ≠ 0 := ne_of_gt (Real.exp_pos b)
  field_simp

theorem scaledC_matches_basis (a b : ℝ) : scaledC a b = basisC a b := by
  have e1 : Real.exp (-a-2*b) = (Real.exp a*(Real.exp b)^2)⁻¹ := by
    rw [show -a-2*b = -(a+(b+b)) by ring, Real.exp_neg, Real.exp_add, Real.exp_add]
    simp [pow_two]
  have e2 : Real.exp (-2*b) = ((Real.exp b)^2)⁻¹ := by
    rw [show -2*b = -(b+b) by ring, Real.exp_neg, Real.exp_add]
    simp [pow_two]
  unfold scaledC basisC
  rw [e1, e2, Real.cosh_eq, Real.cosh_eq]
  simp only [Real.exp_add, Real.exp_neg]
  have ha : Real.exp a ≠ 0 := ne_of_gt (Real.exp_pos a)
  have hb : Real.exp b ≠ 0 := ne_of_gt (Real.exp_pos b)
  field_simp

theorem scaled_bed_condition (b : ℝ) : scaledS (-b) b = 0 := by
  rw [scaledS_matches_basis]
  simp [basisS]

/-! Elementary differential identities underlying the Fourier field. -/

theorem basisS_derivative (a b : ℝ) :
    HasDerivAt (fun x => basisS x b) (basisC a b) a := by
  simpa [basisS, basisC] using
    (((hasDerivAt_id a).add_const b).sinh.div_const (Real.cosh b))

theorem basisC_derivative (a b : ℝ) :
    HasDerivAt (fun x => basisC x b) (basisS a b) a := by
  simpa [basisS, basisC] using
    (((hasDerivAt_id a).add_const b).cosh.div_const (Real.cosh b))

theorem modeS_derivative (j y kd : ℝ) :
    HasDerivAt (fun t => basisS (j*t) (j*kd)) (j*basisC (j*y) (j*kd)) y := by
  simpa [Function.comp_def, mul_comm] using!
    (basisS_derivative (j*y) (j*kd)).comp y ((hasDerivAt_id y).const_mul j)

theorem modeC_derivative (j y kd : ℝ) :
    HasDerivAt (fun t => basisC (j*t) (j*kd)) (j*basisS (j*y) (j*kd)) y := by
  simpa [Function.comp_def, mul_comm] using!
    (basisC_derivative (j*y) (j*kd)).comp y ((hasDerivAt_id y).const_mul j)

theorem phaseCos_derivative (j x : ℝ) :
    HasDerivAt (fun t => Real.cos (j*t)) (-j*Real.sin (j*x)) x := by
  simpa [mul_comm] using ((hasDerivAt_id x).const_mul j).cos

theorem phaseSin_derivative (j x : ℝ) :
    HasDerivAt (fun t => Real.sin (j*t)) (j*Real.cos (j*x)) x := by
  simpa [mul_comm] using ((hasDerivAt_id x).const_mul j).sin

theorem vertical_second_derivative (j y kd : ℝ) :
    deriv (fun t => deriv (fun v => basisS (j*v) (j*kd)) t) y =
    j^2*basisS (j*y) (j*kd) := by
  have h : (fun t => deriv (fun v => basisS (j*v) (j*kd)) t) =
      (fun t => j*basisC (j*t) (j*kd)) := by
    funext t
    exact (modeS_derivative j t kd).deriv
  rw [h, ((modeC_derivative j y kd).const_mul j).deriv]
  ring

theorem horizontal_second_derivative (j x : ℝ) :
    deriv (fun t => deriv (fun v => Real.cos (j*v)) t) x =
    -j^2*Real.cos (j*x) := by
  have h : (fun t => deriv (fun v => Real.cos (j*v)) t) =
      (fun t => -j*Real.sin (j*t)) := by
    funext t
    exact (phaseCos_derivative j t).deriv
  rw [h, ((phaseSin_derivative j x).const_mul (-j)).deriv]
  ring

theorem separated_mode_laplace (j x y kd B : ℝ) :
    B*basisS (j*y) (j*kd)*
      deriv (fun t => deriv (fun v => Real.cos (j*v)) t) x +
    B*Real.cos (j*x)*
      deriv (fun t => deriv (fun v => basisS (j*v) (j*kd)) t) y = 0 := by
  rw [horizontal_second_derivative, vertical_second_derivative]
  ring



def psiReadme (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    B j * basisS ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.cos ((j:ℝ)*X)

def psiPython (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    B j * scaledS ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.cos ((j:ℝ)*X)

def uReadme (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    (j:ℝ) * B j * basisC ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.cos ((j:ℝ)*X)

def uPython (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    (j:ℝ) * B j * scaledC ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.cos ((j:ℝ)*X)

def vReadme (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    (j:ℝ) * B j * basisS ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.sin ((j:ℝ)*X)

def vPython (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :=
  ∑ j ∈ Finset.Icc 1 N,
    (j:ℝ) * B j * scaledS ((j:ℝ)*Y) ((j:ℝ)*kd) * Real.sin ((j:ℝ)*X)

theorem fourier_psi_equiv (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :
    psiPython N B X Y kd = psiReadme N B X Y kd := by
  simp only [psiPython, psiReadme, scaledS_matches_basis]

theorem fourier_u_equiv (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :
    uPython N B X Y kd = uReadme N B X Y kd := by
  simp only [uPython, uReadme, scaledC_matches_basis]

theorem fourier_v_equiv (N : ℕ) (B : ℕ → ℝ) (X Y kd : ℝ) :
    vPython N B X Y kd = vReadme N B X Y kd := by
  simp only [vPython, vReadme, scaledS_matches_basis]

theorem fourier_bed (N : ℕ) (B : ℕ → ℝ) (X kd : ℝ) :
    psiReadme N B X (-kd) kd = 0 ∧ vReadme N B X (-kd) kd = 0 := by
  simp [psiReadme, vReadme, basisS]


structure Parameters where
  N : ℕ
  Hoverd : ℝ
  height : ℝ
  currentDepth : ℝ 
  currentCriterion : ℕ 

def globalReadme (p : Parameters) (z : ℕ → ℝ) : Fin 8 → ℝ := fun i =>
  match i.val with
  | 0 => z 2 - z 1 * p.Hoverd
  | 1 => z 2 - p.height * (z 3)^2
  | 2 => z 4 * z 3 - 2*Real.pi
  | 3 => z 5 + z 7 - z 4
  | 4 => z 1 * (z 6 + z 7 - z 4) - z 8
  | 5 => z (p.currentCriterion+4) - p.currentDepth * Real.sqrt (z 1)
  | 6 => z 10 + z (p.N+10) + 2 * ∑ i ∈ Finset.Icc 1 (p.N-1), z (10+i)
  | _ => z 10 - z (p.N+10) - z 2


def globalPython (p : Parameters) (z : ℕ → ℝ) : Fin 8 → ℝ := fun i =>
  match i.val with
  | 0 => z 2 - z 1 * p.Hoverd
  | 1 => z 2 - p.height * z 3 * z 3
  | 2 => z 4 * z 3 - 2*Real.pi
  | 3 => z 5 + z 7 - z 4
  | 4 => z 1 * (z 6 + z 7 - z 4) - z 8
  | 5 => z (p.currentCriterion+4) - p.currentDepth * Real.sqrt (z 1)
  | 6 => z 10 + z (p.N+10) + ∑ i ∈ Finset.Icc 1 (p.N-1), 2 * z (10+i)
  | _ => z 10 - z (p.N+10) - z 2

theorem global_residuals_equiv (p : Parameters) (z : ℕ → ℝ) :
    globalPython p z = globalReadme p z := by
  funext i
  simp only [globalPython, globalReadme, pow_two, mul_assoc, Finset.mul_sum]

def streamlineReadme (N : ℕ) (B : ℕ → ℝ) (X Y kd U q : ℝ) :=
  psiReadme N B X Y kd - U*Y - q
def streamlinePython (N : ℕ) (B : ℕ → ℝ) (X Y kd U q : ℝ) :=
  psiPython N B X Y kd - q - U*Y

theorem streamline_residual_equiv (N : ℕ) (B : ℕ → ℝ) (X Y kd U q : ℝ) :
    streamlinePython N B X Y kd U q = streamlineReadme N B X Y kd U q := by
  unfold streamlinePython streamlineReadme
  rw [fourier_psi_equiv]
  ring

def bernoulliReadme (N : ℕ) (B : ℕ → ℝ) (X Y kd U r : ℝ) :=
  ((uReadme N B X Y kd-U)^2 + (vReadme N B X Y kd)^2)/2 + Y-r
def bernoulliPython (N : ℕ) (B : ℕ → ℝ) (X Y kd U r : ℝ) :=
  (1/2:ℝ)*((-U+uPython N B X Y kd)^2 +
    vPython N B X Y kd * vPython N B X Y kd)+Y-r

theorem bernoulli_residual_equiv (N : ℕ) (B : ℕ → ℝ) (X Y kd U r : ℝ) :
    bernoulliPython N B X Y kd U r = bernoulliReadme N B X Y kd U r := by
  unfold bernoulliPython bernoulliReadme
  rw [fourier_u_equiv, fourier_v_equiv]
  ring

theorem equation_count (N : ℕ) : 8+(N+1)+(N+1)=2*N+10 := by omega
theorem default_equation_count : 8+(50+1)+(50+1)=110 := by norm_num


theorem eulerian_closure (ce U c : ℝ) (h : ce+U-c=0) : ce=c-U := by linarith

theorem stokes_closure (kd cs U c q : ℝ) (hk : kd ≠ 0)
    (h : kd*(cs+U-c)-q=0) : cs=c-(U*kd-q)/kd := by
  field_simp
  nlinarith [h]

theorem wavelength_conversion (kd d : ℝ) (hk : kd ≠ 0) (hd : d ≠ 0) :
    2*Real.pi/(kd/d) = d*(2*Real.pi/kd) := by field_simp

theorem crest_trough (zc zt k H : ℝ) (hk : k ≠ 0)
    (h : zc-zt-k*H=0) : zc/k-zt/k=H := by
  field_simp
  nlinarith [h]


theorem material_acceleration (u v c ux uy : ℝ) :
    (-c*ux)+u*ux+v*uy=(u-c)*ux+v*uy := by ring

/-! Integral properties: README section 8; Fenton (1990), section 7.
The identities below check algebraic translation, not the derivation of the
integral identities from the continuous free-boundary problem. -/

def potentialCode (lower terminal : ℝ) := lower + terminal^2/16
def potentialAmplitude (lower terminal : ℝ) := lower + (terminal/2)^2/4

theorem potential_matches_surface (lower terminal : ℝ) :
    potentialCode lower terminal = potentialAmplitude lower terminal := by
  unfold potentialCode potentialAmplitude
  ring

theorem terminal_mode_energy : potentialCode 0 1 = 1/16 := by
  norm_num [potentialCode]

def betaAbsolute (c ce r : ℝ) := 2*r-c^2+2*c*ce
def stressCode (c ce kd r impulse K P : ℝ) :=
  4*K-3*P+betaAbsolute c ce r*kd-2*ce*impulse
def powerCode (c ce kd r impulse K P : ℝ) :=
  c*(3*K-2*P)+0.5*betaAbsolute c ce r*(impulse+c*kd)-2*c*ce*impulse
def stressFenton (ce kd beta impulse K P : ℝ) := 4*K-3*P+beta*kd-2*ce*impulse
def powerFenton (c ce kd beta impulse K P : ℝ) :=
  c*(3*K-2*P)+beta*(impulse+c*kd)/2-2*c*ce*impulse

theorem beta_matches_Fenton (c ce r : ℝ) :
    betaAbsolute c ce r = 2*r-c*(c-2*ce) := by unfold betaAbsolute; ring

theorem stress_matches_Fenton (c ce kd r impulse K P : ℝ) :
    stressCode c ce kd r impulse K P =
    stressFenton ce kd (2*r-c*(c-2*ce)) impulse K P := by
  unfold stressCode stressFenton betaAbsolute
  ring

theorem power_matches_Fenton (c ce kd r impulse K P : ℝ) :
    powerCode c ce kd r impulse K P =
    powerFenton c ce kd (2*r-c*(c-2*ce)) impulse K P := by
  unfold powerCode powerFenton betaAbsolute
  ring

theorem impulse_flux_identity (kd c ce U q : ℝ) (h : ce+U-c=0) :
    q+kd*ce = c*kd-(U*kd-q) := by
  have hc : ce = c-U := by linarith
  rw [hc]
  ring

theorem kinetic_flux_form (c ce kd U q : ℝ) :
    (c*(q+kd*ce)-ce*(U*kd-q))/2 =
    (c*(q+kd*ce)+ce*(q-U*kd))/2 := by ring

theorem stokes_equals_mean_transport (c Q d rho I : ℝ)
    (hd : d ≠ 0) (hr : rho ≠ 0) (h : I=rho*(c*d-Q)) :
    I/(rho*d) = c-Q/d := by rw [h]; field_simp

def dimensionalStress (rho d ce beta I K P : ℝ) := 4*K-3*P+rho*beta*d-2*ce*I
def dimensionalPower (rho d c ce beta I K P : ℝ) :=
  c*(3*K-2*P)+beta*(I+rho*c*d)/2-2*c*ce*I

theorem uniform_current_stress (rho d U : ℝ) :
    dimensionalStress rho d U (U^2) (rho*d*U) (rho*d*U^2/2) 0 = rho*d*U^2 := by
  unfold dimensionalStress; ring

theorem uniform_current_power (rho d c U : ℝ) :
    dimensionalPower rho d c U (U^2) (rho*d*U) (rho*d*U^2/2) 0 = rho*d*U^3/2 := by
  unfold dimensionalPower; ring

theorem absolute_orbital_relation (mean orbital : ℝ) :
    (mean^2+orbital)-mean^2=orbital := by ring

theorem bernoulli_coordinate_shift (g d eta u v c r : ℝ) :
    (r+g*d)-g*(d+eta)-((u-c)^2+v^2)/2 =
    r-g*eta-((u-c)^2+v^2)/2 := by ring

theorem momentum_frame_change (u c rho p : ℝ) :
    p+rho*(u-c)^2 = (p+rho*u^2)-2*c*(rho*u)+rho*c^2 := by ring

theorem pressure_zero_from_dynamic (u v c eta r g : ℝ)
    (h : ((u-c)^2+v^2)/2+g*eta-r=0) :
    r-g*eta-((u-c)^2+v^2)/2=0 := by linarith

/-! Dimensional conversion expressed using positive length and velocity
scales, so no floating exponent or square-root identity is assumed. -/

theorem stress_scaling (rho g d kd c ce r I K P : ℝ) (hk : kd ≠ 0) :
    rho*g*d^2*(stressCode c ce kd r I K P/kd^2) =
    4*(rho*g*d^2*K/kd^2)-3*(rho*g*d^2*P/kd^2)+
    rho*g*d^2*betaAbsolute c ce r/kd-2*rho*g*d^2*ce*I/kd^2 := by
  unfold stressCode
  field_simp

theorem power_unit_scaling (rho length speed c ce beta I K P : ℝ) :
    dimensionalPower rho length (speed*c) (speed*ce) (speed^2*beta)
      (rho*length*speed*I) (rho*length*speed^2*K) (rho*length*speed^2*P) =
    rho*length*speed^3*powerFenton c ce 1 beta I K P := by
  unfold dimensionalPower powerFenton
  ring

theorem stress_unit_scaling (rho length speed ce beta I K P : ℝ) :
    dimensionalStress rho length (speed*ce) (speed^2*beta)
      (rho*length*speed*I) (rho*length*speed^2*K) (rho*length*speed^2*P) =
    rho*length*speed^2*stressFenton ce 1 beta I K P := by
  unfold dimensionalStress stressFenton
  ring

theorem wavelength_period_relation (L T : ℝ) (hT : T ≠ 0) :
    (L/T)*T=L := by field_simp

theorem acceleration_scaling (u v c ux vx speed length g : ℝ)
    (hl : length ≠ 0) (hs : speed^2=g*length) :
    (speed*u-speed*c)*(speed/length*ux)+speed*v*(speed/length*vx) =
    g*((u-c)*ux+v*vx) := by
  field_simp
  rw [hs]
  ring

/-! Jacobian algebra. These proofs check the quotient-rule and chain-rule
expressions once elementary derivatives are supplied. The full assembled
matrix is independently checked against high-precision differences. -/

theorem basis_depth_quotient_S (j ep em e : ℝ) (_h : 1+e ≠ 0) :
    ((2*j*em)*(1+e)-(ep-em)*(-2*j*e))/(1+e)^2 =
    2*j*(em+e*ep)/(1+e)^2 := by ring

theorem basis_depth_quotient_C (j ep em e : ℝ) (_h : 1+e ≠ 0) :
    ((-2*j*em)*(1+e)-(ep+em)*(-2*j*e))/(1+e)^2 =
    2*j*(e*ep-em)/(1+e)^2 := by ring

theorem bernoulli_directional_derivative (u U v du dU dv deta dr : ℝ) :
    ((2*(u-U)*(du-dU)+2*v*dv)/2+deta-dr) =
    (u-U)*(du-dU)+v*dv+deta-dr := by ring

theorem streamline_directional_derivative (dpsi U eta dU deta dq : ℝ) :
    dpsi-(dU*eta+U*deta)-dq=dpsi-eta*dU-U*deta-dq := by ring

theorem column_equilibration (a scale dx : ℝ) (hs : scale ≠ 0) :
    (a/scale)*(scale*dx)=a*dx := by field_simp

/-! Necessary acceptance conditions; no assertion of existence, uniqueness
or convergence for every input is made. -/

def accepted (N : ℕ) (it : ℕ) (err crit kd tol : ℝ) (residual : Fin N → ℝ) :=
  it>1 ∧ err<crit*|kd| ∧ ∀ i, |residual i|≤tol

theorem acceptance_requires_residuals (N it : ℕ) (err crit kd tol : ℝ)
    (residual : Fin N → ℝ) (h : accepted N it err crit kd tol residual) :
    ∀ i, |residual i|≤tol := h.2.2

theorem stagnation_is_not_sufficient :
    ¬ accepted 1 2 0 (1/10^10) 1 (1/10^10) (fun _ => 1) := by
  intro h
  have hh := h.2.2 (0 : Fin 1)
  norm_num at hh

theorem postprocess_requires_residual_check (N : ℕ) (tol : ℝ)
    (after : Fin N → ℝ) (i : Fin N) (hi : tol < |after i|) :
    ¬ (∀ j, |after j| ≤ tol) := by
  intro h
  exact (not_le_of_gt hi) (h i)

theorem final_fraction_is_target (height fraction : ℝ) (h : fraction=1) :
    fraction*height=height := by rw [h]; ring

theorem residual_bound_not_exact_equality :
    ∃ x : ℝ, |x|≤(1/10^10) ∧ x≠0 := by
  refine ⟨1/10^11, ?_⟩
  norm_num

theorem power_energy_ratio (F E : ℝ) (h : E≠0) : (F/E)*E=F := by field_simp

/-! This mathematical ratio alone does not establish a nonlinear wave-packet
group speed; the report keeps that interpretation separate. -/

end
end FentonProduction


/- Kernel dependency audit for every theorem. -/
#print axioms FentonProduction.stableS_eq_basisS
#print axioms FentonProduction.stableC_eq_basisC
#print axioms FentonProduction.pySinh_eq
#print axioms FentonProduction.pyCosh_eq
#print axioms FentonProduction.bed_streamline_basis
#print axioms FentonProduction.exponential_branch
#print axioms FentonProduction.scaledS_matches_basis
#print axioms FentonProduction.scaledC_matches_basis
#print axioms FentonProduction.scaled_bed_condition
#print axioms FentonProduction.basisS_derivative
#print axioms FentonProduction.basisC_derivative
#print axioms FentonProduction.modeS_derivative
#print axioms FentonProduction.modeC_derivative
#print axioms FentonProduction.phaseCos_derivative
#print axioms FentonProduction.phaseSin_derivative
#print axioms FentonProduction.vertical_second_derivative
#print axioms FentonProduction.horizontal_second_derivative
#print axioms FentonProduction.separated_mode_laplace
#print axioms FentonProduction.fourier_psi_equiv
#print axioms FentonProduction.fourier_u_equiv
#print axioms FentonProduction.fourier_v_equiv
#print axioms FentonProduction.fourier_bed
#print axioms FentonProduction.global_residuals_equiv
#print axioms FentonProduction.streamline_residual_equiv
#print axioms FentonProduction.bernoulli_residual_equiv
#print axioms FentonProduction.equation_count
#print axioms FentonProduction.default_equation_count
#print axioms FentonProduction.eulerian_closure
#print axioms FentonProduction.stokes_closure
#print axioms FentonProduction.wavelength_conversion
#print axioms FentonProduction.crest_trough
#print axioms FentonProduction.material_acceleration
#print axioms FentonProduction.potential_matches_surface
#print axioms FentonProduction.terminal_mode_energy
#print axioms FentonProduction.beta_matches_Fenton
#print axioms FentonProduction.stress_matches_Fenton
#print axioms FentonProduction.power_matches_Fenton
#print axioms FentonProduction.impulse_flux_identity
#print axioms FentonProduction.kinetic_flux_form
#print axioms FentonProduction.stokes_equals_mean_transport
#print axioms FentonProduction.uniform_current_stress
#print axioms FentonProduction.uniform_current_power
#print axioms FentonProduction.absolute_orbital_relation
#print axioms FentonProduction.bernoulli_coordinate_shift
#print axioms FentonProduction.momentum_frame_change
#print axioms FentonProduction.pressure_zero_from_dynamic
#print axioms FentonProduction.stress_scaling
#print axioms FentonProduction.power_unit_scaling
#print axioms FentonProduction.stress_unit_scaling
#print axioms FentonProduction.wavelength_period_relation
#print axioms FentonProduction.acceleration_scaling
#print axioms FentonProduction.basis_depth_quotient_S
#print axioms FentonProduction.basis_depth_quotient_C
#print axioms FentonProduction.bernoulli_directional_derivative
#print axioms FentonProduction.streamline_directional_derivative
#print axioms FentonProduction.column_equilibration
#print axioms FentonProduction.acceptance_requires_residuals
#print axioms FentonProduction.stagnation_is_not_sufficient
#print axioms FentonProduction.postprocess_requires_residual_check
#print axioms FentonProduction.final_fraction_is_target
#print axioms FentonProduction.residual_bound_not_exact_equality
#print axioms FentonProduction.power_energy_ratio
