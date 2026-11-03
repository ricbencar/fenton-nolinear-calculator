# The Mathematics of Nonlinear Steady Water Waves

*A theoretical and computational treatise on the Rienecker–Fenton Fourier approximation, finite-amplitude dispersion, wave–current interaction, nonlinear free-surface boundary conditions, Newton–SVD collocation, kinematics, wave integrals, numerical verification, and the relationship with Stokes and cnoidal theories*

**Scientific scope:** two-dimensional, progressive, periodic, finite-amplitude gravity waves over a horizontal, impermeable bed; homogeneous, incompressible, inviscid, irrotational flow.  
**Principal numerical formulation:** Rienecker and Fenton (1981); Fenton (1988, 1990, 1999).  
**Reference implementation:** the single-file Rust solver `fenton_gui.rs`, with a numerical core independent of the Win32 graphical interface.  
**Units:** SI unless a nondimensional quantity is explicitly indicated.

***

## Abstract

The accurate description of a finite-amplitude periodic surface gravity wave is a nonlinear free-boundary problem: the fluid velocity field and the position of the water surface must satisfy Laplace's equation in an unknown domain, impermeability at the seabed, conservation of mass, and both the kinematic and dynamic free-surface boundary conditions. Linear Airy theory supplies invaluable physical intuition, but its free-surface linearization omits the very mechanisms that determine finite-amplitude wave shape, the nonlinear wavelength, asymmetry, near-crest kinematics and integral transport properties. Fifth-order Stokes theory and cnoidal theory recover progressively more of these effects in complementary depth regimes, yet their respective perturbation expansions are not uniformly accurate over the entire domain of steady periodic waves.

The Rienecker–Fenton method avoids a wave-height perturbation expansion. Instead, it represents a harmonic stream function by a truncated Fourier basis satisfying the field equation and bottom boundary condition identically; it determines the free surface, Fourier coefficients, wavelength, currents and Bernoulli constant by simultaneously enforcing the nonlinear conditions at selected surface collocation points. This paper derives that formulation from first principles, with particular attention to coordinate systems, current definitions, dimensional scaling, the exact structure of the finite nonlinear equations, analytic Jacobians, stabilized Newton iteration, height continuation and independent numerical checks. It also derives velocity, acceleration, pressure, energy, momentum, radiation stress and energy flux, places Stokes and cnoidal theories in their proper mathematical context, and separates accurate *finite-system evaluation* from the much stronger claim of a rigorous existence theorem for the infinite-dimensional free-boundary problem.

**Keywords:** steady gravity waves; stream function; nonlinear dispersion; Fourier collocation; Rienecker–Fenton; Stokes waves; cnoidal waves; Eulerian current; Stokes mass transport; Newton method; singular-value decomposition; Gauss–Legendre quadrature; finite precision; interval arithmetic; radiation stress; wave power.

***

## 1. Motivation, historical development and the mathematical problem

### 1.1 Why a nonlinear solution is required

Surface gravity waves are simultaneously dispersive and nonlinear. **Dispersion** arises from gravity, inertia, wavelength and finite water depth; **nonlinearity** arises because pressure and velocity conditions must be imposed on the unknown moving free surface and contain products or squares of the flow variables. For modest amplitudes, linear solutions can predict wavelength and subsurface kinematics satisfactorily. As wave amplitude rises, however, crest elevations are no longer the mirror image of trough elevations, harmonic energy is redistributed among spatial modes, the nonlinear relation between period and wavelength changes, and velocities beneath crests can depart significantly from linear predictions.

The classical progression of methods is instructive:

1. **Airy linear theory** linearizes the free-surface conditions about the mean level.
2. **Stokes theory** introduces a Fourier expansion with coefficients expressed as successive powers of a small wave-steepness parameter, allowing finite-amplitude corrections.
3. **Cnoidal theory** expands about the shallow-water regime, leading to elliptic-function waveforms with extended troughs and relatively narrow crests.
4. **Chappelear and Dean Fourier/stream-function methods** compute Fourier coefficients by numerical free-surface matching rather than by a perturbation series.
5. **Rienecker and Fenton (1981)** solve a closed collocation system by Newton's method, enforcing the nonlinear boundary conditions pointwise. Fenton (1988) reformulates the system in wavenumber-scaled variables and simplifies the historical FORTRAN implementation. Fenton's later reviews place this method among numerical treatments of more general wave problems [1–4].

The crucial distinction is that *Fourier approximation* and *Stokes perturbation theory* both employ trigonometric series but assign fundamentally different mathematical roles to their coefficients. In Stokes theory the coefficients are generated analytically order by order in amplitude. In the Fourier approximation method the coefficients are unknowns determined numerically from the fully nonlinear boundary conditions. For a fixed Fourier order, the nonlinear boundary conditions themselves are **not** Taylor-linearized.

### 1.2 Problem definition

Consider one train of identical waves moving at constant speed $c$ in the positive horizontal direction. Let

- $H>0$: crest-to-trough height of a **single regular wave**, not spectral significant wave height $H_{m0}$;
- $T>0$: observed period in a chosen fixed reference frame;
- $d>0$: mean depth from horizontal bed to still-water level;
- $L>0$: wavelength and $k=2\pi/L$ its wavenumber;
- $U_c$: prescribed collinear current, positive in the wave-propagation direction;
- $g$: gravitational acceleration;
- $\rho$: constant fluid density.

The nominal engineering input–output problem is

$$
\boxed{(H,T,d,U_c)\longmapsto\{L,c,\eta,\mathbf u,p,\text{integral properties}\}.}
$$

This representation is **incomplete unless the nature of $U_c$ is identified**. In particular, prescribing the Eulerian mean current and prescribing depth-integrated mass-transport velocity lead to different wave-frame constraints. The principal GUI adopts the Eulerian-current convention. The underlying mathematical system can also be closed using a Stokes/mass-transport current criterion [1–4, 11].

### 1.3 Model assumptions and governing domain

The assumptions are: (i) homogeneous density; (ii) incompressibility; (iii) negligible viscosity outside thin boundary layers; (iv) irrotational interior motion; (v) two-dimensionality; (vi) flat, rigid and impermeable bed; (vii) gravity as the restoring force; (viii) negligible surface tension; (ix) a smooth, non-overturning, periodic surface of permanent form. The fluid domain at time $t$ is

$$
\Omega(t)=\{(x,y):0<x<L,\;-d<y<\eta(x,t)\},
$$

where $y=0$ denotes the still-water level and $y=-d$ the bed. Periodicity is imposed in $x$; the free-surface profile is stationary after a Galilean transformation. These assumptions are materially different from an unsteady, three-dimensional, rotational, breaking-wave or wave–structure interaction model [3, 4, 8].

***

## 2. Reference frames, current conventions and kinematic identities

### 2.1 Laboratory and travelling coordinates

Let $(x,y)$ be a frame fixed to the bed, and introduce the dimensional travelling coordinate

$$
\xi=x-ct.
$$

A wave of permanent form is represented by

$$
\eta(x,t)=\eta(\xi),\qquad \eta(\xi+L)=\eta(\xi).
$$

Introduce the wave phase $\theta=k\xi$. Then

$$
\theta=k(x-ct)=kx-\omega t,\qquad \omega=kc=\frac{2\pi}{T},\qquad c=\frac{L}{T}.
$$

The latter identity relates **apparent phase speed** and period measured in the chosen fixed frame; it is not an independent theoretical prediction of wave speed on an unspecified current. In the wave frame the horizontal and vertical velocities are $(U,V)$, related to fixed-frame velocities $(u,v)$ by

$$
\boxed{U=u-c,\qquad V=v,\qquad u=U+c.}
$$

For any field $f(\xi,y)$,

$$
\left.\frac{\partial f}{\partial t}\right|_{x,y}=-c\frac{\partial f}{\partial x},
$$

which later gives exact relations between local time derivatives and spatial derivatives of the wave field.

### 2.2 The Eulerian current $\bar u_1$

For a steady periodic irrotational wave, the horizontal velocity averaged over one period at a fixed elevation below the trough is independent of the selected elevation, provided the averaging interval remains immersed. With the mean wave-frame horizontal velocity written $-\bar U$,

$$
\bar U=-\frac1L\int_0^L U(\xi,y)\,d\xi,
$$

and its corresponding fixed-frame Eulerian mean velocity is

$$
\boxed{\bar u_1=c-\bar U.}
$$

The sign convention is essential: for a positive-propagating wave the background flow in the wave frame is usually negative. Setting $\bar u_1=0$ selects the reference frame for which $c=\bar U$, traditionally associated with Stokes' first wave-speed convention. It is not an identity true for every laboratory frame [3, §2, pp. 2–3].

### 2.3 Stokes or mass-transport current $\bar u_2$

Define the **positive wave-frame discharge magnitude** $Q$ by

$$
\int_{-d}^{\eta(\xi)}U(\xi,y)\,dy=-Q.
$$

Incompressibility and impermeable boundaries make $Q$ independent of $\xi$. The fixed-frame depth-integrated discharge is

$$
\int_{-d}^{\eta}u\,dy=\int_{-d}^{\eta}(U+c)\,dy=-Q+c(d+\eta).
$$

Averaging over one wavelength with $\langle\eta\rangle=0$ yields $cd-Q$. The mean mass-transport velocity therefore is

$$
\boxed{\bar u_2=c-\frac Qd.}
$$

Consequently the two distinct conventions are related by

$$
\bar u_2-\bar u_1=\bar U-\frac Qd.
$$

The difference is a nonlinear mass-transport effect. A zero-Eulerian-current wave generally **does not** have zero depth-integrated mass transport. Conversely a wave in a closed flume may require $\bar u_2=0$ even though the Eulerian mean current is not zero [1, 3, 11].

### 2.4 Wave-volume-flux parameter $q$

An auxiliary variable reduces the algebraic complexity of the nonlinear constraints:

$$
\boxed{q=\bar U d-Q.}
$$

Equivalently,

$$
Q=\bar U d-q,
\qquad
\bar u_2=c-\bar U+\frac qd=\bar u_1+\frac qd.
$$

The distinction among $Q$, $q$, $\bar u_1$ and $\bar u_2$ prevents common dimensional and sign mistakes. $Q$ and $q$ have units $\mathrm{m^2/s}$ per unit width, while $\bar u_1$, $\bar u_2$, $c$ and $\bar U$ have units $\mathrm{m/s}$.

### 2.5 Why current is indispensable when period is prescribed

A Galilean transformation by an arbitrary uniform speed $V_0$ changes the apparent wave speed and fixed-frame horizontal velocity, but leaves the relative wave profile unchanged:

$$
\tilde c=c+V_0,\qquad \tilde u=u+V_0,
\qquad \tilde U=\tilde u-\tilde c=U.
$$

It follows that specifying $H$, $d$ and $L$ determines the **relative** steady wave solution, but does not determine its speed in a particular externally chosen frame. If instead an apparent period is prescribed, $c=L/T$ introduces an additional constraint that cannot be evaluated without a current or flux convention. In the linear approximation this distinction appears as Doppler shift:

$$
\boxed{\omega=\omega_\mathrm{int}+kU_c,\qquad
\omega_\mathrm{int}^2=gk\tanh(kd).}
$$

The positive intrinsic-frequency branch is $\omega-kU_c=+\sqrt{gk\tanh(kd)}$. An opposing current can admit blocking or multiple-root complications; squaring this equation loses the branch sign. In the finite-amplitude problem the current enters the *entire nonlinear collocation system*, not simply a post-processing correction to linear wavelength [1, 3, 5].

***

## 3. Field equations and free-boundary conditions

### 3.1 Conservation of mass and irrotationality

For constant-density incompressible two-dimensional flow,

$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0.
$$

For irrotational motion,

$$
\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}=0.
$$

A scalar velocity potential $\phi$ exists with $(u,v)=\nabla\phi$, and incompressibility implies

$$
\boxed{\nabla^2\phi=\phi_{xx}+\phi_{yy}=0.}
$$

A stream function $\psi$ exists for the *wave-frame* field:

$$
\boxed{U=\psi_y,\qquad V=-\psi_\xi.}
$$

Its defining derivatives automatically satisfy mass conservation. Irrotationality then requires

$$
\boxed{\psi_{\xi\xi}+\psi_{yy}=0.}
$$

Thus the governing differential equation in the fluid interior is **linear** despite the full physical problem being nonlinear. Nonlinearity is located principally in the boundary conditions imposed on the unknown free surface.

### 3.2 Bed condition

At the plane bed $y=-d$, impermeability gives

$$
V(\xi,-d)=0.
$$

Because the bed is a streamline, the arbitrary additive constant in $\psi$ is chosen so that

$$
\boxed{\psi(\xi,-d)=0.}
$$

For a sloping bed $y=b(x)$ the kinematic condition would instead be $v=u b_x$ in a stationary frame, so the simple harmonic basis developed below would no longer satisfy the bed condition identically. This is one reason a local steady Fourier calculation is not itself a shoaling or refraction model [4, 7, 8].

### 3.3 Free-surface kinematic condition

A material surface can be represented by $F(x,y,t)=y-\eta(x,t)=0$. The exact requirement that fluid particles do not cross it is

$$
\frac{DF}{Dt}=F_t+uF_x+vF_y=0.
$$

Therefore

$$
\boxed{\eta_t+u\eta_x=v\qquad\text{at }y=\eta.}
$$

For a permanent wave $\eta_t=-c\eta_x$,

$$
\boxed{(u-c)\eta_\xi=v\qquad\text{at }y=\eta(\xi).}
$$

Using the stream-function definition, the above states that the derivative of $\psi$ along the surface vanishes:

$$
\frac{d}{d\xi}\psi(\xi,\eta(\xi))=\psi_\xi+\psi_y\eta_\xi=-V+U\eta_\xi=0.
$$

The free surface is thus another streamline. With the discharge convention of §2.3,

$$
\boxed{\psi(\xi,\eta(\xi))=-Q.}
$$

No small-slope approximation has been introduced.

### 3.4 Dynamic free-surface condition and Bernoulli's law

For irrotational inviscid motion, unsteady Bernoulli in the fixed frame is

$$
\phi_t+\frac12(u^2+v^2)+gy+\frac p\rho=C(t).
$$

In the travelling frame it becomes a stationary Bernoulli relationship:

$$
\frac12(U^2+V^2)+gy+\frac p\rho=r,
$$

where $r$ is the Bernoulli constant with origin at still-water level $y=0$. At the free surface, gauge pressure is zero, hence

$$
\boxed{\frac12\left[U^2(\xi,\eta)+V^2(\xi,\eta)\right]+g\eta(\xi)=r.}
$$

Equivalently,

$$
\boxed{\frac12\left[\psi_y^2+\psi_\xi^2\right]_{y=\eta}+g\eta-r=0.}
$$

If vertical position is instead measured from the bed by $z=y+d$, the corresponding Bernoulli constant is $R=r+gd$ and

$$
\boxed{\frac12(U^2+V^2)+gz+\frac p\rho=R,
\qquad R=r+gd.}
$$

A difference of $gd$ between two presentations of the Bernoulli constant is therefore a **datum transformation**, not a physical difference in flow [3, §7; 11, Appendix B.6].

### 3.5 Mean elevation and crest-to-trough height

The mean-water-level datum requires

$$
\boxed{\langle\eta\rangle=\frac1L\int_0^L\eta(\xi)\,d\xi=0.}
$$

For a symmetric periodic wave with crest at $\xi=0$ and trough at $\xi=L/2$,

$$
\boxed{\eta(0)-\eta(L/2)=H.}
$$

The two equations remove the otherwise arbitrary vertical translation and scale the wave to the prescribed height.

***

## 4. Nondimensionalization and the Fourier stream function

### 4.1 Natural wavenumber–gravity scales

A convenient set of scales for the Fourier representation is

$$
\ell_k=\frac1k,\qquad V_k=\sqrt{\frac gk},\qquad
 t_k=\frac1{\sqrt{gk}},\qquad
 \psi_k=\sqrt{\frac g{k^3}}.
$$

Define

$$
X=k(x-ct),\qquad Y=ky,\qquad D=kd,\qquad \zeta(X)=k\eta(X),
$$

and nondimensional variables

$$
C=c\sqrt{\frac kg},\qquad A=\bar U\sqrt{\frac kg},\qquad
\mathcal Q=q\sqrt{\frac{k^3}g},\qquad
\mathcal R=\frac{rk}{g}.
$$

Let the dimensionless stream function be $\Psi=\psi\sqrt{k^3/g}$. Then

$$
\frac{\partial\Psi}{\partial Y}=U\sqrt{\frac kg},\qquad
-\frac{\partial\Psi}{\partial X}=V\sqrt{\frac kg}.
$$

The physical bed is $Y=-D$, the mean surface is $Y=0$, and the actual surface is $Y=\zeta(X)$.

### 4.2 Fourier ansatz and its exact field properties

For a symmetric steady wave, introduce a truncated cosine expansion of order $N$:

$$
\boxed{\Psi_N(X,Y)=-A(Y+D)+\sum_{j=1}^N B_j\,S_j(Y;D)\cos(jX),}
$$

where

$$
\boxed{S_j(Y;D)=\frac{\sinh[j(Y+D)]}{\cosh(jD)}.}
$$

Each term $S_j(Y;D)\cos(jX)$ is harmonic because

$$
\frac{\partial^2}{\partial X^2}\left(S_j\cos jX\right)=-j^2S_j\cos jX,
$$

$$
\frac{\partial^2}{\partial Y^2}\left(S_j\cos jX\right)=+j^2S_j\cos jX.
$$

Their sum is zero. At the bed $Y=-D$, $S_j(-D;D)=0$, and the linear mean-flow term also vanishes. Consequently

$$
\nabla^2\Psi_N=0\quad\text{throughout the fluid},\qquad
\Psi_N(X,-D)=0
$$

**for any finite coefficients $B_j$**. This property is central: neither Laplace's equation nor the impermeable-bed boundary condition requires a numerical approximation or collocation residual. Only the surface and global conditions remain to be solved [1, §2; 2, pp. 358–360].

### 4.3 The companion vertical velocity basis

Define

$$
\boxed{C_j(Y;D)=\frac{\cosh[j(Y+D)]}{\cosh(jD)}.}
$$

The derivative identities are

$$
\frac{\partial S_j}{\partial Y}=jC_j,\qquad
\frac{\partial C_j}{\partial Y}=jS_j.
$$

The wave-frame nondimensional velocities become

$$
\boxed{\widehat U=-A+\sum_{j=1}^N jB_j C_j(Y;D)\cos(jX),}
$$

$$
\boxed{\widehat V=\sum_{j=1}^N jB_j S_j(Y;D)\sin(jX).}
$$

Returning to the physical laboratory frame,

$$
\boxed{u=c-\bar U+\sqrt{\frac gk}\sum_{j=1}^N jB_jC_j(ky;kd)\cos[jk(x-ct)],}
$$

$$
\boxed{v=\sqrt{\frac gk}\sum_{j=1}^N jB_jS_j(ky;kd)\sin[jk(x-ct)].}
$$

The mean of each harmonic over $X\in[0,2\pi]$ is zero, directly recovering the Eulerian-current identity $\bar u_1=c-\bar U$.

### 4.4 Stable scaled-exponential evaluation

Direct computation of $\sinh(jD)$ and $\cosh(jD)$ can overflow for high $jD$, even if their **ratios** remain finite and well behaved. Set

$$
E_j=e^{-2jD},\qquad P_j(Y)=e^{jY},\qquad Q_j(Y)=e^{-jY-2jD},\qquad D_j=1+E_j.
$$

Then the equivalent formulas are

$$
\boxed{S_j=\frac{P_j-Q_j}{D_j},\qquad C_j=\frac{P_j+Q_j}{D_j}.}
$$

The transformations follow by multiplying numerator and denominator of the hyperbolic ratios by $2e^{-jD}$. The large exponentials $e^{jD}$ never need to be constructed. At the seabed,

$$
S_j(-D;D)=0,\qquad
C_j(-D;D)=\frac{2e^{-jD}}{1+e^{-2jD}}=\operatorname{sech}(jD).
$$

For fixed $Y\le0$ and $D\to\infty$,

$$
S_j(Y;D)\to e^{jY},\qquad C_j(Y;D)\to e^{jY}.
$$

The exponential limit is appropriate for deep water, but **the finite-depth code must retain the reflected exponential** if exact impermeability at the finite seabed is required. The normalization by $\cosh(jD)$ and the choice of a reference depth were emphasized by Rienecker and Fenton to improve deep-water numerical behavior [1, pp. 122–123; 2, p. 360; 11, Appendix A].

***

## 5. The complete finite nonlinear collocation system

### 5.1 Why one solves for nodal surface elevations

Introduce $N+1$ equally spaced phase nodes over a symmetric half-wave:

$$
\boxed{X_m=\frac{m\pi}{N},\qquad m=0,1,\ldots,N.}
$$

The wave crest is at $m=0$ and the trough at $m=N$. Treat the surface values $\zeta_m=\zeta(X_m)$ as unknowns, rather than prescribing a shape. This avoids introducing a separate Fourier series of unknown surface coefficients into the nonlinear equations. The numerical surface Fourier coefficients can be reconstructed after solution.

The model enforces the kinematic and dynamic boundary conditions at each of these $N+1$ nodes. Together with eight global conditions, this produces $2N+10$ equations in $2N+10$ unknowns, exactly the dimension of Fenton's 1988 wavenumber-based formulation [2, pp. 358–360].

### 5.2 Full state-vector definition

A 1-indexed vector $z$ is used:

| Index | Nondimensional unknown | Physical interpretation |
|:---|:---|:---|
| $z_1$ | $kd$ | Depth relative to wavelength |
| $z_2$ | $kH$ | Total wave height |
| $z_3$ | $T\sqrt{gk}$ | Period |
| $z_4$ | $c\sqrt{k/g}$ | Apparent celerity |
| $z_5$ | $\bar u_1\sqrt{k/g}$ | Eulerian current |
| $z_6$ | $\bar u_2\sqrt{k/g}$ | Mass-transport current |
| $z_7$ | $\bar U\sqrt{k/g}$ | Positive magnitude of mean opposing wave-frame flow |
| $z_8$ | $q\sqrt{k^3/g}$ | Wave-volume-flux variable |
| $z_9$ | $rk/g$ | Mean-level Bernoulli constant |
| $z_{10+m}$, $0\le m\le N$ | $\zeta_m=k\eta(X_m)$ | Surface collocation ordinates |
| $z_{N+10+j}$, $1\le j\le N$ | $B_j$ | Fourier stream-function coefficients |

The active vector contains

$$
9+(N+1)+N=\boxed{2N+10}
$$

unknowns. The Rust implementation retains an unused element $z[0]$ so the numerical indices agree with the analytical notation. Thus **$N=50$ produces a 110-dimensional nonlinear system**, not a 50-equation system.

### 5.3 Global prescribed-input parameters

For a period-specified wave define

$$
h_r=\frac Hd,\qquad h_s=\frac{H}{gT^2},\qquad
\chi=\frac{U_c}{\sqrt{gd}}.
$$

These symbols are small only when the physical conditions make them small; their definitions do not assume linear waves. During height continuation replace $H$ in $h_r$ and $h_s$ by the current trial height $H_\lambda=\lambda H$, where $0<\lambda\le1$.

### 5.4 The eight global residual equations

The first two residuals enforce physical height and period consistency:

$$
\boxed{F_1=z_2-h_rz_1=0,\qquad F_2=z_2-h_s z_3^2=0.}
$$

Indeed, $h_r z_1=(H/d)(kd)=kH=z_2$, and

$$
h_sz_3^2=\frac{H}{gT^2}(T^2gk)=kH.
$$

The period–celerity identity is

$$
\boxed{F_3=z_4z_3-2\pi=0,}
$$

since $(c\sqrt{k/g})(T\sqrt{gk})=ckT=2\pi$. The current and flux definitions become

$$
\boxed{F_4=z_5+z_7-z_4=0,}
$$

$$
\boxed{F_5=z_1(z_6+z_7-z_4)-z_8=0.}
$$

A current selector $s$ determines which mean current is prescribed:

$$
\boxed{F_6=z_{s+4}-\chi\sqrt{z_1}=0,\qquad
s=\begin{cases}1&\text{Eulerian current},\\2&\text{mass-transport current}.\end{cases}}
$$

To check the nondimensional scaling note that

$$
\chi\sqrt{z_1}=\frac{U_c}{\sqrt{gd}}\sqrt{kd}=U_c\sqrt{\frac kg}.
$$

A discrete mean-level condition, equivalent to trapezoidal integration of the half-wave cosine interpolant, is

$$
\boxed{F_7=\zeta_0+\zeta_N+2\sum_{m=1}^{N-1}\zeta_m=0.}
$$

The crest-to-trough condition is

$$
\boxed{F_8=\zeta_0-\zeta_N-z_2=0.}
$$

The factor of two in the interior-node terms of $F_7$ arises from symmetry and endpoint half weights. It is **not** an arbitrary discrete average; the constant Fourier mode of the interpolating surface is being set to zero.

### 5.5 Surface sums evaluated at a collocation point

For each $m=0,\ldots,N$, set $D=z_1$, $\zeta_m=z_{10+m}$ and $b_j=z_{N+10+j}$. Evaluate

$$
\mathcal P_m=\sum_{j=1}^N b_jS_j(\zeta_m;D)\cos(jX_m),
$$

$$
\mathcal U_m=\sum_{j=1}^N jb_jC_j(\zeta_m;D)\cos(jX_m),
$$

$$
\mathcal V_m=\sum_{j=1}^N jb_jS_j(\zeta_m;D)\sin(jX_m).
$$

Here $\mathcal U_m$ is the **harmonic contribution**, not the total horizontal wave-frame velocity. The total is $\mathcal U_m-z_7$. Keeping this distinction is essential in differentiating Bernoulli's residual correctly.

### 5.6 Surface streamline equations

At the free surface, $\Psi=-Q\sqrt{k^3/g}$. From the Fourier representation,

$$
\Psi(X_m,\zeta_m)=-z_7(\zeta_m+D)+\mathcal P_m.
$$

Using $z_8=z_7D-Q\sqrt{k^3/g}$, the residual simplifies to

$$
\boxed{F_{9+m}=\mathcal P_m-z_8-z_7\zeta_m=0,
\qquad m=0,\ldots,N.}
$$

These are **$N+1$ nonlinear kinematic equations**. Nonlinearity enters because the normalized hyperbolic functions $S_j$ are evaluated at the unknown surface ordinate $\zeta_m$.

### 5.7 Surface Bernoulli equations

The nondimensional dynamic condition is

$$
\frac12\left[(\mathcal U_m-z_7)^2+\mathcal V_m^2\right]+\zeta_m-z_9=0.
$$

Thus

$$
\boxed{F_{N+10+m}=\frac12\left[(\mathcal U_m-z_7)^2+\mathcal V_m^2\right]+\zeta_m-z_9=0,
\quad m=0,\ldots,N.}
$$

The squared velocities provide essential finite-amplitude nonlinear coupling. Combined with the eight global equations and the $N+1$ streamline equations, the entire system is

$$
\boxed{\mathbf F(\mathbf z;H,T,d,U_c,s)=\mathbf0,
\qquad \mathbf F\in\mathbb R^{2N+10}.}
$$

For $N=50$, $8+51+51=110$ equations. The original steady-wave papers present equivalent closed systems in slightly different coordinate conventions [1, §2.1; 2, pp. 358–360; 11, Appendix B.3].

### 5.8 Direct wavelength input is a different closure

A problem specified by $(H,L,d)$ does not need the period relation used above to determine the *relative* wave field. One may fix $kd=2\pi d/L$ and replace the appropriate global equation by

$$
F_{L}=z_1-\frac{2\pi d}{L}=0
$$

or use the equivalent height–wavelength relation $z_2=2\pi H/L$, according to the chosen unknown set. A separately specified current or frame then determines apparent $c$ and $T$. Fenton's historical **FOURIER** program accommodates both the period and wavelength modes; the particular `fenton_gui.rs` entry point studied here exposes the **period-specified** input mode. These interfaces should not be described as identical [2, 11].

***
## 6. Exact differentiation and the analytical Jacobian

### 6.1 Newton's method requires a well-defined Jacobian

Let $M=2N+10$. Newton iteration linearizes $\mathbf F(\mathbf z)=0$ about the current state:

$$
\mathbf F(\mathbf z+\Delta\mathbf z)
=\mathbf F(\mathbf z)+\mathbf J(\mathbf z)\Delta\mathbf z
+O(\|\Delta\mathbf z\|^2),
$$

where

$$
\boxed{J_{ij}=\frac{\partial F_i}{\partial z_j},\qquad
\mathbf J\in\mathbb R^{M\times M}.}
$$

The Newton correction satisfies

$$
\boxed{\mathbf J\Delta\mathbf z=-\mathbf F.}
$$

An analytical Jacobian avoids reevaluating the full nonlinear residual for every finite-difference perturbation and avoids finite-difference truncation and cancellation errors. The original Rienecker–Fenton (1981) paper derived analytical derivatives; Fenton's simplified 1988 historical FORTRAN program used numerical differentiation. The Rust implementation returns to explicit analytical derivatives but combines them with a stabilized singular-value linear solve [1, §2.2; 2, p. 360].

### 6.2 Differentiation of the normalized vertical basis

Write again

$$
P_j=e^{j\zeta},\qquad Q_j=e^{-j\zeta-2jD},\qquad E_j=e^{-2jD},\qquad W_j=1+E_j,
$$

and

$$
S_j=\frac{P_j-Q_j}{W_j},\qquad C_j=\frac{P_j+Q_j}{W_j}.
$$

At fixed $D$,

$$
\boxed{\frac{\partial S_j}{\partial\zeta}=jC_j,
\qquad \frac{\partial C_j}{\partial\zeta}=jS_j.}
$$

At fixed $\zeta$, use

$$
\frac{\partial Q_j}{\partial D}=-2jQ_j,\qquad
\frac{\partial W_j}{\partial D}=-2jE_j.
$$

The quotient rule gives

$$
\frac{\partial S_j}{\partial D}
=\frac{2jQ_jW_j+2jE_j(P_j-Q_j)}{W_j^2}
=\boxed{\frac{2j(Q_j+E_jP_j)}{W_j^2}},
$$

$$
\frac{\partial C_j}{\partial D}
=\frac{-2jQ_jW_j+2jE_j(P_j+Q_j)}{W_j^2}
=\boxed{\frac{2j(E_jP_j-Q_j)}{W_j^2}}.
$$

These formulas are useful because they remain within the same scaled-exponential representation as the residual function; one need not differentiate numerically unstable large hyperbolic functions.

### 6.3 Derivatives of collocation sums

At node $m$, define $a_m=\mathcal U_m-z_7$ and $b_m=\mathcal V_m$. For brevity denote a partial derivative with respect to $D=z_1$ by a subscript $D$, and with respect to local $\zeta_m$ by a subscript $\zeta$. Then

$$
\mathcal P_{m,D}
=\sum_{j=1}^N B_j\frac{\partial S_j}{\partial D}\cos(jX_m),
$$

$$
\mathcal U_{m,D}
=\sum_{j=1}^N jB_j\frac{\partial C_j}{\partial D}\cos(jX_m),
$$

$$
\mathcal V_{m,D}
=\sum_{j=1}^N jB_j\frac{\partial S_j}{\partial D}\sin(jX_m).
$$

The derivatives with respect to the local surface coordinate are

$$
\boxed{\mathcal P_{m,\zeta}=\mathcal U_m,}
$$

$$
\boxed{\mathcal U_{m,\zeta}
=\sum_{j=1}^N j^2 B_j S_j(\zeta_m;D)\cos(jX_m),}
$$

$$
\boxed{\mathcal V_{m,\zeta}
=\sum_{j=1}^N j^2 B_j C_j(\zeta_m;D)\sin(jX_m).}
$$

The identity $\mathcal P_{m,\zeta}=\mathcal U_m$ is particularly valuable: it means that the surface-ordinate derivative of the streamline residual is the **total horizontal wave-frame velocity**, $a_m$.

### 6.4 Jacobian of the eight global constraints

All unlisted derivatives in the following expressions vanish. For $F_1$ and $F_2$,

$$
\frac{\partial F_1}{\partial z_1}=-h_r,\quad
\frac{\partial F_1}{\partial z_2}=1;
\qquad
\frac{\partial F_2}{\partial z_2}=1,\quad
\frac{\partial F_2}{\partial z_3}=-2h_sz_3.
$$

For the phase/period condition,

$$
\frac{\partial F_3}{\partial z_3}=z_4,\qquad
\frac{\partial F_3}{\partial z_4}=z_3.
$$

For the Eulerian identity,

$$
\frac{\partial F_4}{\partial z_4}=-1,\qquad
\frac{\partial F_4}{\partial z_5}=1,\qquad
\frac{\partial F_4}{\partial z_7}=1.
$$

For the mass-flux identity,

$$
\begin{aligned}
F_{5,z_1}&=z_6+z_7-z_4,\\
F_{5,z_4}&=-z_1,\\
F_{5,z_6}&=z_1,\\
F_{5,z_7}&=z_1,\\
F_{5,z_8}&=-1.
\end{aligned}
$$

For the selected mean-current constraint,

$$
\boxed{F_{6,z_1}=-\frac\chi{2\sqrt{z_1}},\qquad
F_{6,z_{s+4}}=1.}
$$

For the mean-level condition,

$$
F_{7,z_{10+m}}=\begin{cases}1,&m=0\text{ or }N,\\2,&1\le m\le N-1,\end{cases}
$$

while

$$
F_{8,z_{10}}=1,\qquad F_{8,z_{N+10}}=-1,\qquad F_{8,z_2}=-1.
$$

These derivatives impose a useful implementation invariant: the global residuals and their derivatives must use the **same** current selector and mean-level weights. A mismatch can create a Newton correction that appears numerically reasonable but does not linearize the actual equations being solved.

### 6.5 Jacobian of the streamline residual

For

$$
K_m=\mathcal P_m-z_8-z_7\zeta_m,
$$

we have

$$
\boxed{\frac{\partial K_m}{\partial z_1}=\mathcal P_{m,D},\quad
\frac{\partial K_m}{\partial z_7}=-\zeta_m,\quad
\frac{\partial K_m}{\partial z_8}=-1.}
$$

The derivative with respect to its own surface ordinate is

$$
\boxed{\frac{\partial K_m}{\partial\zeta_m}=\mathcal U_m-z_7=a_m.}
$$

For Fourier coefficient $B_j$,

$$
\boxed{\frac{\partial K_m}{\partial B_j}
=S_j(\zeta_m;D)\cos(jX_m).}
$$

For $\ell\ne m$, $\partial K_m/\partial\zeta_\ell=0$ because this collocation equation is evaluated only at its own surface node.

### 6.6 Jacobian of Bernoulli's residual

For

$$
B_m=\frac12\bigl(a_m^2+b_m^2\bigr)+\zeta_m-z_9,
$$

application of the chain rule yields

$$
\boxed{\frac{\partial B_m}{\partial z_1}
=a_m\mathcal U_{m,D}+b_m\mathcal V_{m,D},}
$$

$$
\boxed{\frac{\partial B_m}{\partial z_7}=-a_m,
\qquad \frac{\partial B_m}{\partial z_9}=-1,}
$$

and

$$
\boxed{\frac{\partial B_m}{\partial\zeta_m}
=1+a_m\mathcal U_{m,\zeta}+b_m\mathcal V_{m,\zeta}.}
$$

For $B_j$,

$$
\boxed{\frac{\partial B_m}{\partial B_j}
=j\left[a_mC_j(\zeta_m;D)\cos(jX_m)
+b_mS_j(\zeta_m;D)\sin(jX_m)\right].}
$$

The formulas reproduce the structure assembled by the Rust `jacobian()` method. Their consistency with the nonlinear residual is critical. In particular, differentiating $\frac12 a_m^2$ gives $a_m$ times the derivative of $a_m$, not the derivative of a linear-wave approximation.

### 6.7 Conditioning and column equilibration

The Jacobian contains physical variables of different scales, including mean-current, period, Bernoulli and Fourier coefficients. A raw condition number can be distorted by the arbitrary units and relative scales of its columns. The numerical linear solve first computes

$$
s_j=\max_i|J_{ij}|,\qquad
\widetilde J_{ij}=J_{ij}/s_j
$$

for nonzero columns, replacing a zero scale by unity. Define the diagonal matrix $D_s=\operatorname{diag}(s_j)$. Then

$$
\widetilde{\mathbf J}=\mathbf J D_s^{-1},\qquad
\widetilde{\mathbf J}\,\Delta\widetilde{\mathbf z}=-\mathbf F,
\qquad \Delta\mathbf z=D_s^{-1}\Delta\widetilde{\mathbf z}.
$$

This transformation does not change the exact solution of a nonsingular linear system. It changes the numerical geometry of the linear problem to mitigate scaling disparities; it is *not* a guarantee that a nearly singular nonlinear wave has become well posed.

***

## 7. Small-amplitude theory and initialization

### 7.1 Linear Airy wave in finite depth

For infinitesimal wave amplitude, let

$$
\eta(x,t)=a\cos(kx-\omega t),\qquad a=H/2.
$$

A suitable linear velocity potential, for zero uniform current and bed at $y=-d$, is

$$
\phi_1(x,y,t)=\frac{ag}{\omega}
\frac{\cosh[k(y+d)]}{\cosh(kd)}\sin(kx-\omega t).
$$

Differentiation gives

$$
u_1=\frac{agk}{\omega}\frac{\cosh[k(y+d)]}{\cosh(kd)}\cos\theta,
$$

$$
v_1=\frac{agk}{\omega}\frac{\sinh[k(y+d)]}{\cosh(kd)}\sin\theta.
$$

Linearizing the kinematic condition to $\eta_t=v_1(0)$ and the dynamic condition to $\phi_{1,t}+g\eta=0$ yields

$$
\boxed{\omega^2=gk\tanh(kd).}
$$

The phase velocity relative to the zero-current frame is

$$
\boxed{c_{\mathrm{lin}}=\frac\omega k
=\sqrt{\frac gk\tanh(kd)}.}
$$

The limiting regimes are

$$
kd\gg1:\qquad c_{\mathrm{lin}}\sim\sqrt{\frac gk},\qquad
L\sim\frac{gT^2}{2\pi},
$$

$$
kd\ll1:\qquad c_{\mathrm{lin}}\sim\sqrt{gd},\qquad L\sim T\sqrt{gd}.
$$

Linear theory supplies an initial approximation; it does not determine the final solution at finite $H$.

### 7.2 Solving for the initial $kd$ from a known period

Introduce

$$
\sigma=\omega\sqrt{\frac dg}
=\frac{2\pi}{T\sqrt{g/d}}.
$$

Writing $D=kd$, linear dispersion becomes

$$
\boxed{\sigma^2=D\tanh D.}
$$

A convenient Fenton–McKee-type explicit approximation is

$$
\boxed{D_0=\frac{\sigma^2}{\left[\tanh(\sigma^{3/2})\right]^{2/3}}.}
$$

Since $\tanh s=1/\coth s$, this is equivalently

$$
D_0=\sigma^2\left[\coth(\sigma^{3/2})\right]^{2/3}.
$$

For small $\sigma$, the denominator behaves like $\sigma$, giving $D_0\sim\sigma$. For large $\sigma$ the hyperbolic tangent tends to unity and $D_0\sim\sigma^2$. Thus the approximation has the correct shallow- and deep-water asymptotes. Fenton and McKee's approximation is a starting estimate—not a substitute for solving the final nonlinear wave problem [3, §3; 11, Appendix B.1].

### 7.3 First-order initial state

With an initial $D=kd$, define

$$
C_0=\sqrt{\tanh D},\qquad kH=D(H/d).
$$

The zero-current first-order values are

$$
z_1=D,\quad z_2=kH,\quad z_4=C_0,\quad
z_3=\frac{2\pi}{C_0},
$$

$$
z_7=C_0,\quad z_8=0,\quad z_9=\frac12C_0^2,
$$

$$
\zeta_m^{(0)}=\frac{kH}{2}\cos\left(\frac{m\pi}{N}\right),
$$

$$
\boxed{B_1^{(0)}=\frac{kH}{2C_0},\qquad B_j^{(0)}=0\ (j\ge2).}
$$

The current state variable $z_5$ or $z_6$ can be initialized from the selected current criterion. A full nonlinear solution then adjusts the entire state, including $D$, the free surface and every Fourier coefficient. For significant current or steep wave height the initial guess need not satisfy the target equations closely; this motivates height continuation.

### 7.4 Why linearization is a physically useful but mathematically limited seed

At first order the surface is symmetric about the mean level and has a single harmonic. At finite amplitude, higher harmonics of the surface and of velocity are generated by the quadratic terms in Bernoulli's condition and by evaluating the vertical basis at the displaced free surface. Consequently the nonlinear wavelength need not equal the Airy prediction, even if its period and depth are unchanged. Newton iteration would fail more often if asked to jump directly from the linear state to a steep target wave, because it may cross a basin boundary or encounter ill-conditioning before reaching the correct solution branch.

***

## 8. Nonlinear solution, singular values and height continuation

### 8.1 Local Newton correction and globalization

If $\mathbf F$ is differentiable and its Jacobian is nonsingular at a root, undamped Newton iteration can converge quadratically sufficiently near that root. Far from the root it can diverge, and the Jacobian can be poorly conditioned. The practical update therefore uses a step fraction $\alpha$:

$$
\boxed{\mathbf z_{n+1}=\mathbf z_n+\alpha_n\Delta\mathbf z_n,
\qquad 0<\alpha_n\le1.}
$$

A residual merit function is

$$
\Phi(\mathbf z)=\frac12\|\mathbf F(\mathbf z)\|_2^2
=\frac12\sum_{i=1}^{M}F_i(\mathbf z)^2.
$$

The numerical algorithm first tests $\alpha=1$, then $1/2$, $1/4$, and so forth until the candidate is finite, has $kd>0$ and does not increase the merit function. Thus it performs a *residual-nonincreasing backtracking search*. It is similar in purpose to an Armijo safeguard, although the implementation does not require the usual Armijo slope parameter. If no acceptable step is found, the current continuation stage is treated as failed rather than accepting an unstable Newton correction.

### 8.2 Column-equilibrated singular-value decomposition

Let $\widetilde{\mathbf J}$ be the column-scaled Jacobian. Its singular-value decomposition is

$$
\widetilde{\mathbf J}=\mathbf U\Sigma\mathbf V^{\mathsf T},
$$

where $\mathbf U$ and $\mathbf V$ are orthogonal and $\Sigma$ contains nonnegative singular values $\sigma_1,\ldots,\sigma_M$. The least-squares correction with a truncated pseudoinverse is

$$
\boxed{\Delta\widetilde{\mathbf z}
=-\mathbf V\Sigma_{\mathrm{cut}}^+
\mathbf U^{\mathsf T}\mathbf F.}
$$

The application uses a relative singular-value cutoff:

$$
\sigma_j\le10^{-15}\sigma_{\max}\quad\Rightarrow\quad
\frac1{\sigma_j}\ \text{is replaced by zero in the pseudoinverse}.
$$

For retained modes $\sigma_j>10^{-15}\sigma_{\max}$, the reciprocal $1/\sigma_j$ is used. Truncation suppresses the amplification of residual and rounding errors in directions where the Jacobian has insufficient numerical rank. The price is that the correction may solve a regularized linearized problem rather than exactly inverting $\mathbf J$.

### 8.3 One-sided Jacobi SVD

A one-sided Jacobi SVD orthogonalizes matrix columns through plane rotations. For columns $\mathbf a_p,\mathbf a_q$ define

$$
p=\mathbf a_p^{\mathsf T}\mathbf a_p,\qquad
q=\mathbf a_q^{\mathsf T}\mathbf a_q,\qquad
r=\mathbf a_p^{\mathsf T}\mathbf a_q.
$$

For non-negligible $r$, a stable plane rotation can be constructed from

$$
\tau=\frac{q-p}{2r},\qquad
t=\frac{\operatorname{sign}(\tau)}{|\tau|+\sqrt{1+\tau^2}},
$$

$$
\gamma=\frac1{\sqrt{1+t^2}},\qquad
s=t\gamma.
$$

Both working columns and the accumulated right singular vectors are rotated using the orthogonal matrix

$$
\mathbf G=\begin{pmatrix}\gamma&s\\-s&\gamma\end{pmatrix}
$$

with the sign convention chosen consistently with the implementation's column updates. Once orthogonalized, each transformed column norm estimates a singular value. The algorithm employs repeated sweeps, checks relative column orthogonality against a binary64-epsilon-scaled threshold, and fails explicitly if its sweep limit is exceeded. This gives a robust approach to dense matrices of moderate order without forming normal equations $\mathbf J^{\mathsf T}\mathbf J$, which would square the condition number.

### 8.4 Condition-number screening

The spectral condition estimate is

$$
\boxed{\kappa_2(\widetilde{\mathbf J})
=\frac{\sigma_{\max}}{\sigma_{\min}}.}
$$

If the smallest singular value is effectively zero, the condition is treated as infinite. The Rust solver tests a numerical cutoff near

$$
\boxed{\kappa_2\le10^{11}.}
$$

If the condition is too large, it may restart the entire computation at a smaller Fourier order, with an active-order sequence such as

$$
50\rightarrow40\rightarrow32\rightarrow25\rightarrow20\rightarrow16\rightarrow12.
$$

A reduced order decreases matrix size and may improve conditioning, but it also reduces spatial resolution. It is therefore **not** an unconditional accuracy improvement. A reduced-order solution still has to pass the off-grid free-surface and integral checks. Condition number is a sensitivity indicator for the *linearized discrete equations*, not a proof about the physical existence or uniqueness of a steady wave.

### 8.5 Continuation in wave height

Let $H_\mathrm{target}$ be the requested wave height. Introduce a continuation parameter $\lambda\in[0,1]$ and solve the family

$$
\boxed{\mathbf F(\mathbf z;\lambda H_\mathrm{target},T,d,U_c)=\mathbf0.}
$$

The initial small-height state is close to linear theory. Having solved successive stages at $(\lambda_{n-1},\mathbf z_{n-1})$ and $(\lambda_n,\mathbf z_n)$, predict the next state by secant extrapolation:

$$
\boxed{\mathbf z^{\mathrm{pred}}_{n+1}
=\mathbf z_n+
\frac{\lambda_{n+1}-\lambda_n}{\lambda_n-\lambda_{n-1}}
(\mathbf z_n-\mathbf z_{n-1}).}
$$

The extrapolation becomes useful once two accepted states are available. The next Newton solve starts from this predictor. If a stage fails, the preceding accepted state is restored and the continuation increment is reduced; after a successful stage the increment may increase toward a chosen cap. The final accepted state must correspond to $\lambda=1$. A solution at $\lambda=0.8$, for example, is not a solution of the user-requested wave height and must not be reported as if it were [2, p. 361; 4, §3.1; 11, Appendix B.5].

### 8.6 Residual tolerances and convergence are separate criteria

Useful norms are

$$
\|\mathbf F\|_{\infty}=\max_i|F_i|,
\qquad
\|\mathbf F\|_{\mathrm{RMS}}=
\sqrt{\frac1M\sum_{i=1}^{M}F_i^2}.
$$

They are not interchangeable: for $M$ equations,

$$
\|\mathbf F\|_{\mathrm{RMS}}\le\|\mathbf F\|_{\infty}
\le\sqrt M\,\|\mathbf F\|_{\mathrm{RMS}}.
$$

The implementation uses nominal intermediate and final infinity-norm limits of

$$
\boxed{\|\mathbf F\|_\infty\le10^{-8}
\quad\text{(intermediate height stages)},}
$$

$$
\boxed{\|\mathbf F\|_\infty\le10^{-10}
\quad\text{(final requested wave height)}.}
$$

It additionally checks that the Newton surface correction has become small. A small correction without a small residual is not sufficient: the solver may stagnate near a singular or inconsistent linearization. Conversely a small collocation residual does not alone establish accurate behavior **between** the collocation nodes. The independent verification procedures are therefore essential.

### 8.7 Numerical complexity and practical performance

With $M=2N+10$, assembling a dense Jacobian requires $O(N^2)$ basic harmonic operations, while a dense SVD generally requires work scaling approximately like $O(M^3)$. Higher $N$ also increases the cost of dense phase/depth integrations and interval checks. The optimal $N$ is not necessarily the largest available: it is the smallest order giving adequate resolution while remaining comfortably conditioned. The code uses fixed numerical choices to balance these effects; empirical runtime should be measured on the intended CPU rather than inferred solely from an optimization flag.

***

## 9. Reconstruction of the continuous free surface

### 9.1 Cosine transform from the collocation ordinates

The solved surface is known at $X_m=m\pi/N$. The physical elevation for arbitrary phase can be reconstructed by a cosine interpolant. Define the discrete cosine coefficients

$$
\boxed{E_j=\frac2N\sum_{m=0}^N{}''\zeta_m
\cos\left(\frac{jm\pi}{N}\right),\qquad 0\le j\le N,}
$$

where $\sum{}''$ means that the endpoint terms $m=0,N$ each carry weight $1/2$. The full interpolant is

$$
\boxed{\zeta_N(X)=\frac{E_0}{2}
+\sum_{j=1}^{N-1}E_j\cos(jX)
+\frac{E_N}{2}\cos(NX).}
$$

For the still-water-level datum, $E_0=0$. The **last Fourier harmonic has the half weight $E_N/2$**; omitting that factor would distort the energy and high-order spectral diagnostics [11, Appendix B.6]. The dimensional elevation is $\eta_N=\zeta_N/k$ and the absolute depth above the bed is $d+\eta_N$.

### 9.2 Derivatives and curvature of the surface

Because $X=k(x-ct)$,

$$
\eta_x=\frac{d\zeta_N}{dX},
\qquad
\eta_{xx}=k\frac{d^2\zeta_N}{dX^2}.
$$

Thus

$$
\boxed{\eta_x=-\sum_{j=1}^{N-1}jE_j\sin(jX)
-\frac{NE_N}{2}\sin(NX),}
$$

$$
\boxed{\eta_{xx}=-k\left[
\sum_{j=1}^{N-1}j^2E_j\cos(jX)
+\frac{N^2E_N}{2}\cos(NX)\right].}
$$

For a graph-type free surface, curvature is

$$
\mathcal K(X)=\frac{\eta_{xx}}{\left(1+\eta_x^2\right)^{3/2}}.
$$

Higher derivatives amplify higher Fourier modes by increasing powers of $j$. Consequently coefficients that appear numerically small in elevation can materially affect curvature or derivatives near a sharp crest. A small $|E_N|$ is necessary for good resolution but is not sufficient by itself.

### 9.3 Why Fourier convergence may become difficult near the highest wave

For a sufficiently smooth periodic function, Fourier coefficients often decay rapidly. A nearly limiting steep wave develops a progressively narrow, sharp crest and extended trough region. Its Fourier spectrum broadens, so many harmonics may be required to resolve the crest. The same high $N$ can make the Jacobian increasingly ill-conditioned, leading to a practical competition between spectral resolution and floating-point stability. At the exact highest limiting wave the idealized crest may cease to be smooth, so exponential coefficient decay cannot be assumed uniformly as the limit is approached. This observation is central to the 1981 and 1999 treatments [1, 4, 11].

### 9.4 Collocation accuracy versus off-grid accuracy

The solver imposes the nonlinear equations at $X_m$ exactly up to its numerical tolerance. At an arbitrary phase $X$, one must re-evaluate

$$
\mathcal R_K(X)=\mathcal P(X,\zeta_N(X))-z_8-z_7\zeta_N(X),
$$

$$
\mathcal R_B(X)=\frac12\left[(\mathcal U(X,\zeta_N)-z_7)^2
+\mathcal V(X,\zeta_N)^2\right]+\zeta_N(X)-z_9.
$$

The maximum of $|\mathcal R_K|$ and $|\mathcal R_B|$ over an off-grid phase scan detects unresolved oscillations that collocation residuals cannot reveal. A finite scan is not a rigorous supremum over the continuous phase interval; it is a strong numerical diagnostic whose effectiveness depends on scan density and the smoothness of the solution.

***

## 10. Complete velocity field and differential kinematics

### 10.1 Fixed-frame and moving-frame velocities

Once the nonlinear state has been determined, the velocity field is available at any physically admissible position beneath the free surface. Introduce the harmonic sums

$$
\mathcal U_h(X,Y)=\sum_{j=1}^{N}jB_jC_j(Y;D)\cos(jX),
\qquad
\mathcal V_h(X,Y)=\sum_{j=1}^{N}jB_jS_j(Y;D)\sin(jX).
$$

The nondimensional wave-frame velocities are

$$
\widehat U=-z_7+\mathcal U_h,\qquad \widehat V=\mathcal V_h.
$$

Using $\sqrt{g/k}$ as the velocity scale, the dimensional laboratory-frame field is

$$
\boxed{u(X,y)=\bar u_1+\sqrt{\frac gk}\,\mathcal U_h(X,ky),}
$$

$$
\boxed{v(X,y)=\sqrt{\frac gk}\,\mathcal V_h(X,ky).}
$$

The formula for $u$ is not obtained by adding a current to a previously solved zero-current nonlinear wave: $B_j$, $k$, $z_7$ and the surface profile are themselves outputs of the **current-dependent** nonlinear system. The background-current term enters the field reconstruction only *after* those coupled coefficients have been found. The bed condition follows immediately:

$$
\boxed{v(X,-d)=0}
$$

because every $S_j(-D;D)$ vanishes. A result for velocity at the bed is therefore a horizontal velocity only. The mean fixed-frame horizontal velocity at a constant immersed level is

$$
\frac{1}{2\pi}\int_0^{2\pi}u(X,y)\,dX=\bar u_1,
$$

since all nonzero cosine harmonics integrate to zero.

### 10.2 Differential operators and dimensionless derivatives

Because $X=k(x-ct)$ and $Y=ky$,

$$
\frac{\partial}{\partial x}=k\frac{\partial}{\partial X},
\qquad
\frac{\partial}{\partial y}=k\frac{\partial}{\partial Y},
\qquad
\frac{\partial}{\partial t}=-kc\frac{\partial}{\partial X}
$$

for any quantity depending only on the travelling coordinates. Define

$$
A_j=jB_j,\qquad W_j(X,Y)=A_j C_j(Y;D)\cos(jX),
\qquad Z_j(X,Y)=A_j S_j(Y;D)\sin(jX).
$$

The first derivatives of the **physical** velocities are

$$
\boxed{u_x=-\sqrt{gk}\sum_{j=1}^{N}j^2B_jC_j(Y;D)\sin(jX),}
$$

$$
\boxed{u_y=\sqrt{gk}\sum_{j=1}^{N}j^2B_jS_j(Y;D)\cos(jX),}
$$

$$
\boxed{v_x=\sqrt{gk}\sum_{j=1}^{N}j^2B_jS_j(Y;D)\cos(jX),}
$$

$$
\boxed{v_y=\sqrt{gk}\sum_{j=1}^{N}j^2B_jC_j(Y;D)\sin(jX).}
$$

There is an important **sign correction** in the fourth formula: differentiation of $\sin(jX)$ with respect to $Y$ leaves its sign unchanged, whereas $u_x$ differentiates $\cos(jX)$ and receives a minus sign. Consequently

$$
\boxed{u_y=v_x,\qquad u_x+v_y=0.}
$$

These are precisely the irrotationality and incompressibility identities, and they hold identically for every finite Fourier truncation. They constitute useful analytic checks on any code implementing the velocity kernel [1, 2, 10].

### 10.3 Local accelerations

A fixed sensor measures the local derivatives

$$
\boxed{a_{x,\mathrm{local}}=\frac{\partial u}{\partial t}=-cu_x,}
$$

$$
\boxed{a_{y,\mathrm{local}}=\frac{\partial v}{\partial t}=-cv_x.}
$$

Their Fourier expressions are

$$
\boxed{a_{x,\mathrm{local}}=c\sqrt{gk}\sum_{j=1}^{N}j^2B_jC_j\sin(jX),}
$$

$$
\boxed{a_{y,\mathrm{local}}=-c\sqrt{gk}\sum_{j=1}^{N}j^2B_jS_j\cos(jX).}
$$

The distinction between these *local* derivatives and particle-following accelerations is fundamental for hydrodynamic loading calculations. For instance, the separate flow-field output of a historical FOURIER implementation can report local time acceleration, whereas the Rust GUI reports the material horizontal acceleration at its searched extrema. Those two quantities agree only when the advective contribution happens to vanish [10, pp. 3–7].

### 10.4 Material accelerations

The acceleration experienced by a fluid particle is the material derivative of velocity:

$$
\boxed{\frac{Du}{Dt}=u_t+uu_x+vu_y,}
$$

$$
\boxed{\frac{Dv}{Dt}=v_t+uv_x+vv_y.}
$$

For a permanent wave, the explicit time derivatives may be eliminated:

$$
\boxed{a_x=(u-c)u_x+vu_y,}
$$

$$
\boxed{a_y=(u-c)v_x+vv_y.}
$$

Using $v_x=u_y$ and $v_y=-u_x$, the vertical expression can also be written

$$
\boxed{a_y=(u-c)u_y-vu_x.}
$$

These identities avoid finite differencing in time. For a Fourier kernel evaluated at dimensionless $(X,Y)$, define

$$
\widehat u=u\sqrt{\frac kg},\quad
\widehat v=v\sqrt{\frac kg},\quad
\widehat c=c\sqrt{\frac kg},
$$

$$
\widehat u_X=-\sum_{j=1}^{N}j^2B_jC_j\sin(jX),\qquad
\widehat u_Y=\sum_{j=1}^{N}j^2B_jS_j\cos(jX).
$$

Then the dimensionless material accelerations, normalized by $g$, are simply

$$
\boxed{\frac{a_x}{g}=(\widehat u-\widehat c)\widehat u_X+
\widehat v\widehat u_Y,}
$$

$$
\boxed{\frac{a_y}{g}=(\widehat u-\widehat c)\widehat u_Y-
\widehat v\widehat u_X.}
$$

This derivation supplies a strong implementation check: **velocity is scaled by $\sqrt{g/k}$ while material acceleration is scaled by $g$**. No extra factor of $k$ or $1/k$ should remain in the final acceleration after nondimensional differentiation.

### 10.5 A useful vector identity for potential flow

Irrotationality also gives

$$
(\mathbf u\cdot\nabla)\mathbf u
=\nabla\left(\frac12|\mathbf u|^2\right)
-\mathbf u\times(\nabla\times\mathbf u)
=\nabla\left(\frac12|\mathbf u|^2\right).
$$

Thus in a fixed frame

$$
\mathbf a=-c\,\partial_x\mathbf u+\nabla\left(\frac12|\mathbf u|^2\right).
$$

The equivalence of this expression and the direct material-derivative formulas provides an independent analytical consistency check for smooth irrotational waves.

### 10.6 Local and global velocity extrema

The surface is evaluated at $y=\eta(X)$, *not* at the still-water plane $y=0$. Define

$$
U_{s,\max}=\max_{0\le X\le\pi}|u(X,\eta(X))|,
\qquad
U_{b,\max}=\max_{0\le X\le\pi}|u(X,-d)|.
$$

Likewise,

$$
A_{s,\max}=\max_{0\le X\le\pi}|a_x(X,\eta(X))|,
\qquad
A_{b,\max}=\max_{0\le X\le\pi}|a_x(X,-d)|.
$$

A practical numerical search first samples phase and then refines brackets near candidate maxima with a one-dimensional optimizer. The current Rust implementation samples at least $\max(257,8N+1)$ phases and uses golden-section refinement. This is a **high-resolution numerical maximum search**, not a rigorous proof of the global supremum. The value at crest phase is a *signed value* and should not be confused with the maximized absolute magnitude.

***

## 11. Pressure, Bernoulli head and free-surface pressure conditions

### 11.1 Pressure from the wave-frame Bernoulli invariant

The wave-frame Bernoulli relationship is

$$
\frac12[(u-c)^2+v^2]+gy+\frac p\rho=r.
$$

At any point $-d\le y\le\eta(X)$, the gauge pressure is therefore

$$
\boxed{p(X,y)=\rho\left[r-gy-\frac12\left((u(X,y)-c)^2+v(X,y)^2\right)\right].}
$$

With a bed-origin elevation $z=y+d$ and $R=r+gd$,

$$
\boxed{p(X,z)=\rho\left[R-gz-\frac12\left((u-c)^2+v^2\right)\right].}
$$

These formulas are identical, and the latter is convenient for integration over $0\le z\le d+\eta$. Pressure has units of pascals. The combination $R$ has units of $\mathrm{m^2/s^2}$, because it is a specific energy or head multiplied by $g$.

### 11.2 Nondimensional pressure

Using the wave-number scale,

$$
\widehat p=\frac{k p}{\rho g},\qquad \widehat r=\frac{kr}{g}=z_9,
$$

one obtains

$$
\boxed{\widehat p(X,Y)=z_9-Y-\frac12\left[(\mathcal U_h-z_7)^2+\mathcal V_h^2\right].}
$$

Using mean depth as the reference,

$$
\frac{p}{\rho gd}=\frac{z_9}{D}-\frac YD
-\frac{(\mathcal U_h-z_7)^2+\mathcal V_h^2}{2D}.
$$

At the free surface $Y=\zeta(X)$, the dynamic equation demands $p=0$. Thus the interpolated-surface Bernoulli defect may be interpreted directly as a nondimensional **spurious free-surface pressure**:

$$
\boxed{\frac{kp(X,\eta_N(X))}{\rho g}=-\mathcal R_B(X).}
$$

This interpretation makes the continuous off-grid residual test physically transparent. A small collocation defect does not preclude oscillatory pressure errors between nodes if the Fourier surface has insufficient spectral resolution.

### 11.3 Hydrostatic and dynamic components

It is often useful to decompose the pressure into a hydrostatic reference and a dynamic perturbation. Let

$$
p_{\mathrm{hyd}}(y)=\rho g(-y),
$$

which is the resting pressure below still-water level $y=0$. Then

$$
\boxed{p-p_{\mathrm{hyd}}=\rho\left[r-\frac12((u-c)^2+v^2)\right].}
$$

This reference is not unique. It is particularly important not to mistake the wave-frame kinetic term $(u-c)^2$ for the fixed-frame kinetic energy density $u^2$. The distinction matters whenever one reconstructs pressure or calculates the momentum and energy flux in a fixed laboratory frame.

### 11.4 Local and absolute pressure extrema

A pressure value at a given phase and elevation is generally not the maximum pressure at that location. For engineering application one may define

$$
p_{\max}(y)=\max_{X:\,y\le\eta(X)}p(X,y),
\quad
p_{\min}(y)=\min_{X:\,y\le\eta(X)}p(X,y).
$$

The restriction $y\le\eta(X)$ is necessary: points above the instantaneous water surface are in air and cannot be evaluated as wet fluid points. Near the free surface, the set of admissible phases depends on elevation. Consequently a direct phase scan for subsurface pressure must enforce this wet-domain condition.

***

## 12. Volume flux, wave impulse and mass transport

### 12.1 Uniform wave-frame discharge

The wave-frame discharge at any horizontal phase is

$$
\int_{-d}^{\eta(X)}(u-c)\,dy=-Q.
$$

The constancy of this integral can be derived without solving the field equations. Differentiate with respect to $x$:

$$
\frac{d}{dx}\int_{-d}^{\eta}(u-c)\,dy
=\int_{-d}^{\eta}u_x\,dy+(u_s-c)\eta_x.
$$

Since $u_x=-v_y$,

$$
\frac{d}{dx}\int_{-d}^{\eta}(u-c)\,dy
=-v_s+v_b+(u_s-c)\eta_x=0,
$$

using $v_b=0$ and $v_s=(u_s-c)\eta_x$. This is a particularly valuable physical consistency test because it combines continuity, the bed boundary condition and the free-surface kinematic condition.

### 12.2 Laboratory-frame discharge and its mean

At each phase, the fixed-frame discharge is

$$
q_{\mathrm{lab}}(X)=\int_{-d}^{\eta(X)}u\,dy
=-Q+c[d+\eta(X)].
$$

The mean over a wavelength is

$$
\boxed{\overline q_{\mathrm{lab}}=cd-Q=d\bar u_2.}
$$

The wave-frame volume-flux correction $q=\bar U d-Q$ gives an equivalent identity

$$
\boxed{\bar u_2=\bar u_1+\frac qd.}
$$

The distinct definitions have practical consequences. In a closed laboratory tank with no net mass flux, $\bar u_2=0$ but $\bar u_1$ need not vanish. In an open coast with a measured Eulerian current, the solver should impose $\bar u_1=U_c$ and calculate the mass-transport correction, not arbitrarily set both means to the measured value.

### 12.3 Mean horizontal impulse

The mean integrated horizontal momentum (frequently called *impulse* in Fenton's wave-integral notation) per unit horizontal area and unit width is

$$
I=\left\langle\int_{-d}^{\eta}\rho u\,dy\right\rangle.
$$

Using the mean discharge identity,

$$
\boxed{I=\rho(cd-Q)=\rho d\bar u_2.}
$$

The units are $\mathrm{kg\,m^{-1}s^{-1}}$, equivalent to $\mathrm{kg/(m\,s)}$. In the source code the wave-number-scaled equivalent employs the state combination

$$
\boxed{P=z_8+D z_5,}
$$

and

$$
\boxed{I=\rho\sqrt{\frac{g}{k^3}}\,P.}
$$

To verify the identity, $z_8=q\sqrt{k^3/g}$ and $Dz_5=kd\,\bar u_1\sqrt{k/g}=\bar u_1 d\sqrt{k^3/g}$. Therefore

$$
P=(q+\bar u_1d)\sqrt{\frac{k^3}{g}}=(\bar u_2d)\sqrt{\frac{k^3}{g}},
$$

and the dimensional impulse follows.

### 12.4 Mean transport is not the classical small-amplitude particle drift at each depth

The quantity $\bar u_2$ is a **depth-integrated mass-transport velocity**. It should not be identified with the depth-dependent Lagrangian Stokes drift of a fluid particle in the weakly nonlinear approximation. For a narrow-band small-amplitude wave in quiescent finite-depth water, the classical leading-order Stokes drift profile is

$$
\boxed{u_{\mathrm{SD}}(y)=\frac{a^2\omega k}{2}
\frac{\cosh[2k(y+d)]}{\sinh^2(kd)},\qquad a=\frac H2,}
$$

under the standard no-mean-current convention to second order in amplitude. Its vertical structure is not resolved by simply quoting the single scalar $\bar u_2$. Higher-order particle trajectories also differ from first-order sinusoidal orbits and should be obtained by integrating the nonlinear velocity field in time if particle drift is required.

***

## 13. Kinetic energy, potential energy and total energy

### 13.1 Mean potential energy relative to still water

The gravitational potential energy of a vertical water column is obtained by integrating $\rho gy$ from the bed to the free surface, subtracting the energy of a resting column of constant depth. At phase $X$ the difference is

$$
\Delta V(X)=\rho g\left[\int_{-d}^{\eta(X)}y\,dy
-\int_{-d}^0 y\,dy\right]
=\frac12\rho g\eta(X)^2.
$$

Averaging over one wavelength therefore gives

$$
\boxed{V=\frac{\rho g}{2}\langle\eta^2\rangle
=\frac{\rho g}{4\pi}\int_0^{2\pi}\eta(X)^2\,dX.}
$$

The physical units are $\mathrm{J/m^2}$ for the wave-energy quantity per unit horizontal sea-surface area. This is positive even though $\langle\eta\rangle=0$.

### 13.2 Parseval's identity and the last Fourier harmonic

For

$$
\zeta_N(X)=\sum_{j=1}^{N-1}E_j\cos(jX)+\frac{E_N}{2}\cos(NX),
$$

orthogonality gives

$$
\left\langle\zeta_N^2\right\rangle
=\frac12\sum_{j=1}^{N-1}E_j^2+\frac18 E_N^2.
$$

Because $\eta=\zeta/k$,

$$
\boxed{V=\frac{\rho g}{k^2}\left[
\frac14\sum_{j=1}^{N-1}E_j^2+\frac1{16}E_N^2\right].}
$$

The last harmonic receives the $1/16$ coefficient rather than $1/4$. The reason is purely the discrete cosine interpolation convention, not a physical suppression of the last mode. This exact finite-Fourier formula is a useful alternative to surface-elevation quadrature and is implemented by the solver [2, 11].

### 13.3 Mean kinetic energy in the fixed frame

The dimensional kinetic energy per unit horizontal area is

$$
\boxed{T_K=\left\langle\int_{-d}^{\eta(X)}
\frac12\rho\left[u(X,y)^2+v(X,y)^2\right]dy\right\rangle.}
$$

The formula is in the laboratory frame and includes any kinetic contribution associated with the Eulerian background current. By writing $u=U+c$, one can relate it to the moving-frame momentum and discharge invariants. For a steady wave in Fenton's convention, the integrated identity is

$$
\boxed{T_K=\frac12\left(cI-\rho\bar u_1 Q\right).}
$$

This relation has correct dimensions: both $cI$ and $\rho\bar u_1Q$ have units $\mathrm{J/m^2}$. The right side depends on which laboratory frame is chosen, as any kinetic energy naturally does. It is not the same as the kinetic energy of the perturbation with the background current removed.

### 13.4 Kinetic energy in terms of the finite Fourier state

Using the nondimensional values

$$
C=z_4,\quad E=z_5,\quad A=z_7,\quad D=z_1,\quad q^*=z_8,
$$

and defining the dimensionless impulse-like combination

$$
P=q^*+DE,
$$

as well as the flux combination

$$
\mathcal Q=DA-q^*,
$$

the energy in **wave-number units** is

$$
\boxed{\widehat T_K=\frac12(CP-E\mathcal Q),}
$$

where

$$
\boxed{T_K=\frac{\rho g}{k^2}\widehat T_K.}
$$

To establish this, use $I=\rho\sqrt{g/k^3}P$, $Q=\sqrt{g/k^3}\mathcal Q$, $c=C\sqrt{g/k}$ and $\bar u_1=E\sqrt{g/k}$ in the dimensional energy identity. The resulting common factor is $\rho g/k^2$. This provides a precise correspondence between Fenton's integral invariant and the solver's `kinetic_energy` field.

### 13.5 Total energy and a note on background current

The total mean energy reported by the wave solver is

$$
\boxed{E=T_K+V.}
$$

Because $T_K$ is calculated with laboratory-frame absolute velocities, $E$ in a current-bearing case includes the current-dependent contribution implied by that convention. When comparing different current strengths, a physically meaningful **wave-excess energy** may require subtracting a consistently defined uniform-flow reference. One must not blindly interpret every current-bearing $E$ as the wave perturbation energy of a quiescent sea. Definitions must be matched before computing transmission or energy-flux coefficients.

### 13.6 Linear small-amplitude limit

For $\eta=a\cos(kx-\omega t)$ with $a=H/2$, the potential energy is

$$
V=\frac12\rho g\left\langle a^2\cos^2\theta\right\rangle
=\frac14\rho ga^2=\frac1{16}\rho gH^2.
$$

With zero background current, linear wave theory also gives

$$
T_K=V+O(H^3),
$$

and hence the leading-order energy

$$
\boxed{E=\frac18\rho gH^2+O(H^3).}
$$

This is the familiar regular-wave energy relationship. The absence of a background current and the convention for measuring wave energy are needed for its direct application.

***

## 14. Bed velocities, radiation stress and momentum flux

### 14.1 Bed orbital velocity and its mean square

At the bed $Y=-D$, the harmonic vertical-velocity factors vanish, and the horizontal fixed-frame velocity is

$$
\boxed{u_b(X)=\bar u_1+\sqrt{\frac gk}
\sum_{j=1}^{N}jB_j\operatorname{sech}(jD)\cos(jX).}
$$

Define the *orbital* bed velocity as the oscillatory component obtained after subtracting the solved Eulerian current:

$$
\widetilde u_b(X)=u_b(X)-\bar u_1.
$$

Its mean is zero. Orthogonality of the cosine harmonics implies

$$
\frac{1}{2\pi}\int_0^{2\pi}\cos(jX)\cos(lX)\,dX=
\begin{cases}1/2,&j=l\ge1,\\0,&j\ne l.\end{cases}
$$

Therefore

$$
\boxed{\left\langle\widetilde u_b^2\right\rangle
=\frac{g}{2k}\sum_{j=1}^{N}
\left[jB_j\operatorname{sech}(jD)\right]^2.}
$$

The RMS orbital velocity is

$$
\boxed{u_{b,\mathrm{orb,RMS}}=
\sqrt{\left\langle\widetilde u_b^2\right\rangle}.}
$$

This closed finite-sum formula avoids hundreds of evaluations of an oscillatory depth-dependent velocity field. It is both faster and, for the retained Fourier model, more accurate than phase sampling. It is one of the nine engineering quantities independently re-evaluated by the interval checker.

### 14.2 Absolute mean-square bed velocity

The **absolute** bed velocity is $u_b=\bar u_1+\widetilde u_b$, hence

$$
\left\langle u_b^2\right\rangle
=\bar u_1^2+2\bar u_1\left\langle\widetilde u_b\right\rangle
+\left\langle\widetilde u_b^2\right\rangle.
$$

Since the orbital average vanishes,

$$
\boxed{\left\langle u_b^2\right\rangle=
\bar u_1^2+\left\langle\widetilde u_b^2\right\rangle.}
$$

In Fenton's wave-integral identities the same mean square can also be recovered from the Bernoulli constant and wave speed:

$$
\boxed{\left\langle u_b^2\right\rangle
=2r-c(c-2\bar u_1).}
$$

**Derivation.** At the bed, $v=0$ and $y=-d$; the Bernoulli equation gives

$$
\frac12(u_b-c)^2-gd+\frac{p_b}{\rho}=r.
$$

This equation by itself contains the bed pressure $p_b$, so the mean-square identity is **not** obtained merely by setting the bed pressure to zero. It follows from the steady-wave integral identities and their relation to the Fourier state. Equivalence of the Bernoulli/invariant expression and the Fourier-orthogonality expression is a valuable independent check. One must not confuse the two means in the presence of current [3, §7].

### 14.3 Radiation stress $S_{xx}$

The radiation-stress component in the propagation direction is defined by a wavelength-average of the horizontal momentum flux after removing the hydrostatic still-water contribution:

$$
\boxed{S_{xx}=\left\langle\int_{-d}^{\eta(X)}
\left[p(X,y)+\rho u(X,y)^2\right]dy\right\rangle
-\frac12\rho gd^2.}
$$

For the present fully nonlinear, steady and irrotational model, the integral can be rewritten as

$$
\boxed{S_{xx}=4T_K-3V+
\rho d\left\langle u_b^2\right\rangle-2\bar u_1 I.}
$$

All terms have units of force per unit crest length, $\mathrm{N/m}$. The expression emphasizes why the *absolute* mean-square bed velocity, rather than the orbital mean square alone, is required in the general wave-current formula. In zero current the terms involving $\bar u_1$ vanish, simplifying the expression [3, §7, p. 16].

### 14.4 Radiation stress in the Fourier state

Using the same state combinations as §13, let

$$
\widehat V=\frac{k^2V}{\rho g},\qquad
\widehat T_K=\frac{k^2T_K}{\rho g},
$$

and define

$$
\boxed{\beta=2z_9-z_4^2+2z_4z_5.}
$$

The physical meaning is

$$
\boxed{\beta=\frac{k}{g}\left\langle u_b^2\right\rangle.}
$$

Then

$$
\boxed{\widehat S_{xx}=4\widehat T_K-3\widehat V
+\beta D-2z_5P,}
$$

$$
\boxed{S_{xx}=\frac{\rho g}{k^2}\widehat S_{xx}.}
$$

Every term in $\widehat S_{xx}$ is dimensionless. The formula is useful as a regression target because it is assembled directly from a small number of state entries and a single Parseval surface sum.

### 14.5 Momentum flux in the moving frame

The phase-independent momentum-flux invariant of a steady wave is a different quantity from $S_{xx}$:

$$
\boxed{S=\int_{-d}^{\eta(X)}
\left[p(X,y)+\rho(u(X,y)-c)^2\right]dy.}
$$

For an exact steady solution this quantity is independent of $X$. To see why, integrate the steady horizontal Euler momentum equation over a vertical cross section, accounting for the fact that pressure is zero on the actual sloping surface. Its horizontal flux derivative vanishes, leaving a constant flux along the wave.

Expanding the quadratic $(u-c)^2$ and using the discharge identity yields the relationship

$$
\boxed{S=S_{xx}-2cI+\rho d c^2+\frac12\rho gd^2.}
$$

The difference between $S$ and $S_{xx}$ includes the wave-frame transformation and the hydrostatic reference. Therefore a table entry called *Momentum flux* must not be interpreted as another label for *Radiation stress*. The Rust program also evaluates the momentum-flux defining integral at a representative phase, whereas the radiation stress uses a wave-averaged invariant. The independence of $S$ from phase provides an additional cross-check when the finite Fourier equations are satisfied accurately.

### 14.6 Linear limit of radiation stress

For a small-amplitude monochromatic wave with no mean current, denote the linear intrinsic phase speed by $c_0=\omega/k$ and the linear group speed by $c_{g0}$. Let

$$
n=\frac{c_{g0}}{c_0}=\frac{1}{2}\left(1+\frac{2kd}{\sinh(2kd)}\right).
$$

Then to leading order the classical radiation-stress relation is

$$
\boxed{S_{xx}=E\left(2n-\frac12\right)+O(H^3).}
$$

In deep water $n\to1/2$, so $S_{xx}\to E/2$. In shallow water $n\to1$, so $S_{xx}\to3E/2$. These limiting identities are valuable for verification of the small-amplitude Fourier calculations [3, 4].

***

## 15. Wave power, energy flux and energy-transport velocity

### 15.1 Exact mean energy flux in the laboratory frame

The instantaneous horizontal flux of mechanical energy per unit crest width is the horizontal velocity multiplied by total pressure-plus-mechanical-energy density. In the mean-level coordinate system,

$$
\mathscr F(X,y)=
\left[p+\rho gy+\frac12\rho(u^2+v^2)\right]u.
$$

The mean energy flux (wave power) is

$$
\boxed{F=\left\langle\int_{-d}^{\eta(X)}\mathscr F(X,y)\,dy\right\rangle.}
$$

The units are $\mathrm{W/m}$. A term associated with the choice of gravitational datum may be added to a general current-bearing energy flux; physical comparisons require a consistent datum and reference flow. The formula above uses the same still-water-level datum as $r$ and the post-processing adopted here.

### 15.2 Nonlinear integral identity for power

Using the wave-frame Bernoulli relationship, mass transport and the identities for $T_K$, $V$ and the bed mean square, Fenton's steady-wave result is

$$
\boxed{F=c(3T_K-2V)+\frac12\left\langle u_b^2\right\rangle
(I+\rho cd)-2c\bar u_1I.}
$$

The dimensions of the terms are consistent:

$$
[cT_K]=[\langle u_b^2\rangle I]=\mathrm{kg\,s^{-3}}=\mathrm{W/m}.
$$

This formula includes the general current-related terms. A zero-current formula cannot be transferred to a wave with a nonzero prescribed current by merely changing $c$.

### 15.3 Expression in the numerical state

With $P=z_8+Dz_5$, $C=z_4$, $\beta=2z_9-C^2+2Cz_5$ and the dimensionless kinetic and potential energies already defined, the dimensionless power is

$$
\boxed{\widehat F=C(3\widehat T_K-2\widehat V)
+\frac12\beta(P+CD)-2Cz_5P.}
$$

The dimensional scaling is

$$
\boxed{F=\rho\left(\frac{g}{k}\right)^{3/2}\frac1k\,\widehat F
=\frac{\rho g^{3/2}}{k^{5/2}}\widehat F.}
$$

Equivalently, in depth-based nondimensional form,

$$
\boxed{\widehat F_d=\frac{F}{\rho g^{3/2}d^{5/2}}
=\frac{\widehat F}{(kd)^{5/2}}.}
$$

This is the quantity commonly printed in the nineteenth row of a Fenton `SOLUTION.RES`-style table.

### 15.4 Energy-flux-to-energy ratio

The program also reports

$$
\boxed{C_E=\frac FE,}
$$

with units $\mathrm{m/s}$ provided $E\ne0$. It is often named an *energy-flux velocity*. **This number should not automatically be called the linear group velocity**. For finite-amplitude waves, energy and flux depend on nonlinear wave shape, current convention and frame; with a strong background current the ratio may be strongly affected by non-wave transport terms. For engineering reporting, $F/E$ is the precise mathematical definition being used.

### 15.5 Linear group velocity as a separate concept

For linear waves with no background current,

$$
\omega_0(k)=\sqrt{gk\tanh(kd)}.
$$

The phase speed is

$$
\boxed{c_0=\frac{\omega_0}{k}=\sqrt{\frac gk\tanh(kd)}.}
$$

Differentiation with respect to wavenumber gives

$$
\begin{aligned}
\frac{d\omega_0}{dk}
&=\frac{g\left[\tanh(kd)+kd\operatorname{sech}^2(kd)\right]}{2\omega_0}\\
&=\frac{c_0}{2}\left[1+\frac{2kd}{\sinh(2kd)}\right].
\end{aligned}
$$

Thus

$$
\boxed{c_{g0}=\frac{c_0}{2}
\left[1+\frac{2kd}{\sinh(2kd)}\right].}
$$

For small-amplitude waves without current,

$$
F=c_{g0}E+O(H^3).
$$

The limiting cases are $c_{g0}/c_0\to1$ as $kd\to0$ and $c_{g0}/c_0\to1/2$ as $kd\to\infty$. In a weak uniform current the linear laboratory-frame group speed gains a $U_c$ Doppler term, but a fully nonlinear steady-wave problem requires more careful treatment of the conserved energy and action when comparing wave packets [4, 8].

***

## 16. Wave limiting height, steepness and nonlinearity indicators

### 16.1 Independent dimensionless measures

Three useful dimensionless groups are

$$
\boxed{\text{relative depth}\;\mu=\frac dL,\qquad
\text{wave steepness}\;s=\frac HL,\qquad
\text{relative height}\;h_r=\frac Hd.}
$$

They are not independent: $s=\mu h_r$. The wavenumber-depth parameter is

$$
D=kd=2\pi\frac dL=2\pi\mu.
$$

A wave with $H/d$ large may be strongly nonlinear despite modest apparent steepness $H/L$ if the wavelength is long. Conversely, a short deep-water wave may be steep although $H/d$ is small.

### 16.2 Ursell number

The Ursell number measures the relative importance of nonlinearity and shallowness:

$$
\boxed{Ur=\frac{HL^2}{d^3}=
\frac{H/d}{(d/L)^2}.}
$$

In terms of the wavenumber-scaled variables,

$$
\boxed{Ur=4\pi^2\frac{H/d}{(kd)^2}.}
$$

Fenton also discusses the Stokes–Ursell number

$$
\boxed{SU=\frac{kH/2}{(kd)^3}=
\frac{Ur}{8\pi^2}.}
$$

These parameters are **regime indicators**, not universal acceptance thresholds: the actual accuracy of an approximate theory depends on wave height, depth and the specific quantity to be predicted. The fully nonlinear Fourier method does not require $Ur$ to be small, but very high Ursell numbers can necessitate many Fourier modes and can drive numerical ill-conditioning [3, §6; 11, §3].

### 16.3 Miche-type height screening

A commonly used empirical wave-height indicator is

$$
\boxed{H_{\mathrm{Miche}}=0.142L\tanh(kd),}
$$

with ratio

$$
\boxed{\mathcal B_M=\frac{H}{H_{\mathrm{Miche}}}.}
$$

This is a **screening calculation**, not a unique determination of the exact nonlinear limiting wave at arbitrary depth, current and geometry. It does not incorporate wave breaking as a dissipative dynamical process. A computed wave satisfying the inviscid equations can still be physically unstable or susceptible to breaking under perturbations.

### 16.4 Fenton's empirical highest-steady-wave relation

A well-known fit to numerically determined highest steady waves, quoted by Fenton (1990, Eq. 32) and in his `FOURIER` program instructions, uses the dimensionless wavelength $\lambda=L/d$:

$$
\boxed{
\frac{H_m}{d}=
\frac{0.141063\lambda+0.0095721\lambda^2+0.0077829\lambda^3}
{1+0.0788340\lambda+0.0317567\lambda^2+0.0093407\lambda^3}.}
$$

This is a **rational fit**, not a rigorous mathematical theorem for the existence of a periodic solution. It may be evaluated stably for $\lambda>1$ by dividing numerator and denominator by $\lambda^3$:

$$
\boxed{
\frac{H_m}{d}=
\frac{0.0077829+0.0095721/\lambda+0.141063/\lambda^2}
{0.0093407+0.0317567/\lambda+0.0788340/\lambda^2+1/\lambda^3}.}
$$

At the long-wave limit,

$$
\lim_{\lambda\to\infty}\frac{H_m}{d}
=\frac{0.0077829}{0.0093407}\approx0.8332.
$$

The ratio $H/H_m$ provides a useful theoretical proximity indicator. The approximation $H\lesssim0.98H_m$ is sometimes used as a practical spectral-solution guideline, but the 98% factor is not a general breaking threshold. Near-limiting waves have sharp crests and broad Fourier spectra; a numerical method may lose spectral accuracy **before** the nominal theoretical limit [3, §6; 11, §3].

### 16.5 Crest stagnation

At the crest of a wave propagating in the positive direction, define the laboratory-frame crest velocity $u_c=u(0,\eta(0))$. The crest stagnation margin is

$$
\boxed{\Delta u_{\mathrm{crest}}=c-u_c.}
$$

In the travelling wave frame, the corresponding horizontal velocity is $u_c-c=-\Delta u_{\mathrm{crest}}$. As the limiting wave is approached, the crest may tend toward stagnation in the moving frame. A small positive margin is therefore a useful diagnostic, but numerical ill-conditioning, finite harmonic resolution and the physical stability of the wave must also be considered.

### 16.6 Theoretical height versus practical breaking

The theory of a highest steady wave concerns a mathematical family of periodic gravity-wave solutions of an ideal fluid. Field and laboratory breaking also involve perturbations, local slopes, turbulence, bed changes, wave groups and three-dimensional effects that are absent from the steady inviscid model. A simple $H/d=0.55$ engineering screen, sometimes associated with observed shallow-water breaking, is neither an exact theorem nor universally appropriate to every bathymetry and wave/current condition. Thus the program's numerical acceptance, the Fenton highest-wave ratio, Miche screening and observed breaking criteria answer **different questions** and should be reported separately.

***

## 17. Numerical verification: residuals, spectral resolution and integral identities

### 17.1 A hierarchy of questions, not one tolerance

An accurate scientific calculation must answer several distinct questions:

1. **Equation solve:** do the discrete unknowns satisfy the $2N+10$ nonlinear collocation equations?
2. **Surface resolution:** does the reconstructed surface satisfy the free-surface conditions *between* collocation nodes?
3. **Derived quantities:** are energy, momentum and power calculated consistently with independent defining expressions?
4. **Arithmetic reliability:** do the reported finite-precision values lie within conservative enclosures of their defining formulas?
5. **Model suitability:** is the idealized steady, smooth, irrotational wave an appropriate representation of the actual physical case?

A positive answer to the first question does not automatically imply positive answers to the other four. In particular, a small residual is not a proof of small error in the state when the Jacobian is nearly singular. Nor does a finite system's numerical accuracy establish that the truncated Fourier representation is a uniformly accurate approximation of the exact infinite-dimensional water-wave solution.

### 17.2 Collocation residual norm

Let $M=2N+10$. The maximum norm is

$$
\boxed{R_\infty=\max_{1\le i\le M}|F_i(\mathbf z)|.}
$$

The root-mean-square residual is

$$
\boxed{R_2=\sqrt{\frac1M\sum_{i=1}^MF_i(\mathbf z)^2}.}
$$

These norms are related by

$$
\boxed{R_2\le R_\infty\le\sqrt M\,R_2.}
$$

For the same numerical residual vector, an RMS criterion and a maximum criterion are not interchangeable unless the factor $\sqrt M$ is accounted for. The Rust solver uses a maximum residual limit $10^{-8}$ at intermediate continuation stages and $10^{-10}$ at the final requested wave height. It also requires sufficiently small Newton corrections.

### 17.3 Re-evaluation after dimensional normalization

The nonlinear solve uses dimensionless state variables, while engineering outputs are dimensional. After the converged value $D=kd$ has been obtained, calculate

$$
\boxed{k=\frac Dd,\qquad L=\frac{2\pi}{k},\qquad c=\frac LT.}
$$

The period and celerity variables can then be reconstructed from physical inputs:

$$
z_3^{\mathrm{check}}=T\sqrt{gk},\qquad
z_4^{\mathrm{check}}=c\sqrt{\frac kg}.
$$

For the Eulerian-current case,

$$
z_5^{\mathrm{check}}=U_c\sqrt{\frac kg},\qquad
z_7^{\mathrm{check}}=z_4^{\mathrm{check}}-z_5^{\mathrm{check}}.
$$

Because normalization introduces additional rounded floating-point operations, the complete nonlinear residual is evaluated once more at the reconstructed state. If the final residual is no longer within tolerance, the computation is rejected. A printed result must not be assumed consistent merely because Newton iteration previously converged.

### 17.4 Independent surface-interpolation test

At arbitrary phase $X$, calculate the reconstructed elevation $\zeta_N(X)$, the stream-function Fourier sum $\mathcal P(X)$ and harmonic velocities $\mathcal U_h(X)$, $\mathcal V_h(X)$. Define the off-grid residuals

$$
\boxed{R_K(X)=\mathcal P(X)-z_8-z_7\zeta_N(X),}
$$

$$
\boxed{R_B(X)=\frac12\left[(\mathcal U_h-z_7)^2+\mathcal V_h^2\right]
+\zeta_N(X)-z_9.}
$$

A discrete phase-scan diagnostic is

$$
\boxed{R_{\mathrm{surf}}=\max_{X\in\mathcal G}
\max\left(|R_K(X)|,|R_B(X)|\right),}
$$

where $\mathcal G$ is a fine, independent phase grid. The Rust code uses at least $\max(2048,32N)$ phases across a wavelength and a nominal acceptance threshold of $10^{-6}$. At the collocation points $X_m$, both residuals are small by construction; the off-grid scan is intended to reveal the missing information **between** those nodes. Scanning does not establish a rigorous bound on all $X\in[0,2\pi]$ unless supplemented by derivative/Lipschitz or interval arguments.

### 17.5 Defining-integral verification

The analytic spectral formulas for $I,T_K,V,S_{xx},F$ may be compared with direct quadrature of their defining integrals. Let

$$
\mathcal I_A=\frac1{2\pi}\int_0^{2\pi}\int_{-d}^{\eta(X)}
\mathscr F_A(X,y)\,dy\,dX,
$$

where $A$ denotes the selected integral property and $\mathscr F_A$ its dimensional integrand. Introduce a Gauss–Legendre coordinate $\xi\in[-1,1]$:

$$
y(X,\xi)=-d+\frac{d+\eta(X)}{2}(1+\xi),
$$

$$
\boxed{dy=\frac{d+\eta(X)}{2}\,d\xi.}
$$

For $n_y$ quadrature nodes $\xi_i$ and weights $w_i$,

$$
\int_{-d}^{\eta(X)}f(X,y)\,dy
\approx\frac{d+\eta(X)}{2}\sum_{i=1}^{n_y}
 w_i f\left(X,-d+\frac{d+\eta(X)}{2}(1+\xi_i)\right).
$$

A periodic phase trapezoidal rule then approximates the outer integral:

$$
\frac1{2\pi}\int_0^{2\pi}g(X)\,dX
\approx\frac1{n_x}\sum_{m=0}^{n_x-1}
 g\left(\frac{2\pi m}{n_x}\right).
$$

For smooth periodic functions this rule can converge very rapidly. For a near-limiting, poorly resolved crest its error can be much larger, so increasing $n_x$ and $n_y$ is essential for a meaningful check.

### 17.6 Coarse/fine quadrature comparison

The production implementation performs independent integrations using two discretizations:

$$
(n_x,n_y)_{\mathrm{coarse}}
=\bigl(\max(256,8N),64\bigr),
$$

$$
(n_x,n_y)_{\mathrm{fine}}
=\bigl(\max(512,16N),96\bigr).
$$

For a reported quantity $q_{\mathrm{spectral}}$, a fine integral value $q_f$ and a coarse result $q_c$, define

$$
\boxed{e_{\mathrm{formula}}=|q_{\mathrm{spectral}}-q_f|,
\qquad e_{\mathrm{quadrature}}=|q_c-q_f|.}
$$

The numerical budget for a dimensional reference quantity is

$$
\boxed{\varepsilon(q)=10^{-7}+10^{-8}|q|.}
$$

Both differences must fit the chosen budget. The absolute term is in the *specific SI units of the compared quantity*; it is not a dimensionless universal precision theorem. A coarse/fine agreement may fail to detect a common systematic bias shared by both discretizations, so the additional analytic identity provides an important independent test.

### 17.7 Spectral momentum and energy-flux identities

At the stored state, reconstruct the independent bed mean square

$$
\beta_b=\left\langle u_b^2\right\rangle
=\bar u_1^2+\frac g{2k}\sum_{j=1}^{N}
(jB_j\operatorname{sech}(jD))^2.
$$

Then compute the alternative expressions

$$
\boxed{S_{xx}^{(\mathrm{ind})}
=4T_K-3V+\rho\beta_b d-2\bar u_1I,}
$$

$$
\boxed{F^{(\mathrm{ind})}
=c(3T_K-2V)+\frac12\beta_b(I+\rho cd)-2c\bar u_1I.}
$$

Compare these with the spectral state formulas for radiation stress and wave power. The two routes use the same physical theory but partly different combinations of the numerical state, so discrepancies can detect an inconsistency in the Fourier coefficients or post-processing conventions.

### 17.8 Residual and output *ratios* as acceptance diagnostics

A normalized discrepancy or ratio expresses an error relative to its allowance:

$$
\boxed{r_q=\frac{|q_{\mathrm{reported}}-q_{\mathrm{reference}}|}
{10^{-7}+10^{-8}|q_{\mathrm{reference}}|}.}
$$

Then $r_q\le1$ meets that check. The numerical code stores several metrics of this type, including spectral comparison ratios, coarse/fine quadrature ratios and finite-formula enclosure ratios. Such ratios are **diagnostic compliance measures** and should not be interpreted as frequentist confidence probabilities or statements such as “99.9% certain.” A ratio less than one indicates that the particular tested mathematical discrepancy falls inside its chosen tolerance, no more and no less.

***

## 18. Finite-precision arithmetic, interval enclosures and the meaning of verification

### 18.1 Finite precision is a separate source of uncertainty

An ideal exact arithmetic evaluation of the finite Fourier residual is distinct from the result of evaluating that residual on a binary computer. Normal IEEE-754 binary64 arithmetic has 53 bits of significand precision. Its machine epsilon is

$$
\boxed{\epsilon_{\mathrm{mach}}=2^{-52}\approx2.220446049250313\times10^{-16}.}
$$

For a normal, correctly rounded elementary operation, the standard first-order error model is

$$
\operatorname{fl}(a\circ b)=(a\circ b)(1+\delta),
\qquad |\delta|\lesssim\epsilon_{\mathrm{mach}},
$$

under the usual assumptions of no overflow, troublesome underflow, or exceptional rounding. This heuristic model does not, by itself, rigorously bound an arbitrarily long sequence of nonlinear operations.

### 18.2 Interval values and inclusion

Represent a real quantity by an interval

$$
\boxed{A=[a_-,a_+],\qquad a_-\le a_+.}
$$

An interval expression is useful if it is *inclusion isotonic*: when every input belongs to its interval, every valid exact result is contained in the resulting output interval. For addition and subtraction,

$$
A+B=[a_-+b_-,a_++b_+],
$$

$$
A-B=[a_--b_+,a_+-b_-],
$$

with endpoint rounding directed outward. For multiplication,

$$
\boxed{AB=\left[
\min\{a_-b_-,a_-b_+,a_+b_-,a_+b_+\},
\max\{a_-b_-,a_-b_+,a_+b_-,a_+b_+\}
\right].}
$$

Division is admissible only if the denominator interval excludes zero. If $0\notin B$,

$$
\boxed{A/B=A\cdot[1/b_+,1/b_-].}
$$

The reciprocal endpoints must again be enlarged outward. An interval implementation that silently divides by an interval crossing zero is mathematically invalid for these checks.

### 18.3 Adjacent representable numbers

Let $\operatorname{nextdown}(x)$ and $\operatorname{nextup}(x)$ denote the adjacent binary64 values immediately below and above a finite $x$. The Rust implementation evaluates an ordinary floating-point endpoint operation, then enlarges the result by one representable number in the outward direction. In the absence of exceptional arithmetic behavior, this supplies conservative enclosures for the elementary rounded operations used by the checker. Exact zero and signed zero require special handling because their adjacent values are subnormal.

One additional unit in the last place is not, by itself, a rigorous error bound for all complex transcendental libraries. For functions such as sine, cosine and exponential, the implementation therefore uses explicitly bounded Taylor constructions instead of relying solely on ordinary system-library `sin`, `cos` or `exp` results. The scope and validity of those constructions must still be checked against the arithmetic assumptions of the compiler and hardware.

### 18.4 Exponential enclosure by range reduction

For $x\ge0$, choose an integer $m$ such that

$$
r=\frac{x}{2^m},\qquad 0\le r\le\frac18.
$$

Use the truncated Taylor polynomial

$$
P_n(r)=\sum_{j=0}^n\frac{r^j}{j!},
$$

and the remainder inequality

$$
\boxed{0\le e^r-P_n(r)\le\frac{e^r r^{n+1}}{(n+1)!}.}
$$

Since $r\le1/8$, one may choose a convenient elementary upper bound for $e^r$ and use outward-rounded interval arithmetic to evaluate the remainder. Recover $e^x$ by $m$ squarings and $e^{-x}$ by reciprocal, with the required interval/domain checks. The Rust checker uses a 32-term Taylor evaluation together with a remainder enlargement and repeated squaring.

### 18.5 Trigonometric enclosure at collocation nodes

For the collocation phases

$$
X_m=\frac{m\pi}{N},
$$

one first encloses $\pi$ by a narrow interval around its binary64 representation. Symmetry folds phases into an interval such as $[-\pi,\pi]$, after which Taylor series give

$$
\sin x=\sum_{j=0}^{n}(-1)^j\frac{x^{2j+1}}{(2j+1)!}+R_{\sin},
$$

$$
\cos x=\sum_{j=0}^{n}(-1)^j\frac{x^{2j}}{(2j)!}+R_{\cos}.
$$

The omitted terms are enclosed by remainder bounds valid on the reduced argument interval. Because the Fourier residual contains products of trigonometric factors and exponentials, these interval evaluations need to be internally consistent with the algebra used by the residual checker.

### 18.6 Square-root enclosure and domain control

For $A=[a_-,a_+]$ with $a_-\ge0$,

$$
\sqrt A=[\sqrt{a_-},\sqrt{a_+}],
$$

again with outward-rounded endpoints. One can verify that a proposed lower endpoint squared does not exceed $a_-$ and that a proposed upper endpoint squared is not less than $a_+$. If the bracket cannot be verified within its finite adjustment budget, the checker should fail rather than return a spurious certificate. Inputs that cause a negative radicand interval must likewise be rejected.

### 18.7 Finite residual enclosure

Let $\mathcal F_i(\mathbf z)$ be an interval evaluation of the $i$th nonlinear equation at the stored finite vector. For an interval $A=[a_-,a_+]$ define

$$
\operatorname{mag}(A)=\max(|a_-|,|a_+|).
$$

The checker forms

$$
\boxed{R_{\mathrm{interval}}=
\max_{1\le i\le M}\operatorname{mag}\bigl(\mathcal F_i(\mathbf z)\bigr).}
$$

The nominal acceptance test is $R_{\mathrm{interval}}\le10^{-10}$, with the tolerance rounded conservatively. This bounds the evaluation of the **finite system at the stored state**, subject to the interval implementation's arithmetic correctness; it is not a bound on the distance between the stored state and an exact root.

### 18.8 Independent enclosures for selected output formulas

For a reported scalar $q_r$ and interval enclosure $Q=[q_-,q_+]$, define the worst-case discrepancy

$$
\boxed{e_+(q_r,Q)=\max(|q_r-q_-|,|q_r-q_+|).}
$$

The smallest possible magnitude in that interval is

$$
m(Q)=\begin{cases}
0,&q_-\le0\le q_+,\\
\min(|q_-|,|q_+|),&\text{otherwise}.
\end{cases}
$$

A conservative output budget is

$$
\boxed{b_-(Q)=\operatorname{down}
\left(10^{-7}+10^{-8}m(Q)\right).}
$$

The output test accepts only if

$$
\boxed{e_+(q_r,Q)\le b_-(Q).}
$$

The implementation includes nine outputs:

| Output name | Physical quantity | Dimensional units |
|:---|:---|:---|
| `L` | Nonlinear wavelength | m |
| `c` | Fixed-frame phase speed | m/s |
| `Impulse` | Integrated mean horizontal momentum | kg/(m s) |
| `KineticEnergy` | Mean kinetic energy | J/m² |
| `PotentialEnergy` | Mean potential energy | J/m² |
| `EnergyDensity` | Sum of kinetic and potential energy | J/m² |
| `Sxx` | Horizontal radiation stress | N/m |
| `Power` | Mean energy flux | W/m |
| `MeanSquareBedVelocity` | Mean square **orbital** bed velocity | m²/s² |

These nine do **not** include every quantity shown in the engineering report, and their finite-formula enclosures should not be described as certification of all kinematic and breaking indicators.

### 18.9 Why the certificate does not prove an exact root exists

Suppose $\mathbf z_*\in\mathbb R^M$ is a numerical state with a small residual. To bound an exact nearby root $\mathbf z^\dagger$ one additionally needs information about invertibility and variation of the Jacobian in a region. For example, a Newton–Kantorovich or interval-Newton argument might use a verified bound on

$$
\|\mathbf J(\mathbf z)^{-1}\|,
$$

and a Lipschitz or derivative enclosure for $\mathbf J$ throughout a neighborhood. A small $\|\mathbf F(\mathbf z_*)\|$ alone does not supply such information. The Rust certificate accordingly records that **root existence has not been proved**. A further gap remains between a root of the $M$-equation truncated Fourier problem and the exact infinite-dimensional water-wave boundary-value problem.

### 18.10 Precision-preserving compilation

The numerical core uses Rust `f64` throughout, including its interval endpoints. When compiling for a target CPU, the developer should preserve normal IEEE floating-point semantics and avoid enabling transformations that freely reassociate expressions. Hardware fused multiply-add can produce a different rounding outcome than separately evaluated multiplication and addition; that may matter to interval algorithms built around a specific operation sequence. The production GNU command therefore uses `-C target-feature=-fma` with `-C panic=unwind`. CPU-specific optimization, link-time optimization and single code-generation-unit compilation may improve throughput, but accuracy remains an empirical and analytical validation question, not a property established solely by the optimizer switches.

***

## 19. Stokes perturbation theory and the finite-amplitude dispersion correction

### 19.1 Mathematical role of Stokes theory

Stokes theory and the fully nonlinear Fourier approximation share trigonometric basis functions, but differ in the treatment of the coefficients. A Stokes expansion assumes an amplitude parameter $\epsilon$ sufficiently small that the unknown surface, velocity potential, celerity and integral quantities admit ordered power series. In the wave-number scaling used by Fenton (1985), choose

$$
\boxed{\epsilon=\frac{kH}{2},\qquad D=kd.}
$$

The theory is constructed by substituting formal expansions into the exact kinematic and dynamic free-surface conditions and matching coefficients of powers of $\epsilon$ at each order. The nonlinear Fourier method, by contrast, solves numerically for the coefficients of the *truncated* basis without assuming they form a power series in $\epsilon$. This distinction explains why the Fourier approach generally remains useful beyond the comfortable validity range of a low-order Stokes approximation [5].

### 19.2 Surface expansion to fifth order

The fifth-order Stokes surface in Fenton's physically meaningful height parameter has the structure

$$
\begin{aligned}
\zeta(X)={}&\epsilon\cos X+\epsilon^2 B_{22}\cos(2X)\\
&+\epsilon^3B_{31}[\cos X-\cos(3X)]\\
&+\epsilon^4[B_{42}\cos(2X)+B_{44}\cos(4X)]\\
&+\epsilon^5[-(B_{53}+B_{55})\cos X
+B_{53}\cos(3X)+B_{55}\cos(5X)]\\
&+O(\epsilon^6).
\end{aligned}
$$

Here $B_{ij}=B_{ij}(D)$ are **Stokes perturbation coefficients**, not the $B_j$ solved directly by the Rust Fourier collocation model. The notational similarity should not obscure the different roles of the coefficients. By construction, the mean level is zero, and the odd-harmonic combinations give the prescribed crest-to-trough height. The series explicitly displays the emergence of second and higher spatial harmonics as amplitude increases.

Introduce Fenton's auxiliary depth function

$$
\boxed{S=\operatorname{sech}(2D).}
$$

Two representative coefficients are

$$
\boxed{B_{22}=\frac{\coth D\,(1+2S)}{2(1-S)},}
$$

$$
\boxed{B_{31}=-\frac{3(1+3S+3S^2+2S^3)}{8(1-S)^3}.}
$$

The remaining $B$ coefficients, together with the mean-flow, flux and Bernoulli coefficients, are collected in Appendix A. Their dependence on $(1-S)^{-p}$ explains why fifth-order Stokes expansions become increasingly difficult to use as $D\to0$, where $S\to1$.

### 19.3 Mean wave-frame speed

Fenton expresses the mean wave-frame velocity magnitude as

$$
\boxed{\bar U\sqrt{\frac kg}
=C_0+\epsilon^2C_2+\epsilon^4C_4+O(\epsilon^6),}
$$

where

$$
\boxed{C_0=\sqrt{\tanh D},}
$$

$$
\boxed{C_2=\sqrt{\tanh D}\,
\frac{2+7S^2}{4(1-S)^2},}
$$

$$
\boxed{C_4=\sqrt{\tanh D}\,
\frac{4+32S-116S^2-400S^3-71S^4+146S^5}
{32(1-S)^5}.}
$$

These coefficients have the correct linear limit $C_0=\sqrt{\tanh(kd)}$. Fenton's formula is expressed for a given $k$; if $T$ rather than $L$ is specified, one must solve a *nonlinear equation for $k$* involving the current convention. Assuming $c=\bar U$ silently sets the Eulerian mean current to zero.

### 19.4 Wavelength from a period and a prescribed Eulerian current

The Eulerian-current identity is

$$
\boxed{c=\bar u_1+\bar U.}
$$

With $c=2\pi/(kT)$, the fifth-order Stokes wavenumber relation becomes

$$
\boxed{
\frac{2\pi}{T\sqrt{gk}}
-\bar u_1\sqrt{\frac kg}
-C_0(D)-\epsilon^2C_2(D)-\epsilon^4C_4(D)=0,
}
$$

where $D=kd$ and $\epsilon=kH/2$ are functions of the unknown $k$. Even with an explicit series for the coefficient functions, determining wavelength from period is thus a **scalar nonlinear root problem**. Because current affects the measured period, omitting $\bar u_1$ is not a benign high-order approximation.

### 19.5 Stokes mass-flux and Bernoulli expressions

The fifth-order moving-frame volume-flux expansion is

$$
\boxed{Q\sqrt{\frac{k^3}{g}}
=C_0D+\epsilon^2(C_2D+D_2)
+\epsilon^4(C_4D+D_4)+O(\epsilon^6).}
$$

The auxiliary coefficients begin with

$$
\boxed{D_2=-\frac12\sqrt{\coth D},}
$$

$$
\boxed{D_4=\sqrt{\coth D}\,
\frac{2+4S+S^2+2S^3}{8(1-S)^3}.}
$$

Since $q=\bar Ud-Q$, it follows that

$$
\boxed{q\sqrt{\frac{k^3}{g}}=-\epsilon^2D_2-\epsilon^4D_4+O(\epsilon^6).}
$$

Similarly, the bed-origin Bernoulli constant is

$$
\boxed{\frac{Rk}{g}=D+\frac12C_0^2
+\epsilon^2E_2+\epsilon^4E_4+O(\epsilon^6),}
$$

with

$$
\boxed{E_2=\tanh D\,
\frac{2+2S+5S^2}{4(1-S)^2},}
$$

$$
\boxed{E_4=\tanh D\,
\frac{8+12S-152S^2-308S^3-42S^4+77S^5}{32(1-S)^5}.}
$$

The mean-level constant is $r=R-gd$ and consequently $rk/g=Rk/g-D$.

### 19.6 Breakdown of the perturbation expansion in shallow water

As $D\to0$,

$$
\tanh D\sim D,\qquad S=\operatorname{sech}(2D)\sim1-2D^2.
$$

The nominal steepness parameter $\epsilon=kH/2$ can therefore be small while higher Stokes coefficients become large. The effective shallow-water expansion parameter is proportional to

$$
\boxed{\frac{\epsilon}{D^3}=
\frac{kH/2}{(kd)^3}=\frac{Ur}{8\pi^2}.}
$$

This quantity must also be sufficiently small if successive Stokes terms are to decrease. The Fourier approximation avoids this particular *perturbation-series* limitation, although sufficiently long and high waves still generate broad harmonic spectra and difficult numerical conditioning. The selection between Stokes and Fourier methods should therefore depend on the desired accuracy and the hydrodynamic outputs of interest, not merely on the fact that both use Fourier harmonics [3, 5, 11].

***

## 20. Cnoidal and solitary-wave theory as the complementary long-wave limit

### 20.1 Why long waves require a different asymptotic organization

For finite-amplitude waves whose wavelength is large compared with depth, $D=kd\ll1$. Long-wave theory treats the horizontal variation as slowly varying rather than treating the wave height as the sole small parameter. The leading nonlinear/dispersive model is the Korteweg–de Vries (KdV) equation. Cnoidal theory is a periodic travelling-wave solution of that class of asymptotic equations, expressed using Jacobi elliptic functions. Fenton's higher-order cnoidal treatments systematically improve that representation [3, 6].

### 20.2 KdV equation in a standard shallow-water scaling

For a wave on quiescent water of mean depth $d$, let $c_0=\sqrt{gd}$. A conventional leading-order KdV equation for the surface displacement is

$$
\boxed{\eta_t+c_0\eta_x+
\frac{3c_0}{2d}\eta\eta_x+
\frac{c_0d^2}{6}\eta_{xxx}=0.}
$$

The first term is unsteady evolution, the second is linear long-wave propagation, the third is nonlinear steepening, and the fourth is weak dispersion. The balance of the third and fourth terms permits a coherent permanent wave instead of the shock-like steepening predicted by nondispersive shallow-water theory. This KdV equation is a **long-wave asymptotic model** and is not the exact potential-flow boundary-value problem solved by the Fourier method.

### 20.3 Derivation of the cnoidal waveform

Seek a travelling solution $\eta(x,t)=f(\xi)$ with $\xi=x-ct$. Substitution gives

$$
(c_0-c)f'+\frac{3c_0}{2d}ff'
+\frac{c_0d^2}{6}f'''=0.
$$

After one integration,

$$
\boxed{(c_0-c)f+\frac{3c_0}{4d}f^2
+\frac{c_0d^2}{6}f''=C_*,}
$$

where $C_*$ is a constant. Introduce the Jacobi elliptic function $\operatorname{cn}(s\mid m)$, with parameter $0<m<1$, and the trial waveform

$$
\boxed{f(\xi)=\eta_t+H\operatorname{cn}^2(\kappa\xi\mid m).}
$$

Here $\eta_t$ is the trough elevation relative to still water. The identity

$$
\boxed{\frac{d^2}{ds^2}\operatorname{cn}^2(s\mid m)
=2(1-m)+(8m-4)\operatorname{cn}^2(s\mid m)
-6m\operatorname{cn}^4(s\mid m)}
$$

allows the integrated equation to be separated into independent powers of $\operatorname{cn}^2$. Matching the coefficient of the fourth-power term yields

$$
\boxed{\kappa^2=\frac{3H}{4md^3}.}
$$

Matching the second-power term gives the corresponding leading-order celerity relation

$$
\boxed{c=c_0+\frac{3c_0\eta_t}{2d}
+\frac{c_0H(2m-1)}{2md}.}
$$

These relations are the leading KdV/cnoidal approximation, **not** the final celerity produced by the nonlinear Fourier solver.

### 20.4 Periodicity and complete elliptic integrals

The square of the Jacobi elliptic function $\operatorname{cn}(s\mid m)$ is periodic with real period $2K(m)$, where

$$
\boxed{K(m)=\int_0^{\pi/2}
\frac{d\varphi}{\sqrt{1-m\sin^2\varphi}}}
$$

is the complete elliptic integral of the first kind. Because $\kappa L=2K(m)$,

$$
\boxed{L=\frac{2K(m)}{\kappa}
=4K(m)\sqrt{\frac{md^3}{3H}}.}
$$

This expression connects the elliptic parameter $m$ to the wave geometry. It displays explicitly the long-wave scaling $L\propto d^{3/2}H^{-1/2}$, modulated by the elliptic integral and its parameter.

### 20.5 Determination of the trough level

Let the complete elliptic integral of the second kind be

$$
\boxed{E(m)=\int_0^{\pi/2}
\sqrt{1-m\sin^2\varphi}\,d\varphi.}
$$

The average of $\operatorname{cn}^2$ over one real period is

$$
\boxed{\left\langle\operatorname{cn}^2\right\rangle
=\frac{E(m)-(1-m)K(m)}{mK(m)}.}
$$

Imposing $\langle\eta\rangle=0$ gives

$$
\boxed{\eta_t=-H\,
\frac{E(m)-(1-m)K(m)}{mK(m)}.}
$$

Thus the crest elevation is $\eta_c=\eta_t+H$. This demonstrates analytically why long nonlinear waves have **sharper, taller crests and flatter, extended troughs** than a symmetric linear cosine wave. The initial periodic cnoidal waveform is not sinusoidal even though it is a smooth travelling wave.

### 20.6 Solitary-wave limit

As $m\to1^-$,

$$
K(m)\to\infty,
\qquad
\operatorname{cn}(s\mid m)\to\operatorname{sech}(s).
$$

Therefore $L\to\infty$, the periodic cnoidal train tends to an isolated elevation pulse, and the leading-order solitary profile becomes

$$
\boxed{\eta(\xi)\simeq H\operatorname{sech}^2
\left(\sqrt{\frac{3H}{4d^3}}\,\xi\right).}
$$

The corresponding leading-order solitary-wave speed is

$$
\boxed{c\simeq\sqrt{gd}\left(1+\frac{H}{2d}\right).}
$$

This is an asymptotic limit; it is not valid for arbitrary solitary-wave height. The need for a different representation near the solitary-wave limit explains why a Fourier approximation of a *finite-length periodic* wave becomes costly when very long flat trough regions coexist with a highly localized crest.

### 20.7 Higher-order cnoidal theories

Fenton developed higher-order cnoidal expressions using a perturbation parameter associated with *shallowness*, not merely with the relative height $H/d$. The solution can be organized in schematic form as

$$
\frac{\eta}{h}=1+\sum_{n=1}^{P}\alpha^{2n}Y_n(\theta;m),
\qquad
\frac{u_b}{Q/h}=1+\sum_{n=1}^{P}\alpha^{2n}U_n(\theta;m),
$$

where $h$ is the trough depth, $\theta=\alpha\xi/h$, and each $Y_n$, $U_n$ is a polynomial in $\operatorname{cn}^2(\theta\mid m)$ of increasing degree. Higher-order coefficient matching enforces the long-wave-ordered potential-flow equations. Fenton's later treatment also develops approximations for elliptic integrals and functions that converge rapidly as $m\to1$, precisely where straightforward classical series can become inefficient [6].

### 20.8 Domains of usefulness

Stokes and cnoidal theory are not fundamentally different categories of waves; they are different *expansion methods applied to the same steady gravity-wave family*. The Stokes series is organized about small steepness in relatively deeper water, while the cnoidal series is organized about small shallowness in longer waves. Their domains overlap. Fourier collocation supplies a high-accuracy numerical reference over much of both regions, provided its Fourier spectrum is resolved and the numerical problem remains well conditioned [3, 4, 6, 11].

***

## 21. The Rienecker–Fenton method within the wider family of numerical wave methods

### 21.1 Historical Fourier approximation and stream-function methods

The essential feature of a spectral approximation is to use a global basis capable of representing smooth periodic spatial variation to high accuracy. The basis is chosen so that the interior differential equations and some boundary conditions are exactly satisfied. Remaining unknown coefficients and surface ordinates are determined from the nonlinear free-boundary conditions. This family includes the earlier approaches of Chappelear and Dean and the explicit collocation formulation by Rienecker and Fenton [1, 4].

### 21.2 Why solving the nonlinear boundary conditions directly is powerful

For a Fourier basis of order $N$, increasing $N$ does not change the perturbation order: the model continues to impose the **full quadratic Bernoulli condition** and exact kinematic streamline condition on the finite representation. The principal approximation is spectral truncation and pointwise collocation. This is different from adding, say, the next Stokes term in powers of $H$ and means that an accurate numerical Fourier solution can be a reference against which analytical approximations are evaluated.

### 21.3 Original Fenton program versus the Rust numerical system

Fenton's 1988 FORTRAN program adopted numerically differentiated Jacobian entries to reduce programming complexity. The Rust solver whose mathematics is presented here instead implements an analytical dense Jacobian, a column-equilibrated Jacobi SVD, residual-decreasing line search, adaptive height continuation, conditioning-based order selection, independent integral checks and interval arithmetic. The **physical equations remain in the same Rienecker–Fenton tradition**, but these numerical algorithms are implementation decisions and should not be attributed wholesale to the 1981 or 1988 source programs. Likewise the Fenton 1988 program includes deep-water and wavelength-specified modes, while the present GUI provides finite-depth, period-specified inputs [1, 2, 11].

### 21.4 Why the 1982 Fourier time-stepping paper is a different problem

Fenton and Rienecker (1982) studied **unsteady** fully nonlinear waves, including solitary-wave interactions, using Fourier approximation for horizontal variation and temporal advancement of the boundary conditions. In that problem the surface and Fourier coefficients evolve in time. The steady travelling-wave problem considered here transforms to a stationary boundary-value problem and contains no physical time-stepping integration in its solver. Its *height continuation* is a numerical homotopy in wave amplitude, not the time evolution of a wave train. The 1982 paper supplies valuable background on spectral differentiation and the limitations introduced by high-frequency content, but it should not be cited as the exact time-integration algorithm of the Rust GUI [9].

### 21.5 Extension to variable bathymetry

A steady Fourier wave at fixed depth can be a useful **local** approximation along a gently varying seabed, and a sequence of such waves may form a shoaling calculation when conservation constraints are imposed. It is not a substitute for a spatial propagation model when reflections, strong slopes, rapid nonlinear evolution, group interactions or breaking are important. Later Fenton papers on boundary-integral and local polynomial approximation methods treat more general, time-dependent domains [7, 8]. Those approaches solve a different mathematical problem and have different accuracy and computational costs.

### 21.6 Why finite elements, BEM and other potential-flow models remain useful

The Fourier method assumes a horizontal bed and spatial periodicity. Boundary element methods can address irregular boundaries and the effects of structures; finite-element methods can also accommodate complex geometry; high-order spectral and local polynomial methods can evolve nonlinear waves over more complicated topography. In each case, the underlying irrotational model still relies on Laplace's equation with nonlinear surface conditions, but the discretization and boundary information differ. For marine engineering, the correct tool depends on whether one needs the properties of a single **incoming regular wave** or the actual wave transformation and forces near a physical structure.

***

## 22. Worked problems and quantitative interpretation

### 22.1 From dimensional inputs to the nondimensional state

Consider the period-specified regular wave

$$
H=0.50\ \mathrm m,\qquad T=8.0\ \mathrm s,\qquad
 d=12.0\ \mathrm m,\qquad U_c=0.
$$

Using $g=9.80665\ \mathrm{m/s^2}$, calculate the dimensionless input groups

$$
\boxed{\frac Hd=\frac{0.5}{12}=0.0416666667,}
$$

$$
\boxed{\tau_d=T\sqrt{\frac gd}
=8\sqrt{\frac{9.80665}{12}}\approx7.2320213864,}
$$

$$
\boxed{h_s=\frac{H}{gT^2}
=\frac{0.5}{9.80665\times64}\approx7.96653\times10^{-4}.}
$$

The Eulerian current Froude number is zero:

$$
\chi=\frac{U_c}{\sqrt{gd}}=0.
$$

No wavelength is assumed at this stage. It is an unknown determined by the full nonlinear equations.

### 22.2 Linear reference and nonlinear wavelength

For zero current, linear theory solves

$$
\left(\frac{2\pi}{8}\right)^2=9.80665\,k_0\tanh(12k_0),
$$

with the positive root approximately

$$
k_0\approx0.0828549562\ \mathrm{rad/m},
\qquad L_0=\frac{2\pi}{k_0}\approx75.8335481\ \mathrm m.
$$

The regression value stored in the Rust solver's numerical tests for this case is

$$
L_{\mathrm{Fourier,ref}}\approx75.88257\ \mathrm m.
$$

The difference is

$$
\boxed{L_{\mathrm{Fourier,ref}}-L_0\approx0.04902\ \mathrm m,}
$$

which is a small finite-amplitude correction in this relatively mild wave. These numerical values are an *illustrative benchmark from the supplied source and its reference constants*; no new independent Rust execution is asserted by their presentation here.

The corresponding nonlinear wavenumber and depth are approximately

$$
k\approx\frac{2\pi}{75.88257}=0.08280143\ \mathrm{rad/m},
\qquad D=kd\approx0.99361716.
$$

The relative wavelength, steepness and Ursell number are

$$
\frac Ld\approx6.32355,\qquad
\frac HL\approx0.00658913,\qquad
Ur=\frac{HL^2}{d^3}\approx1.66614.
$$

This example lies in a regime where a small-amplitude reference provides a close first estimate, while the numerical Fourier solution also supplies crest/trough asymmetry, velocities and integrals.

### 22.3 A stronger nonlinear wave with collinear current

Consider the input used in the GUI test configuration:

$$
H=3.0\ \mathrm m,\qquad T=9.0\ \mathrm s,
\qquad d=5.0\ \mathrm m,\qquad U_c=+1.0\ \mathrm{m/s}.
$$

The basic dimensionless quantities are

$$
\boxed{\frac Hd=0.6,\qquad
T\sqrt{\frac gd}\approx12.6042742750,
\qquad \frac{U_c}{\sqrt{gd}}\approx0.1428086981.}
$$

The reference value embedded in the supplied Rust regression tests is

$$
\boxed{L_{\mathrm{Fourier,ref}}\approx78.82722\ \mathrm m.}
$$

It implies

$$
\boxed{k\approx0.0797083183\ \mathrm{rad/m},
\quad kd\approx0.3985415912,
\quad c\approx8.7585802351\ \mathrm{m/s}.}
$$

The computed steepness and Ursell number from this reference wavelength are

$$
\boxed{\frac HL\approx0.03805792,
\qquad Ur\approx149.12954.}
$$

For comparison, the positive-frequency linear current-modified dispersion relation gives

$$
(\omega-k_0U_c)^2=gk_0\tanh(k_0d),
\qquad \omega=\frac{2\pi}{9},
$$

with $L_{0,\mathrm{current}}\approx70.03375$ m. The nonlinear and linear values are materially different, indicating that an Airy-only wavelength would omit important finite-amplitude effects at this height/depth ratio.

The supplied Rust regression also records a reference wave power of approximately $73.50878$ kW/m for the $U_c=1$ case. Its complete acceptance, including surface and integral checks, must still be established by executing the actual build on the intended machine.

### 22.4 The current convention is part of the wave specification

At $H=3$ m, $T=9$ s and $d=5$ m, the source documentation lists approximately

$$
L(U_c=0)=68.70676\ \mathrm m,
\qquad L(U_c=+1)=78.82722\ \mathrm m,
$$

for the prescribed **Eulerian** current. The difference is approximately $10.12$ m, emphasizing the role of current in the finite-amplitude wavelength. These numbers should not be reused for a Stokes-mass-transport current of the same numeric magnitude, because that would be a differently constrained physical problem.

### 22.5 A practical dimensional consistency audit

For any accepted wave, the following independent dimensional identities should be evaluated:

$$
\boxed{kL=2\pi,\qquad cT=L,\qquad kd=\frac{2\pi d}{L},}
$$

$$
\boxed{\bar u_1=c-\bar U,\qquad
\bar u_2=c-\frac Qd,\qquad q=\bar U d-Q,}
$$

$$
\boxed{E=T_K+V,\qquad I=\rho(cd-Q),}
$$

and the independent forms of $F$ and $S_{xx}$ given in §§14–15. An algebraically and numerically consistent report should satisfy these identities within tolerances appropriate to the dimensional scales, not only show five matching printed decimals.

### 22.6 Interpretation and limitations of a benchmark value

A numerical test case is valuable as a reproducibility target, not as a substitute for physical validation. A regression test that confirms $L$ to five decimals and power to a quoted tolerance does not establish accuracy for every other wave height or current, and it does not replace numerical surface/interpolation and integral checks. The high $H/d$ of the second example is specifically a reason to examine the theoretical limiting-wave ratio and physical breaking considerations rather than assume a regular steady wave must be realized in the field.

***

## 23. Wave–current interaction, dispersion branches and blocking

### 23.1 Intrinsic frequency and observed frequency

For linear waves on a uniform current $U_c$, the positive intrinsic-frequency branch satisfies

$$
\boxed{\sigma=\omega-kU_c>0,\qquad
\sigma^2=gk\tanh(kd).}
$$

The measured or absolute frequency is $\omega=2\pi/T$, whereas $\sigma$ is the frequency in a frame moving with the uniform current. The relation is a Doppler shift. Squaring the dispersion equation without retaining the branch sign can introduce spurious mathematical solutions that violate $\sigma>0$.

### 23.2 Group velocity and blocking

The linear absolute frequency is

$$
\omega(k)=kU_c+\sqrt{gk\tanh(kd)}.
$$

At constant depth and uniform current,

$$
\boxed{\frac{d\omega}{dk}=U_c+c_{g0}(k).}
$$

A blocking point on an adverse current is associated with vanishing laboratory-frame group velocity,

$$
\boxed{U_c+c_{g0}=0,}
$$

within the assumptions of linear theory. The corresponding dispersion equation can develop repeated roots or branches that terminate, so root-finding must track the physically appropriate branch. The finite-amplitude problem can have additional nonlinear changes in the branch structure and cannot be treated as a universal algebraic Doppler correction.

### 23.3 Why an adverse-current calculation can fail

An input with strong adverse current may be physically impossible as a single steady progressive-wave solution at the prescribed laboratory period; alternatively it may be mathematically possible but numerically hard to reach from the chosen initial branch. A non-converged Newton–SVD solve, by itself, does not distinguish these cases. A robust study may inspect a continuation sequence in current or wavelength, monitor the intrinsic frequency and compare with wave-blocking criteria and independent high-resolution solvers.

### 23.4 Galilean transformation and invariants

A uniform change of laboratory reference velocity by a constant $U_0$ transforms

$$
u'=u-U_0,\qquad c'=c-U_0,
$$

but leaves the relative wave-frame field $u-c$ unchanged. Therefore the shape of an ideal inviscid steady wave at fixed $(H,L,d)$ is Galilean-invariant, while apparent period, Eulerian current, laboratory kinetic energy and energy flux are frame-dependent. This provides a useful conceptual check on current conventions: the physical relative flow cannot depend on the arbitrary observer speed, yet the reported speed $c=L/T$ certainly can.

### 23.5 Practical interpretation of currents

A single scalar collinear current in an irrotational potential-flow model is a *uniform mean-current convention*, not a resolved vertical shear profile. Actual currents measured at one instrument elevation may exhibit vertical variation and directional spreading. Using such a measured speed as $U_c$ assumes the idealized current convention is representative for the wave problem. In strong shear, vorticity is nonzero and the Fourier irrotational formulation does not describe the complete wave–current interaction.

***

## 24. Shoaling, propagation over variable depth and the limits of local steady solutions

### 24.1 Slowly varying bathymetry as a sequence of local steady waves

A frequently used engineering approximation for waves crossing gently varying bathymetry assumes that each location can be represented by a locally periodic nonlinear wave of mean depth $d(x)$ and wave parameters that change slowly in the horizontal direction. The *local* steady wave problem then supplies $L(x)$, the velocity profile and integral properties. The spatial evolution must be governed by additional conservation or balance equations.

### 24.2 Conservation constraints

For an idealized stationary medium without wave dissipation or reflection, a non-current case is often modelled with conserved period and energy flux:

$$
\boxed{T(x)=T_0,\qquad F(x)\approx F_0.}
$$

The mass-flux condition or return-current convention must also be specified consistently. With waves on currents, energy flux in the fixed frame is not always the appropriate conserved wave quantity by itself; wave-action conservation and mean-flow coupling are generally required. These relations are **propagation modelling assumptions beyond the stand-alone Fourier solver**, not extra constraints automatically imposed by `fenton_gui.rs`.

### 24.3 Linear shoaling coefficient

In the small-amplitude, no-current limit with conserved energy flux, $E=\rho gH^2/8$ and $F=c_{g0}E$. If $H_0$ and $c_{g0,0}$ denote offshore values, then

$$
\frac18\rho gH(x)^2c_g(x)
=\frac18\rho gH_0^2c_{g0,0}.
$$

Hence the conventional linear shoaling coefficient is

$$
\boxed{K_s=\frac{H(x)}{H_0}
=\sqrt{\frac{c_{g0,0}}{c_g(x)}}.}
$$

This relationship follows from linear energy conservation. It is not exact for a finite-amplitude wave with current, partial reflection, dissipation, strong bathymetric slopes or wave breaking.

### 24.4 Nonlinear shoaling using the Fourier solver

A nonlinear local-steady shoaling procedure can, in principle, determine the unknown $H(x)$ by solving simultaneously for the local Fourier wave and an outer energy-flux condition

$$
\boxed{\mathcal G(H,d,T,\text{current constraint})
=F(H,T,d,U_c)-F_0=0.}
$$

The outer iteration must handle a current or mass-transport closure at every depth and may require momentum/mean-flow balances as well. Near the highest wave or a wave-blocking point, the mapping may cease to be monotone or a smooth steady solution may be physically inadmissible. Numerical continuation in the spatial coordinate then becomes as important as continuation in wave height.

### 24.5 When a propagation model is necessary

Local steady theory is inadequate for strong reflection, diffraction, multimodal wave interactions, horizontal current shear, bar-trough bathymetry that changes within a wavelength, breaking, run-up at structures and rapidly evolving nonlinear groups. The 1982 unsteady Fourier paper and the 1993/1996 propagation papers address classes of such problems through temporal spectral, boundary-integral or local polynomial methods. Their use of the same potential-flow foundations should not be mistaken for the same computational algorithm [7–9].

***

## 25. Engineering use: wave kinematics and forces on slender members

### 25.1 Pointwise kinematics from a nonlinear incident-wave field

For a vertical circular pile or another sufficiently slender member, the steady nonlinear Fourier solution supplies the incident horizontal velocity $u(x,y,t)$ and acceleration $a_x(x,y,t)$. In a current-bearing case, decide explicitly whether the drag force uses total relative water–member velocity or an orbital-only component. The appropriate acceleration term must be chosen for the selected Morison-type formulation and body-motion convention.

### 25.2 Classical Morison-type decomposition

For a fixed, slender, vertical circular cylinder of diameter $D_p$, a conventional horizontal force per unit member length is

$$
\boxed{f_x(y,t)=\frac12\rho C_DD_p\,u|u|
+\rho C_M\frac{\pi D_p^2}{4}\,a_x.}
$$

Here $C_D$ is a drag coefficient and $C_M$ is an inertia coefficient. The formula's validity depends on assumptions about member slenderness, diffraction, flow separation, roughness and the interpretation of acceleration. If a chosen engineering standard defines the inertia term using a different relative or local acceleration, that definition takes precedence. The nonlinear Fourier program itself **does not calculate these Morison coefficients**; it only supplies incident-flow kinematics. Westberg (1984, pp. 3–9) discusses accelerations and an extension of Fenton-style kinematics to horizontal forces and moments on vertical cylinders [10].

### 25.3 Depth integration of force and base moment

Let $z=y+d$ denote height above the seabed. The integrated horizontal loading at time $t$ is

$$
\boxed{F_x(t)=\int_0^{d+\eta(t)}f_x(z,t)\,dz.}
$$

The corresponding overturning moment about the seabed is

$$
\boxed{M_b(t)=\int_0^{d+\eta(t)}z\,f_x(z,t)\,dz.}
$$

For a member that is not submerged throughout the entire wave period, the wet integration limits must be truncated to the actual submerged geometry. Peak integrated force and peak integrated moment generally occur at different phases; integrating the local peak envelope over depth would not reproduce either true simultaneous time-dependent maximum.

### 25.4 Current, acceleration convention and marine growth

The drag term is sensitive to total relative velocity, including mean current and structural motion. A member's hydrodynamic diameter and surface roughness may differ from its nominal steel dimensions because of marine growth. The coefficients $C_D$ and $C_M$ vary with Reynolds and Keulegan–Carpenter numbers and may be affected by local geometry and wave steepness. These uncertainties are outside the ideal potential-flow solver, so high mathematical accuracy in wave kinematics should not be equated with equally high accuracy in empirical structural force coefficients.

### 25.5 Pressure-driven loads versus Morison loads

For large or closely spaced structures, the incident field is altered by diffraction and radiation, and potential-flow pressures must be solved around the body rather than simply inserted into an isolated-pile Morison expression. The steady-wave pressure from §11 is the pressure in the **undisturbed incoming wave field**. It is not a wave-impact pressure distribution on a wall and does not account for impulsive breaking loads, air entrapment or nonpotential turbulence.

***

## 26. Numerical software architecture and reproducibility of the mathematical model

### 26.1 Mapping of mathematical objects to solver components

The single-file Rust implementation separates several mathematical responsibilities into internal modules, without requiring third-party Rust crates:

| Mathematical responsibility | Rust component | Principal function or object |
|:---|:---|:---|
| Physical inputs, nondimensionalization and active unknowns | `solver` | `FentonStreamFunction::new` |
| Periodic harmonic evaluation | `solver` | `init_trig_tables`, `surface_keta`, `point` |
| Complete nonlinear equations | `solver` | `eqns_for_state_into` |
| Analytical derivatives | `solver` | `jacobian` |
| Newton corrections and line search | `solver` | `newton_step` |
| Adaptive height progression | `solver` | `solve_internal` |
| Stable Jacobian inversion | `linalg` | `svd_solve` |
| Dimensional integral invariants | `solver` | `calc_integral_props` |
| Pressure/velocity-based independent integrals | `solver` | `defining_integrals` |
| Residual and output enclosures | `verified` / `solver` | `Interval`, `verify_finite_system` |
| Output residual and accuracy gate | `solver` | `verify_accuracy` |
| Engineering report and Fenton-style tables | `report` | `generate_output`, `print_solution_flat` |
| Windows graphical interface | `windows_gui` | Win32 event loop and worker thread |

The mathematical core works with arrays of $f64$ values and does not depend on the Win32 graphical system. The GUI is a front end; it cannot make a rejected numerical state acceptable merely through presentation.

### 26.2 The data path from four inputs to engineering results

The computational dependencies can be represented schematically as

$$
(H,T,d,U_c)\longrightarrow(h_r,h_s,\chi)
\longrightarrow\mathbf z^{(0)}
\longrightarrow\mathbf z^{(*)}
\longrightarrow\{k,L,c,\eta,u,v,p\}
\longrightarrow\{T_K,V,I,S_{xx},F,\ldots\}.
$$

At the same time, an independent verification path is evaluated:

$$
\mathbf z^{(*)}\longrightarrow
\{\mathbf F,\ R_{\mathrm{surf}},\ \mathcal I_A,\ \mathcal F_i\}
\longrightarrow\text{accept/reject}.
$$

The report should be generated only from accepted states. The no-current and imposed-current problems are independent nonlinear solves; the current-bearing case must not be obtained by algebraically adjusting a zero-current report.

### 26.3 Zero-current and nonzero-current problems can be solved concurrently

If a nonzero current is specified, the program can run the no-current comparison and the current-bearing solution on distinct threads, since each case has its own state, Jacobian, Fourier coefficients and verification record. This is *task-level parallelism*, not a change in the mathematical equations, the floating-point precision, or the sequence of operations within a single Newton/SVD solve. Any speed increase depends on available CPU resources and workload balance.

### 26.4 Numerical tests must be separated from formatted output comparisons

Formatting tests verify table widths, Unicode scientific symbols, labels and report wrapping. Numerical tests instead verify physical identities, residual limits, reference wavelengths, power and verification flags. Passing a report-width test does not prove that a numerical solution has converged, and passing a numerical benchmark does not prove that the report retains all digits required for further calculations. Values should be retained as `f64` until formatting; fixed decimal places are a presentation decision.

### 26.5 A recommended regression set

A useful production test matrix should cover independent corners of the model's applicability:

| Test category | Purpose | Expected property |
|:---|:---|:---|
| Mild zero-current wave | Linear consistency | $L$ close to the Airy value; low higher harmonics |
| Moderate nonlinear wave | Finite-amplitude physics | Nonlinear correction and asymmetric free surface |
| Following Eulerian current | Doppler/current closure | $\bar u_1=U_c$ and physical branch selection |
| Adverse current | Root/branch and conditioning | Clear acceptance or rejection, no false positive |
| Shallow long wave | Spectral resolution | Several modes and off-grid surface checks may be needed |
| Short high wave | Conditioning and limiting behavior | No acceptance of under-resolved solution |
| Near-limiting wave | Spectral-tail and crest diagnostics | Independent physical and numerical warnings |
| Invalid inputs | Domain control | Rejected before Newton iteration |
| Independent integral identities | Post-processing | Errors within declared dimensional budgets |
| Interval input/operation domain | Verification integrity | Invalid interval operation causes rejection |

The test matrix is not a claim about a particular code's universal validity. It specifies meaningful mathematical properties to assess in each environment.

***

## 27. Numerical failure modes and defensible scientific interpretation

### 27.1 Failure of local convergence

A Newton iteration may fail because the initial guess is outside the basin of attraction, because the nonlinear solution branch folds, because the Jacobian is nearly singular, or because the input combination does not admit a smooth steady solution. Adaptive height continuation mitigates initial-guess difficulty but cannot eliminate these fundamental limitations.

### 27.2 Small residual with large state error

Near a singular Jacobian, the linear sensitivity estimate is

$$
\|\Delta\mathbf z\|\lesssim
\|\mathbf J^{-1}\|\,\|\mathbf F\|.
$$

Even a residual as small as $10^{-10}$ can correspond to an appreciable state uncertainty when $\|\mathbf J^{-1}\|$ is large. This motivates column scaling, singular-value diagnostics and conservative rejection thresholds. A condition number computed from a rounded finite Jacobian is an **estimate**, not a verified condition-number bound.

### 27.3 False spectral convergence

A discrete Fourier series can exactly interpolate the values at $N+1$ nodes while representing the continuous curve poorly between them, especially for a sharp crest with unresolved harmonics. A residual test only at the collocation points can therefore appear excellent while the reconstructed surface exhibits oscillations or excessive Bernoulli defects off grid. High-density surface scanning, coefficient-decay inspection and independent quadrature are complementary checks.

### 27.4 Incorrect wave branch

An algebraic nonlinear system may admit multiple solutions, including states whose dominant harmonic structure does not correspond to the intended fundamental wave. Fenton noted that for some very long waves a method could converge to a wave with one third of the intended wavelength. Tracking a continuous branch from a small-amplitude seed and inspecting Fourier harmonics reduces that risk [4, 11].

### 27.5 Floating-point reproducibility and scientific accuracy

Two CPU/compiler combinations may return slightly different last digits because of order of operations, library transcendentals, vectorization and arithmetic contraction. These differences are not necessarily a loss of *engineering* accuracy. However, when the acceptance criterion depends on a bound as small as $10^{-10}$, small numerical changes may legitimately alter a pass/fail decision near a threshold. That sensitivity should be treated transparently rather than concealed by looser convergence tolerances.

### 27.6 Model inadequacy despite numerical success

A wave may be a highly accurate solution of the smooth irrotational steady equations yet be inappropriate for a physical application because breaking, wind forcing, turbulent dissipation, wave grouping, strongly varying current or irregular bathymetry dominate. The mathematics certifies neither the input environment nor the full-scale engineering consequence. Such limitations should accompany any use of the reported velocities, accelerations or loads in critical design.

### 27.7 Differences between a numerical certificate and a mathematical theorem

The most useful logical distinction is

$$
\boxed{
\text{verified arithmetic at a finite state}
\not\Longrightarrow
\text{existence of a nearby exact finite-system root}
\not\Longrightarrow
\text{rigorous continuum-wave error bound}.
}
$$

A full exact-root theorem would require bounds on the inverse Jacobian and nonlinear remainder. A continuum theorem would additionally need a controlled spectral truncation estimate. A physical design validation would require confirmation that the mathematical model is appropriate to measured boundary conditions and the actual load mechanism. The implemented checks are significant evidence of **numerical internal consistency** within their stated scope, without conflating these logically distinct levels.

***

## 28. Conclusions

The mathematics of a finite-amplitude steady surface gravity wave is a nonlinear free-boundary problem in an otherwise harmonic, incompressible and irrotational flow. The wave field, free surface, phase speed, mean current and integral properties form a coupled system; wavelength cannot generally be inferred from height, depth and measured period without a consistent current convention.

The Rienecker–Fenton Fourier approximation makes efficient use of the analytical structure of the problem. Every basis function satisfies Laplace's equation and the impermeable bed condition; Newton collocation determines the nonlinear free surface and harmonic coefficients simultaneously. Analytical Jacobians, scaled SVD corrections, safeguarded height continuation and conditioning-aware order selection provide practical numerical robustness. Fourier interpolation and integral verification are necessary because a small residual at the chosen nodes does not establish adequate spectral resolution over the entire wave.

The post-processing identities for velocity, material acceleration, Bernoulli pressure, volume flux, kinetic and potential energy, wave impulse, radiation stress and energy flux follow from the same accepted nonlinear state and its reference-frame conventions. Independent integral evaluations and conservative finite-precision enclosures offer important internal checks. They support numerical confidence but do not constitute proofs of an exact root, spectral convergence for the continuous equations, or the physical stability of a particular high wave.

Stokes and cnoidal theories arise as complementary perturbation organizations of the same physical wave family, most naturally suited to relatively shorter and longer waves, respectively. Their approximate formulas retain considerable explanatory value, whereas numerical Fourier collocation is well suited as a reference when nonlinear wave shape, detailed kinematics and integral quantities are required over a broad range of physically admissible conditions.

For engineering practice the essential requirement is to keep all of the following **simultaneously** consistent: the prescribed regular-wave height and period, the definition of current, the chosen frame, the dimensional scaling, the Fourier resolution, the numerical acceptance budgets and the physical limits of the inviscid steady-wave model.

***

## Appendix A. Fifth-order Stokes coefficients in Fenton's wave-height formulation

### A.1 Purpose, conventions and relationship to the numerical solver

For completeness, this appendix assembles the coefficient functions that appear in Fenton's fifth-order Stokes theory [5, Table 1; 3, Table 1]. They are provided for analytical comparison, asymptotic interpretation and independent checking of mild-wave numerical solutions. **They are not the coefficient-generation algorithm used by the fully nonlinear Rust Fourier solver.**

Define

$$
D=kd,\qquad \epsilon=\frac{kH}{2},\qquad
S=\operatorname{sech}(2D),\qquad
C_0=\sqrt{\tanh D}.
$$

All coefficients below depend on $D$, either explicitly or through $S$. It is important that $S$ is **$\operatorname{sech}(2D)$**, not $\operatorname{sech}^2(D)$. Confusing these two functions produces seriously incorrect coefficients even at intermediate depths.

### A.2 Stream-function and potential coefficients $A_{ij}$

The dimensionless coefficients of the harmonics in the fifth-order potential expansion are

$$
\boxed{A_{11}=\frac1{\sinh D},}
$$

$$
\boxed{A_{22}=\frac{3S^2}{2(1-S)^2},}
$$

$$
\boxed{A_{31}=
\frac{-4-20S+10S^2-13S^3}
{8\sinh D\,(1-S)^3},}
$$

$$
\boxed{A_{33}=
\frac{-2S^2+11S^3}
{8\sinh D\,(1-S)^3},}
$$

$$
\boxed{A_{42}=
\frac{12S-14S^2-264S^3-45S^4-13S^5}
{24(1-S)^5},}
$$

$$
\boxed{A_{44}=
\frac{10S^3-174S^4+291S^5+278S^6}
{48(3+2S)(1-S)^5},}
$$

$$
\boxed{A_{51}=
\frac{-1184+32S+13232S^2+21712S^3+20940S^4+12554S^5
-500S^6-3341S^7-670S^8}
{64\sinh D\,(3+2S)(4+S)(1-S)^6},}
$$

$$
\boxed{A_{53}=
\frac{4S+105S^2+198S^3-1376S^4-1302S^5-117S^6+58S^7}
{32\sinh D\,(3+2S)(1-S)^6},}
$$

$$
\boxed{A_{55}=
\frac{-6S^3+272S^4-1552S^5+852S^6+2029S^7+430S^8}
{64\sinh D\,(3+2S)(4+S)(1-S)^6}.}
$$

An equivalent potential representation for the laboratory-frame oscillatory velocity field, using the bed-origin vertical coordinate $z=y+d$, is

$$
\boxed{
\phi(x,z,t)=\bar u_1x+
C_0\sqrt{\frac{g}{k^3}}
\sum_{i=1}^{5}\epsilon^i
\sum_{j\in\mathcal J_i}A_{ij}\cosh(jkz)\sin(jX)
+\phi_0(t)+O(\epsilon^6),
}
$$

where $X=k(x-ct)$, $\mathcal J_i$ contains the harmonics $j$ of parity compatible with order $i$, and $\phi_0(t)$ is an arbitrary gauge function of time. The formula is schematic as to the temporal gauge but explicit about the spatial potential harmonics. Differentiating with respect to $x$ and $z$ produces the corresponding Stokes velocity components. No factor $\cosh(jD)$ is inserted in addition to the supplied $A_{ij}$, because their coefficients already include the required depth dependence.

### A.3 Free-surface coefficients $B_{ij}$

The harmonic coefficients of the surface-height series in §19.2 are

$$
\boxed{B_{11}=1,}
$$

$$
\boxed{B_{22}=\frac{\coth D\,(1+2S)}{2(1-S)},}
$$

$$
\boxed{B_{31}=
-\frac{3(1+3S+3S^2+2S^3)}{8(1-S)^3},
\qquad B_{33}=-B_{31},}
$$

$$
\boxed{B_{42}=
\coth D\,
\frac{6-26S-182S^2-204S^3-25S^4+26S^5}
{6(3+2S)(1-S)^4},}
$$

$$
\boxed{B_{44}=
\coth D\,
\frac{24+92S+122S^2+66S^3+67S^4+34S^5}
{24(3+2S)(1-S)^4},}
$$

$$
\boxed{B_{53}=
\frac{9(132+17S-2216S^2-5897S^3-6292S^4-2687S^5
+194S^6+467S^7+82S^8)}
{128(3+2S)(4+S)(1-S)^6},}
$$

$$
\boxed{B_{55}=
\frac{5(300+1579S+3176S^2+2949S^3+1188S^4+675S^5
+1326S^6+827S^7+130S^8)}
{384(3+2S)(4+S)(1-S)^6},}
$$

$$
\boxed{B_{51}=-(B_{53}+B_{55}).}
$$

The identities $B_{33}=-B_{31}$ and $B_{51}=-(B_{53}+B_{55})$ help maintain the prescribed crest-to-trough height in this particular wave-height-based perturbation series. They are excellent algebraic regression checks for any independent fifth-order Stokes implementation.

### A.4 Mean-speed, discharge and Bernoulli coefficients

The mean-speed coefficient set is

$$
\boxed{C_0=\sqrt{\tanh D},}
$$

$$
\boxed{C_2=\sqrt{\tanh D}\,
\frac{2+7S^2}{4(1-S)^2},}
$$

$$
\boxed{C_4=\sqrt{\tanh D}\,
\frac{4+32S-116S^2-400S^3-71S^4+146S^5}
{32(1-S)^5}.}
$$

The discharge-offset functions are

$$
\boxed{D_2=-\frac12\sqrt{\coth D},}
$$

$$
\boxed{D_4=\sqrt{\coth D}\,
\frac{2+4S+S^2+2S^3}{8(1-S)^3}.}
$$

The Bernoulli functions are

$$
\boxed{E_2=\tanh D\,
\frac{2+2S+5S^2}{4(1-S)^2},}
$$

$$
\boxed{E_4=\tanh D\,
\frac{8+12S-152S^2-308S^3-42S^4+77S^5}
{32(1-S)^5}.}
$$

In the strict $\epsilon\to0$ limit, the higher-order terms vanish and only $C_0$ contributes to the linear dispersion. As $D\to0$, many denominators shrink and the Stokes asymptotic expansion loses uniform validity. These formulas are appropriately used as *asymptotic analytical checks* and should not be substituted indiscriminately for the numerical Fourier model under very high or long waves.

### A.5 Deep-water coefficients and the classic small-amplitude limit

When $D\to\infty$, $\tanh D\to1$, $S\to0$, and

$$
\boxed{C_0\to1,\qquad C_2\to\frac12,\qquad C_4\to\frac18.}
$$

Consequently

$$
\bar U\sqrt{\frac kg}
=1+\frac12\epsilon^2+\frac18\epsilon^4+O(\epsilon^6)
$$

in the corresponding deep-water fifth-order convention. These are coefficients in a chosen amplitude and mean-current convention; comparisons with other published Stokes expansions require the same definition of the expansion parameter and wave speed.

***

## Appendix B. Algebraic consistency identities for code verification

### B.1 Identity checks that do not require a full nonlinear solve

The following exact mathematical statements hold independently of the values of the Fourier coefficients, provided $D>0$:

$$
S_j(-D;D)=0,
\qquad C_j(-D;D)=\operatorname{sech}(jD),
$$

$$
\frac{\partial S_j}{\partial Y}=jC_j,
\qquad
\frac{\partial C_j}{\partial Y}=jS_j,
$$

$$
\left(\partial_{XX}+\partial_{YY}\right)
[S_j(Y;D)\cos(jX)]=0.
$$

They are fundamental unit tests for the harmonic kernel. If these fail, numerical refinement or better Newton tolerances cannot repair the physics.

### B.2 Self-consistency of the global state

For any accepted state,

$$
\boxed{z_1=kd,\quad z_2=kH,\quad z_3=T\sqrt{gk},
\quad z_4=c\sqrt{k/g}.}
$$

Therefore

$$
\boxed{\frac{z_2}{z_1}=\frac Hd,
\quad z_3z_4=2\pi,
\quad\frac{2\pi}{z_1}=\frac Ld.}
$$

The mean-current identity requires

$$
\boxed{z_4=z_5+z_7.}
$$

The wave-flux identity requires

$$
\boxed{z_8=z_1(z_6+z_7-z_4).}
$$

These statements provide quick diagnostics before calculating velocities or integrals.

### B.3 Free-surface conservation checks

The crest-trough height and mean-level conditions can be checked without a pointwise flow evaluation:

$$
\boxed{\zeta_0-\zeta_N=kH,}
$$

$$
\boxed{\zeta_0+\zeta_N
+2\sum_{m=1}^{N-1}\zeta_m=0.}
$$

The discrete cosine transform should then give a zero constant surface coefficient, up to arithmetic rounding.

### B.4 Fourier orthogonality and Parseval tests

For positive integers $j,l$,

$$
\frac1\pi\int_0^{2\pi}
\cos(jX)\cos(lX)\,dX=\delta_{jl},
$$

and similarly for sine-sine products. The cosine-sine cross integrals vanish. The reconstructed potential-energy identity and the bed-velocity Parseval identity in §§13–14 should agree with sufficiently fine direct integration. These tests isolate Fourier normalization errors, especially half weights on the last surface harmonic.

### B.5 Differential consistency of the velocity field

Evaluate the analytically reconstructed derivatives at many interior points:

$$
\boxed{u_x+v_y=0,\qquad u_y-v_x=0.}
$$

Since the Fourier representation satisfies both equations identically, their numerical errors should be driven chiefly by floating-point evaluation and not by an $O(H)$ model approximation.

### B.6 Bernoulli pressure check

At an accepted collocation point, the pressure should vanish:

$$
\boxed{p(X_m,\eta_m)\approx0.}
$$

For off-grid phases the residual pressure is a direct measure of the interpolation's Bernoulli defect; see §11.2. A pressure formula that uses $u^2$ rather than $(u-c)^2$ in the wave-frame Bernoulli equation will fail this test whenever $c\ne0$.

### B.7 Discharge invariance check

Choose arbitrary phases $X_a,X_b$ and compute

$$
Q(X)=-\int_{-d}^{\eta(X)}(u-c)\,dy.
$$

Then, for a sufficiently resolved numerical quadrature,

$$
\boxed{Q(X_a)-Q(X_b)\approx0.}
$$

This tests continuity and free-surface kinematics independently of the algebraic state expressions for $Q$.

### B.8 Dimensional analysis as an error detector

Each expression should be dimensionally consistent before its numerical value is compared. For example,

$$
[k]=\mathrm{m^{-1}},\quad [\psi]=\mathrm{m^2/s},
\quad [Q]=[q]=\mathrm{m^2/s},
$$

$$
[I]=\mathrm{kg/(m\,s)},\quad [E]=[T_K]=[V]=\mathrm{J/m^2},
$$

$$
[S_{xx}]=[S]=\mathrm{N/m},\quad
[F]=\mathrm{W/m},\quad [C_E]=\mathrm{m/s}.
$$

Dimensional checks often reveal a missing $k$, $d$, $g$, or $\rho$ factor more quickly than a high-precision comparison in one numerical case.

***

## Appendix C. Comprehensive notation and nondimensional scales

### C.1 Physical variables and constants

| Symbol | Definition | SI unit |
|:---|:---|:---|
| $x$ | Horizontal coordinate fixed to seabed | m |
| $y$ | Vertical coordinate, zero at still water | m |
| $z=y+d$ | Height above the seabed | m |
| $\xi=x-ct$ | Dimensional travelling coordinate | m |
| $t$ | Time | s |
| $H$ | Crest-to-trough height of a regular wave | m |
| $T$ | Apparent regular-wave period | s |
| $L$ | Wavelength | m |
| $d$ | Mean still-water depth | m |
| $g$ | Gravitational acceleration | m/s² |
| $\rho$ | Fluid density | kg/m³ |
| $k=2\pi/L$ | Wavenumber | rad/m |
| $\omega=2\pi/T$ | Angular frequency | rad/s |
| $c=L/T$ | Apparent phase speed | m/s |
| $\eta$ | Free-surface elevation relative to still water | m |
| $\eta_c,\eta_t$ | Crest and trough elevations | m |
| $u,v$ | Fixed-frame fluid velocities | m/s |
| $U=u-c,V=v$ | Wave-frame fluid velocities | m/s |
| $\phi$ | Velocity potential | m²/s |
| $\psi$ | Wave-frame stream function | m²/s |
| $p$ | Gauge pressure | Pa |
| $r$ | Bernoulli constant, mean-level datum | m²/s² |
| $R=r+gd$ | Bernoulli constant, bed datum | m²/s² |
| $\bar U$ | Positive magnitude of mean opposing wave-frame flow | m/s |
| $\bar u_1$ | Eulerian mean current | m/s |
| $\bar u_2$ | Depth-integrated mass-transport current | m/s |
| $Q$ | Positive magnitude of wave-frame volume discharge | m²/s |
| $q=\bar U d-Q$ | Wave-flux correction | m²/s |
| $\theta=k(x-ct)$ | Dimensionless phase | rad |
| $N$ | Fourier order | – |

### C.2 Nondimensional variables

| Symbol | Nondimensional definition | Interpretation |
|:---|:---|:---|
| $X$ | $k(x-ct)$ | Wave phase |
| $Y$ | $ky$ | Vertical phase coordinate |
| $D$ | $kd$ | Relative depth |
| $\zeta$ | $k\eta$ | Surface elevation |
| $\Psi$ | $\psi\sqrt{k^3/g}$ | Wave-frame stream function |
| $\widehat u$ | $u\sqrt{k/g}$ | Laboratory velocity |
| $\widehat v$ | $v\sqrt{k/g}$ | Vertical velocity |
| $z_1$ | $kd$ | Active depth |
| $z_2$ | $kH$ | Active height |
| $z_3$ | $T\sqrt{gk}$ | Period variable |
| $z_4$ | $c\sqrt{k/g}$ | Phase-speed variable |
| $z_5$ | $\bar u_1\sqrt{k/g}$ | Eulerian-current variable |
| $z_6$ | $\bar u_2\sqrt{k/g}$ | Mass-transport-current variable |
| $z_7$ | $\bar U\sqrt{k/g}$ | Mean wave-frame velocity magnitude |
| $z_8$ | $q\sqrt{k^3/g}$ | Flux correction |
| $z_9$ | $rk/g$ | Bernoulli offset |
| $z_{10+m}$ | $k\eta(X_m)$ | Nodal surface ordinate |
| $z_{N+10+j}$ | $B_j$ | Harmonic stream-function coefficient |
| $h_r$ | $H/d$ | Height-depth ratio |
| $h_s$ | $H/(gT^2)$ | Period-form height variable |
| $\chi$ | $U_c/\sqrt{gd}$ | Prescribed current Froude number |
| $\mu$ | $d/L$ | Relative depth |
| $s$ | $H/L$ | Wave steepness |
| $Ur$ | $HL^2/d^3$ | Ursell number |
| $SU$ | $Ur/(8\pi^2)$ | Stokes–Ursell parameter |

### C.3 Flux and energy scales

Dimensionalizing an equation from wave-number-based units uses the scales

$$
\boxed{[\mathrm{length}]=\frac1k,
\quad [\mathrm{velocity}]=\sqrt{\frac gk},
\quad [\mathrm{time}]=\frac1{\sqrt{gk}}.}
$$

The derived scales are

$$
\boxed{[\mathrm{flux}]=\sqrt{\frac g{k^3}},
\qquad [\mathrm{Bernoulli}]=\frac gk,}
$$

$$
\boxed{[\mathrm{energy/area}]=[\mathrm{stress}]
=\frac{\rho g}{k^2},}
$$

$$
\boxed{[\mathrm{impulse}]=\rho\sqrt{\frac g{k^3}},
\qquad [\mathrm{power/width}]=\frac{\rho g^{3/2}}{k^{5/2}}.}
$$

An equivalent water-depth-based scale set is

$$
\boxed{[\mathrm{velocity}]_d=\sqrt{gd},
\quad [\mathrm{flux}]_d=\sqrt{gd^3},
\quad [\mathrm{energy/area}]_d=\rho gd^2,}
$$

$$
\boxed{[\mathrm{impulse}]_d=\rho\sqrt{gd^3},
\qquad [\mathrm{power/width}]_d=\rho g^{3/2}d^{5/2}.}
$$

In particular, a dimensionless energy calculated in wave-number units $\widehat E_k$ converts to its depth-based value by

$$
\widehat E_d=\frac{\widehat E_k}{(kd)^2},
$$

and a dimensionless wave power $\widehat F_k$ converts by

$$
\widehat F_d=\frac{\widehat F_k}{(kd)^{5/2}}.
$$

### C.4 The nineteen quantities in a Fenton-style `SOLUTION.RES` table

The Rust report displays dimensional values alongside corresponding depth-scaled nondimensional values. The mathematical form of each row is:

| No. | Quantity | Nondimensional output |
|---:|:---|:---|
| 1 | Water depth $d$ | $d/d=1$ |
| 2 | Wavelength $L$ | $L/d$ |
| 3 | Wave height $H$ | $H/d$ |
| 4 | Wave period $T$ | $T\sqrt{g/d}$ |
| 5 | Wave speed $c$ | $c/\sqrt{gd}$ |
| 6 | Eulerian current $\bar u_1$ | $\bar u_1/\sqrt{gd}$ |
| 7 | Stokes current $\bar u_2$ | $\bar u_2/\sqrt{gd}$ |
| 8 | Mean fluid speed $\bar U$ | $\bar U/\sqrt{gd}$ |
| 9 | Wave volume flux $q$ | $q/\sqrt{gd^3}$ |
| 10 | Reduced Bernoulli constant $r$ | $r/(gd)$ |
| 11 | Volume flux $Q$ | $Q/\sqrt{gd^3}$ |
| 12 | Bernoulli constant $R$ | $R/(gd)$ |
| 13 | Momentum flux $S$ | $S/(\rho gd^2)$ |
| 14 | Impulse $I$ | $I/(\rho\sqrt{gd^3})$ |
| 15 | Kinetic energy $T_K$ | $T_K/(\rho gd^2)$ |
| 16 | Potential energy $V$ | $V/(\rho gd^2)$ |
| 17 | Absolute mean square of bed velocity $\langle u_b^2\rangle$ | $\langle u_b^2\rangle/(gd)$ |
| 18 | Radiation stress $S_{xx}$ | $S_{xx}/(\rho gd^2)$ |
| 19 | Wave power $F$ | $F/(\rho g^{3/2}d^{5/2})$ |

The **absolute** bed velocity used in row 17 and in the invariant for radiation stress is not the same as the **orbital** mean square, which subtracts $\bar u_1$. The latter is a separate physical/engineering diagnostic. This distinction is important when interpreting the relative-current case.

***

## Appendix D. Source-oriented guide to the supplied scientific literature

The mathematical literature has overlapping but distinguishable scopes. The following guide identifies the role of each primary source and the most useful parts for independent verification.

| Source | Location in source | Use in this article |
|:---|:---|:---|
| Rienecker & Fenton (1981) [1] | Journal pages 119–125, especially §2 | Original finite Fourier collocation, streamline/Bernoulli equations, Newton derivatives, current closures |
| Fenton (1988) [2] | Journal pages 358–362, equation and table sections | Wavenumber-based state $z_i$, $2N+10$ equations, period versus wavelength input, height continuation, integral output variables |
| Fenton (1990) [3] | §§2–7, especially §7 | Distinction between Eulerian and Stokes current, Stokes/cnoidal comparison, highest-wave fit, corrected wave integral relations |
| Fenton (1985) [5] | Journal pages 218–225 and Tables 1–2 | Fifth-order Stokes formulation in $\epsilon=kH/2$, explicit coefficient functions and current-aware dispersion |
| Fenton (1999) [4] | Chapter §§2–3 and periodic-wave subsection | Mathematical interpretation of fully nonlinear Fourier approximation; relationship with more general numerical propagation methods |
| Fenton's program instructions [11] | §§3–4 and Appendices A–B | Physical applicability, height limits, input conventions, output scaling, Fourier interpolation, operational interpretation |
| Fenton & Rienecker (1982) [9] | §§2–3, especially Fourier representation and time advancement | Extension of Fourier approximations to unsteady nonlinear wave interactions; distinguished from steady solver |
| Westberg (1984) [10] | Research-report pp. 3–9 | Derived material accelerations, pile force and moment engineering interpretations |
| Fenton (1993) [7] | Boundary-integral formulation and results | Shoaling simulation and numerical limitations of a time-dependent boundary integral method |
| Fenton & Kennedy (1996) [8] | Proceedings pp. 1130–1143 | Local polynomial and boundary-integral methods for nonlinear propagation over varying bed topography |
| Fenton's cnoidal theory chapter [6] | Theory, numerical cnoidal theory, applicability | Long-wave expansions, elliptic functions and solitary-wave connection |

***

## References

1. **Rienecker, M. M., and Fenton, J. D. (1981).** “A Fourier approximation method for steady water waves.” *Journal of Fluid Mechanics*, **104**, 119–137. Foundational numerical Fourier/stream-function collocation formulation, nonlinear boundary conditions and current-aware steady-wave calculation. [Author-hosted paper](https://johndfenton.com/Papers/Rienecker81%2BFenton-A-Fourier-approximation-method-for-steady-water-waves.pdf).

2. **Fenton, J. D. (1988).** “The numerical solution of steady water wave problems.” *Computers & Geosciences*, **14**(3), 357–368. Wavenumber-based nondimensional state, period/wavelength and current closures, Newton solution and engineering outputs. [Author-hosted paper](https://johndfenton.com/Papers/Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf).

3. **Fenton, J. D. (1990).** “Nonlinear wave theories.” In B. Le Méhauté and D. M. Hanes (eds.), *The Sea, Volume 9: Ocean Engineering Science*, Part A, John Wiley & Sons, New York, pp. 3–25. Comparative Stokes, cnoidal and numerical Fourier theory; current definitions; wave integrals. [Author-hosted chapter](https://johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf).

4. **Fenton, J. D. (1999).** “Numerical methods for nonlinear waves.” In P. L.-F. Liu (ed.), *Advances in Coastal and Ocean Engineering*, Vol. 5, World Scientific, Singapore, pp. 241–324. Review of fully nonlinear numerical methods, including steady-wave Fourier approximation and unsteady wave propagation. [Author-hosted chapter](https://johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf).

5. **Fenton, J. D. (1985).** “A fifth-order Stokes theory for steady waves.” *Journal of Waterway, Port, Coastal, and Ocean Engineering*, **111**(2), 216–234. Explicit high-order wave-height expansion and current-aware application. [Author-hosted paper](https://johndfenton.com/Papers/Fenton85d-A-fifth-order-Stokes-theory-for-steady-waves.pdf).

6. **Fenton, J. D. (1998/1999).** “The cnoidal theory of water waves.” Chapter 2 of *Developments in Offshore Engineering*, J. B. Herbich (ed.), Gulf Publishing, Houston; the distributed chapter PDF identifies the volume as 1998. Full cnoidal-theory treatment and long-wave numerical developments. [Author-hosted chapter](https://johndfenton.com/Papers/Fenton99Cnoidal-The-cnoidal-theory-of-water-waves.pdf).

7. **Fenton, J. D. (1993).** “Simulating wave shoaling with boundary integral equations.” Conference-paper manuscript distributed by the author. Development and limitations of a boundary-integral time-evolution method. [Author-hosted paper](https://johndfenton.com/Papers/Fenton93-Simulating-wave-shoaling-with-boundary-integral-equations.pdf).

8. **Fenton, J. D., and Kennedy, A. B. (1996).** “Fast methods for computing the shoaling of nonlinear waves.” *Proceedings of the 25th International Conference on Coastal Engineering*, Vol. 1, Orlando, pp. 1130–1143. Local polynomial approximation and boundary-integral methods for nonlinear propagation over varying topography. [Author-hosted paper](https://johndfenton.com/Papers/Fenton96%2BKennedy-Fast-methods-for-computing-the-shoaling-of-nonlinear-waves.pdf).

9. **Fenton, J. D., and Rienecker, M. M. (1982).** “A Fourier method for solving nonlinear water-wave problems: Application to solitary-wave interactions.” *Journal of Fluid Mechanics*, **118**, 411–443. Time-dependent Fourier spectral model, distinct from the periodic steady-state boundary-value problem. [Author-hosted paper](https://johndfenton.com/Papers/Fenton82c%2BRienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf).

10. **Westberg, R. J., Jr. (1984).** *Application of a Fourier Approximation Method for the Solution of Steady Wave Problems to the Microcomputer.* CE 299 research report, Department of Civil Engineering, University of California, Berkeley, 7 December 1984. Physical wave kinematics, accelerations, Morison-style force and moment extensions and computational performance. (User-supplied research report.)

11. **Fenton, J. D. (2022, distributed instructions).** *Use of the Programs FOURIER, CNOIDAL and STOKES for Steady Waves.* Program instructions dated 31 October 2022. Documentation of Fourier data, output scalings, highest-wave criteria, coefficient conventions and reference cases. [Author-hosted instructions](https://johndfenton.com/Steady-waves/Instructions.pdf).

12. **Fenton, J. D., and McKee, W. D. (1990).** “On calculating the lengths of water waves.” *Coastal Engineering*, **14**, 499–513. Approximate explicit finite-depth linear-dispersion wavelength relations used for initialization and comparison. [Author-hosted paper](https://johndfenton.com/Papers/Fenton90c%2BMcKee-On-calculating-the-lengths-of-water-waves.pdf).

13. **Dean, R. G. (1965).** “Stream function representation of nonlinear ocean waves.” *Journal of Geophysical Research*, **70**, 4561–4572. Foundational stream-function numerical method preceding Rienecker–Fenton.

14. **Chappelear, J. E. (1961).** “Direct numerical calculation of wave properties.” *Journal of Geophysical Research*, **66**, 501–508. Early Fourier approach to nonlinear periodic waves.

15. **Korteweg, D. J., and de Vries, G. (1895).** “On the change of form of long waves advancing in a rectangular canal, and on a new type of long stationary waves.” *Philosophical Magazine*, Fifth Series, **39**, 422–443. Classical equation and periodic cnoidal-wave solutions.

16. **Morison, J. R., O'Brien, M. P., Johnson, J. W., and Schaaf, S. A. (1950).** “The force exerted by surface waves on piles.” *Petroleum Transactions, AIME*, **189**, 149–154. Classical drag-plus-inertia approach for slender members; referenced here as an engineering application beyond the potential-flow wave solver.

***