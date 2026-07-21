# Fenton Nonlinear Wave Suite

## Theoretical formulation, solver architecture, equations, and usage

This repository computes steady, two-dimensional, periodic finite-amplitude gravity waves in finite water depth, with optional collinear current. The reference calculation is a stream-function / Fourier collocation method in the Rienecker-Fenton and Fenton numerical-wave tradition: the Fourier basis satisfies the field equation and the bed condition, while Newton iteration solves the nonlinear free-surface boundary conditions and global wave-current constraints.

The central engineering calculation is the nonlinear mapping

$$
L=F(H,T,d,U_c),
$$

where $H$ is wave height, $T$ is the apparent period in the fixed frame, $d$ is still-water depth, $U_c$ is the imposed collinear current under the selected current convention, and $L$ is the dynamically consistent wavelength. In Fenton's formulation, a steady wave train is fundamentally defined by three length scales, $H$, $d$ and $L$. A period can replace wavelength only when the current or wave-speed convention is also specified, because the observed period is Doppler-shifted by current. The program therefore solves wavelength, wave speed, free-surface profile, Fourier coefficients, wave-frame flux, current variables, Bernoulli constant and integral quantities as one coupled nonlinear problem.

| Theoretical point | Implementation consequence |
| --- | --- |
| The steady-wave model is inviscid, incompressible, irrotational and two-dimensional over a horizontal bed. | The field equation is Laplace's equation; the nonlinear difficulty is concentrated at the free surface. |
| The Fourier terms are chosen to satisfy Laplace's equation and the bed condition exactly. | The unknowns are solved from free-surface and global constraints, not by fitting a precomputed profile. |
| If $T$ is specified instead of $L$, current must also be specified or assumed. | The solver contains explicit Eulerian-current and Stokes/mass-transport current closures. |
| The method is not a Stokes or cnoidal perturbation expansion. | Increasing $N$ increases numerical resolution rather than changing an analytical order of wave-height expansion. |
| The wavelength-only and surrogate layers are secondary tools. | They should not be used when velocities, pressures, integral quantities or near-limiting behavior are required. |

***
## 1. Repository components

| Group | Files | Role in the toolchain |
| --- | --- | --- |
| Main GUI solver | `fenton_gui.py` | Python GUI nonlinear solver and report writer. |
| Native GUI solver | `fenton_gui.cpp` | Windows C++ GUI implementation. |
| Wavelength API | `function.py` | Callable wavelength function `L(H,T,d,Uc)`. |
| Reference-style solver | `fourier.cpp` | Fenton-style standalone output files. |
| Tables | `tables.py` | Batch wavelength-table generation. |
| Pade surrogate | `pade.py`, `pade.bas` | Rational approximation model. |
| Neural surrogate | `neural.py`, `neural.bas` | Tanh-network approximation model. |
| Genetic model | `genetic.py`, `genetic.bas` | Explicit symbolic-regression wavelength formula. |
| Effective-variable model | `effective.bas` | Effective-period/effective-depth wavelength approximation. |
| Nomograms | `nomogram_plots.py`, `nomogen_plots.py` | Printable design-chart generation. |
| Excel/VBA | `*.bas`, `wavelenght.xlsm` | Spreadsheet implementation layer. |

The hierarchy is intentional. The stream-function / Fourier solver is the physical reference model. The API, Excel modules, surrogate formulas and nomograms are operational derivatives intended for faster wavelength estimation, tabulation, spreadsheet use or graphical design work.

***

## 2. Physical problem and travelling-wave formulation

The physical model is the classical steady progressive wave problem used in Fenton's papers: a homogeneous inviscid fluid, irrotational motion, incompressibility, a flat impermeable bed and a single periodic wave train. The approximation is local in the engineering sense: it is appropriate where the bed can be treated as horizontal over one wavelength and where the wave is non-breaking and changes slowly compared with the periodic motion.

A fixed frame $(x,y)$ is used with $x$ in the propagation direction, $y=0$ at the still-water level, $y=-d$ at the bed, and $y=\tilde{\eta}(x,t)$ at the free surface. A wave-frame coordinate is introduced by

$$
X_{\mathrm{dim}}=x-ct,
$$

where $c$ is the wave celerity in the fixed frame. In this frame the surface is steady:

$$
\tilde\eta(x,t)=\eta(X_{\mathrm{dim}}).
$$

The period, wavelength, wave number and celerity are related by

$$
c=\frac{L}{T}, \qquad k=\frac{2\pi}{L}, \qquad \omega=\frac{2\pi}{T}.
$$

In the program outputs, the reported celerity is the value obtained from the final solved wavelength and the specified period, so the printed geometry remains consistent with this kinematic identity.

Periodicity requires that the surface and the full flow field repeat every wavelength:

$$
\eta(X_{\mathrm{dim}}+L)=\eta(X_{\mathrm{dim}}).
$$

The specified wave height is the crest-to-trough difference,

$$
H=\eta_{\text{crest}}-\eta_{\text{trough}},
$$

which becomes, in the nondimensional collocation equations, a difference between the crest and trough ordinates.

| Input | Meaning | Typical units |
| --- | --- | --- |
| $H$ | Crest-to-trough wave height | m |
| $T$ | Apparent wave period in the fixed frame | s |
| $d$ | Still-water depth | m |
| $U_c$ | Collinear current under the selected convention | m/s |
| $N$ | Fourier order / number of half-wave intervals | dimensionless |

The current $U_c$ is not a harmless correction added after solving a no-current wave. In a steady-wave theory the period observed in the fixed frame is affected by current, and the fixed-frame horizontal velocities require a current convention. The program therefore separates Eulerian current, Stokes/mass-transport current, mean wave-frame velocity and volume flux.

***

## 3. Governing equations before discretization

### 3.1 Potential-flow equations

The physical velocity components in the fixed frame are obtained from a velocity potential $\phi$:

$$
u=\frac{\partial\phi}{\partial x}, \qquad v=\frac{\partial\phi}{\partial y}.
$$

Incompressibility requires zero divergence,

$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0,
$$

and irrotationality requires zero vorticity,

$$
\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}=0.
$$

Together they imply that the potential satisfies Laplace's equation throughout the unknown water domain:

$$
\nabla^2\phi=\frac{\partial^2\phi}{\partial x^2}+\frac{\partial^2\phi}{\partial y^2}=0,
$$

with the vertical range

$$
-d\le y\le \eta(X_{\mathrm{dim}}).
$$

The bed is impermeable, so the normal velocity there is zero:

$$
v(x,-d,t)=0.
$$

### 3.2 Free-surface kinematic condition

A fluid particle on the free surface remains on the free surface. In the fixed frame this condition is

$$
\frac{\partial\eta}{\partial t}+u\frac{\partial\eta}{\partial x}=v \qquad \text{at } y=\eta(x,t).
$$

In the wave frame, $U=u-c$ and $V=v$, and the free surface is steady. The same condition becomes

$$
U\eta_x=V \qquad \text{on } y=\eta(X_{\mathrm{dim}}),
$$

or equivalently

$$
(u-c)\eta_x=v \qquad \text{on } y=\eta(X_{\mathrm{dim}}).
$$

This is why the stream-function formulation is efficient: in the moving frame the free surface is a streamline, so the kinematic boundary condition can be enforced by making the stream function constant on the unknown surface.

### 3.3 Free-surface dynamic condition

The pressure is atmospheric and therefore constant on the free surface. Bernoulli's equation in the fixed frame is

$$
\frac{1}{2}(u^2+v^2)+gy+\frac{p}{\rho}+\frac{\partial\phi}{\partial t}=B(t).
$$

In the frame moving with the wave, the dynamic condition on the free surface can be written as

$$
\frac{1}{2}\left[(u-c)^2+v^2\right]+g\eta=R,
$$

where the arbitrary datum in the gravitational term is absorbed into the Bernoulli constant. This is the main nonlinear boundary condition: velocities are squared and evaluated on the unknown free-surface elevation. The nonlinear solver must therefore find the surface and the velocity field simultaneously.

### 3.4 Linear dispersion as initialization only

For infinitesimal-amplitude waves without current, linear Airy theory gives

$$
\omega^2=gk\tanh(kd).
$$

In wavelength form,

$$
\left(\frac{2\pi}{T}\right)^2=g\frac{2\pi}{L}\tanh\left(\frac{2\pi d}{L}\right).
$$

The familiar limiting estimates are

$$
L_0\approx\frac{gT^2}{2\pi}
$$

for deep water and

$$
c\approx\sqrt{gd}, \qquad L\approx T\sqrt{gd}
$$

for shallow water. With a collinear current, the positive intrinsic-frequency branch is

$$
\omega-kU_c=\sqrt{gk\tanh(kd)},
$$

or equivalently

$$
\left(\frac{2\pi}{T}-kU_c\right)^2=gk\tanh(kd),
$$

with branch selection required in adverse-current or near-blocking cases. These formulas are useful for initialization, checking and wavelength-only baselines. They are not the final finite-amplitude theory used by the full solver.

***

## 4. Stream-function formulation

### 4.1 Definition and boundary interpretation

For two-dimensional incompressible flow in the wave frame, the stream function $\psi$ is defined by

$$
U=\frac{\partial\psi}{\partial y}, \qquad V=-\frac{\partial\psi}{\partial x},
$$

where $U=u-c$ and $V=v$. Continuity is then satisfied identically. Irrotationality becomes

$$
\nabla^2\psi=0.
$$

With the constant in $\psi$ chosen conveniently, the impermeable bed and the moving-frame free surface are streamlines:

$$
\begin{array}{c}
\displaystyle \psi(x,-d)=0,\\[3pt]
\displaystyle \psi(x,\eta(x))=-Q,
\end{array}
$$

where $Q$ is the volume flux per unit width in the moving frame. The sign convention follows Fenton's steady-frame formulation: in the wave frame the mean flow is usually in the negative direction under a wave propagating in the positive fixed-frame direction.

### 4.2 Nondimensional wave-frame coordinates

The numerical formulation uses the wave phase and the mean-level vertical coordinate,

$$
X=k(x-ct), \qquad Y=ky, \qquad k=\frac{2\pi}{L}.
$$

One wavelength corresponds to

$$
X\in[0,2\pi],
$$

and the bed is located at

$$
Y=-kd.
$$

The nondimensional free surface is written as

$$
Y=\zeta(X)=k\eta(X).
$$

Define the nondimensional stream function

$$
\Psi=\psi\sqrt{\frac{k^3}{g}}.
$$

The nondimensional wave-frame velocities are

$$
\hat U=U\sqrt{\frac{k}{g}}=\frac{\partial\Psi}{\partial Y}, \qquad
\hat V=V\sqrt{\frac{k}{g}}=-\frac{\partial\Psi}{\partial X}.
$$

The Cauchy-Riemann relations between the nondimensional potential and stream function are

$$
\Phi_X=\Psi_Y, \qquad \Phi_Y=-\Psi_X.
$$

The field equation remains Laplace's equation in nondimensional form:

$$
\Psi_{XX}+\Psi_{YY}=0.
$$

The bed and free-surface streamline conditions are

$$
\Psi(X,-kd)=0,
$$

and

$$
\Psi(X,\zeta(X))=-Q\sqrt{\frac{k^3}{g}}.
$$

The dynamic boundary condition becomes

$$
\frac{1}{2}\left(\Psi_X^2+\Psi_Y^2\right)+\zeta(X)=\frac{Rk}{g},
$$

or, if the Bernoulli offset $r=R$ is measured relative to the selected vertical datum,

$$
\frac{1}{2}\left(\Psi_X^2+\Psi_Y^2\right)+\zeta(X)-\frac{rk}{g}=0.
$$

The deep-water limiting case is represented by

$$
kd\to\infty,
$$

with the finite-depth hyperbolic functions replaced by their exponential limits.

***

## 5. Fourier representation and collocation

### 5.1 Harmonic basis satisfying the field equation

The stream-function / Fourier method is efficient because each vertical basis function is harmonic and satisfies the bottom condition before the nonlinear equations are assembled. A finite-depth nondimensional representation consistent with Fenton's 1988 formulation is

$$
\Psi(X,Y)=-\bar U^\ast(Y+kd)+\sum_{j=1}^{N}B_j\frac{\sinh[j(Y+kd)]}{\cosh(jkd)}\cos(jX),
$$

where

$$
\bar U^\ast=\bar U\sqrt{\frac{k}{g}}.
$$

At the bed $Y=-kd$, both $Y+kd$ and the hyperbolic sine terms vanish, so

$$
\Psi(X,-kd)=0.
$$

Differentiating the series gives the nondimensional wave-frame velocity components:

$$
\begin{array}{c}
\displaystyle \hat U(X,Y)=-\bar U^\ast+
\sum_{j=1}^{N}jB_j\frac{\cosh[j(Y+kd)]}{\cosh(jkd)}\cos(jX),\\[8pt]
\displaystyle \hat V(X,Y)=
\sum_{j=1}^{N}jB_j\frac{\sinh[j(Y+kd)]}{\cosh(jkd)}\sin(jX).
\end{array}
$$

For numerical stability, the hyperbolic ratio is evaluated in the equivalent form

$$
\frac{\sinh[j(Y+kd)]}{\cosh(jkd)}
=\sinh(jY)+\cosh(jY)\tanh(jkd).
$$

The reusable finite-depth function is therefore

$$
S_j(Y)=\sinh(jY)+\cosh(jY)\tanh(jkd).
$$

Similarly,

$$
C_j(Y)=\cosh(jY)+\sinh(jY)\tanh(jkd).
$$

For large $jkd$, the code avoids overflow by using

$$
\tanh(jkd)\to1,
$$

which gives the stable limiting behavior

$$
S_j(Y)\to e^{jY}, \qquad C_j(Y)\to e^{jY}.
$$

### 5.2 Free-surface representation

A symmetric periodic wave can be represented by cosine harmonics. For the surface itself one may write

$$
\eta(x)=\sum_{n=0}^{N}a_n\cos(nkx),
$$

or in nondimensional phase form,

$$
\zeta(X)=\sum_{n=0}^{N}a_n\cos(nX).
$$

The implementation follows the Rienecker-Fenton/Fenton numerical approach more directly: it solves for the point ordinates of the free surface at collocation nodes rather than using the surface Fourier coefficients as the primary unknowns. A standard half-wave distribution is

$$
x_m=\frac{m\pi}{Nk}, \qquad m=0,1,\dots,N,
$$

or in nondimensional form

$$
X_m=\frac{m\pi}{N}, \qquad m=0,1,\dots,N.
$$

At every node,

$$
\zeta_m=\zeta(X_m)=k\eta(X_m).
$$

Using the flux convention

$$
q^\ast=\bar U^\ast kd-Q\sqrt{\frac{k^3}{g}},
$$

the kinematic free-surface condition can be written as the residual

$$
R_m^{(\psi)}=
\sum_{j=1}^{N}B_jS_j(\zeta_m)\cos(jX_m)-\bar U^\ast\zeta_m-q^\ast=0.
$$

The dynamic condition is enforced at the same nodes:

$$
R_m^{(B)}=
\frac{1}{2}\left[\hat U_m^2+\hat V_m^2\right]+\zeta_m-r^\ast=0,
$$

where

$$
r^\ast=\frac{rk}{g},
$$

and

$$
\begin{array}{c}
\displaystyle \hat U_m=-\bar U^\ast+\sum_{j=1}^{N}jB_jC_j(\zeta_m)\cos(jX_m),\\[8pt]
\displaystyle \hat V_m=\sum_{j=1}^{N}jB_jS_j(\zeta_m)\sin(jX_m).
\end{array}
$$

This separation is important when current is prescribed because the same relative wave form can correspond to different fixed-frame velocities and flux/current conventions.

***

## 6. Dimensionless variables and `SOLUTION.RES` interpretation

The solver reports quantities using two nondimensional systems: one based on wave number $k$, and one based on water depth $d$. The $k$-based form is natural for Fourier theory, while the $d$-based form is convenient for engineering interpretation and comparison between depths.

The core nondimensional engineering groups are

$$
\hat H=\frac{H}{d}, \qquad \hat L=\frac{L}{d}, \qquad \hat c=\frac{c}{\sqrt{gd}}, \qquad \hat U_c=\frac{U_c}{\sqrt{gd}}.
$$

Depth and steepness are interpreted using

$$
\frac{d}{L}, \qquad \frac{H}{L}, \qquad kd=\frac{2\pi d}{L}.
$$

The two common dimensionless periods are

$$
\tau_d=T\sqrt{\frac{g}{d}}, \qquad \tau_k=T\sqrt{gk}=\tau_d\sqrt{kd}.
$$

The associated depth-based frequency variables are

$$
\sigma_d=\frac{2\pi}{\tau_d}, \qquad \alpha_d=\sigma_d^2.
$$

The Ursell number, current Froude number and Doppler parameter used for regime interpretation are

$$
Ur=\frac{HL^2}{d^3},
$$

$$
Fr_c=\frac{U_c}{\sqrt{gd}},
$$

and

$$
D=\frac{U_cT}{L}.
$$

The intrinsic wave-current frequency is

$$
\omega-kU_c=\omega_{\text{intrinsic}}(k).
$$

These nondimensional groups explain why the finite-amplitude wavelength is generally not the linear Airy wavelength:

$$
L\ne L_{\text{Linear-Airy}}.
$$

### 6.1 Dimensionless-variable table

| Quantity | $k$-based nondimensional value | $d$-based nondimensional value |
| --- | --- | --- |
| Water depth | $kd$ | $d/d=1$ |
| Wave length | $k\lambda=2\pi$ | $\lambda/d$ |
| Wave height | $kH$ | $H/d$ |
| Wave period | $T\sqrt{gk}$ | $T\sqrt{g/d}$ |
| Wave speed | $c\sqrt{k/g}$ | $c/\sqrt{gd}$ |
| Eulerian current | $\bar u_1\sqrt{k/g}$ | $\bar u_1/\sqrt{gd}$ |
| Stokes current / mass-transport current | $\bar u_2\sqrt{k/g}$ | $\bar u_2/\sqrt{gd}$ |
| Mean wave-frame speed | $\bar U\sqrt{k/g}$ | $\bar U/\sqrt{gd}$ |
| Wave volume flux | $q\sqrt{k^3/g}$ | $q/\sqrt{gd^3}$ |
| Bernoulli offset | $rk/g$ | $r/gd$ |
| Volume flux | $Q\sqrt{k^3/g}$ | $Q/\sqrt{gd^3}$ |
| Bernoulli constant | $Rk/g$ | $R/gd$ |
| Momentum flux | $Sk^2/(\rho g)$ | $S/(\rho gd^2)$ |
| Impulse | $I\sqrt{k^3}/(\rho\sqrt{g})$ | $I/(\rho\sqrt{gd^3})$ |
| Kinetic energy | $T_k k^2/(\rho g)$ | $T_k/(\rho gd^2)$ |
| Potential energy | $V_p k^2/(\rho g)$ | $V_p/(\rho gd^2)$ |
| Bed velocity square | $u_b^2 k/g$ | $u_b^2/(gd)$ |
| Radiation stress | $S_{xx}k^2/(\rho g)$ | $S_{xx}/(\rho gd^2)$ |
| Wave power | $Fk^{5/2}/(\rho g^{3/2})$ | $F/(\rho g^{3/2}d^{5/2})$ |
| Fourier coefficients | $B_j$ | $E_j$ |

The current definitions used by Fenton can be summarized as

$$
\bar u_1=c-\bar U,
$$

and

$$
\bar u_2=c-\frac{Q}{d}.
$$

The flux variable used in the nonlinear equations is

$$
q=\bar U d-Q,
$$

and the Bernoulli offset is commonly written as

$$
r=R-gd
$$

when the Bernoulli constant is shifted from a bed-origin formulation to a mean-level formulation.

***

## 7. The `z[i]` state vector and solver description

### 7.1 Physical meaning of the `z[i]` equations

The `z[i]` equations are not empirical curve-fitting equations. They are the algebraic, nondimensional form of the steady nonlinear free-boundary problem. Solving them means finding a single state vector that simultaneously satisfies:

1. the prescribed wave height and water depth,
2. the specified wavelength or period relation,
3. the selected Eulerian-current or Stokes-current closure,
4. the mean-level convention,
5. the crest-to-trough wave-height condition,
6. the free-surface streamline condition at every collocation node,
7. Bernoulli's equation at every collocation node, and
8. the Fourier representation of a harmonic velocity field.

The unknown vector can be understood as

$$
\mathbf z=[\text{global scalars}\mid\text{surface ordinates}\mid\text{Fourier coefficients}].
$$

The primary physical correspondence in the Fenton-style 1-based state-vector convention is

$$
\begin{array}{c}
\displaystyle z_1\leftrightarrow kd, \qquad z_2\leftrightarrow kH, \qquad z_3\leftrightarrow T\sqrt{gk},\\[5pt]
\displaystyle z_4\leftrightarrow c\sqrt{\frac{k}{g}}, \qquad z_5\leftrightarrow \bar u_1\sqrt{\frac{k}{g}}, \qquad z_6\leftrightarrow \bar u_2\sqrt{\frac{k}{g}},\\[5pt]
\displaystyle z_7\leftrightarrow \bar U\sqrt{\frac{k}{g}}, \qquad z_8\leftrightarrow q\sqrt{\frac{k^3}{g}}, \qquad z_9\leftrightarrow \frac{rk}{g}.
\end{array}
$$

| Block | Indices | Physical/numerical content |
| --- | --- | --- |
| Global scalars | `z[1]...z[9]` | $kd$, $kH$, period/celerity variables, current variables, flux variable and Bernoulli offset. |
| Free surface | `z[10]...z[N+10]` | Nodal ordinates $k\eta_m$ from crest to trough. |
| Fourier modes | `z[N+11]...z[2N+10]` | Harmonic stream-function coefficients $B_j$. |

For the default $N=50$, this indexing corresponds to 110 active 1-based state entries in the Fenton-style indexing convention used by the code. Some ports pack only the independent optimization variables into a shorter zero-based vector; the physical meaning remains the same.

### 7.2 Current definitions inside the state vector

The program explicitly distinguishes Eulerian current, Stokes/mass-transport current, wave-frame mean velocity and flux. The Eulerian-current variable is connected to celerity and mean wave-frame velocity by

$$
\bar u_1=c-\bar U.
$$

The Stokes-current or mass-transport current is connected to volume flux by

$$
\bar u_2=c-\frac{Q}{d}.
$$

The wave volume flux convention is

$$
q=\bar U d-Q.
$$

This is not a cosmetic distinction. Prescribing Eulerian current and prescribing Stokes current close different mathematical problems. The residual system therefore selects the active current equation through a current-criterion index.

### 7.3 Residual equations

The residual vector is the numerical form of the complete nonlinear wave problem. At convergence, every component of

$$
\mathbf F(\mathbf z)=\mathbf 0
$$

must vanish. These equations are not optional diagnostics: they are the constraints that define the solved wave. The first eight residuals impose the global wave, current, flux, mean-level and height constraints. The next two blocks are collocation equations imposed at every free-surface node.

For Fourier order $N$, the half-wave grid contains $N+1$ collocation nodes,

$$
X_m=\frac{m\pi}{N}, \qquad m=0,1,\ldots,N.
$$

The complete residual system therefore contains

$$
8+(N+1)+(N+1)=2N+10
$$

equations. For the default value $N=50$, this gives 110 nonlinear algebraic equations.

| Residual group | Number of equations | Unknowns constrained | Physical purpose |
| --- | ---: | --- | --- |
| Global scalar equations | 8 | $z_1$ to $z_9$ and current selector | Enforce depth, height, period or wavelength, celerity, current, flux, datum and crest-trough height. |
| Streamline collocation equations | $N+1$ | Surface ordinates and Fourier coefficients | Force every free-surface node to lie on one moving-frame streamline. |
| Bernoulli collocation equations | $N+1$ | Surface ordinates, velocities and Bernoulli offset | Force constant total head along the unknown free surface. |
| Complete system | $2N+10$ | Full state vector $\mathbf z$ | Defines the nonlinear free-boundary wave solution. |

The residuals evaluated by the implementation are

$$
\begin{array}{rcl}
\displaystyle r_1&=&\displaystyle z_2-z_1(H/d),\\[3pt]
\displaystyle r_2&=&\displaystyle z_2-H_s z_3^2,\\[3pt]
\displaystyle r_3&=&\displaystyle z_4z_3-2\pi,\\[3pt]
\displaystyle r_4&=&\displaystyle z_5+z_7-z_4,\\[3pt]
\displaystyle r_5&=&\displaystyle z_1(z_6+z_7-z_4)-z_8,\\[3pt]
\displaystyle r_6&=&\displaystyle z_{c+4}-U_c\sqrt{z_1},\\[3pt]
\displaystyle r_7&=&\displaystyle z_{10}+z_{N+10}+2\sum_{i=1}^{N-1}z_{10+i},\\[3pt]
\displaystyle r_8&=&\displaystyle z_{10}-z_{N+10}-z_2,\\[3pt]
\displaystyle r_{9+m}&=&\displaystyle \psi_m-z_8-z_7z_{10+m},\\[3pt]
\displaystyle r_{N+10+m}&=&\displaystyle \frac{1}{2}\left[(-z_7+u_m)^2+v_m^2\right]+z_{10+m}-z_9,
\end{array}
$$

where $m=0,1,\ldots,N$ in the two collocation blocks.

#### 7.3.1 Notation used in the residual equations

| Symbol | Meaning in the residual system | Practical interpretation |
| --- | --- | --- |
| $m$ | Collocation-node index, $m=0,1,\ldots,N$ | Moves from crest to trough over the symmetric half-wave. |
| $H/d$ | User-specified relative height | Dimensional input height divided by still-water depth. |
| $H_s$ | Continuation-step value of the period-form height parameter | For a period-specified problem, it represents the current step of $H/(gT^2)$, so $H_s z_3^2=kH$ at convergence for that step. |
| $z_{10+m}$ | Free-surface ordinate at node $m$ | Nondimensional surface elevation $k\eta_m$ solved directly by Newton iteration. |
| $\psi_m$ | Fourier stream-function contribution evaluated at node $m$ | The sum $\sum B_jS_j(\zeta_m)\cos(jX_m)$, excluding the mean-flow term. |
| $u_m$, $v_m$ | Fourier velocity contributions at node $m$ | Used with $-z_7$ to form the total wave-frame velocities in Bernoulli's equation. |
| $z_7$ | Mean wave-frame velocity variable | Appears in both streamline and dynamic boundary conditions. |
| $z_8$ | Flux-related wave-frame constant | Couples the streamline level to the volume-flux convention. |
| $z_9$ | Bernoulli offset | The scalar total-head value, after datum shift, that must be the same at every surface node. |
| $z_{c+4}$ | Current-selected state variable | Selects the active current definition. The selector $c$ is not the wave celerity. |

#### 7.3.2 Global residual interpretation

The residual equations are listed above in the same order in which they are solved. The table below explains the role of each global residual without changing the formula notation.

| Residual | Constraint imposed | What zero residual means |
| --- | --- | --- |
| $r_1$ | Height-depth consistency | Since $z_1=kd$ and $z_2=kH$, the condition $r_1=0$ imposes $kH=kd(H/d)$. The solved nondimensional geometry therefore respects the specified physical $H$ and $d$. |
| $r_2$ | Period-form height closure | The current continuation value of $H/(gT^2)$ is consistent with $kH$ and $T\sqrt{gk}$. This is the period-specified counterpart of directly prescribing $H/L$. |
| $r_3$ | Period/celerity closure | Since $z_3=T\sqrt{gk}$ and $z_4=c\sqrt{k/g}$, $r_3=0$ imposes $ckT=2\pi$, equivalent to $c=L/T$. |
| $r_4$ | Eulerian-current identity | The state variables satisfy $\bar u_1=c-\bar U$. |
| $r_5$ | Stokes-current and flux identity | The mass-transport current, mean-flow variable and flux-related constant remain mutually consistent. The flux convention cannot drift independently of the velocity field. |
| $r_6$ | Prescribed-current equation | The user-specified current is imposed on the selected current variable. The selector $c$ in $z_{c+4}$ chooses the active current criterion and is not the wave celerity. |
| $r_7$ | Mean-level convention | The discrete free surface has the required zero-mean datum. This removes the otherwise arbitrary vertical translation of the surface. |
| $r_8$ | Crest-trough wave height | The difference between crest ordinate and trough ordinate equals the solved nondimensional height $kH$. The solver cannot converge to a nearby wave with the wrong height. |

#### 7.3.3 Collocation residual interpretation

The two collocation residual families are applied at every node $m=0,1,\ldots,N$.

| Residual family | Number of equations | Boundary condition represented | What zero residual means |
| --- | ---: | --- | --- |
| $r_{9+m}$ | $N+1$ | Kinematic free-surface condition | The unknown free surface is a streamline in the moving frame. No fluid crosses the solved surface. |
| $r_{N+10+m}$ | $N+1$ | Dynamic free-surface condition | Bernoulli head is constant along the whole solved surface, so pressure is atmospheric at every free-surface node. |

The collocation residuals are the part of the system that makes the solution fully nonlinear. The surface ordinates are unknown, and the velocities are evaluated on that unknown surface. Newton iteration therefore updates the surface profile, the Fourier coefficients, the celerity/current variables and the Bernoulli offset simultaneously.

| Solver consequence | Reason |
| --- | --- |
| The wavelength is not obtained from linear dispersion alone. | The period closure is solved together with finite-amplitude free-surface and Bernoulli constraints. |
| The current affects the nonlinear solution, not only the final reported celerity. | Current variables enter the global residuals and the wave-frame velocity in the boundary conditions. |
| The Fourier coefficients are not fitted after the surface is known. | They are solved inside the same nonlinear system that determines the surface ordinates. |
| The mean water level is fixed by the residual system. | The datum equation $r_7=0$ removes vertical translation freedom. |
| The final wave height is enforced directly. | The crest-trough equation $r_8=0$ prevents convergence to a nearby wave with the wrong height. |

Once $z_1=kd$ has converged, the wavelength follows directly:

$$
kd=\frac{2\pi d}{\lambda}.
$$

Thus

$$
\frac{\lambda}{d}=\frac{2\pi}{kd}, \qquad L=\lambda, \qquad c=\frac{L}{T}.
$$

The residual variable $z_4$ enforces the same period-celerity closure during Newton iteration; final reports and API diagnostics should still derive the dimensional celerity from the final $L/T$ relation.

### 7.4 Newton iteration and stabilized linear solve

The residual system is written compactly as

$$
\mathbf F(\mathbf z)=0.
$$

Newton iteration linearizes the residuals about the current iterate:

$$
\begin{array}{c}
\displaystyle \mathbf J(\mathbf z^{(m)})\,\Delta\mathbf z^{(m)}=-\mathbf F(\mathbf z^{(m)}),\\[5pt]
\displaystyle \mathbf z^{(m+1)}=\mathbf z^{(m)}+\Delta\mathbf z^{(m)}.
\end{array}
$$

The Jacobian entries are

$$
\mathbf J_{ij}=\frac{\partial F_i}{\partial z_j}.
$$

Fenton's 1988 program evaluates the derivatives numerically, which simplifies coding when different choices are allowed for finite/deep water, wavelength/period input and Eulerian/Stokes current specification. Some ports may use analytic or mixed analytic/numerical derivatives, but the Newton correction has the same interpretation.

The most sensitive entries involve derivatives of the streamline and Bernoulli residuals with respect to Fourier coefficients, surface ordinates, $kd$, celerity variables and the Bernoulli offset. Typical derivative families are

$$
\frac{\partial R_m^{(\psi)}}{\partial B_n}, \qquad \frac{\partial R_m^{(\psi)}}{\partial \eta_j}, \qquad \frac{\partial R_m^{(\psi)}}{\partial k},
$$

and

$$
\frac{\partial R_m^{(B)}}{\partial B_n}, \qquad \frac{\partial R_m^{(B)}}{\partial \eta_j}, \qquad \frac{\partial R_m^{(B)}}{\partial R}, \qquad \frac{\partial R_m^{(B)}}{\partial c}.
$$

For a dimensional Bernoulli residual,

$$
R_m^{(B)}=\frac{1}{2}\left[(u_m-c)^2+v_m^2\right]+g\eta_m-R,
$$

a representative chain-rule derivative is

$$
\frac{\partial R_m^{(B)}}{\partial z_j}=(u_m-c)\frac{\partial u_m}{\partial z_j}+v_m\frac{\partial v_m}{\partial z_j}+g\frac{\partial\eta_m}{\partial z_j}-\frac{\partial R}{\partial z_j}.
$$

Because near-limiting waves and adverse-current cases can make the Jacobian ill-conditioned, the solvers use stabilized linear algebra. When singular value decomposition is used,

$$
\mathbf J=\mathbf U\Sigma\mathbf V^T,
$$

and the pseudoinverse Newton correction is

$$
\Delta\mathbf z=-\mathbf V\Sigma^{-1}\mathbf U^T\mathbf F.
$$

Small singular values may be truncated or damped. In addition, continuation in wave height solves a sequence of easier problems before the final target wave height is reached. This follows Fenton's practical strategy for high or long waves: a direct jump from a linear initial solution may fail, while a sequence of increasing-height solves can converge robustly.

***

## 8. Post-processing and engineering outputs

After convergence, the nondimensional state is converted back to dimensional engineering quantities. The fundamental conversions are

$$
k=\frac{2\pi}{L}, \qquad c=\frac{L}{T}.
$$

This post-processing convention avoids treating wavelength and celerity as independent reported quantities in period-specified runs.

The moving-frame volume flux $Q$ and the wave flux variable $q$ are related by

$$
q=\bar U d-Q.
$$

The mean free-surface level is

$$
\bar\eta=\frac{1}{L}\int_0^L\eta(x)\,dx.
$$

In the usual mean-level convention,

$$
\bar\eta=0.
$$

The energy calculation separates potential and kinetic components. The total energy is

$$
E=E_p+E_k,
$$

with potential energy relative to the still-water datum

$$
E_p=\frac{\rho g}{L}\int_0^L\frac{\eta^2(x)}{2}\,dx,
$$

and kinetic energy

$$
E_k=\frac{\rho}{L}\int_0^L\int_{-d}^{\eta(x)}\frac{u^2+v^2}{2}\,dy\,dx.
$$

In the GUI report, the main integral quantities are evaluated with the following dimensional expressions:

$$
\begin{array}{rcl}
\displaystyle I&=&\displaystyle \rho(cd-Q),\\[3pt]
\displaystyle T_k&=&\displaystyle \frac{1}{2}(cI-Q\rho U_c),\\[3pt]
\displaystyle V_p&=&\displaystyle \frac{1}{2}\rho g\langle\eta^2\rangle,\\[3pt]
\displaystyle E&=&\displaystyle T_k+V_p,\\[3pt]
\displaystyle F&=&\displaystyle c(3T_k-2V_p)+\frac{1}{2}\langle u_b^2\rangle(I+\rho cd)+\frac{1}{2}\rho Q U_c^2,\\[3pt]
\displaystyle S_{xx}&=&\displaystyle 4T_k-3V_p+\rho\langle u_b^2\rangle d+2\rho I U_c,\\[3pt]
\displaystyle U_{\text{drift}}&=&\displaystyle \frac{I}{\rho d},\\[3pt]
\displaystyle C_g&=&\displaystyle \frac{F}{E}.
\end{array}
$$

These quantities should be interpreted with the same current convention used in the solve. The Eulerian-current and Stokes-current closures lead to different fixed-frame flux and drift interpretations even when the relative wave form is similar.

The mean square bed velocity diagnostic is evaluated after subtracting the imposed Eulerian current:

$$
u_{b2}=\left\langle\left(u_{\text{bed}}(t)-\bar u_1\right)^2\right\rangle.
$$

The stability diagnostics include the Miche limiting-height estimate,

$$
H_{\text{Miche}}=0.142L\tanh(kd),
$$

and the saturation ratio,

$$
\text{saturation}=\frac{H}{H_{\text{Miche}}}.
$$

The Miche estimate is a practical screening formula rather than the exact limiting wave curve. Near-limiting cases should also be judged by convergence behavior, crest sharpness, Fourier-coefficient decay and the physical plausibility of the branch selected in current.

| Output group | Examples | Use |
| --- | --- | --- |
| Geometry | $L$, $k$, $c$ | Defines the solved wave. |
| Free surface | $\eta_m$, Fourier amplitudes | Crest/trough shape and asymmetry. |
| Kinematics | $u$, $v$, acceleration | Loads, velocity checks and flow fields. |
| Integrals | $I$, $E$, $F$, $S_{xx}$ | Momentum, energy, power and radiation-stress diagnostics. |
| Stability | $Ur$, $H_{\text{Miche}}$, saturation | Practical validity screening. |

***

## 9. Wavelength-only kernels and surrogate formulas

The repository includes approximate layers for speed, spreadsheet deployment and design charts. These formulas are placed here because they are computational derivatives of the same physical mapping $L=F(H,T,d,U_c)$.

### 9.1 Current-modified linear baseline

For wavelength-only calculations, a common baseline solves

$$
(H,T,d,U_c)\mapsto L.
$$

With a collinear current, the baseline residual in $k$ is

$$
F(k)=\left(\frac{2\pi}{T}-kU_c\right)^2-gk\tanh(kd).
$$

Newton iteration is then

$$
k_{n+1}=k_n-\frac{F(k_n)}{F'(k_n)},
$$

where

$$
F'(k)=-2U_c\left(\frac{2\pi}{T}-kU_c\right)-g\tanh(kd)-gkd\frac{1}{\cosh^2(kd)}.
$$

The wavelength is recovered from

$$
L=\frac{2\pi}{k}.
$$

Without current, the same baseline can be expressed as the fixed-point wavelength equation

$$
L_{\text{lin}}=\frac{gT^2}{2\pi}\tanh\!\left(\frac{2\pi d}{L_{\text{lin}}}\right).
$$

### 9.2 Common surrogate structure

Most reduced-order models correct a physically meaningful baseline:

$$
L=L_{\text{base}}C(\mathbf x).
$$

Typical nondimensional features include

$$
\frac{H}{L}, \qquad \frac{d}{L}, \qquad kd, \qquad \frac{U_c}{c}, \qquad \frac{U_c}{\sqrt{gd}}.
$$

The implemented engineering feature set includes

$$
\begin{array}{c}
\displaystyle \text{WaveSteepness}=\frac{H}{L_{\text{lin}}},\\
\displaystyle \text{RelativeDepth}=\frac{d}{L_{\text{lin}}},\\
\displaystyle \text{DopplerFactor}=\frac{U_cT}{L_{\text{lin}}},\\
\displaystyle \text{UrsellNumber}=\frac{HL_{\text{lin}}^2}{d^3},\\
\displaystyle \text{CurrentFroude}=\frac{U_c}{\sqrt{gd}}.
\end{array}
$$

For symbolic-regression and neural workflows, additional logarithmic and current features may be used:

$$
\begin{array}{c}
\displaystyle x_0=\ln\!\left(\frac{d}{L}\right), \qquad x_1=\ln\!\left(\frac{H}{L}\right), \qquad x_2=\ln\!\left(\frac{H}{d}\right),\\
\displaystyle x_3=\ln(Ur), \qquad Ur=\frac{HL^2}{d^3},\\
\displaystyle x_4=\frac{U_c}{\sqrt{gd}}, \qquad x_5=\frac{U_cT}{L}, \qquad x_6=\frac{U_c}{C_0},
\end{array}
$$

with the deep-water scale quantities

$$
C_0=\frac{gT}{2\pi}, \qquad L_0=\frac{gT^2}{2\pi},
$$

and the additional features

$$
x_7=\ln\!\left(\frac{H}{L_0}\right), \qquad x_8=\ln\!\left(T\sqrt{\frac{g}{d}}\right).
$$

### 9.3 Pade/rational surrogate

The Pade surrogate represents a correction with a rational function,

$$
P(x)=\frac{a_0+a_1x+a_2x^2+\cdots+a_mx^m}{1+b_1x+b_2x^2+\cdots+b_nx^n}.
$$

The depth-regime variable is usually

$$
\mu=\frac{d}{L}.
$$

The regime split used for training and interpretation is

$$
\begin{array}{c}
\displaystyle \mu<0.05 \quad \text{(shallow)},\\
\displaystyle 0.05\le\mu<0.5 \quad \text{(intermediate)},\\
\displaystyle \mu\ge0.5 \quad \text{(deep)}.
\end{array}
$$

The resulting Pade wavelength expression is

$$
L=L_{\text{base}}C_{\text{Pade}}(\mathbf x),
$$

or explicitly

$$
L=L_{\text{base}}\frac{p_0+p_1x_1+p_2x_2+\cdots}{1+q_1x_1+q_2x_2+\cdots}.
$$

### 9.4 Neural-network surrogate

The neural surrogate is a global tanh-network approximation,

$$
\hat y(\mathbf x)=b_o+\sum_{j=1}^{M}w_j^{(o)}\tanh\!\left(b_j+\sum_{i=1}^{p}w_{ji}^{(h)}x_i\right).
$$

The inputs may be standardized as

$$
\tilde x_i=\frac{x_i-\mu_i}{\sigma_i},
$$

or min-max scaled as

$$
\tilde x_i=\frac{x_i-x_i^{\min}}{x_i^{\max}-x_i^{\min}}.
$$

Training minimizes mean-square loss,

$$
\mathcal L=\frac{1}{N}\sum_{n=1}^{N}(y_n-\hat y_n)^2,
$$

and prediction quality is commonly summarized by

$$
\mathrm{MAPE}=\frac{100}{N}\sum_{n=1}^{N}\left|\frac{y_n-\hat y_n}{y_n}\right|.
$$

### 9.5 Genetic symbolic-regression formula

The genetic symbolic-regression model is an explicit correction over a baseline,

$$
L\approx L_{\text{base}}C_{\text{GEP}}(\mathbf x),
$$

or, more generally,

$$
L\approx G(\mathbf x).
$$

The documented explicit expression is

$$
L_{\text{genetic}}=L_{\text{lin}}\left[\exp\!\left(\frac{345U_cT}{509L_{\text{lin}}}\right)+\frac{1052\pi U_c}{1052\pi U_c+435gT}-\frac{5\tanh\!\left(\frac{U_cT}{L_{\text{lin}}}\ln\!\left(T\sqrt{\frac{g}{d}}\right)+\frac{2117}{961}\right)}{67\ln\!\left(\frac{H}{d}\right)}\right].
$$

For readability and implementation checking, it can be decomposed as

$$
L_{\text{genetic}}=L_{\text{lin}}M,
$$

where

$$
\begin{array}{c}
\displaystyle M=t_1+t_2+t_3,\\
\displaystyle t_1=\exp\!\left(\frac{345U_cT}{509L_{\text{lin}}}\right), \qquad t_2=\frac{1052\pi U_c}{1052\pi U_c+435gT},\\
\displaystyle t_3=-\frac{5\tanh\!\left(\frac{U_cT}{L_{\text{lin}}}\ln\!\left(T\sqrt{\frac{g}{d}}\right)+\frac{2117}{961}\right)}{67\ln\!\left(\frac{H}{d}\right)}.
\end{array}
$$

### 9.6 Effective transformed-variable formula

The effective-variable model transforms the input period and depth, then solves a linear dispersion relation with transformed variables. The transformed period and depth are

$$
\begin{array}{c}
\displaystyle T_{\text{eff}}=T-\frac{93}{943}\left(U_c^2-\frac{H}{1819\cdot708}\right)+\frac{500}{739}U_c,\\
\displaystyle d_{\text{eff}}=d+\frac{H^{3/4}}{d^{1/4}}+\frac{986}{873}U_c\tanh(H).
\end{array}
$$

Equivalently,

$$
T_{\text{eff}}=T+\Delta_T, \qquad d_{\text{eff}}=d+\Delta_d,
$$

with

$$
\begin{array}{c}
\displaystyle \Delta_T=\left(U_c^2-\frac{H}{1819\cdot708}\right)\left(-\frac{93}{943}\right)+\frac{500}{739}U_c,\\
\displaystyle \Delta_d=\sqrt{\frac{H}{\sqrt{d/H}}}+\frac{986}{873}U_c\tanh(H)=\frac{H^{3/4}}{d^{1/4}}+\frac{986}{873}U_c\tanh(H).
\end{array}
$$

The algebraic identity in the depth transformation is

$$
\sqrt{\frac{H}{\sqrt{d/H}}}=\frac{H^{3/4}}{d^{1/4}}.
$$

The final effective wavelength is obtained from

$$
L_{\text{effective}}=\mathcal L_{\text{lin}}(T_{\text{eff}},d_{\text{eff}}),
$$

where

$$
\left(\frac{2\pi}{T_{\text{eff}}}\right)^2=gk\tanh(kd_{\text{eff}}), \qquad L_{\text{effective}}=\frac{2\pi}{k}.
$$

The transformed period must remain nonzero:

$$
T_{\text{eff}}\ne0.
$$

***

## 10. Program features

| Feature | Description |
| --- | --- |
| Fully nonlinear stream-function solve | Solves surface shape, velocity field, wavelength and integral quantities. |
| Current-aware formulation | Separates Eulerian, Stokes, mean-flow and flux conventions. |
| Continuation | Ramps wave height to improve nonlinear convergence. |
| Fourier order | Default $N=50$ for high-resolution free-surface representation. |
| GUI output | Reports geometry, kinematics, integrals, stability and classification. |
| Wavelength API | Provides callable `L(H,T,d,Uc)` for scripts and tables. |
| Surrogates | Provides Pade, neural, symbolic-regression and effective-variable approximations. |
| Excel/VBA support | Allows spreadsheet deployment of linear, nonlinear and surrogate methods. |
| Nomograms | Generates design charts for rapid engineering use. |

***

## 11. Installation and use

### 11.1 Python GUI

```bash
python -m venv build_env
build_env\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -U "numpy>=1.20"
python -m pip install -U "numba>=0.57"
python fenton_gui.py
```

Numba is optional. If Numba is unavailable, the solver falls back to the pure-NumPy path.

### 11.2 Build the C++ GUI executable

```bash
g++ fenton_gui.cpp -o fenton_gui.exe -O3 -std=c++20 ^
   -march=native -flto=auto -fopenmp -static -static-libgcc ^
   -static-libstdc++ -mwindows -pthread ^
   -lgdi32 -luser32 -lkernel32 -lcomctl32
```

For reproducible no-current and current-case comparisons, keep the same CPU family, compiler family and optimization flags.

### 11.3 Build a Python standalone executable

```bash
pip install pyinstaller
pyinstaller --onefile --noconsole fenton_gui.py
```

The executable is written to the `dist` directory.

### 11.4 Wavelength-only use

Use `function.py` when only the wavelength is required by another program or table. The intended callable form is:

```python
L(H, T, d, Uc)
```

***

## 12. Validity and numerical cautions

| Condition | Practical implication |
| --- | --- |
| Non-breaking steady waves | The formulation is for a periodic wave of permanent form; breaking, strong turbulence and air entrainment are outside the model. |
| Irrotational flow | Strong shear, vorticity, boundary-layer dominance or rapidly varying currents require a different model. |
| Horizontal-bed assumption | The steady solver is a local horizontal-bed theory; shoaling over variable bathymetry requires an additional propagation/shoaling model or a local approximation. |
| Very long waves | Fourier coefficients become broad-banded as the solitary-wave limit is approached; increase $N$ and inspect coefficient decay. |
| Near highest waves | Sharp crests degrade convergence and conditioning; use continuation and inspect the free-surface shape. |
| Opposing current | Wave blocking, multiple branches or no physical branch may occur; the chosen root must be checked. |
| Current criterion | Do not mix Eulerian-current and Stokes/mass-transport current definitions. |
| Surrogate use | Do not extrapolate outside the training/calibration range. Surrogates are wavelength-estimation tools, not replacements for the nonlinear kinematic solution. |
| Units | Use coherent SI inputs throughout. |

Fenton's theory comparisons give the following practical hierarchy. Fifth-order Stokes theory is useful for shorter waves, cnoidal theory for sufficiently long waves, with the Ursell number often used as a regime guide. The Fourier approximation method is the preferred reference method when high accuracy is required over a broad range of depths and heights, especially for velocities, pressures, integral quantities or waves close to the boundary of the perturbation theories.

The nonlinear solver is therefore the preferred method when free-surface shape, velocities, accelerations, energy, radiation stress or near-breaking behavior matter. The wavelength-only kernels, symbolic formulas, Pade model, neural model and effective-variable model are appropriate for fast wavelength estimates only when their calibration range is respected.

***

## 13. Dedication

This work is dedicated to **Dr John D. Fenton** (1945-), Australian civil engineer, hydraulic engineer, and applied mathematician, in recognition of his distinguished and continuing scientific contributions to mathematical hydraulics, computational fluid mechanics, and nonlinear water-wave theory. Further information, publications, technical resources, and software by Dr Fenton are available on his official website: [https://johndfenton.com](https://johndfenton.com).

***

## 14. References

1. **Fenton, J. D. (1999).** *Numerical methods for nonlinear waves.* In P. L.-F. Liu (ed.), *Advances in Coastal and Ocean Engineering*, Vol. 5, World Scientific, Singapore, pp. 241–324.  
   **Relevance:** Review of numerical methods for nonlinear waves, including the Fourier approximation method for steadily propagating periodic waves and more general propagation methods.  
   **URL:** [Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf](https://johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf)

2. **Fenton, J. D. (1999).** *The cnoidal theory of water waves.* In J. B. Herbich (ed.), *Developments in Offshore Engineering*, Gulf Publishing, Houston, pp. 275–337.  
   **Relevance:** Modern presentation of cnoidal theory for long finite-depth waves, including practical application, accuracy limits and numerical cnoidal ideas.  
   **URL:** [Fenton99Cnoidal-The-cnoidal-theory-of-water-waves.pdf](https://johndfenton.com/Papers/Fenton99Cnoidal-The-cnoidal-theory-of-water-waves.pdf)

3. **Fenton, J. D. & Kennedy, A. B. (1996).** *Fast methods for computing the shoaling of nonlinear waves.* In *Proceedings of the 25th International Conference on Coastal Engineering*, Vol. 1, Orlando, pp. 1130–1143.  
   **Relevance:** Nonlinear wave propagation and shoaling over varying bathymetry.  
   **URL:** [Fenton96+Kennedy-Fast-methods-for-computing-the-shoaling-of-nonlinear-waves.pdf](https://johndfenton.com/Papers/Fenton96%2BKennedy-Fast-methods-for-computing-the-shoaling-of-nonlinear-waves.pdf)

4. **Fenton, J. D. (1995).** *A numerical cnoidal theory for steady water waves.* In *Proceedings of the 12th Australasian Coastal and Ocean Engineering Conference*, Melbourne, pp. 175–180.  
   **Relevance:** Numerical cnoidal formulation for finite-depth long-wave calculations.  
   **URL:** [Fenton95-A-numerical-cnoidal-theory-for-steady-water-waves.pdf](https://johndfenton.com/Papers/Fenton95-A-numerical-cnoidal-theory-for-steady-water-waves.pdf)

5. **Townsend, M. & Fenton, J. D. (1995).** *Numerical comparisons of wave analysis methods.* In *Proceedings of the 12th Australasian Coastal and Ocean Engineering Conference*, Melbourne, pp. 169–173.  
   **Relevance:** Wave analysis and inversion from pressure measurements; conditioning and method comparison.  
   **URL:** [Townsend95+Fenton-Numerical-comparisons-of-wave-analysis-methods.pdf](https://johndfenton.com/Papers/Townsend95%2BFenton-Numerical-comparisons-of-wave-analysis-methods.pdf)

6. **Kennedy, A. B. & Fenton, J. D. (1995).** *Simulation of the propagation of surface gravity waves using local polynomial approximation.* In *Proceedings of the 12th Australasian Coastal and Ocean Engineering Conference*, Melbourne, pp. 287–292.  
   **Relevance:** Nonlinear wave propagation and shoaling over varying bathymetry.  
   **URL:** [Kennedy95+Fenton-Simulation-of-the-propagation-of-surface-gravity-waves-using-local-polynomial-approximation.pdf](https://johndfenton.com/Papers/Kennedy95%2BFenton-Simulation-of-the-propagation-of-surface-gravity-waves-using-local-polynomial-approximation.pdf)

7. **Fenton, J. D. (1993).** *Simulating wave shoaling with boundary integral equations.* In *Proceedings of the 11th Australasian Conference on Coastal and Ocean Engineering*, Townsville, pp. 71–76.  
   **Relevance:** Boundary-integral formulation for nonlinear wave transformation, with singularity subtraction.  
   **URL:** [Fenton93-Simulating-wave-shoaling-with-boundary-integral-equations.pdf](https://johndfenton.com/Papers/Fenton93-Simulating-wave-shoaling-with-boundary-integral-equations.pdf)

8. **Fenton, J. D. (1990).** *Nonlinear wave theories.* In B. Le Méhauté & D. M. Hanes (eds.), *The Sea: Ocean Engineering Science, Part A*, Vol. 9, Wiley, New York, pp. 3–25.  
   **Relevance:** Key source for current definitions, steady-wave governing equations, Stokes/cnoidal regime guidance, Fourier approximation methods and integral properties.  
   **URL:** [Fenton90b-Nonlinear-wave-theories.pdf](https://johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf)

9. **Drennan, W. M., Fenton, J. D. & Donelan, M. A. (1990).** *Numerical simulation of nonlinear wave groups.* In *Proceedings of the 11th Annual Conference of the Canadian Applied Mathematics Society*, Halifax.  
   **Relevance:** Relevant to nonlinear and free-surface wave modelling.  
   **URL:** [Drennan90-Numerical-simulation-of-nonlinear-wave-groups.pdf](https://johndfenton.com/Papers/Drennan90-Numerical-simulation-of-nonlinear-wave-groups.pdf)

10. **Fenton, J. D. & McKee, W. D. (1990).** *On calculating the lengths of water waves.* *Coastal Engineering*, **14**, 499–513.  
    **Relevance:** Wavelength determination for nonlinear waves in finite depth.  
    **URL:** [Fenton90c+McKee-On-calculating-the-lengths-of-water-waves.pdf](https://johndfenton.com/Papers/Fenton90c%2BMcKee-On-calculating-the-lengths-of-water-waves.pdf)

11. **Fenton, J. D. (1988).** *The numerical solution of steady water wave problems.* *Computers and Geosciences*, **14**, 357–368.  
    **Relevance:** Core practical program for steady periodic nonlinear waves; finite/deep water, wavelength/period input, Eulerian or Stokes current, Newton solution with numerical Jacobian.  
    **URL:** [Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf](https://johndfenton.com/Papers/Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf)

12. **Fenton, J. D. (1986).** *Polynomial approximation and water waves.* In *Proceedings of the 20th International Conference on Coastal Engineering*, Vol. 1, Taipei, pp. 193–207.  
    **Relevance:** Relevant to nonlinear and free-surface wave modelling.  
    **URL:** [Fenton86-Polynomial-approximation-and-water-waves.pdf](https://johndfenton.com/Papers/Fenton86-Polynomial-approximation-and-water-waves.pdf)

13. **Fenton, J. D. (1985).** *A fifth-order Stokes theory for steady waves.* *Journal of Waterway, Port, Coastal, and Ocean Engineering*, **111**, 216–234.  
    **Relevance:** Closed-form fifth-order Stokes expansion for steady waves in finite depth.  
    **URL:** [Fenton85d-A-fifth-order-Stokes-theory-for-steady-waves.pdf](https://johndfenton.com/Papers/Fenton85d-A-fifth-order-Stokes-theory-for-steady-waves.pdf)

14. **Fenton, J. D. (1983).** *On the application of steady wave theories.* In *Proceedings of the 6th Australasian Conference on Coastal and Ocean Engineering*, Christchurch, pp. 65–70.  
    **Relevance:** Guidance on applicability and limits of steady-wave theories in engineering.  
    **URL:** [Fenton83-On-the-application-of-steady-wave-theories.pdf](https://johndfenton.com/Papers/Fenton83-On-the-application-of-steady-wave-theories.pdf)

15. **Fenton, J. D. & Rienecker, M. M. (1982).** *A Fourier method for solving nonlinear water wave problems.* *Journal of Fluid Mechanics*, **118**, 411–443.  
    **Relevance:** Time-dependent Fourier method for full nonlinear water-wave equations, used for solitary-wave interactions; relevant background but not the main steady-wave collocation solver.  
    **URL:** [Fenton82c+Rienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf](https://johndfenton.com/Papers/Fenton82c%2BRienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf)

16. **Schwartz, L. W. & Fenton, J. D. (1982).** *Strongly-nonlinear waves.* In M. Van Dyke, J. V. Wehausen & J. L. Lumley (eds.), *Annual Review of Fluid Mechanics*, **14**, 39–60.  
    **Relevance:** Fundamental properties and approximations for strongly nonlinear wave motion.  
    **URL:** [Schwartz82-Strongly-nonlinear-waves.pdf](https://johndfenton.com/Papers/Schwartz82-Strongly-nonlinear-waves.pdf)

17. **Rienecker, M. M. & Fenton, J. D. (1981).** *A Fourier approximation method for steady water waves.* *Journal of Fluid Mechanics*, **104**, 119–137.  
    **Relevance:** Foundational Fourier approximation method for steadily progressing periodic waves; finite Fourier series, collocation of nonlinear free-surface conditions and Newton solution.  
    **URL:** [Rienecker81+Fenton-A-Fourier-approximation-method-for-steady-water-waves.pdf](https://johndfenton.com/Papers/Rienecker81%2BFenton-A-Fourier-approximation-method-for-steady-water-waves.pdf)

18. **Fenton, J. D. & Rienecker, M. M. (1980).** *Accurate numerical solutions for nonlinear waves.* In *Proceedings of the 17th International Conference on Coastal Engineering*, Sydney, pp. 50–69.  
    **Relevance:** Benchmark accurate numerical solutions for nonlinear wave profiles and kinematics.  
    **URL:** [Fenton80+Rienecker-Accurate-numerical-solutions-for-nonlinear-waves.pdf](https://johndfenton.com/Papers/Fenton80%2BRienecker-Accurate-numerical-solutions-for-nonlinear-waves.pdf)

19. **Fenton, J. D. (1979).** *A high-order cnoidal wave theory.* *Journal of Fluid Mechanics*, **94**, 129–161.  
    **Relevance:** Cnoidal-wave theory for finite-depth long waves and associated numerical formulations.  
    **URL:** [Fenton79-A-high-order-cnoidal-wave-theory.pdf](https://johndfenton.com/Papers/Fenton79-A-high-order-cnoidal-wave-theory.pdf)

20. **Fenton, J. D. & Mills, D. A. (1976).** *Shoaling waves: numerical solution of exact equations.* In D. G. Provis & R. Radok (eds.), *Proceedings of the IUTAM Symposium on Waves on Water of Variable Depth*, Canberra, Springer-Verlag, pp. 93–100.  
    **Relevance:** Nonlinear wave propagation and shoaling over varying bathymetry.  
    **URL:** [Fenton76+Mills-Shoaling-waves-Numerical-solution-of-exact-equations.pdf](https://johndfenton.com/Papers/Fenton76%2BMills-Shoaling-waves-Numerical-solution-of-exact-equations.pdf)

21. **Fenton, J. D. (1972).** *A ninth-order solution for the solitary wave.* *Journal of Fluid Mechanics*, **53**, 257–271.  
    **Relevance:** High-order analytic and series solution for solitary waves.  
    **URL:** [Fenton72-A-ninth-order-solution-for-the-solitary-wave.pdf](https://johndfenton.com/Papers/Fenton72-A-ninth-order-solution-for-the-solitary-wave.pdf)
