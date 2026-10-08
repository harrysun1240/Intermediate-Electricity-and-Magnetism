# Chapter 3: Special Techniques (The Electrostatic Potential)

**Contents**

1. [Why Laplace's equation](#1-why-laplaces-equation)
2. [Two key properties of Laplace's equation](#2-two-key-properties-of-laplaces-equation)
3. [Method of relaxation](#3-method-of-relaxation)
4. [Uniqueness theorems](#4-uniqueness-theorems)
5. [The method of images](#5-the-method-of-images)
6. [Separation of variables (Cartesian coordinates)](#6-separation-of-variables-cartesian-coordinates)
7. [Separation of variables (spherical coordinates)](#7-separation-of-variables-spherical-coordinates)
8. [Separation of variables (cylindrical coordinates)](#8-separation-of-variables-cylindrical-coordinates)
9. [Multipole expansion](#9-multipole-expansion)

---

## 1. Why Laplace's equation

In most real problems it is easier to solve a differential equation for $V$ than to do the integral for $V$ directly.

The potential of a charge distribution $\rho$ is:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int \frac{\rho(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\, d^3r'
```

This has two problems:

- The integral is often hard to evaluate.
- With conductors, we usually do not know $\rho$. The charge moves around on the conductor until it reaches equilibrium, and we don't know where it ends up.

**Poisson's equation.** Combine Gauss's law, $\nabla\cdot\vec{E} = \rho/\varepsilon_0$, with $\vec{E} = -\nabla V$. Then $\nabla\cdot(-\nabla V) = \rho/\varepsilon_0$, which gives:

```math
\nabla^2 V = -\frac{\rho}{\varepsilon_0}
```

We solve this together with suitable boundary conditions.

**Laplace's equation.** Very often we want $V$ in a region with no charge ($\rho = 0$). Then Poisson's equation becomes Laplace's equation:

```math
\nabla^2 V = \frac{\partial^2 V}{\partial x^2} + \frac{\partial^2 V}{\partial y^2} + \frac{\partial^2 V}{\partial z^2} = 0
```

## 2. Two key properties of Laplace's equation

### Property 1: the mean value property

If $V$ satisfies Laplace's equation, then $V$ at any point equals the average of $V$ over any sphere centered on that point (as long as the sphere contains no charge):

```math
V(\vec{r}) = \frac{1}{4\pi R^2}\oint_{\text{sphere}} V\, da
```

**Proof.** First find the average potential over a sphere caused by a single point charge $q$ outside it. Then use superposition.

**Step 1: set up coordinates.** Put the origin at the center of the sphere (radius $R$). Put $q$ on the $z$-axis at distance $z > R$, so the charge sits at $\vec{r}_q = z\,\hat{z}$. A point on the sphere is:

```math
\vec{r}_s = R\sin\theta\cos\phi\,\hat{x} + R\sin\theta\sin\phi\,\hat{y} + R\cos\theta\,\hat{z}
```

**Step 2: distance from $q$ to a point on the sphere.**

```math
|\vec{r}_s - \vec{r}_q| = \sqrt{R^2\sin^2\theta + (R\cos\theta - z)^2} = \sqrt{z^2 + R^2 - 2zR\cos\theta}
```

So the potential on the sphere is:

```math
V(R,\theta,\phi) = \frac{q}{4\pi\varepsilon_0\sqrt{z^2 + R^2 - 2zR\cos\theta}}
```

**Step 3: average over the sphere.** The area element on the sphere is $da = R^2\sin\theta\,d\theta\,d\phi$. Nothing depends on $\phi$, so the $\phi$ integral just gives $2\pi$:

```math
V_{\text{avg}} = \frac{1}{4\pi R^2}\int_0^{2\pi}\!d\phi\int_0^{\pi}\!d\theta\, R^2\sin\theta\, V = \frac{2\pi R^2}{4\pi R^2}\,\frac{q}{4\pi\varepsilon_0}\int_0^{\pi}\frac{\sin\theta\, d\theta}{\sqrt{z^2+R^2-2zR\cos\theta}}
```

**Step 4: substitute.** Let $u = z^2 + R^2 - 2zR\cos\theta$. Then $du = 2zR\sin\theta\,d\theta$, so $\sin\theta\,d\theta = du/(2zR)$. The limits change as follows: $\theta = 0$ gives $u = (z-R)^2$, and $\theta = \pi$ gives $u = (z+R)^2$.

```math
V_{\text{avg}} = \frac{q}{8\pi\varepsilon_0}\cdot\frac{1}{2zR}\int_{(z-R)^2}^{(z+R)^2} u^{-1/2}\, du = \frac{q}{8\pi\varepsilon_0}\cdot\frac{1}{2zR}\Big[2u^{1/2}\Big]_{(z-R)^2}^{(z+R)^2}
```

**Step 5: evaluate.** Because $z > R$, the square root of $(z-R)^2$ is $z - R$ (not $R - z$). So:

```math
V_{\text{avg}} = \frac{q}{8\pi\varepsilon_0}\cdot\frac{2}{2zR}\Big[(z+R) - (z-R)\Big] = \frac{q}{8\pi\varepsilon_0}\cdot\frac{2R}{zR} = \frac{q}{4\pi\varepsilon_0 z}
```

This is exactly the potential of $q$ at the center of the sphere.

**Step 6: superposition.** Any charge distribution outside the sphere is a sum of point charges. Each one gives "average = value at the center," so the sum does too. No charge is allowed inside the sphere, because Laplace's equation needs $\rho = 0$ there. This proves the property.

### Property 2: no local maxima or minima

If $V$ satisfies Laplace's equation in a region, $V$ has no local maximum or minimum inside that region. The extreme values can only happen on the boundary.

**Proof.** Suppose $V$ had a local maximum at point $\vec{r}$. Then we could draw a small sphere around $\vec{r}$ on which every value of $V$ is smaller than $V(\vec{r})$. The average over that sphere would then be smaller than $V(\vec{r})$. This contradicts Property 1. The same argument, with "larger" instead of "smaller," rules out a minimum.

## 3. Method of relaxation

The mean value property gives a simple numerical way to solve Laplace's equation.

1. Start with the known values of $V$ on the boundary, and reasonable guesses for $V$ on a grid of interior points.
2. Sweep through the grid. Replace $V$ at each point with the average of its nearest neighbors.
3. Repeat until the numbers barely change.

**Why this works.** On a 2D grid with spacing $h$, the second derivatives can be approximated by finite differences:

```math
\nabla^2 V \approx \frac{V(x+h,y) + V(x-h,y) + V(x,y+h) + V(x,y-h) - 4V(x,y)}{h^2}
```

Setting this to zero gives $V(x,y) =$ the average of its four neighbors. So the grid version of Laplace's equation is exactly the averaging rule used in step 2.

## 4. Uniqueness theorems

Laplace's equation alone does not fix $V$. We also need boundary conditions. "Suitable" means strong enough to fix the answer, but not so strong that they contradict each other. The uniqueness theorems tell us which conditions are enough.

### First uniqueness theorem

The solution to $\nabla^2 V = 0$ in a volume $\mathcal{V}$ is unique if $V$ is specified everywhere on the boundary surface.

- The volume may contain cavities, as long as $V$ is known on their surfaces too.
- The outer surface may be at infinity.
- Specifying $V$ on the surface is called a **Dirichlet boundary condition**.

**Proof.**

1. Suppose there were two solutions, $V_1$ and $V_2$, with $\nabla^2 V_1 = 0$ and $\nabla^2 V_2 = 0$.
2. Define $V_3 = V_1 - V_2$. Since $\nabla^2$ is linear, $\nabla^2 V_3 = \nabla^2 V_1 - \nabla^2 V_2 = 0$. So $V_3$ also satisfies Laplace's equation.
3. On the boundary, $V_1$ and $V_2$ both equal the given values. So $V_3 = 0$ everywhere on the boundary.
4. By Property 2, $V_3$ has no maximum or minimum inside the volume. Its largest and smallest values are on the boundary, and both are $0$.
5. So $V_3 = 0$ everywhere inside, which means $V_1 = V_2$.

Note: this shows a solution is unique, but not that one exists. Proving existence is harder (it has been done).

### Extension to Poisson's equation

Uniqueness still holds if there is charge inside ($\rho \neq 0$), as long as $\rho$ is known.

1. Take two solutions: $\nabla^2 V_1 = -\rho/\varepsilon_0$ and $\nabla^2 V_2 = -\rho/\varepsilon_0$.
2. Then $V_3 = V_1 - V_2$ gives $\nabla^2 V_3 = -\rho/\varepsilon_0 - (-\rho/\varepsilon_0) = 0$.
3. $V_3$ satisfies Laplace's equation and is $0$ on the boundary, so by the argument above $V_3 = 0$ and $V_1 = V_2$.

**Result:** $V$ is uniquely determined in a volume if (a) $\rho$ is known throughout the region and (b) $V$ is known on all boundaries.

### Second uniqueness theorem

In a volume $\mathcal{V}$ where $\rho$ is known, surrounded by conductors with known total charge $Q_i$ on each one, the electric field $\vec{E}$ is uniquely determined. The outer boundary can be a conductor or infinity.

**Proof.**

**Step 1: two candidate fields.** Suppose $\vec{E}_1$ and $\vec{E}_2$ both work. In the volume:

```math
\nabla\cdot\vec{E}_1 = \frac{\rho}{\varepsilon_0}, \qquad \nabla\cdot\vec{E}_2 = \frac{\rho}{\varepsilon_0}
```

**Step 2: Gauss's law on each surface.** For a surface just enclosing the $i$-th conductor, both fields give the same charge:

```math
\oint_{S_i}\vec{E}_1\cdot d\vec{a} = \frac{Q_i}{\varepsilon_0}, \qquad \oint_{S_i}\vec{E}_2\cdot d\vec{a} = \frac{Q_i}{\varepsilon_0}
```

For the outer boundary, both enclose the total charge $Q_{\text{tot}} = \sum_i Q_i + \int\rho\, d\tau$:

```math
\oint_{\text{outer}}\vec{E}_1\cdot d\vec{a} = \frac{Q_{\text{tot}}}{\varepsilon_0}, \qquad \oint_{\text{outer}}\vec{E}_2\cdot d\vec{a} = \frac{Q_{\text{tot}}}{\varepsilon_0}
```

**Step 3: the difference field.** Define $\vec{E}_3 = \vec{E}_1 - \vec{E}_2$. Subtracting the equations above:

```math
\nabla\cdot\vec{E}_3 = 0 \text{ in } \mathcal{V}, \qquad \oint_{S_i}\vec{E}_3\cdot d\vec{a} = 0 \text{ over each boundary surface}
```

**Step 4: a product rule.** Both fields are electrostatic, so we can write $\vec{E}_3 = -\nabla V_3$ with $V_3 = V_1 - V_2$. Using the product rule $\nabla\cdot(f\vec{A}) = f(\nabla\cdot\vec{A}) + \vec{A}\cdot\nabla f$:

```math
\nabla\cdot(V_3\vec{E}_3) = V_3(\nabla\cdot\vec{E}_3) + \vec{E}_3\cdot\nabla V_3 = 0 - E_3^2 = -E_3^2
```

**Step 5: divergence theorem.** Integrate over the volume:

```math
\int_{\mathcal{V}}\nabla\cdot(V_3\vec{E}_3)\,d\tau = \oint_S V_3\vec{E}_3\cdot d\vec{a} = -\int_{\mathcal{V}} E_3^2\, d\tau
```

**Step 6: the surface term is zero.** Each conductor is an equipotential for both solutions. So $V_3 = V_1 - V_2$ is a constant on each conducting surface (but not necessarily the same constant on different ones). It can come out of each surface integral:

```math
\oint_S V_3\vec{E}_3\cdot d\vec{a} = \sum_i V_3^{(i)}\oint_{S_i}\vec{E}_3\cdot d\vec{a} = 0
```

Each integral is zero by Step 3. (If the outer boundary is at infinity, $V_3 \to 0$ there, so that part is also zero.)

**Step 7: conclude.** So $\int E_3^2\, d\tau = 0$. Since $E_3^2 \geq 0$ everywhere, the only way the integral can vanish is $\vec{E}_3 = 0$ everywhere in the volume. So $\vec{E}_1 = \vec{E}_2$.

## 5. The method of images

The trick: find any setup whose potential matches the boundary conditions in the region we care about. By uniqueness, that potential is the answer there.

### Point charge above a grounded conducting plane

A charge $q$ sits at height $d$ above an infinite grounded conducting plane (the $xy$-plane). What is $V$ for $z > 0$?

The answer is **not** just $q/(4\pi\varepsilon_0|\vec{r} - d\hat{z}|)$. The charge $q$ induces a surface charge $\sigma(x,y)$ on the plane, and that also contributes to $V$.

We need to solve Poisson's equation for $z > 0$, with a single point charge at $d\hat{z}$, subject to:

- (a) $V = 0$ when $z = 0$ (the plane is grounded)
- (b) $V \to 0$ far from the charge (for $z > 0$)

**The image solution.** Remove the plane. Put the real charge $+q$ at $d\hat{z}$ and an image charge $-q$ at $-d\hat{z}$:

```math
V(\vec{r}) = \frac{q}{4\pi\varepsilon_0}\left(\frac{1}{|\vec{r}-d\hat{z}|} - \frac{1}{|\vec{r}+d\hat{z}|}\right) = \frac{q}{4\pi\varepsilon_0}\left(\frac{1}{\sqrt{x^2+y^2+(z-d)^2}} - \frac{1}{\sqrt{x^2+y^2+(z+d)^2}}\right)
```

Check the conditions:

- At $z = 0$, the two distances are equal, so $V = 0$.
- Far away ($x^2 + y^2 + z^2 \gg d^2$), both terms go to $0$, so $V \to 0$.
- The only charge in the region $z > 0$ is $q$ at $d\hat{z}$. (The image charge is at $z < 0$, outside the region.)

So this $V$ fits every requirement of the original problem in the region $z > 0$ (but not for $z < 0$, where the real $V$ is $0$). By uniqueness, it is **the** solution for $z > 0$.

### Induced surface charge

Just above a conductor, the field is $\vec{E} = (\sigma/\varepsilon_0)\,\hat{n}$. Here $\hat{n} = \hat{z}$ and $\vec{E} = -\nabla V$, so:

```math
\sigma = -\varepsilon_0 \left.\frac{\partial V}{\partial z}\right|_{z\to 0}
```

Take the derivative:

```math
\frac{\partial V}{\partial z} = \frac{q}{4\pi\varepsilon_0}\left(\frac{-(z-d)}{[x^2+y^2+(z-d)^2]^{3/2}} + \frac{(z+d)}{[x^2+y^2+(z+d)^2]^{3/2}}\right)
```

At $z = 0$ the two terms are each $d/(x^2 + y^2 + d^2)^{3/2}$, so they add:

```math
\left.\frac{\partial V}{\partial z}\right|_{z=0} = \frac{qd}{2\pi\varepsilon_0(x^2+y^2+d^2)^{3/2}} \quad\Rightarrow\quad \sigma(x,y) = \frac{-qd}{2\pi(x^2+y^2+d^2)^{3/2}}
```

**Total induced charge.** Use polar coordinates on the plane: $s^2 = x^2 + y^2$, $da = s\,ds\,d\phi$.

```math
q_{\text{ind}} = \int_0^{2\pi}\!d\phi\int_0^{\infty}\! s\,ds\,\frac{-qd}{2\pi(s^2+d^2)^{3/2}} = -qd\int_0^{\infty}\frac{s\,ds}{(s^2+d^2)^{3/2}}
```

Let $u = s^2 + d^2$, so $du = 2s\,ds$, and $u$ runs from $d^2$ to $\infty$:

```math
q_{\text{ind}} = -qd\cdot\frac{1}{2}\int_{d^2}^{\infty}u^{-3/2}\,du = -qd\cdot\frac{1}{2}\Big[-2u^{-1/2}\Big]_{d^2}^{\infty} = -qd\cdot\frac{1}{2}\cdot\frac{2}{d} = -q
```

The plane carries a total induced charge of exactly $-q$, as expected.

### Force and energy

The charge $q$ is attracted toward the plane. The force is the same as the force from the image charge $-q$ at distance $2d$:

```math
\vec{F} = \frac{1}{4\pi\varepsilon_0}\frac{q(-q)}{(2d)^2}\,\hat{z} = -\frac{1}{4\pi\varepsilon_0}\frac{q^2}{4d^2}\,\hat{z}
```

**Energy** = the work we must do to bring $q$ in from infinity to height $d$. At height $z$ the force on $q$ is $\vec{F}(z) = -\frac{q^2}{4\pi\varepsilon_0 (4z^2)}\,\hat{z}$. We must push with $-\vec{F}$, so:

```math
W = -\int_{\infty}^{d}\vec{F}\cdot d\vec{l} = \int_{\infty}^{d}\frac{q^2}{4\pi\varepsilon_0}\frac{dz}{4z^2} = \frac{q^2}{4\pi\varepsilon_0}\cdot\frac{1}{4}\Big[-\frac{1}{z}\Big]_{\infty}^{d} = -\frac{q^2}{4\pi\varepsilon_0(4d)}
```

This is **half** the energy of a real $`q`$ and $`-q`$ pair at separation $`2d`$ (which is $`-q^2/(4\pi\varepsilon_0 \cdot 2d)`$). The reason: the real field exists only for $`z > 0`$ ($`\vec{E} = 0`$ inside the conductor), so only half the field energy of the two-charge system is actually there.

### Other image problems

- Any fixed charge distribution near a grounded conducting plane: use its mirror image with the opposite charge density.
- Conducting spheres and cylinders: we must figure out where to put the image charges. They can never be placed in the region where we are calculating $V$.

### Point charge near a grounded conducting sphere

A charge $q$ sits on the $z$-axis at distance $a$ from the center of a grounded conducting sphere of radius $R$ ($a > R$). Find $V$ outside the sphere. (Inside, $V = 0$.)

Put an image charge $q'$ on the $z$-axis inside the sphere at $b\hat{z}$, with $b < R$. The potential is:

```math
V(r,\theta) = \frac{1}{4\pi\varepsilon_0}\left(\frac{q}{\sqrt{r^2+a^2-2ar\cos\theta}} + \frac{q'}{\sqrt{r^2+b^2-2br\cos\theta}}\right)
```

Can we choose $q'$ and $b$ so that $V(R,\theta) = 0$ for every $\theta$?

Setting $V = 0$ at $r = R$ and squaring gives:

```math
\frac{R^2+b^2-2bR\cos\theta}{R^2+a^2-2aR\cos\theta} = \left(\frac{q'}{q}\right)^2 = \text{constant}
```

For this ratio to be the same for every $\theta$, the top must be a fixed multiple of the bottom. So the constant terms and the $\cos\theta$ terms must be in the same ratio:

```math
\frac{R^2+b^2}{R^2+a^2} = \frac{bR}{aR} = \frac{b}{a}
```

Cross-multiply: $a(R^2 + b^2) = b(R^2 + a^2)$. Rearranging gives $(a - b)(R^2 - ab) = 0$. Since $b \neq a$:

```math
b = \frac{R^2}{a} = \left(\frac{R}{a}\right)R \;<\; R
```

The constant ratio is then $b/a = R^2/a^2$, so $(q'/q)^2 = R^2/a^2$. The image must have the opposite sign to cancel $q$ on the sphere:

```math
q' = -\frac{R}{a}\,q
```

## 6. Separation of variables (Cartesian coordinates)

The idea: build the solution as a sum of simple pieces. Each piece solves Laplace's equation and meets some of the boundary conditions. Then pick the combination that meets the rest.

### Orthogonal and complete functions

- A set of functions $f_n(y)$ is **complete** if any other function $f(y)$ can be written as a linear combination of them: $f(y) = \sum_{n=1}^{\infty} C_n f_n(y)$.
- Example: $\sin(n\pi y/a)$ is complete on $0 \le y \le a$. (Proving completeness for a given set is hard.)
- A set is **orthogonal** if the integral of the product of any two different members is zero:

```math
\int_0^a f_n(y)\,f_m(y)\,dy = 0 \quad\text{for } n \neq m
```

**Goal:** find a convenient complete, orthogonal set to express the solution. Ideally each function solves Laplace's equation and also satisfies some boundary conditions. Then the problem reduces to finding the linear combination that satisfies the remaining ones.

### Separating the equation

Convenient basis functions are products, such as $X(x)Y(y)Z(z)$ or $R(r)\Theta(\theta)\Phi(\phi)$. Try $V(x,y,z) = X(x)Y(y)Z(z)$:

```math
0 = \nabla^2 V = X''YZ + XY''Z + XYZ''
```

Divide by $V = XYZ$ (careful where $V = 0$):

```math
\frac{X''}{X} + \frac{Y''}{Y} + \frac{Z''}{Z} = 0
```

The first term depends only on $x$, the second only on $y$, the third only on $z$. The sum is zero for all $x$, $y$, $z$. Changing $x$ alone can't change the other two terms, so the first term must be a constant. The same goes for the others:

```math
\frac{X''}{X} = C_1, \qquad \frac{Y''}{Y} = C_2, \qquad \frac{Z''}{Z} = C_3, \qquad C_1 + C_2 + C_3 = 0
```

Each is an ordinary differential equation. For $X'' = C_1 X$:

| Case | Solution |
| --- | --- |
| $C_1 > 0$ | $Ae^{\sqrt{C_1}\,x} + Be^{-\sqrt{C_1}\,x}$ (growing/decaying) |
| $C_1 = 0$ | $Ax + B$ |
| $C_1 < 0$ | $A\cos(\sqrt{-C_1}\,x) + B\sin(\sqrt{-C_1}\,x)$ (oscillating) |

The same holds for $Y$ and $Z$. The boundary conditions decide which case to use in each direction. Try to satisfy as many as possible with each basis function.

### Example: rectangular metal pipe

An infinitely long rectangular metal pipe runs along $x$. Its cross-section is $0 \le y \le a$ and $0 \le z \le b$. The pipe walls are grounded, and the end at $x = 0$ is held at potential $V_0(y,z)$. Find $V$ inside.

**Boundary conditions:**

| Where | Condition |
| --- | --- |
| Bottom, $y = 0$ | $V = 0$ |
| Top, $y = a$ | $V = 0$ |
| Back, $z = 0$ | $V = 0$ |
| Front, $z = b$ | $V = 0$ |
| Far end, $x \to \infty$ | $V \to 0$ |
| Near end, $x = 0$ | $V = V_0(y,z)$ |

**Choose the cases.** $V$ must vanish at two points in $y$ (and in $z$). An exponential or linear function can't do that unless it is zero everywhere. So use oscillating functions in $y$ and $z$, and exponentials in $x$. Set $C_2 = -k^2$, $C_3 = -l^2$, $C_1 = k^2 + l^2$:

```math
X(x) = A e^{\sqrt{k^2+l^2}\,x} + B e^{-\sqrt{k^2+l^2}\,x}, \quad Y(y) = C\sin(ky) + D\cos(ky), \quad Z(z) = E\sin(lz) + F\cos(lz)
```

**Apply the boundary conditions one at a time:**

- $V \to 0$ as $x \to \infty$ ⇒ $A = 0$
- $V = 0$ at $y = 0$ ⇒ $D = 0$
- $V = 0$ at $z = 0$ ⇒ $F = 0$
- $V = 0$ at $y = a$ ⇒ $\sin(ka) = 0$ ⇒ $k = n\pi/a$, $n = 1, 2, 3, \ldots$
- $V = 0$ at $z = b$ ⇒ $\sin(lb) = 0$ ⇒ $l = m\pi/b$, $m = 1, 2, 3, \ldots$

So each basis function is:

```math
V_{nm}(x,y,z) = C_{nm}\, e^{-\pi\sqrt{(n/a)^2+(m/b)^2}\,x}\,\sin\!\left(\frac{n\pi y}{a}\right)\sin\!\left(\frac{m\pi z}{b}\right)
```

Each one satisfies Laplace's equation inside the pipe and every boundary condition except the one at $x = 0$.

**Combine.** The full solution is a sum of these:

```math
V(x,y,z) = \sum_{n=1}^{\infty}\sum_{m=1}^{\infty} V_{nm}(x,y,z)
```

**Fit the last condition.** At $x = 0$:

```math
V_0(y,z) = \sum_{n=1}^{\infty}\sum_{m=1}^{\infty} C_{nm}\sin\!\left(\frac{n\pi y}{a}\right)\sin\!\left(\frac{m\pi z}{b}\right)
```

To pick out one coefficient, multiply both sides by $\sin(n'\pi y/a)\sin(m'\pi z/b)$ and integrate over the cross-section. Use the orthogonality relation:

```math
\int_0^L \sin\!\left(\frac{n\pi u}{L}\right)\sin\!\left(\frac{n'\pi u}{L}\right)du = \frac{L}{2}\,\delta_{nn'}
```

On the right, only the $n = n'$, $m = m'$ term survives, giving $\frac{a}{2}\cdot\frac{b}{2}\,C_{n'm'}$. So:

```math
C_{nm} = \frac{4}{ab}\int_0^a\!dy\int_0^b\!dz\; V_0(y,z)\sin\!\left(\frac{n\pi y}{a}\right)\sin\!\left(\frac{m\pi z}{b}\right)
```

## 7. Separation of variables (spherical coordinates)

With azimuthal symmetry, every solution of Laplace's equation is a sum of $r^l$ and $r^{-(l+1)}$ terms times Legendre polynomials.

### Setting up

Laplace's equation in spherical coordinates:

```math
\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial V}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial V}{\partial\theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2 V}{\partial\phi^2} = 0
```

Try $V = R(r)P(\theta)Q(\phi)$. Substitute, divide by $RPQ$, and multiply by $r^2$:

```math
\frac{1}{R}\frac{d}{dr}\left(r^2\frac{dR}{dr}\right) + \frac{1}{P\sin\theta}\frac{d}{d\theta}\left(\sin\theta\frac{dP}{d\theta}\right) + \frac{1}{Q\sin^2\theta}\frac{d^2Q}{d\phi^2} = 0
```

**Limit to azimuthal symmetry** ($V$ does not depend on $\phi$). The last term drops. The first term depends only on $r$ and the second only on $\theta$, so each must be a constant. Choose the constant to be $l(l+1)$:

```math
\frac{1}{R}\frac{d}{dr}\left(r^2\frac{dR}{dr}\right) = l(l+1), \qquad \frac{1}{P\sin\theta}\frac{d}{d\theta}\left(\sin\theta\frac{dP}{d\theta}\right) = -l(l+1)
```

(Why $l(l+1)$? Trying $R = r^\alpha$ in the radial equation gives $\alpha(\alpha+1) =$ constant. Writing the constant as $l(l+1)$ makes the two roots simply $\alpha = l$ and $\alpha = -(l+1)$.)

### Radial solution

```math
R(r) = A r^l + \frac{B}{r^{l+1}}
```

**Check:** $dR/dr = Al\,r^{l-1} - (l+1)B/r^{l+2}$. Multiply by $r^2$:

```math
r^2\frac{dR}{dr} = Al\,r^{l+1} - \frac{(l+1)B}{r^l} \;\Rightarrow\; \frac{d}{dr}\left(r^2\frac{dR}{dr}\right) = Al(l+1)r^l + \frac{l(l+1)B}{r^{l+1}} = l(l+1)R \;\checkmark
```

### Angular solution: Legendre polynomials

The solutions are $P(\theta) = P_l(\cos\theta)$, $l = 0, 1, 2, 3, \ldots$ They are given by the Rodrigues formula:

```math
P_l(x) = \frac{1}{2^l\, l!}\left(\frac{d}{dx}\right)^l (x^2-1)^l
```

| $l$ | $P_l(x)$ |
| --- | --- |
| $0$ | $1$ |
| $1$ | $x$ |
| $2$ | $\frac{1}{2}(3x^2 - 1)$ |
| $3$ | $\frac{1}{2}(5x^3 - 3x)$ |

Properties:

- $P_l(x)$ is a polynomial of order $l$. It has only even powers of $x$ if $l$ is even, and only odd powers if $l$ is odd.
- $P_l(1) = 1$.
- They are complete and orthogonal on $-1 \le x \le 1$:

```math
\int_{-1}^{1} P_l(x)\,P_{l'}(x)\,dx = \frac{2}{2l+1}\,\delta_{ll'}
```

The angular equation has other solutions, but they blow up at $\theta = 0$ or $\theta = \pi$, so they are usually ruled out on physical grounds. For example, the second solution for $l = 0$ is $\ln(\tan(\theta/2))$. This excludes the second solution for $l = 0, 1, 2, \ldots$ and both solutions for non-integer $l$.

### General solution (azimuthal symmetry)

```math
V(r,\theta) = \sum_{l=0}^{\infty}\left(A_l r^l + \frac{B_l}{r^{l+1}}\right)P_l(\cos\theta)
```

### Example 1: inside a sphere with $V_0(\theta)$ on its surface

Find $V$ inside a hollow sphere of radius $R$ with $V_0(\theta)$ given on the surface.

- $B_l = 0$ for every $l$, or $V$ would blow up at the origin.
- We need: $\sum_l A_l R^l P_l(\cos\theta) = V_0(\theta)$.

**Fourier trick.** Multiply both sides by $P_{l'}(\cos\theta)\sin\theta$ and integrate from $0$ to $\pi$:

```math
\sum_{l=0}^{\infty} A_l R^l \int_0^{\pi} P_l(\cos\theta)P_{l'}(\cos\theta)\sin\theta\,d\theta = \int_0^{\pi} V_0(\theta)P_{l'}(\cos\theta)\sin\theta\,d\theta
```

Let $x = \cos\theta$, so $dx = -\sin\theta\,d\theta$. The limits $\theta = 0 \to \pi$ become $x = 1 \to -1$, and the minus sign flips them back to $-1 \to 1$. By orthogonality only $l = l'$ survives on the left, giving $\frac{2}{2l'+1}A_{l'}R^{l'}$. So:

```math
A_l = \frac{2l+1}{2R^l}\int_0^{\pi} V_0(\theta)\,P_l(\cos\theta)\sin\theta\,d\theta
```

### Example 2: outside a sphere with $V_0(\theta)$ on its surface

```math
B_l = \frac{(2l+1)R^{l+1}}{2}\int_0^{\pi} V_0(\theta)\,P_l(\cos\theta)\sin\theta\,d\theta
```

### Example 3: uncharged metal sphere in a uniform field

An uncharged metal sphere of radius $R$ is centered at the origin in a uniform field $\vec{E}_0 = E_0\,\hat{z}$. Find $V$ outside.

**Fix the constants by symmetry.**

- The sphere is an equipotential. Set $V_{\text{sphere}} = 0$.
- Far away, $\vec{E} \to E_0\,\hat{z}$, so $V \to -E_0 z + C$ for some constant $C$.
- Reflect the whole setup in the $xy$-plane and flip the sign of every charge. You get back the same setup. That operation reverses the in-plane part of $\vec{E}$ on the $xy$-plane itself, so that part must be zero. So $\vec{E}$ points along $\hat{z}$ at every point of the $xy$-plane.
- Then $\int\vec{E}\cdot d\vec{l} = 0$ along any path in the $xy$-plane, so the plane is an equipotential.
- The plane touches the sphere (on its equator), where $V = 0$. So $V = 0$ on the whole $xy$-plane, and $C = 0$.

**Boundary conditions:** $V = 0$ at $r = R$, and $V \to -E_0 r\cos\theta$ for $r \gg R$.

**Solve.** At $r = R$: $A_l R^l + B_l/R^{l+1} = 0$, so $B_l = -A_l R^{2l+1}$. Then:

```math
V(r,\theta) = \sum_{l=0}^{\infty} A_l\left(r^l - \frac{R^{2l+1}}{r^{l+1}}\right)P_l(\cos\theta)
```

For large $r$ this becomes $\sum_l A_l r^l P_l(\cos\theta)$. Match it to $-E_0 r\cos\theta$ term by term, using $P_1(\cos\theta) = \cos\theta$: $A_1 = -E_0$, and $A_l = 0$ for $l \neq 1$.

```math
V(r,\theta) = -E_0\left(r - \frac{R^3}{r^2}\right)\cos\theta
```

The first term is the applied field. The second is from the charge induced on the sphere.

**Induced charge.** $\partial V/\partial r = -E_0(1 + 2R^3/r^3)\cos\theta$. At $r = R$ this is $-3E_0\cos\theta$, so:

```math
\sigma(\theta) = -\varepsilon_0\left.\frac{\partial V}{\partial r}\right|_{r=R} = 3\varepsilon_0 E_0\cos\theta
```

It is positive on the northern half and negative on the southern half.

### Example 4: surface charge glued on a spherical shell

A surface charge $\sigma_0(\theta)$ is glued onto a spherical shell of radius $R$. Find $V$ inside and outside.

Inside (finite at the origin) and outside (zero at infinity):

```math
V_{\text{in}}(r,\theta) = \sum_{l=0}^{\infty} A_l r^l P_l(\cos\theta), \qquad V_{\text{out}}(r,\theta) = \sum_{l=0}^{\infty}\frac{B_l}{r^{l+1}}P_l(\cos\theta)
```

**Condition 1: $V$ is continuous at $r = R$.** Match the coefficient of each $P_l$:

```math
A_l R^l = \frac{B_l}{R^{l+1}} \;\Rightarrow\; B_l = A_l R^{2l+1}
```

**Condition 2: the jump in $`\vec{E}`$.** Across a surface charge, $`\vec{E}_{\text{out}} - \vec{E}_{\text{in}} = (\sigma/\varepsilon_0)\,\hat{r}`$. With $`\vec{E} = -\nabla V`$:

```math
-\left.\frac{\partial V_{\text{out}}}{\partial r}\right|_{R} + \left.\frac{\partial V_{\text{in}}}{\partial r}\right|_{R} = \frac{\sigma_0(\theta)}{\varepsilon_0}
```

Take the derivatives and use $B_l = A_l R^{2l+1}$:

```math
\sum_{l=0}^{\infty}\left[\frac{(l+1)B_l}{R^{l+2}} + lA_lR^{l-1}\right]P_l = \sum_{l=0}^{\infty}\big[(l+1) + l\big]A_lR^{l-1}P_l = \sum_{l=0}^{\infty}(2l+1)R^{l-1}A_l\,P_l(\cos\theta) = \frac{\sigma_0(\theta)}{\varepsilon_0}
```

Multiply by $P_{l'}(\cos\theta)\sin\theta$ and integrate. The orthogonality integral gives $\frac{2}{2l+1}\delta_{ll'}$, which cancels the $(2l+1)$:

```math
A_l = \frac{1}{2\varepsilon_0 R^{l-1}}\int_0^{\pi}\sigma_0(\theta)\,P_l(\cos\theta)\sin\theta\,d\theta
```

**Case 1: $\sigma_0(\theta) = k\cos\theta = kP_1(\cos\theta)$.** Only $l = 1$ survives:

```math
\int_0^{\pi}\sigma_0(\theta)\,P_l(\cos\theta)\sin\theta\,d\theta = k\int_{-1}^{1}P_l(u)\,P_1(u)\,du = k\cdot\frac{2}{2(1)+1}\,\delta_{l1} = \frac{2k}{3}\,\delta_{l1}
```

So:

```math
A_l = \frac{k}{2\varepsilon_0}\cdot\frac{2}{3}\,\delta_{l1} = \frac{k}{3\varepsilon_0}\,\delta_{l1} \;\Rightarrow\; V_{\text{in}} = \frac{k}{3\varepsilon_0}\,r\cos\theta, \quad V_{\text{out}} = \frac{kR^3}{3\varepsilon_0}\frac{\cos\theta}{r^2}
```

**Case 2: $\sigma_0(\theta) = \sigma_0$ (constant) $= \sigma_0 P_0(\cos\theta)$.** Only $l = 0$ survives:

```math
\int_0^{\pi}\sigma_0\,P_l(\cos\theta)\sin\theta\,d\theta = \sigma_0\int_{-1}^{1}P_l(u)\,P_0(u)\,du = \sigma_0\cdot\frac{2}{2(0)+1}\,\delta_{l0} = 2\sigma_0\,\delta_{l0}
```

So:

```math
A_l = \frac{2\sigma_0}{2\varepsilon_0 R^{-1}}\,\delta_{l0} = \frac{\sigma_0 R}{\varepsilon_0}\,\delta_{l0} \;\Rightarrow\; V_{\text{in}} = \frac{\sigma_0 R}{\varepsilon_0}, \quad V_{\text{out}} = \frac{\sigma_0 R^2}{\varepsilon_0 r}
```

Check: the total charge is $Q = 4\pi R^2\sigma_0$, so $V_{\text{out}} = Q/(4\pi\varepsilon_0 r)$, the potential of a point charge, as expected.

## 8. Separation of variables (cylindrical coordinates)

In cylindrical coordinates the radial part gives Bessel functions. When nothing depends on $z$, it reduces to powers of $s$ and $\ln s$.

### Separating the equation

Laplace's equation in cylindrical coordinates $(s, \phi, z)$:

```math
\frac{\partial^2 V}{\partial s^2} + \frac{1}{s}\frac{\partial V}{\partial s} + \frac{1}{s^2}\frac{\partial^2 V}{\partial\phi^2} + \frac{\partial^2 V}{\partial z^2} = 0
```

Try $V = S(s)Q(\phi)Z(z)$ and divide by $V$:

```math
\underbrace{\frac{1}{S}\frac{d^2S}{ds^2} + \frac{1}{sS}\frac{dS}{ds} + \frac{1}{Qs^2}\frac{d^2Q}{d\phi^2}}_{\text{function of } s,\,\phi \text{ only}} + \underbrace{\frac{1}{Z}\frac{d^2Z}{dz^2}}_{\text{function of } z \text{ only}} = 0
```

So the $z$-part is a constant. Call it $+k^2$:

```math
\frac{1}{Z}\frac{d^2Z}{dz^2} = k^2, \qquad \frac{1}{S}\frac{d^2S}{ds^2} + \frac{1}{sS}\frac{dS}{ds} + \frac{1}{Qs^2}\frac{d^2Q}{d\phi^2} = -k^2
```

Multiply the second equation by $s^2$. Now one group depends only on $s$ and one only on $\phi$:

```math
\underbrace{\frac{s^2}{S}\frac{d^2S}{ds^2} + \frac{s}{S}\frac{dS}{ds} + k^2s^2}_{s \text{ only}} + \underbrace{\frac{1}{Q}\frac{d^2Q}{d\phi^2}}_{\phi \text{ only}} = 0
```

Set the $\phi$-part equal to $-\nu^2$.

**Summary of the three equations:**

```math
\frac{d^2Z}{dz^2} - k^2 Z = 0, \qquad \frac{d^2Q}{d\phi^2} + \nu^2 Q = 0, \qquad \frac{d^2S}{ds^2} + \frac{1}{s}\frac{dS}{ds} + \left(k^2 - \frac{\nu^2}{s^2}\right)S = 0
```

### Solutions

```math
Z(z) = A_1 e^{kz} + A_2 e^{-kz} \quad (k \neq 0), \qquad Q(\phi) = B_1 e^{i\nu\phi} + B_2 e^{-i\nu\phi} = B_1'\cos(\nu\phi) + B_2'\sin(\nu\phi)
```

$k$ can be anything so far. But $\nu$ must be an integer so that $Q$ is single-valued: going once around, $\phi = 0$ and $\phi = 2\pi$ are the same point, so $Q(0) = Q(2\pi)$.

**Radial equation.** Define $\rho = ks$. Then $d/ds = k\,d/d\rho$, and dividing through by $k^2$:

```math
\frac{d^2S}{d\rho^2} + \frac{1}{\rho}\frac{dS}{d\rho} + \left(1 - \frac{\nu^2}{\rho^2}\right)S = 0
```

This is **Bessel's equation**. For $k \neq 0$ its general solution is $A\,J_\nu(ks) + B\,N_\nu(ks)$, where $J_\nu$ are Bessel functions of the first kind and $N_\nu$ are of the second kind.

### Simplification: $z$-independent problems ($k = 0$)

$Q(\phi)$ has the same form as above, but the radial equation simplifies to:

```math
\frac{d^2S}{ds^2} + \frac{1}{s}\frac{dS}{ds} - \frac{\nu^2}{s^2}S = 0
```

**$\nu \neq 0$.** Try $S = s^\alpha$:

```math
\alpha(\alpha-1)s^{\alpha-2} + \alpha s^{\alpha-2} - \nu^2 s^{\alpha-2} = 0 \;\Rightarrow\; \alpha^2 - \nu^2 = 0 \;\Rightarrow\; \alpha = \pm\nu
```

This gives two independent solutions: $S = A s^\nu + B s^{-\nu}$.

**$\nu = 0$.** Now $s^\nu$ and $s^{-\nu}$ are the same (a constant), so we need a second solution. The equation can be written as:

```math
\frac{d^2S}{ds^2} + \frac{1}{s}\frac{dS}{ds} = \frac{1}{s}\frac{d}{ds}\left(s\frac{dS}{ds}\right) = 0 \;\Rightarrow\; s\frac{dS}{ds} = C \;\Rightarrow\; S = C\ln s + D
```

So the second solution is $\ln s$.

For the angular equation with $\nu = 0$: $Q'' = 0$, so $Q = \alpha_0\phi + \alpha_1$. If $\phi$ covers the full range $0$ to $2\pi$, we must drop $\alpha_0\phi$, since it would make $Q(2\pi) \neq Q(0)$.

### General solution ($z$-independent, full $\phi$ range)

```math
V(s,\phi) = A_0 + B_0\ln s + \sum_{m=1}^{\infty}\Big[s^m(a_m\cos m\phi + b_m\sin m\phi) + s^{-m}(c_m\cos m\phi + d_m\sin m\phi)\Big]
```

## 9. Multipole expansion

Far from a charge distribution, $V$ can be written as a series in $1/r$: the monopole term falls as $1/r$, the dipole as $1/r^2$, and so on.

### Expanding $1/|\vec{r} - \vec{r}\,'|$

We start from the usual integral, but now the field point $\vec{r}$ is far away compared to the size of the charge distribution:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\,d^3r'
```

Let $\theta'$ be the angle between $\vec{r}$ and $\vec{r}\,'$, so $\vec{r}\,'\cdot\hat{r} = r'\cos\theta'$. By the law of cosines:

```math
|\vec{r}-\vec{r}\,'| = \left(r^2 + r'^2 - 2rr'\cos\theta'\right)^{1/2} = r\left(1 + \frac{r'}{r}\Big(\frac{r'}{r} - 2\cos\theta'\Big)\right)^{1/2}
```

If the largest $r'$ is much smaller than $r$, we can use a Taylor series. Define:

```math
\epsilon = \frac{r'}{r}\left(\frac{r'}{r} - 2\cos\theta'\right)
```

The binomial series gives:

```math
\frac{1}{|\vec{r}-\vec{r}\,'|} = \frac{1}{r}(1+\epsilon)^{-1/2} = \frac{1}{r}\left(1 - \frac{1}{2}\epsilon + \frac{3}{8}\epsilon^2 - \frac{5}{16}\epsilon^3 + \cdots\right)
```

**Collect powers of $r'/r$:**

- Order $(r'/r)^1$: only $-\epsilon/2$ contributes, giving $-\frac{1}{2}(-2\cos\theta') = \cos\theta'$.
- Order $(r'/r)^2$: $-\epsilon/2$ gives $-\frac{1}{2}$, and $\frac{3}{8}\epsilon^2$ gives $\frac{3}{8}(4\cos^2\theta')$. Total: $\frac{1}{2}(3\cos^2\theta' - 1)$.
- Order $(r'/r)^3$: $\frac{3}{8}\epsilon^2$ gives $\frac{3}{8}(2)(-2\cos\theta') = -\frac{3}{2}\cos\theta'$, and $-\frac{5}{16}\epsilon^3$ gives $-\frac{5}{16}(-8\cos^3\theta') = \frac{5}{2}\cos^3\theta'$. Total: $\frac{1}{2}(5\cos^3\theta' - 3\cos\theta')$.

```math
\frac{1}{|\vec{r}-\vec{r}\,'|} = \frac{1}{r}\left\{1 + \frac{r'}{r}\cos\theta' + \left(\frac{r'}{r}\right)^2\frac{3\cos^2\theta'-1}{2} + \left(\frac{r'}{r}\right)^3\frac{5\cos^3\theta'-3\cos\theta'}{2} + \cdots\right\}
```

The coefficients are exactly the Legendre polynomials $P_0$, $P_1$, $P_2$, $P_3$. This holds to all orders:

```math
\frac{1}{|\vec{r}-\vec{r}\,'|} = \frac{1}{r}\sum_{n=0}^{\infty}\left(\frac{r'}{r}\right)^n P_n(\cos\theta')
```

### The multipole expansion of $V$

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\sum_{n=0}^{\infty}\frac{1}{r^{n+1}}\int (r')^n P_n(\cos\theta')\,\rho(\vec{r}\,')\,d^3r'
```

| $n$ | Name | Falls off as | Simplest example |
| --- | --- | --- | --- |
| $0$ | Monopole | $1/r$ | A single charge. Coefficient $= \int\rho\,d^3r' = Q$, the total charge |
| $1$ | Dipole | $1/r^2$ | $+q$ and $-q$ side by side. Coefficient $= \int r'\cos\theta'\,\rho\,d^3r'$, the dipole moment |
| $2$ | Quadrupole | $1/r^3$ | Four charges on a square, alternating $+$ and $-$ |
| $3$ | Octopole | $1/r^4$ | Eight charges on a cube, alternating $+$ and $-$ |

### Monopole and dipole terms

The monopole term dominates at large $r$ unless the total charge is zero:

```math
V_{\text{mon}}(\vec{r}) = \frac{Q}{4\pi\varepsilon_0 r}, \qquad Q = \int\rho(\vec{r}\,')\,d^3r'
```

For a neutral object ($Q = 0$), the dipole term dominates. Since $r'\cos\theta' = \vec{r}\,'\cdot\hat{r}$:

```math
V_{\text{dip}}(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\frac{1}{r^2}\int r'\cos\theta'\,\rho(\vec{r}\,')\,d^3r' = \frac{1}{4\pi\varepsilon_0}\frac{\vec{p}\cdot\hat{r}}{r^2}
```

where the **dipole moment** is:

```math
\vec{p} \equiv \int \vec{r}\,'\,\rho(\vec{r}\,')\,d^3r'
```

- For a set of point charges: $\vec{p} = \sum_i q_i\,\vec{r}_i\,'$.
- For $`+q`$ and $`-q`$ separated by $`\vec{d}`$ (pointing from $`-q`$ to $`+q`$): $`\vec{p} = q\vec{r}_+ - q\vec{r}_- = q(\vec{r}_+ - \vec{r}_-) = q\vec{d}`$.
- Dipole moments are vectors and add as vectors.

### Choice of origin

The total charge $Q$ does not depend on where we put the origin. The dipole moment can, because $\vec{r}\,'$ does.

Move the origin by a vector $\vec{a}$. A point that was at $\vec{r}\,'$ is now at $\bar{\vec{r}}\,' = \vec{r}\,' - \vec{a}$. Then:

```math
\bar{\vec{p}} = \int \bar{\vec{r}}\,'\,\rho\,d^3r' = \int(\vec{r}\,'-\vec{a})\,\rho(\vec{r}\,')\,d^3r' = \vec{p} - \vec{a}\,Q
```

So if the total charge is zero, the dipole moment is the same for every choice of origin.

### Electric field of a dipole

Put the dipole at the origin, pointing along $z$ ($\vec{p} = p\,\hat{z}$). Then:

```math
V_{\text{dip}}(r,\theta) = \frac{p\cos\theta}{4\pi\varepsilon_0 r^2}
```

Take $\vec{E} = -\nabla V$ in spherical coordinates:

```math
E_r = -\frac{\partial V}{\partial r} = \frac{2p\cos\theta}{4\pi\varepsilon_0 r^3}, \qquad E_\theta = -\frac{1}{r}\frac{\partial V}{\partial\theta} = \frac{p\sin\theta}{4\pi\varepsilon_0 r^3}, \qquad E_\phi = -\frac{1}{r\sin\theta}\frac{\partial V}{\partial\phi} = 0
```

```math
\vec{E}_{\text{dip}}(r,\theta) = \frac{p}{4\pi\varepsilon_0 r^3}\left(2\cos\theta\,\hat{r} + \sin\theta\,\hat{\theta}\right)
```

**Coordinate-free form:**

```math
\vec{E}_{\text{dip}}(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\frac{1}{r^3}\Big[3(\vec{p}\cdot\hat{r})\,\hat{r} - \vec{p}\Big]
```

*Proof:* write $\hat{z} = \cos\theta\,\hat{r} - \sin\theta\,\hat{\theta}$, so $\vec{p} = p\cos\theta\,\hat{r} - p\sin\theta\,\hat{\theta}$, and $\vec{p}\cdot\hat{r} = p\cos\theta$. Then $3(\vec{p}\cdot\hat{r})\hat{r} - \vec{p} = 3p\cos\theta\,\hat{r} - p\cos\theta\,\hat{r} + p\sin\theta\,\hat{\theta} = p(2\cos\theta\,\hat{r} + \sin\theta\,\hat{\theta})$, which matches the formula above.
