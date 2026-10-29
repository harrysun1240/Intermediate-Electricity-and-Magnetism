# Chapter 4: Electric Fields in Matter

**Contents**

1. [Two kinds of materials](#1-two-kinds-of-materials)
2. [Induced dipoles (stretching)](#2-induced-dipoles-stretching)
3. [Permanent dipoles (rotating)](#3-permanent-dipoles-rotating)
4. [Polarization](#4-polarization)
5. [Field of a polarized object: bound charge](#5-field-of-a-polarized-object-bound-charge)
6. [Field inside a dielectric](#6-field-inside-a-dielectric)
7. [Electric displacement](#7-electric-displacement)
8. [Linear dielectrics](#8-linear-dielectrics)
9. [Example: metal sphere in a dielectric shell](#9-example-metal-sphere-in-a-dielectric-shell)
10. [Example: capacitor filled with a dielectric](#10-example-capacitor-filled-with-a-dielectric)
11. [Boundary value problems with linear dielectrics](#11-boundary-value-problems-with-linear-dielectrics)
12. [Example: dielectric sphere in a uniform field](#12-example-dielectric-sphere-in-a-uniform-field)
13. [Example: point charge above a dielectric half-space](#13-example-point-charge-above-a-dielectric-half-space)
14. [Energy in dielectric systems](#14-energy-in-dielectric-systems)
15. [Forces on dielectrics](#15-forces-on-dielectrics)

---

## 1. Two kinds of materials

There are two broad classes of materials.

- **Conductors** have an essentially unlimited supply of charges that are free to move through the whole material (one or two electrons per atom can roam).
- **Dielectrics** (insulators) have every charge attached to a specific atom or molecule. The charges can only move a little, inside their own atom or molecule.

Even though charges in a dielectric can't travel, an electric field can still distort how the charge is arranged. There are two main ways this happens:

- **Stretching:** the field pulls the positive and negative charge of an atom apart, creating an *induced* dipole moment.
- **Rotating:** molecules that already have a *permanent* dipole moment get turned to line up with the field.

## 2. Induced dipoles (stretching)

An applied field pushes the nucleus one way and the electron cloud the other way. This stretches the atom and gives it a dipole moment. If $`\vec{E}`$ is not too strong, the dipole moment is proportional to the field:

```math
\vec{p} = \alpha\vec{E}
```

The constant $`\alpha`$ is the **atomic polarizability**. Some measured values of $`\alpha/(4\pi\varepsilon_0)`$, in units of $`10^{-30}\ \text{m}^3`$:

| H | He | Li | Be | C | Ne | Na | Ar | K | Cs |
|---|----|----|----|---|----|----|----|---|----|
| 0.667 | 0.205 | 24.3 | 5.6 | 1.76 | 0.396 | 24.1 | 1.64 | 43.4 | 59.6 |

Why the units are a volume: $`\alpha = p/E`$ has units $`(\text{C}\cdot\text{m})/(\text{V/m}) = \text{C}\cdot\text{m}^2/\text{V}`$, and $`4\pi\varepsilon_0`$ has units $`\text{C}/(\text{V}\cdot\text{m})`$. The ratio is $`\text{m}^3`$. Roughly, $`\alpha/(4\pi\varepsilon_0)`$ is the size of the atom cubed.

**Molecules.** A molecule is usually easier to stretch along some directions than others. Then $`\vec{p}`$ need not point along $`\vec{E}`$, and $`\alpha`$ becomes a $`3\times 3`$ matrix, the **polarizability tensor**:

```math
\vec{p} = \alpha\vec{E} \quad\Longrightarrow\quad
\begin{pmatrix} p_x \\ p_y \\ p_z \end{pmatrix} =
\begin{pmatrix} \alpha_{xx} & \alpha_{xy} & \alpha_{xz} \\ \alpha_{yx} & \alpha_{yy} & \alpha_{yz} \\ \alpha_{zx} & \alpha_{zy} & \alpha_{zz} \end{pmatrix}
\begin{pmatrix} E_x \\ E_y \\ E_z \end{pmatrix}
```

## 3. Permanent dipoles (rotating)

Many molecules have a dipole moment even with no field applied. A standard example is water, $`\text{H}_2\text{O}`$:

```math
|\vec{p}| = 6.1\times 10^{-30}\ \text{C}\cdot\text{m}
```

This is large for a molecule. The two hydrogens sit at an angle of about $`105^\circ`$. The electrons cluster near the oxygen, so the oxygen end is negative and the hydrogen end is positive. The dipole moment points from the oxygen toward the hydrogens.

### Force on a charge distribution in a field

Put the origin inside a small charge distribution $`\rho(\vec{r})`$. The total force is:

```math
\vec{F} = \int \rho(\vec{r})\,\vec{E}(\vec{r})\,d^3r
```

The distribution is small, so Taylor expand the field around the origin:

```math
\vec{E}(\vec{r}) = \vec{E}(0) + (\vec{r}\cdot\nabla)\vec{E}(0) + \dots
```

Plug this in. Use $`Q = \int\rho\,d^3r`$ and $`\vec{p} = \int\rho\,\vec{r}\,d^3r`$. Since $`(\vec{r}\cdot\nabla)`$ is linear in $`\vec{r}`$, the integral of $`\rho\,(\vec{r}\cdot\nabla)`$ becomes $`(\vec{p}\cdot\nabla)`$:

```math
\vec{F} = \Big(\int\rho\,d^3r\Big)\vec{E}(0) + \int\rho\,(\vec{r}\cdot\nabla)\vec{E}(0)\,d^3r + \dots = Q\vec{E}(0) + (\vec{p}\cdot\nabla)\vec{E}(0) + \dots
```

The second term is the dipole term. For a neutral molecule ($`Q = 0`$), the force comes from how the field *changes* in space. In a uniform field, a neutral dipole feels no net force.

### Torque on a charge distribution in a field

```math
\vec{N} = \int \vec{r}\times\rho(\vec{r})\vec{E}(\vec{r})\,d^3r = \int\rho(\vec{r})\,\vec{r}\times\big(\vec{E}(0) + \dots\big)\,d^3r = \Big(\int\rho\,\vec{r}\,d^3r\Big)\times\vec{E}(0) + \dots = \vec{p}\times\vec{E}(0) + \dots
```

The torque has size $`pE\sin\theta`$ (with $`\theta`$ the angle between $`\vec{p}`$ and $`\vec{E}`$) and is zero when they are parallel. So the torque tends to rotate $`\vec{p}`$ to line up with $`\vec{E}`$.

## 4. Polarization

Stretching (induced dipoles) and rotating (permanent dipoles) give the same basic result: lots of little dipoles that tend to point along $`\vec{E}`$. For permanent dipoles, the lining up is never complete, because random thermal motion keeps knocking the molecules around.

The material is now **polarized**. We describe this with:

```math
\vec{P}(\vec{r}) \equiv \text{dipole moment per unit volume}
```

$`\vec{P}`$ is called the **polarization**. It is a smoothed-out (macroscopic) quantity. Section 6 explains what that means.

## 5. Field of a polarized object: bound charge

Given $`\vec{P}`$, what field does the object make? First take a point $`\vec{r}`$ **outside** the dielectric.

**Step 1: one point dipole.** A point dipole $`\vec{p}`$ at location $`\vec{r}_p`$ has potential:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\,\frac{\vec{p}\cdot(\vec{r}-\vec{r}_p)}{|\vec{r}-\vec{r}_p|^3}
```

(This is the same as (unit vector from dipole to $`\vec{r}`$) $`\cdot\,\vec{p}`$ divided by distance squared.)

**Step 2: superposition.** Treat the dielectric as a continuous spread of point dipoles. The small volume $`d^3r'`$ at $`\vec{r}\,'`$ holds dipole moment $`\vec{P}(\vec{r}\,')\,d^3r'`$. Add them all up over the volume $`\mathcal{V}`$:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\frac{\vec{P}(\vec{r}\,')\cdot(\vec{r}-\vec{r}\,')}{|\vec{r}-\vec{r}\,'|^3}\,d^3r'
```

**Step 3: a gradient trick.** Taking the gradient with respect to the *primed* coordinates:

```math
\nabla'\left(\frac{1}{|\vec{r}-\vec{r}\,'|}\right) = \frac{\vec{r}-\vec{r}\,'}{|\vec{r}-\vec{r}\,'|^3}
```

(Check one component: $`\partial/\partial x'`$ of $`[(x-x')^2+(y-y')^2+(z-z')^2]^{-1/2}`$ is $`(x-x')/|\vec{r}-\vec{r}\,'|^3`$.) So:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\vec{P}(\vec{r}\,')\cdot\nabla'\left(\frac{1}{|\vec{r}-\vec{r}\,'|}\right)d^3r'
```

**Step 4: integrate by parts.** The product rule $`\nabla'\cdot(f\vec{A}) = f\,(\nabla'\cdot\vec{A}) + \vec{A}\cdot\nabla' f`$, with $`f = 1/|\vec{r}-\vec{r}\,'|`$ and $`\vec{A} = \vec{P}`$, gives:

```math
\vec{P}\cdot\nabla'\left(\frac{1}{|\vec{r}-\vec{r}\,'|}\right) = \nabla'\cdot\left(\frac{\vec{P}}{|\vec{r}-\vec{r}\,'|}\right) - \frac{\nabla'\cdot\vec{P}}{|\vec{r}-\vec{r}\,'|}
```

So:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\nabla'\cdot\left(\frac{\vec{P}(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\right)d^3r' - \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\frac{\nabla'\cdot\vec{P}(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\,d^3r'
```

**Step 5: divergence theorem** on the first integral (surface $`S`$ of the dielectric):

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\oint_S\frac{\vec{P}(\vec{r}\,')\cdot d\vec{a}\,'}{|\vec{r}-\vec{r}\,'|} - \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\frac{\nabla'\cdot\vec{P}(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\,d^3r'
```

**Step 6: read off the bound charges.** The first term looks like the potential of a surface charge, and the second like the potential of a volume charge. Define:

```math
\sigma_b \equiv \vec{P}\cdot\hat{n} \qquad (\hat{n} = \text{outward normal to the surface}), \qquad\qquad \rho_b \equiv -\nabla\cdot\vec{P}
```

Then, with $`d\vec{a}\,' = \hat{n}\,da'`$:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\mathcal{V}}\frac{\rho_b(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\,d^3r' + \frac{1}{4\pi\varepsilon_0}\oint_S\frac{\sigma_b(\vec{r}\,')}{|\vec{r}-\vec{r}\,'|}\,da'
```

- $`\rho_b = -\nabla\cdot\vec{P}`$ is the **bound volume charge density**.
- $`\sigma_b = \vec{P}\cdot\hat{n}`$ is the **bound surface charge density**.

So a polarized object makes the same field as these bound charges. They are real piles of charge, not just a math trick: at the surface the ends of the dipoles are left uncovered, and where $`\vec{P}`$ is not uniform the dipole ends don't cancel inside.

## 6. Field inside a dielectric

What about $`V(\vec{r})`$ for $`\vec{r}`$ **inside** the dielectric?

The result above used the field of a point ("pure") dipole. That formula is right far from a real physical dipole, but not near it. The far fields of a point dipole and a real dipole are the same. The near fields are different (for a real dipole made of two separated charges, the field between the charges points from $`+`$ to $`-`$, opposite to $`\vec{p}`$). Inside a dielectric we are right next to the dipoles, so we need more care.

**The macroscopic field.** The true microscopic field swings wildly from atom to atom, and it isn't useful. We only care about the **macroscopic** field: the field averaged over regions that are large compared to atoms but small compared to the sizes we care about.

### Setup

Define the macroscopic field at $`\vec{r}`$ as the average over a small sphere of radius $`R`$ centered on $`\vec{r}`$. $`R`$ is large compared to atomic scales but small compared to the scales of interest. Split the charges into the ones outside the sphere and the ones inside:

```math
\vec{E}(\vec{r}) = \vec{E}^{\,\text{avg}}_{\text{far}}(\vec{r}) + \vec{E}^{\,\text{avg}}_{\text{near}}(\vec{r})
```

- $`\vec{E}^{\,\text{avg}}_{\text{far}}`$: average field from charges *outside* the sphere (the rest of the dielectric).
- $`\vec{E}^{\,\text{avg}}_{\text{near}}`$: average field from charges *inside* the sphere.

### Two important facts

**(a)** The field averaged over the volume of a sphere, produced by all charges *outside* the sphere, equals the field those charges produce at the *center* of the sphere.

**(b)** The field averaged over the volume of a sphere of radius $`R`$, produced by all charges *inside* the sphere, is:

```math
\vec{E}^{\,\text{avg}} = -\frac{1}{4\pi\varepsilon_0}\frac{\vec{p}}{R^3}
```

where $`\vec{p}`$ is the total dipole moment of the charges inside, using the center of the sphere as the origin.

### Proof of (a) and (b)

Find the average field over the sphere from a single point charge $`q`$. Then use superposition.

**Step 1: set up coordinates.** Put the origin at the center of the sphere, and put $`q`$ on the $`z`$-axis at $`z_q\hat{z}`$. The field at a point $`\vec{r}\,'`$ in the sphere is $`\frac{q}{4\pi\varepsilon_0}\frac{\vec{r}\,'-z_q\hat{z}}{|\vec{r}\,'-z_q\hat{z}|^3}`$. Average it over the sphere's volume $`\frac{4}{3}\pi R^3`$:

```math
\vec{E}_{\text{avg}} = \frac{1}{\frac{4}{3}\pi R^3}\int_0^{2\pi}\!d\phi'\int_0^{\pi}\!d\theta'\sin\theta'\int_0^R\!dr'\,r'^2\;\frac{q}{4\pi\varepsilon_0}\,\frac{\vec{r}\,'-z_q\hat{z}}{|\vec{r}\,'-z_q\hat{z}|^3}
```

A point in the sphere is $`\vec{r}\,' = r'\sin\theta'\cos\phi'\,\hat{x} + r'\sin\theta'\sin\phi'\,\hat{y} + r'\cos\theta'\,\hat{z}`$, and its distance from $`q`$ is:

```math
|\vec{r}\,'-z_q\hat{z}| = \sqrt{r'^2\sin^2\theta' + (r'\cos\theta'-z_q)^2} = \sqrt{z_q^2 + r'^2 - 2z_q r'\cos\theta'}
```

**Step 2: use symmetry.** Everything is symmetric under rotations about the $`z`$-axis, so the $`x`$ and $`y`$ parts average to zero and $`\vec{E}_{\text{avg}}`$ points along $`\hat{z}`$. Only the $`z`$-component $`r'\cos\theta'-z_q`$ survives, and the $`\phi'`$ integral gives $`2\pi`$:

```math
\vec{E}_{\text{avg}} = \frac{\hat{z}}{\frac{4}{3}\pi R^3}\,(2\pi)\int_0^{\pi}\!d\theta'\sin\theta'\int_0^R\!dr'\,r'^2\;\frac{q}{4\pi\varepsilon_0}\,\frac{r'\cos\theta'-z_q}{(z_q^2+r'^2-2z_q r'\cos\theta')^{3/2}}
```

**Step 3: substitute for the angle integral.** Let $`u = z_q^2 + r'^2 - 2z_q r'\cos\theta'`$. Then:

```math
du = 2z_q r'\sin\theta'\,d\theta', \qquad r'\cos\theta' - z_q = \frac{r'^2 - z_q^2 - u}{2z_q}
```

The limits: $`\theta' = 0`$ gives $`u = (z_q - r')^2`$, and $`\theta' = \pi`$ gives $`u = (z_q + r')^2`$. The prefactor is $`\frac{3}{4\pi R^3}\cdot 2\pi = \frac{3}{2R^3}`$, and the substitution brings in $`\frac{1}{2z_q r'}\cdot\frac{1}{2z_q}`$. Together:

```math
\vec{E}_{\text{avg}} = \frac{q\hat{z}}{4\pi\varepsilon_0}\left(\frac{3}{8R^3 z_q^2}\right)\int_0^R\!dr'\,r'\int_{(z_q-r')^2}^{(z_q+r')^2}\!du\;u^{-3/2}\left(r'^2 - z_q^2 - u\right)
```

**Step 4: do the $`u`$ integral.**

```math
\int u^{-3/2}\left(r'^2 - z_q^2 - u\right)du = -2u^{-1/2}\left(r'^2 - z_q^2\right) - 2u^{1/2} = -2\,\frac{r'^2 - z_q^2 + u}{\sqrt{u}}
```

At the upper limit, $`\sqrt{u} = z_q + r'`$:

```math
\frac{r'^2 - z_q^2 + (z_q+r')^2}{z_q + r'} = \frac{2r'^2 + 2z_q r'}{z_q + r'} = 2r'
```

At the lower limit, $`\sqrt{u} = |z_q - r'|`$ (the square root is always positive):

```math
\frac{r'^2 - z_q^2 + (z_q-r')^2}{|z_q - r'|} = \frac{2r'(r' - z_q)}{|z_q - r'|} = \begin{cases} -2r' & z_q > r' \\ +2r' & z_q < r' \end{cases}
```

So the $`u`$ integral is $`-2\,[\,2r' - (\text{lower})\,]`$:

```math
\int_{(z_q-r')^2}^{(z_q+r')^2} = \begin{cases} -2(2r' + 2r') = -8r' & z_q > r' \\ -2(2r' - 2r') = 0 & z_q < r' \end{cases}
```

**Step 5, case (a): $`q`$ outside the sphere** ($`z_q > R`$, so $`z_q > r'`$ for every $`r'`$ in the sphere):

```math
\vec{E}_{\text{avg}} = \frac{q\hat{z}}{4\pi\varepsilon_0}\left(\frac{3}{8R^3 z_q^2}\right)(-8)\int_0^R r'^2\,dr' = \frac{q\hat{z}}{4\pi\varepsilon_0}\left(\frac{-3}{R^3 z_q^2}\right)\left(\frac{R^3}{3}\right) = \frac{q}{4\pi\varepsilon_0}\left(\frac{-\hat{z}}{z_q^2}\right)
```

This is exactly the field of $`q`$ at the center of the sphere.

**Step 6, case (b): $`q`$ inside the sphere** ($`z_q < R`$). Only shells with $`r' < z_q`$ contribute (the rest give $`0`$):

```math
\vec{E}_{\text{avg}} = \frac{q\hat{z}}{4\pi\varepsilon_0}\left(\frac{-3}{R^3 z_q^2}\right)\int_0^{z_q} r'^2\,dr' = \frac{q\hat{z}}{4\pi\varepsilon_0}\left(\frac{-3}{R^3 z_q^2}\right)\left(\frac{z_q^3}{3}\right) = -\frac{1}{4\pi\varepsilon_0}\frac{q z_q\hat{z}}{R^3} = -\frac{1}{4\pi\varepsilon_0}\frac{\vec{p}}{R^3}
```

where $`\vec{p} = q z_q\hat{z}`$ is the dipole moment of $`q`$ about the center.

**Step 7: superposition.** We proved both facts for one point charge. Any charge distribution is a sum of point charges, and both results are linear in $`q`$, so they hold for any distribution.

### Back to the macroscopic field

**Far part.** By fact (a), $`\vec{E}^{\,\text{avg}}_{\text{far}}(\vec{r})`$ is just the field *at* $`\vec{r}`$ from the charges outside the sphere. Those dipoles are far enough from $`\vec{r}`$ that the point-dipole formula is fine:

```math
V^{\,\text{avg}}_{\text{far}}(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\text{outside sphere}}\frac{\vec{P}(\vec{r}\,')\cdot(\vec{r}-\vec{r}\,')}{|\vec{r}-\vec{r}\,'|^3}\,d^3r'
```

**Near part.** The dipoles inside the sphere are too close to $`\vec{r}`$ for the far-field approximation. But by fact (b), their average field is:

```math
\vec{E}^{\,\text{avg}}_{\text{near}}(\vec{r}) = -\frac{1}{4\pi\varepsilon_0}\frac{\vec{p}}{R^3} \qquad (\text{sphere center as origin})
```

The sphere is small enough that $`\vec{P}`$ doesn't change much across it ($`\vec{P}`$ is itself macroscopically averaged). So the total dipole moment inside is $`\vec{p} = \frac{4}{3}\pi R^3\vec{P}`$, and:

```math
\vec{E}^{\,\text{avg}}_{\text{near}} = -\frac{1}{4\pi\varepsilon_0}\cdot\frac{\frac{4}{3}\pi R^3\vec{P}}{R^3} = -\frac{1}{3\varepsilon_0}\vec{P}
```

**Key observation.** This is exactly the field *inside* a sphere of uniform polarization, computed with the point-dipole (bound charge) formula of Section 5.

*Proof.* Let $`\vec{P} = P\hat{z}`$ be uniform. Then $`\rho_b = -\nabla\cdot\vec{P} = 0`$, and $`\sigma_b = \vec{P}\cdot\hat{n} = P\cos\theta`$. We already solved for the potential of a spherical surface charge $`k\cos\theta`$ in Chapter 3 (Example 4 in Section 7). With $`k = P`$:

```math
V(r,\theta) = \begin{cases} \dfrac{P r}{3\varepsilon_0}\cos\theta & r \leq R \\[2ex] \dfrac{P R^3}{3\varepsilon_0 r^2}\cos\theta & r \geq R \end{cases}
```

Inside the sphere $`r\cos\theta = z`$, so $`V = Pz/(3\varepsilon_0)`$ and the field is uniform:

```math
\vec{E} = -\nabla V = -\frac{P}{3\varepsilon_0}\hat{z} = -\frac{1}{3\varepsilon_0}\vec{P}
```

This matches $`\vec{E}^{\,\text{avg}}_{\text{near}}`$. (Outside, the potential is a pure dipole with $`\vec{p} = \frac{4}{3}\pi R^3\vec{P}`$: for a uniformly polarized sphere the total charge is zero, and the quadrupole and higher moments are zero too, so only the dipole term is left.)

Since the bound charge formula came from the point-dipole integral, this means the near part is *also* given by the point-dipole integral, now over the inside of the sphere:

```math
V^{\,\text{avg}}_{\text{near}}(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\text{inside sphere}}\frac{\vec{P}(\vec{r}\,')\cdot(\vec{r}-\vec{r}\,')}{|\vec{r}-\vec{r}\,'|^3}\,d^3r'
```

### Summary of this section

Add the far and near parts. The macroscopic field comes from the potential:

```math
V(\vec{r}) = \frac{1}{4\pi\varepsilon_0}\int_{\text{volume of dielectric}}\frac{\vec{P}(\vec{r}\,')\cdot(\vec{r}-\vec{r}\,')}{|\vec{r}-\vec{r}\,'|^3}\,d^3r'
```

This is the same formula as Section 5, now valid for $`\vec{r}`$ inside the dielectric too. (The integral still converges when $`\vec{r}`$ is inside: the integrand goes like $`1/|\vec{r}-\vec{r}\,'|^2`$, while the volume element near $`\vec{r}`$ goes like $`|\vec{r}-\vec{r}\,'|^2`$.)

- This justifies the earlier calculation, so the bound charges $`\rho_b`$ and $`\sigma_b`$ give the macroscopic field inside the dielectric as well.
- Whatever the actual microscopic charge distribution is, we can replace it with a smooth distribution of point dipoles, as long as we only want the macroscopic field.
- We used a sphere in the argument, but the conclusion does not depend on the shape used.

## 7. Electric displacement

Inside a polarized dielectric there is bound charge: $`\rho_b = -\nabla\cdot\vec{P}`$ inside, and $`\sigma_b = \vec{P}\cdot\hat{n}`$ on the surface. Call every other charge **free charge** $`\rho_f`$. Free charge is anything that is *not* the result of polarizing a dielectric (for example, charge on a conductor, or ions embedded in the material). The total charge density is:

```math
\rho(\vec{r}) = \rho_b(\vec{r}) + \rho_f(\vec{r})
```

**Gauss's law** for the total field $`\vec{E}`$:

```math
\varepsilon_0\nabla\cdot\vec{E} = \rho = \rho_b + \rho_f = -\nabla\cdot\vec{P} + \rho_f
```

Move $`\nabla\cdot\vec{P}`$ to the left:

```math
\nabla\cdot\left(\varepsilon_0\vec{E} + \vec{P}\right) = \rho_f
```

The **electric displacement** is defined by:

```math
\vec{D} \equiv \varepsilon_0\vec{E} + \vec{P}
```

Then Gauss's law becomes:

```math
\nabla\cdot\vec{D} = \rho_f \qquad\text{or, in integral form,}\qquad \oint\vec{D}\cdot d\vec{a} = Q_{f,\text{enc}}
```

This is useful because we usually know (or control) the free charge, but not the bound charge.

**Caution: $`\vec{D}`$ is not just another $`\vec{E}`$.** In electrostatics $`\nabla\times\vec{E} = 0`$, so:

```math
\nabla\times\vec{D} = \varepsilon_0\nabla\times\vec{E} + \nabla\times\vec{P} = \nabla\times\vec{P}
```

This is not necessarily zero. So $`\vec{D}`$ looks similar to $`\vec{E}`$, but there are subtle differences (for example, there is no "Coulomb's law" for $`\vec{D}`$, and it is not always the gradient of a potential).

**How to use it.** When there is enough symmetry (spherical, cylindrical, or planar), find $`\vec{D}`$ from the free charge using Gauss's law. Then find $`\vec{E}`$ and $`\vec{P}`$ from it if needed.

### Boundary conditions

Take a thin pillbox of area $`A`$ and thickness $`h`$ that straddles the boundary, and apply $`\oint\vec{D}\cdot d\vec{a} = Q_{f,\text{enc}}`$ to it. Let $`h \to 0`$ so the box hugs the boundary, and look at each piece:

- Top face (outward normal points up, into the "above" side): flux $`D^\perp_{\text{above}}A`$.
- Bottom face (outward normal points down): flux $`-D^\perp_{\text{below}}A`$.
- Side walls: their area shrinks to zero with $`h`$, so their flux goes to $`0`$.
- Enclosed free charge: any volume charge gives $`\rho_f A h \to 0`$, so only the surface charge is left, $`Q_{f,\text{enc}} = \sigma_f A`$.

So $`D^\perp_{\text{above}}A - D^\perp_{\text{below}}A = \sigma_f A`$. Divide by $`A`$:

```math
D^\perp_{\text{above}} - D^\perp_{\text{below}} = \sigma_f
```

From $`\nabla\times(\vec{D} - \vec{P}) = \varepsilon_0\nabla\times\vec{E} = 0`$: a thin rectangular loop straddling the boundary gives $`\oint(\vec{D}-\vec{P})\cdot d\vec{l} = 0`$, so:

```math
(\vec{D}-\vec{P})^\parallel_{\text{above}}\,l - (\vec{D}-\vec{P})^\parallel_{\text{below}}\,l = 0 \quad\Longrightarrow\quad \vec{D}^\parallel_{\text{above}} - \vec{P}^\parallel_{\text{above}} = \vec{D}^\parallel_{\text{below}} - \vec{P}^\parallel_{\text{below}}
```

(The two short sides of the loop shrink to zero, and $`l`$ is the length of the long sides.) Rearranging:

```math
\vec{D}^\parallel_{\text{above}} - \vec{D}^\parallel_{\text{below}} = \vec{P}^\parallel_{\text{above}} - \vec{P}^\parallel_{\text{below}}
```

## 8. Linear dielectrics

For many materials, as long as $`\vec{E}`$ is not too strong, the polarization is proportional to the field:

```math
\vec{P} = \varepsilon_0\chi_e\vec{E}
```

- $`\chi_e`$ is the **electric susceptibility** of the medium. It is a dimensionless number.
- Materials that satisfy this are called **linear dielectrics**.
- $`\vec{E}`$ here is the *total* field, including the part made by the polarization itself.

In a linear medium:

```math
\vec{D} = \varepsilon_0\vec{E} + \vec{P} = \varepsilon_0(1+\chi_e)\vec{E}
```

Define the **permittivity** of the material $`\varepsilon \equiv \varepsilon_0(1+\chi_e)`$. Then:

```math
\vec{D} = \varepsilon\vec{E}
```

The **relative permittivity**, also called the **dielectric constant**, is:

```math
\varepsilon_r \equiv \frac{\varepsilon}{\varepsilon_0} = 1 + \chi_e
```

| Material | Dielectric constant $`\varepsilon_r`$ |
|----------|------|
| Vacuum   | 1 |
| Helium   | 1.000065 |
| Dry air  | 1.00054 |
| Benzene  | 2.28 |
| Diamond  | 5.7 |
| Silicon  | 11.8 |
| Water    | 80 |
| Ice      | 99 |

**Curl of $`\vec{P}`$ in linear dielectrics.** Inside a *homogeneous* material ($`\chi_e`$ the same everywhere), $`\nabla\times\vec{P} = \varepsilon_0\chi_e\nabla\times\vec{E} = 0`$, so $`\nabla\times\vec{D} = 0`$ too. But at a boundary between different materials, $`\vec{P}`$ jumps, and a loop that crosses the boundary (say between $`\varepsilon_1`$ and $`\varepsilon_2`$) can have $`\oint\vec{P}\cdot d\vec{l} \neq 0`$. So in general you can't find $`\vec{D}`$ from the free charge alone as if it were a vacuum field.

**Crystals.** In a crystal, it is often easier to polarize the material in some directions than others. Then $`\chi_e`$ becomes the **susceptibility tensor**:

```math
\begin{pmatrix} P_x \\ P_y \\ P_z \end{pmatrix} = \varepsilon_0
\begin{pmatrix} \chi_{e,xx} & \chi_{e,xy} & \chi_{e,xz} \\ \chi_{e,yx} & \chi_{e,yy} & \chi_{e,yz} \\ \chi_{e,zx} & \chi_{e,zy} & \chi_{e,zz} \end{pmatrix}
\begin{pmatrix} E_x \\ E_y \\ E_z \end{pmatrix}
```

## 9. Example: metal sphere in a dielectric shell

**Problem.** A metal sphere of radius $`a`$ carries charge $`Q`$. It is surrounded, out to radius $`b`$, by a linear dielectric of permittivity $`\varepsilon`$. Find the potential at the center (relative to infinity), the polarization, and the bound charge.

### Find D with Gauss's law

By spherical symmetry, $`\vec{D}(\vec{r}) = D(r)\,\hat{r}`$. Use a concentric Gaussian sphere of radius $`r`$:

```math
\oint\vec{D}\cdot d\vec{a} = 4\pi r^2 D(r)
```

The free charge enclosed is just the charge on the metal:

```math
Q_{f,\text{enc}} = \begin{cases} 0 & r < a \\ Q & r > a \end{cases}
```

So:

```math
\vec{D}(\vec{r}) = \begin{cases} 0 & r < a \\[1ex] \dfrac{Q}{4\pi r^2}\,\hat{r} & r > a \end{cases}
```

Notice that $`\vec{D}`$ doesn't care whether the dielectric is there or not.

### Find E

$`\vec{P}`$ is hidden in $`\varepsilon`$. In the dielectric, $`\vec{P} = \varepsilon_0\chi_e\vec{E}`$, so $`\vec{D} = \varepsilon_0\vec{E} + \vec{P} = \varepsilon_0(1+\chi_e)\vec{E} = \varepsilon\vec{E}`$. Outside, $`\vec{P} = 0`$, so $`\vec{D} = \varepsilon_0\vec{E}`$. So divide by $`\varepsilon`$ inside the dielectric and by $`\varepsilon_0`$ outside it:

```math
\vec{E}(\vec{r}) = \begin{cases} 0 & r < a \\[1ex] \dfrac{Q}{4\pi\varepsilon r^2}\,\hat{r} & a < r < b \\[2ex] \dfrac{Q}{4\pi\varepsilon_0 r^2}\,\hat{r} & r > b \end{cases}
```

### Potential at the center

```math
V(0) - V(\infty) = -\int_\infty^0\vec{E}\cdot d\vec{l} = \int_0^\infty\vec{E}\cdot d\vec{l} = \int_0^a (0)\,dr + \int_a^b\frac{Q}{4\pi\varepsilon}\frac{dr}{r^2} + \int_b^\infty\frac{Q}{4\pi\varepsilon_0}\frac{dr}{r^2}
```

```math
= \frac{Q}{4\pi\varepsilon}\left[-\frac{1}{r}\right]_a^b + \frac{Q}{4\pi\varepsilon_0}\left[-\frac{1}{r}\right]_b^\infty = \frac{Q}{4\pi\varepsilon}\left(\frac{1}{a} - \frac{1}{b}\right) + \frac{Q}{4\pi\varepsilon_0}\frac{1}{b}
```

```math
V(0) = \frac{Q}{4\pi}\left(\frac{1}{\varepsilon_0 b} + \frac{1}{\varepsilon a} - \frac{1}{\varepsilon b}\right)
```

### Polarization

In the dielectric ($`a < r < b`$):

```math
\vec{P} = \varepsilon_0\chi_e\vec{E} = \varepsilon_0\chi_e\,\frac{Q}{4\pi\varepsilon r^2}\,\hat{r}
```

Simplify the constant. From $`\varepsilon/\varepsilon_0 = 1 + \chi_e`$, we get $`\chi_e = \varepsilon/\varepsilon_0 - 1`$, so:

```math
\frac{\varepsilon_0}{\varepsilon}\chi_e = \frac{\varepsilon_0}{\varepsilon}\left(\frac{\varepsilon}{\varepsilon_0} - 1\right) = 1 - \frac{\varepsilon_0}{\varepsilon} = \frac{\varepsilon - \varepsilon_0}{\varepsilon}
```

```math
\vec{P} = \left(1 - \frac{\varepsilon_0}{\varepsilon}\right)\frac{Q\,\hat{r}}{4\pi r^2}
```

### Bound charge

Volume bound charge (divergence in spherical coordinates):

```math
\rho_b = -\nabla\cdot\vec{P} = -\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2 P_r\right) = -\frac{1}{r^2}\frac{\partial}{\partial r}\left[\frac{Q}{4\pi}\left(1 - \frac{\varepsilon_0}{\varepsilon}\right)\right] = 0
```

Surface bound charge. On the outer surface the outward normal is $`\hat{n} = \hat{r}`$. On the inner surface the outward normal (out of the dielectric) is $`\hat{n} = -\hat{r}`$:

```math
\sigma_b(r = b) = \vec{P}\cdot\hat{r} = \left(1 - \frac{\varepsilon_0}{\varepsilon}\right)\frac{Q}{4\pi b^2}, \qquad \sigma_b(r = a) = -\vec{P}\cdot\hat{r} = -\left(1 - \frac{\varepsilon_0}{\varepsilon}\right)\frac{Q}{4\pi a^2}
```

Since $`\varepsilon > \varepsilon_0`$, the factor $`1 - \varepsilon_0/\varepsilon`$ is positive. So the signs make sense: negative bound charge collects next to the positive metal sphere, and positive bound charge at the outer edge.

**Check.** The total bound charge is zero ($`4\pi a^2\sigma_b(a) + 4\pi b^2\sigma_b(b) = 0`$). And for $`a < r < b`$, the total enclosed charge is $`Q + 4\pi a^2\sigma_b(a) = Q\,\varepsilon_0/\varepsilon`$. Then ordinary Gauss's law gives $`E = Q/(4\pi\varepsilon r^2)`$, matching the result above. The inner bound charge partly screens $`Q`$.

## 10. Example: capacitor filled with a dielectric

**Problem.** A parallel-plate capacitor (plate area $`A`$, separation $`d`$) is filled with a material of dielectric constant $`\varepsilon_r`$. What does this do to the capacitance?

Put $`+Q`$ on the lower plate and $`-Q`$ on the upper plate, so the free surface charge is $`\sigma = Q/A`$. Gauss's law for $`\vec{D}`$ (a pillbox around the lower plate) gives $`\vec{D} = \sigma\hat{z}`$ between the plates, just like in vacuum. Then $`\vec{E} = \vec{D}/\varepsilon`$:

```math
\vec{E} = \frac{\sigma}{\varepsilon_0\varepsilon_r}\,\hat{z} = \frac{Q}{A\varepsilon_r\varepsilon_0}\,\hat{z} \qquad (\text{pointing from the lower plate to the upper plate})
```

The potential difference is $`V = |\vec{E}|\,d`$, so:

```math
C = \left|\frac{Q}{V}\right| = \frac{Q}{Qd/(A\varepsilon_r\varepsilon_0)} = \frac{A}{d}\,\varepsilon_r\varepsilon_0
```

```math
C = \varepsilon_r C_{\text{vac}}
```

The capacitance goes up by a factor of $`\varepsilon_r`$ compared to vacuum.

## 11. Boundary value problems with linear dielectrics

**Bound charge follows free charge.** In a homogeneous linear dielectric, $`\vec{P} = \varepsilon_0\chi_e\vec{E} = (\varepsilon_0\chi_e/\varepsilon)\vec{D}`$, so:

```math
\rho_b = -\nabla\cdot\vec{P} = -\nabla\cdot\left(\frac{\varepsilon_0\chi_e}{\varepsilon}\vec{D}\right) = -\frac{\varepsilon_0\chi_e}{\varepsilon}\,\rho_f = -\left(\frac{\chi_e}{1+\chi_e}\right)\rho_f
```

In other words, $`\rho_b \propto \rho_f`$.

- If no free charge is embedded in a linear dielectric, then $`\rho_b = 0`$. The only bound charge is the surface charge $`\sigma_b`$.
- Then $`\rho = 0`$ inside the dielectric, so the potential there obeys **Laplace's equation**. We can use the separable expansions from Chapter 3.

**Boundary conditions in terms of free charge.** Start from $`D^\perp_{\text{above}} - D^\perp_{\text{below}} = \sigma_f`$ and use $`D = \varepsilon E`$ and $`E^\perp = -\partial V/\partial n`$ (with $`n`$ the coordinate normal to the surface, pointing from below to above):

```math
\varepsilon_{\text{above}}E^\perp_{\text{above}} - \varepsilon_{\text{below}}E^\perp_{\text{below}} = \sigma_f
\quad\Longrightarrow\quad
\varepsilon_{\text{above}}\frac{\partial V_{\text{above}}}{\partial n} - \varepsilon_{\text{below}}\frac{\partial V_{\text{below}}}{\partial n} = -\sigma_f
```

The potential is also continuous across the boundary:

```math
V_{\text{above}} = V_{\text{below}}
```

## 12. Example: dielectric sphere in a uniform field

**Problem.** A sphere of radius $`R`$, made of a homogeneous linear dielectric, is placed in a uniform applied field $`\vec{E}_0 = E_0\hat{z}`$. Find the potential everywhere.

### General solution

The problem has azimuthal symmetry, and there is no free charge, so $`V`$ satisfies Laplace's equation inside and outside. Using the Chapter 3 solutions:

```math
V_{\text{in}}(r,\theta) = \sum_{l=0}^{\infty}A_l r^l P_l(\cos\theta)
```

```math
V_{\text{out}}(r,\theta) = \sum_{l=0}^{\infty}\frac{B_l}{r^{l+1}}P_l(\cos\theta) - E_0 r\cos\theta
```

Inside we drop the $`1/r^{l+1}`$ terms (they blow up at $`r = 0`$). Outside, the general solution also has growing $`C_l r^l P_l(\cos\theta)`$ terms. Far from the sphere, $`V`$ must go back to the potential of the applied field, $`V = -E_0 z = -E_0 r\cos\theta`$ (condition (c) below). So the only growing term allowed is $`l = 1`$ with $`C_1 = -E_0`$, which gives the $`-E_0 r\cos\theta`$ term. All other $`C_l`$ are zero.

### Boundary conditions

- **(a)** $`V_{\text{in}} = V_{\text{out}}`$ at $`r = R`$ (continuity).
- **(b)** $`\varepsilon\,\dfrac{\partial V_{\text{in}}}{\partial r} = \varepsilon_0\,\dfrac{\partial V_{\text{out}}}{\partial r}`$ at $`r = R`$ (no free surface charge, so $`D^\perp`$ is continuous).
- **(c)** $`V_{\text{out}} \to -E_0 r\cos\theta`$ as $`r \to \infty`$. Far away, $`V \to -E_0 z + \text{constant} = -E_0 r\cos\theta + \text{constant}`$, and we set the constant to zero.

### Condition (a): continuity

```math
\sum_{l=0}^{\infty}A_l R^l P_l(\cos\theta) = -E_0 R\cos\theta + \sum_{l=0}^{\infty}\frac{B_l}{R^{l+1}}P_l(\cos\theta)
```

The Legendre polynomials are independent, so the coefficient of each one must match. Remember $`\cos\theta = P_1(\cos\theta)`$:

```math
A_l R^l = \frac{B_l}{R^{l+1}} \quad (l \neq 1), \qquad A_1 R = -E_0 R + \frac{B_1}{R^2}
```

### Condition (b): the normal derivative

```math
\sum_{l=0}^{\infty}\varepsilon\, l A_l R^{l-1}P_l(\cos\theta) = \varepsilon_0\left(-E_0\cos\theta - \sum_{l=0}^{\infty}(l+1)\frac{B_l}{R^{l+2}}P_l(\cos\theta)\right)
```

Divide by $`\varepsilon_0`$ and use $`\varepsilon_r = \varepsilon/\varepsilon_0`$. Matching coefficients:

```math
\varepsilon_r\, l A_l R^{l-1} = -\frac{(l+1)B_l}{R^{l+2}} \quad (l \neq 1), \qquad \varepsilon_r A_1 = -E_0 - \frac{2B_1}{R^3}
```

### Solve

**All $`l \neq 1`$ terms vanish.** For $`l = 0`$, condition (b) says $`0 = -B_0/R^2`$, so $`B_0 = 0`$, and then (a) gives $`A_0 = 0`$. For $`l \geq 2`$, (a) and (b) give two expressions for $`A_l R^l`$:

```math
A_l R^l = \frac{B_l}{R^{l+1}}, \qquad A_l R^l = -\frac{(l+1)B_l}{\varepsilon_r\, l\, R^{l+1}}
```

Setting them equal: $`B_l = -\frac{l+1}{\varepsilon_r l}B_l`$. The factor $`-\frac{l+1}{\varepsilon_r l}`$ is negative, so it can't equal $`1`$. The only option is $`B_l = 0`$, and then $`A_l = 0`$.

**The $`l = 1`$ terms.** From (a) and (b):

```math
A_1 = -E_0 + \frac{B_1}{R^3} = \frac{1}{\varepsilon_r}\left(-E_0 - \frac{2B_1}{R^3}\right)
```

Multiply through by $`\varepsilon_r R^3`$ and collect the $`B_1`$ terms:

```math
\varepsilon_r B_1 + 2B_1 = (\varepsilon_r - 1)E_0 R^3 \quad\Longrightarrow\quad B_1 = \left(\frac{\varepsilon_r - 1}{\varepsilon_r + 2}\right)E_0 R^3
```

```math
A_1 = -E_0 + \frac{\varepsilon_r - 1}{\varepsilon_r + 2}E_0 = \frac{-(\varepsilon_r + 2) + (\varepsilon_r - 1)}{\varepsilon_r + 2}E_0 = -\left(\frac{3}{\varepsilon_r + 2}\right)E_0
```

### Result

```math
V_{\text{in}}(r,\theta) = -\frac{3E_0}{\varepsilon_r + 2}\,r\cos\theta, \qquad V_{\text{out}}(r,\theta) = -E_0 r\cos\theta + \frac{\varepsilon_r - 1}{\varepsilon_r + 2}\,\frac{E_0 R^3}{r^2}\cos\theta
```

**The field inside the dielectric is uniform.** With $`r\cos\theta = z`$:

```math
\vec{E}_{\text{in}} = -\nabla V_{\text{in}} = -\nabla\left(-\frac{3E_0 z}{\varepsilon_r + 2}\right) = \frac{3E_0}{\varepsilon_r + 2}\,\hat{z} = \left(\frac{3}{\varepsilon_r + 2}\right)\vec{E}_0
```

Since $`\varepsilon_r > 1`$, the field inside is weaker than $`\vec{E}_0`$. Outside, the field lines bend toward the sphere. The extra outside term is the potential of a dipole sitting at the center.

**Check.** As $`\varepsilon_r \to \infty`$, $`B_1 \to E_0 R^3`$ and $`\vec{E}_{\text{in}} \to 0`$. This is the uncharged metal sphere in a uniform field from Chapter 3 (Example 3 in Section 7).

## 13. Example: point charge above a dielectric half-space

**Problem.** The entire region below the plane $`z = 0`$ is filled with a uniform linear dielectric of susceptibility $`\chi_e`$. Find the force on a point charge $`q`$ a distance $`d`$ above the origin.

### Bound surface charge

There is no free charge inside the dielectric ($`\rho_f = 0`$), so $`\rho_b = 0`$ (Section 11). The only bound charge is on the surface. The outward normal from the dielectric is $`\hat{n} = \hat{z}`$, so:

```math
\sigma_b = \vec{P}\cdot\hat{n} = P_z = \varepsilon_0\chi_e E_z
```

where $`E_z`$ is the *total* field just inside the dielectric (just below $`z = 0`$). It has two parts.

**Part 1: from $`q`$.** The field of the point charge is:

```math
\vec{E}_q = \frac{q}{4\pi\varepsilon_0}\,\frac{x\hat{x} + y\hat{y} + (z-d)\hat{z}}{\left(x^2 + y^2 + (z-d)^2\right)^{3/2}}
```

At $`z = 0`$, its $`z`$-component is:

```math
-\frac{q d}{4\pi\varepsilon_0\left(x^2+y^2+d^2\right)^{3/2}}
```

**Part 2: from $`\sigma_b`$ itself.** Right next to the sheet, a surface charge acts like an infinite plane, with field $`\sigma_b/(2\varepsilon_0)`$ pointing away from it. Just below the sheet, "away" is $`-\hat{z}`$, so this part is $`-\sigma_b/(2\varepsilon_0)`$. (The rest of the sheet only adds field parallel to the plane at this point.)

**Solve for $`\sigma_b`$.**

```math
\sigma_b = \varepsilon_0\chi_e\left(-\frac{q d}{4\pi\varepsilon_0\left(x^2+y^2+d^2\right)^{3/2}} - \frac{\sigma_b}{2\varepsilon_0}\right)
```

Move the $`\sigma_b`$ term to the left: $`\sigma_b\left(1 + \frac{\chi_e}{2}\right) = -\frac{\chi_e\, q d}{4\pi\left(x^2+y^2+d^2\right)^{3/2}}`$. So:

```math
\sigma_b = -\frac{1}{2\pi}\left(\frac{\chi_e}{\chi_e + 2}\right)\frac{q d}{\left(x^2+y^2+d^2\right)^{3/2}}
```

Now that we know $`q`$ and $`\sigma_b`$, we could get $`\vec{E}`$ directly from Coulomb's law by integrating. That is not easy.

### Method of images

**Total bound charge.** Use $`s^2 = x^2 + y^2`$ and integrate over the plane:

```math
q_b = 2\pi\int_0^\infty s\,ds\;\sigma_b(s) = -\left(\frac{\chi_e}{\chi_e+2}\right)q d\int_0^\infty\frac{s\,ds}{\left(s^2+d^2\right)^{3/2}}
```

Let $`v = s^2 + d^2`$, so $`dv = 2s\,ds`$:

```math
q_b = -\left(\frac{\chi_e}{\chi_e+2}\right)\frac{q d}{2}\int_{d^2}^\infty v^{-3/2}\,dv = -\left(\frac{\chi_e}{\chi_e+2}\right)\frac{q d}{2}\Big[-2v^{-1/2}\Big]_{d^2}^\infty = -\left(\frac{\chi_e}{\chi_e+2}\right)\frac{q d}{2}\cdot\frac{2}{d}
```

```math
q_b = -\left(\frac{\chi_e}{\chi_e+2}\right)q
```

**Guess the image setup.**

- For the region $`z > 0`$: the real charge $`q`$ at $`(0,0,d)`$ plus an image charge $`q_b`$ at $`(0,0,-d)`$.
- For the region $`z < 0`$: a single charge $`q + q_b`$ at $`(0,0,d)`$.

```math
V = \begin{cases} \dfrac{1}{4\pi\varepsilon_0}\left(\dfrac{q}{\sqrt{x^2+y^2+(z-d)^2}} + \dfrac{q_b}{\sqrt{x^2+y^2+(z+d)^2}}\right) & z > 0 \\[3ex] \dfrac{1}{4\pi\varepsilon_0}\,\dfrac{q + q_b}{\sqrt{x^2+y^2+(z-d)^2}} & z < 0 \end{cases}
```

**Check that it works.** This function:

- satisfies Poisson's equation with a point charge $`q`$ at $`(0,0,d)`$ in $`z > 0`$, and Laplace's equation in $`z < 0`$ (the charges it uses there are outside that region);
- goes to zero at infinity;
- is continuous at $`z = 0`$ (both pieces equal $`(q+q_b)/(4\pi\varepsilon_0\sqrt{s^2+d^2})`$ there);
- has a jump in its normal derivative at $`z = 0`$ that matches $`\sigma_b`$.

By uniqueness, it is the answer.

**Another way to find $`q_b`$.** Assume the setup above with unknown image charges: $`q`$ and $`q'`$ at $`(0,0,-d)`$ for $`z > 0`$, and $`q''`$ at $`(0,0,d)`$ for $`z < 0`$. Continuity of $`V`$ at $`z = 0`$ needs $`q + q' = q''`$. Next use continuity of $`D^\perp`$ (no free surface charge). At $`z = 0`$:

```math
\left.\frac{\partial V_{\text{below}}}{\partial z}\right|_{z=0} = \frac{(q + q')\,d}{4\pi\varepsilon_0\left(s^2+d^2\right)^{3/2}}, \qquad \left.\frac{\partial V_{\text{above}}}{\partial z}\right|_{z=0} = \frac{(q - q')\,d}{4\pi\varepsilon_0\left(s^2+d^2\right)^{3/2}}
```

$`D^\perp`$ continuous means $`\varepsilon\,\partial V_{\text{below}}/\partial z = \varepsilon_0\,\partial V_{\text{above}}/\partial z`$, so $`\varepsilon_r(q + q') = q - q'`$:

```math
q' = \left(\frac{1 - \varepsilon_r}{1 + \varepsilon_r}\right)q = \left(\frac{-\chi_e}{2 + \chi_e}\right)q = q_b
```

This matches the total bound charge found above.

### Force on q

The force on $`q`$ is the force from the image charge:

```math
\vec{F} = \frac{q\,q_b}{4\pi\varepsilon_0(2d)^2}\,\hat{z} = -\frac{1}{4\pi\varepsilon_0}\left(\frac{\chi_e}{\chi_e+2}\right)\frac{q^2}{4d^2}\,\hat{z}
```

The force points toward the dielectric (attractive). As $`\chi_e \to \infty`$ this becomes $`-\frac{1}{4\pi\varepsilon_0}\frac{q^2}{4d^2}\hat{z}`$, the result for a grounded conducting plane.

## 14. Energy in dielectric systems

**Idea.** Bring in the free charge a little at a time, letting the dielectric respond as we go.

**Step 1: work for one small step.** The work to bring in an extra free charge $`\Delta\rho_f`$ is:

```math
\Delta W = \int(\Delta\rho_f)\,V\,d^3r
```

where $`V`$ is the potential from all the free (and bound) charge already in place.

**Step 2: write it with D.** Since $`\nabla\cdot\vec{D} = \rho_f`$, adding $`\Delta\rho_f`$ changes $`\vec{D}`$ by $`\Delta\vec{D}`$ with:

```math
\Delta\rho_f = \nabla\cdot(\Delta\vec{D}) \quad\Longrightarrow\quad \Delta W = \int\big(\nabla\cdot(\Delta\vec{D})\big)\,V\,d^3r
```

**Step 3: product rule.** Use $`\nabla\cdot(\Delta\vec{D}\,V) = (\nabla\cdot\Delta\vec{D})\,V + \Delta\vec{D}\cdot\nabla V`$, and $`\nabla V = -\vec{E}`$:

```math
\Delta W = \int\nabla\cdot(\Delta\vec{D}\,V)\,d^3r + \int\Delta\vec{D}\cdot\vec{E}\,d^3r
```

**Step 4: drop the surface term.** The divergence theorem turns the first integral into a surface integral. If we integrate over all space, the surface is at infinity and the integral vanishes. So:

```math
\Delta W = \int\Delta\vec{D}\cdot\vec{E}\,d^3r \qquad (\text{any material})
```

**Step 5: linear dielectrics.** If $`\vec{D} = \varepsilon\vec{E}`$ with $`\varepsilon`$ fixed, then $`\Delta\vec{D} = \varepsilon\,\Delta\vec{E}`$, and:

```math
\Delta(\vec{D}\cdot\vec{E}) = \Delta\vec{D}\cdot\vec{E} + \vec{D}\cdot\Delta\vec{E} = \Delta\vec{D}\cdot\vec{E} + \varepsilon\vec{E}\cdot\Delta\vec{E} = \Delta\vec{D}\cdot\vec{E} + \vec{E}\cdot\Delta\vec{D} = 2\,\Delta\vec{D}\cdot\vec{E}
```

So $`\Delta\vec{D}\cdot\vec{E} = \frac{1}{2}\Delta(\vec{D}\cdot\vec{E})`$, and:

```math
\Delta W = \Delta\left(\frac{1}{2}\int\vec{D}\cdot\vec{E}\,d^3r\right)
```

**Step 6: add up all the steps.** Starting from nothing, the total work is:

```math
W = \frac{1}{2}\int\vec{D}\cdot\vec{E}\,d^3r
```

### What this energy includes

$`W`$ is the work to bring in and assemble the free charges *while letting the dielectric respond*. It includes:

- **(a)** the electrostatic energy of the free charges,
- **(b)** the electrostatic energy of the bound charges,
- **(c)** the work done stretching and rotating the dielectric's atoms and molecules.

Parts (b) and (c) are equal and opposite, so they cancel. To see this, split the energy of all the charges, $`\frac{\varepsilon_0}{2}\int E^2\,d^3r = \frac{1}{2}\int\rho V\,d^3r`$, into a free part and a bound part:

```math
\frac{\varepsilon_0}{2}\int E^2\,d^3r = \underbrace{\frac{1}{2}\int\rho_f V\,d^3r}_{\text{free share}} + \underbrace{\frac{1}{2}\int\rho_b V\,d^3r}_{\text{bound share}}
```

- The bound share is $`\frac{1}{2}\int(-\nabla\cdot\vec{P})V\,d^3r = \frac{1}{2}\int\vec{P}\cdot\nabla V\,d^3r = -\frac{1}{2}\int\vec{P}\cdot\vec{E}\,d^3r`$ (integrating by parts, surface term at infinity dropped).
- The work stored in the stretched molecules is $`+\frac{1}{2}\int\vec{P}\cdot\vec{E}\,d^3r`$. Each linearly responding molecule stores $`\frac{1}{2}\vec{p}\cdot\vec{E}`$, like a stretched spring: $`\int_0^p (p'/\alpha)\,dp' = p^2/(2\alpha) = \frac{1}{2}pE`$.

These two cancel, leaving $`W = \frac{1}{2}\int\rho_f V\,d^3r`$, which equals $`\frac{1}{2}\int\vec{D}\cdot\vec{E}\,d^3r`$ (integrate $`-\vec{D}\cdot\nabla V`$ by parts and use $`\nabla\cdot\vec{D} = \rho_f`$).

**Note:** $`W = \frac{1}{2}\int\vec{D}\cdot\vec{E}\,d^3r`$ is **not** the work to assemble both free and bound charges as plain charges. That would be $`\frac{\varepsilon_0}{2}\int E^2\,d^3r`$, which ignores (c).

## 15. Forces on dielectrics

Dielectrics are pulled into an electric field, just like conductors. The force often comes from the **fringing field** at the edge of a capacitor. For example, a dielectric slab partly inside a parallel-plate capacitor gets pulled in. Computing this directly from the fringing field is hard, so we use energy instead.

**Setup.** Plates of length $`l`$, width $`w`$, separation $`d`$. The slab has been pulled out a distance $`x`$, so a length $`l - x`$ is still inside. Let $`W`$ be the energy of the capacitor.

**Force from energy.** Pull the slab out by $`dx`$. The energy changes by the work we do, $`dW = F_{\text{app}}\,dx`$, where $`F_{\text{app}}`$ is the applied force that balances the electrical force. So:

```math
F_{\text{el}} = -F_{\text{app}} = -\frac{dW}{dx}
```

**Use fixed charge.** With the capacitor isolated, $`Q`$ stays fixed and $`W = Q^2/(2C)`$:

```math
F_{\text{el}} = -\frac{d}{dx}\left(\frac{Q^2}{2C}\right) = \frac{Q^2}{2C^2}\frac{dC}{dx} = \frac{1}{2}V^2\frac{dC}{dx}
```

**Capacitance.** The two parts act like capacitors in parallel: the filled part has length $`l - x`$ and the empty part has length $`x`$:

```math
C = \frac{\varepsilon_0 w}{d}\Big[\varepsilon_r(l - x) + x\Big] = \frac{\varepsilon_0 w}{d}\left(\varepsilon_r l - \chi_e x\right) \qquad (\varepsilon_r = 1 + \chi_e)
```

So $`dC/dx = -\chi_e\varepsilon_0 w/d`$, and:

```math
F_{\text{el}} = \frac{1}{2}V^2\left(-\frac{\chi_e\varepsilon_0 w}{d}\right) = -\frac{\varepsilon_0\chi_e w}{2d}V^2
```

The minus sign means the force points toward smaller $`x`$: it pulls the slab back into the capacitor.

(If the plates were held at fixed $`V`$ by a battery instead, you would also have to count the work done by the battery. Done correctly, you get the same force.)
