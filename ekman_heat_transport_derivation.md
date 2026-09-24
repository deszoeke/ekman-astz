# Ekman transport of heat in the ocean and atmosphere

## Summary

Heat transports by balanced Ekman transports are calculated. These apply to balanced Ekman flows on scales greater than the Rossby radius. The Ekman (mass) transport is equal and opposite in the atmosphere and ocean, regardless of density. 

The atmospheric moist static energy gradients, by including latent heat of water vapor, are expected to be of similar amplitude to those in the ocean, despite its ~4x greater specific enthalpy. 

The atmosphere may experience larger temporary temperature variations (~1 K) than the ocean, but persistent gradients in the ocean are probably stronger owing to its lower velocity and smaller Rossby radius of deformation. 


## Turbulent fluxes and bulk formulae

Turbulent momentum flux (vector) and its scalar magnitude:
$$
\overline{w'\bold{u}'} \;\left[=u_*^{\,2}\right] = C_D U \bold{u}
$$

Turbulent heat flux (scalar):
$$
\overline{w'T'} = u_* T_* = C_H U (T_0 - T_a)\,
$$

Enthalpy flux:
$$
\rho c_{pa} \overline{w'T'} = \rho c_{pa} u_* T_* = \rho c_{pa} C_H U (T_0 - T_a)\,
$$

Scalar wind speed:
$$
U = \lVert \bold{u} \rVert = (\bold{u}\cdot\bold{u})^{1/2}
$$

Wind stress vector:
$$
\bold{\tau} = \rho\,\overline{w'\bold{u}'} = C_D U \bold{u}
$$

## Ekman mass transport

The vertically integrated Ekman transport $\bold{V}_E$ in the ocean, due to the downward stress $\bold{\tau}$ from the atmosphere to the ocean is
$$
\bold{V}_E = -\hat{\bold{k}}\times\bold{\tau}/f
$$

units
$
\left[ \bold{V}_E = \text{kg}\;\text{m}\,\text{s}^{-1}\,\text{m}^{-2} = \text{kg}\,\text{m}^{-1}\,\text{s}^{-1} \right]
$

The atmospheric Ekman transport $-\bold{V}_E$ is equal and opposite the ocean Ekman transport.

The Ekman mass convergence is the curl of the wind stress,
$
-\nabla\cdot\bold{V}_E = -\nabla\times\left(\bold{\tau}/f\right).
$

The curl of the Ekman transport is (minus) the divergence of the wind stress,
$
\nabla\times\bold{V}_E = -\nabla\cdot\left(\bold{\tau}/f\right).
$

The advection of a scalar $H$ by the Ekman transport is
$
-\bold{V}_E\cdot\nabla H.
$

**Dimensional check:**
$$
[\tau] = \text{kg}\,\text{m}^2\,\text{s}^{-2}\,\text{m}^{-3} \;(=\text{kg}\,\text{m}^{-1}\text{s}^{-2})
$$
$$
[\tau/f] = \text{kg}\,\text{m}^2\,\text{s}^{-1}\,\text{m}^{-3} \;(=\text{kg}\,\text{m}^{-1}\text{s}^{-1})
$$
$$
[\nabla\times(\tau/f)] = \text{kg}\,\text{m}\,\text{s}^{-1}\,\text{m}^{-3} \;(=\text{kg}\,\text{m}^{-2}\text{s}^{-1})
$$

## Ekman transport of a scalar

The divergence of the Ekman flux of $H$ is evaluated by the divergence chain rule,
$$
\nabla\cdot\bold{E} = \nabla\cdot(H\bold{V}_E)
= H(\nabla\cdot\bold{V}_E) + \bold{V}_E \cdot \nabla H.
$$

Then convergence of transport is achieved by Ekman pumping and Ekman advection of horizontal gradients of scalar $H$,
$$
-\nabla\cdot\bold{E} = -H\nabla\cdot\bold{V}_E - \bold{V}_E \cdot \nabla H
$$
$$
= \underbrace{-H\,\nabla\times
(\bold{\tau}/f)
}_{\text{Ekman pumping}} \;\underbrace{-\bold{V}_E \cdot \nabla H }_{\text{Ek. adv. of }H}.
$$

More vector derivatives yield the curl of the $H$ transport,
$$
\nabla \times (H\bold{V}_E)
= H \nabla \times \bold{V}_E 
 +\nabla H \times \bold{V}_E
$$
Use $\nabla \times \bold{V}_E = -\nabla \cdot \bold{\tau}/f$ for the first term. For the second,
$$
\nabla H \times \bold{V}_E = \nabla H \times (-\hat{\bold{k}}\, \times \bold{\tau}/f)
$$
$$
= -\bold{\tau}/f \cdot \nabla H.
$$

Now the curl of the enthalpy transport,
$$
\nabla \times \bold{E} = \nabla \times (H\bold{V}_E)
= -H (\nabla \cdot \bold{\tau}/f) - \bold{\tau}/f \cdot \nabla H,
$$
is the sum of $H$ times the wind stress divergence and advection of $H$ by the stress.

$H$ could be any Ekman-layer averaged scalar. An innovation is to sum the ocean and atmospheric Ekman layer scalar transports. The sum of the ocean and atmospheric Ekman mass transport is zero, but it transports the ocean-atmosphere difference of a scalar. We next take $H$ to be the difference in specific enthalpy between the ocean and atmosphere.

## Horizontal enthalpy transports

Specific enthalpies of the atmosphere and ocean Ekman layers are
$$
h_a = c_{pa}T_a, \qquad h_o = c_{po}T_o.
$$

The horizontal enthalpy transports by the Ekman transport in the atmosphere and ocean are,
$$
\bold{E}_a = -\bold{V}_E\,T_a c_{pa} = -\bold{V}_E\,h_a
$$
$$
\bold{E}_o = \bold{V}_E\,T_o c_{po} = \bold{V}_E\,h_o.
$$

The sum over the ocean and atmosphere Ekman layers is
$$
\bold{E} = \bold{V}_E\,(h_o - h_a) = \bold{V}_E\,H,
$$
with the ocean-atmosphere enthalpy difference $H \equiv h_o - h_a$.

**Note:** $H = c_{po}T_o - c_{pa}T_a$ 
- is **not** vertically integrated (unlike mixed-layer MSE).
- is alike, but **not proportional** to $T_o - T_a$ in the surface heat flux.

Different specific heats
$$
c_{po} = 3900\ \text{J}\,\text{kg}^{-1}\text{K}^{-1}
$$
$$
c_{pa} = 1000\ \text{J}\,\text{kg}^{-1}\text{K}^{-1}
$$
make $H$ different than the ocean-air temperature difference. The enthalpy difference can be decomposed into its main contribution from the ocean, and a smaller contribution from the ocean-atmosphere temperature difference:
$$
H = c_o T_o - c_a T_a = c_o \big[\underbrace{\left(1-a\right)T_o}_{\text{mainly }(3/4)\text{ ocean}} + \underbrace{a(T_o-T_a)}_{\substack{\text{small contrib. from}\\ \Delta T \propto \text{sfc. flux}}}\big],$$
with $a=c_{pa}/c_{po}\approx 1/3.9$.

## Moist static energy
The 3.9 times larger specific enthalpy of the ocean makes ocean Ekman temperature transport dominate that in the atmosphere. Yet adding latent heat of vaporization of atmospheric specific humidity anomalies, the atmosphere contributes nearly much to the Ekman layer moist static energy as the ocean temperature anomalies. The latent heat of vaporization of water is
$
L = 2.4\times10^6 \text{ J kg}^{-1}.
$
For a typical tropical specific humidity $q=15\times 10^{-3}$ and $d\ln q_s/dT = 0.07 \text{ K}^{-1}$, the gradient $q$ scales with the temperature as
$$
\Delta q/\Delta T = q (d\ln q_s/d T)
    \approx 1 \text{ g kg}^{-1} \text{ K}^{-1}.
$$
Then
$$
\Delta m_a = c_{pa} \Delta T + L \Delta q 
= [c_{pa} + L(\Delta T/\Delta q)] \Delta T
\approx 3.4c_{pa} \Delta T
$$
The effective specific enthalpy $3.4 c_{pa}$ is nearly as large as $c_{po} = 3.9 c_{pa}$.
If the atmosphere and ocean Ekman layers have the similar temperature gradients, then adding the gradient of latent heating of vapor in the atmospheric layer approximately makes up for the larger specific enthalpy of liquid water.

