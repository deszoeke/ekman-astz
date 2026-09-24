# Ekman flux in the air-sea transition layer

## 1. Setting: two boundary layers sharing one stress

The air-sea transition layer couples an atmospheric boundary layer (ABL) above
the interface to an oceanic mixed layer (OML) below it. Both layers are turbulent,
rotating, and forced by the same interfacial stress $\boldsymbol{\tau}_0$ (the wind
stress on the ocean, equal and opposite to the drag the ocean exerts on the
atmosphere). Away from the thin molecular sublayers at the interface itself, the
horizontal momentum balance in each fluid is a balance between the Coriolis force,
the pressure gradient (geostrophic) force, and the divergence of turbulent
(Reynolds) stress:
$$
f \hat{z} \times (\mathbf{u} - \mathbf{u}_g) = \frac{1}{\rho}\frac{\partial \boldsymbol{\tau}}{\partial z},
$$
where $\mathbf{u}_g$ is the geostrophic wind (atmosphere) or the geostrophic
current (ocean), $f$ is the Coriolis parameter, $\rho$ the fluid density, and
$\boldsymbol{\tau}(z)$ the turbulent stress profile carrying the surface stress
$\boldsymbol{\tau}_0$ into the fluid interior. This is the Ekman balance, and the
layer over which it is dynamically active on each side of the interface is
the Ekman layer.

## 2. Eddy-viscosity closure and the Ekman spiral

Closing the stress with a (possibly depth-varying) eddy viscosity,
$\boldsymbol{\tau} = \rho K \partial \mathbf{u}/\partial z$, and taking $K$
constant gives the classical Ekman equation
$$
f \hat{z} \times (\mathbf{u} - \mathbf{u}_g) = K \frac{\partial^2 \mathbf{u}}{\partial z^2}.
$$
With a no-slip/matching condition at the surface and boundedness far from it,
the solution is the Ekman spiral: the ageostrophic wind (or current) rotates
with depth and decays over the Ekman depth
$$
\delta_E = \sqrt{\frac{2K}{|f|}},
$$
turning $45^\circ$ from the surface stress direction at the interface itself
and further with depth, vanishing at $z \sim \delta_E$. The ocean and
atmosphere spirals turn in the same rotational sense (with $f$) but the ocean's
is far shallower ($\delta_E \sim 10$–$50$ m for typical trade-wind mixing) than
the atmosphere's ($\delta_E \sim 500$–$1000$ m for the ABL), because $K_{ocean}
\ll K_{atmos}$.

## 3. Ekman transport (the depth-integrated flux)

The dynamically important quantity for the transition layer is usually not the
detailed spiral but the depth-integrated ageostrophic transport, obtained by
integrating the momentum balance from the surface to a depth below which the
turbulent stress vanishes:
$$
\mathbf{M}_E = \int (\mathbf{u} - \mathbf{u}_g)\, dz
= \frac{1}{\rho f} \hat{z} \times \boldsymbol{\tau}_0 .
$$
This is exact regardless of the eddy-viscosity profile — it depends only on
the surface stress and $f$, not on $K(z)$. In the ocean, $\mathbf{M}_E$ is the
Ekman transport, directed $90^\circ$ to the right of the wind stress in the
Northern Hemisphere (left in the Southern Hemisphere). In the atmosphere the
same relation gives the depth-integrated cross-isobar mass flux in the
boundary layer, directed $90^\circ$ to the left of the surface stress (i.e.
toward low pressure), since the sign of $f\hat z \times$ flips the roles of
stress and transport between the two fluids sitting on either side of the
interface.

## 4. Surface stress closure

The coupling flux itself is set by the bulk aerodynamic formula
$$
\boldsymbol{\tau}_0 = \rho_{air} C_D |\mathbf{U}_{10} - \mathbf{u}_s|\,(\mathbf{U}_{10} - \mathbf{u}_s),
$$
with $\mathbf{U}_{10}$ the 10 m wind, $\mathbf{u}_s$ the surface current, and
$C_D$ the drag coefficient (itself dependent on wind speed and stability). This
$\boldsymbol{\tau}_0$ is the single number that both Ekman layers respond to:
it sets the oceanic Ekman transport $\mathbf{M}_E = \hat z \times
\boldsymbol{\tau}_0/(\rho_0 f)$ and, through the ABL momentum balance, the
cross-isobaric component of the low-level wind that ventilates the shallow
cumulus layer above the trade-wind ABL.

## 5. Ekman pumping and its relevance to the transition layer

Where the wind stress curl varies horizontally, the divergence of the Ekman
transport drives a vertical velocity at the base of the Ekman layer,
$$
w_E = \hat{z}\cdot \nabla \times \left(\frac{\boldsymbol{\tau}_0}{\rho_0 f}\right)
= \frac{1}{\rho_0}\nabla \times \left(\frac{\boldsymbol{\tau}_0}{f}\right)\cdot \hat z .
$$
In the ocean this is Ekman pumping/suction: convergent Ekman transport
($w_E<0$ pumping down) deepens the mixed layer and suppresses entrainment of
cold thermocline water; divergent transport ($w_E>0$, suction) shoals the
mixed layer and can sharpen SST gradients. Because SST anomalies feed back on
$C_D$ and on boundary-layer stability, spatial structure in the wind-stress
curl (e.g. across SST fronts, or downstream of orographic wakes) couples
directly to mesoscale SST and shallow-cloud patterns through this pumping
term — the same mechanism that appears in the mesoscale SST–wind coupling
discussed elsewhere in this project (see `mesoscale-vs-desiccation.docx`).

## 6. Caveats for the trade-wind transition layer

- **Non-constant $K$.** Both the trade-wind ABL and the ocean mixed layer have
  strongly depth-varying and often convective (not purely shear-driven)
  turbulence, so the classical spiral shape is rarely observed; the transport
  relation $\mathbf{M}_E = \hat z \times \boldsymbol{\tau}_0/(\rho f)$ is more
  robust than the spiral itself because it does not depend on the closure.
- **Finite-depth/mixed-layer effects.** When the turbulent layer is shallow
  and well mixed (common in the oceanic case, and in the convective trade-wind
  ABL under shallow cumulus), the "slab" limit of the Ekman balance is more
  appropriate than the infinite-depth spiral: the whole mixed layer moves
  quasi-uniformly and $\mathbf{M}_E \approx h(\mathbf{u}-\mathbf{u}_g)$ with
  $h$ the mixed-layer depth.
- **Two-way coupling.** Unlike the textbook problem where $\boldsymbol{\tau}_0$
  is prescribed, in the coupled air-sea system the ocean's Ekman-driven SST
  response feeds back on stability and $C_D$, which changes $\boldsymbol{\tau}_0$
  and the atmospheric Ekman transport in turn — relevant wherever this model
  couples surface flux and boundary-layer schemes.

## 7. Summary relations

| Quantity | Relation |
|---|---|
| Ekman balance | $f\hat z\times(\mathbf u-\mathbf u_g) = \rho^{-1}\partial_z\boldsymbol\tau$ |
| Ekman depth | $\delta_E = \sqrt{2K/\vert f\vert}$ |
| Ekman transport | $\mathbf M_E = \hat z \times \boldsymbol\tau_0 /(\rho f)$ |
| Surface stress | $\boldsymbol\tau_0 = \rho_{air} C_D\vert\Delta\mathbf U\vert\Delta\mathbf U$ |
| Ekman pumping | $w_E = \rho_0^{-1}\,\hat z\cdot\nabla\times(\boldsymbol\tau_0/f)$ |
