# Ekman transport of heat in the ocean and atmosphere

Simon de Szoeke

(copyright 2026, MIT License)

## Summary

Heat transports by balanced Ekman transports are calculated. These apply to balanced Ekman flows on scales greater than the Rossby radius. The Ekman (mass) transport is equal and opposite in the atmosphere and ocean, regardless of density. 

The atmospheric moist static energy gradients, by including latent heat of water vapor, are expected to be of similar amplitude to those in the ocean, despite its ~4x greater specific enthalpy. 

The atmosphere may experience larger temporary temperature variations (~1 K) than the ocean, but persistent gradients in the ocean are probably stronger owing to its lower velocity and smaller Rossby radius of deformation. 

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

