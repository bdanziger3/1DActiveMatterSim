# Dissertaion Meeting Notes: Active Solids

## Blake Denziger
### Patrick Piezonka and Till Welker



# February 28 Meeting

- Density plots: position density and polar order (pick a radius and calculate average spin)
- Eventually do it as an animation
- Can plot polarity as color on a 1d hist with height as position density
- Implement orientation correlation and mean squared density



# 16 May

- For no interactions: normalize with `fliprate = v0 = 1`. Boxwidth is the tunable param
- Can keep track of only wrapped positions and can then work out the unwrapped ones by checking for larger than normal jumps between frames


Next steps:
- random initial conditions (maybe offset by non-multiple of dx to get them off lattice)
- change `interactionfliprate` to work like a rate and not a probability

- seperate simulation and physical parameter

- plot the density polarization over time
- look at larger particle numbers and larger box sizes