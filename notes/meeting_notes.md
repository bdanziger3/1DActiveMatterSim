# Dissertaion Meeting Notes: Active Matter

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



# 6 June

- To Do
    - Run Larger Simulations and Systematically vary parameters
    - Look at orientation correlation function

- Could try putting a y offset for particle animations to see when many are on top of each other

- Good idea to keep thermal flip rate always at 1
- Vary the interaction strength
- Look at rate of global polarity shift (all pointing left to all pointing right). See how this varies with interaction strength

- Can try running a few big simulations overnight (10^4 particles, )

- Options:
    1. Density (rho = N/L)
    1. Flipping Rate (R)
    1. System Size (L)

    Keep two constant and vary 1 of these and see what happens to the global reversal time (average time between flips of global polarity)

Good idea for constant parameters:

dt: = 10^-4
total_time: 1,000
N: 10,000

interaction flipping rate: Run from ~0-100

See if I have to worry about fluctuations around 0 (if 1-2 particles flipping causing global spin to change)


to analyze:
- Mean Orientation Time --> Reversal Time
- Orientation correlation function ( should look like decaying exponential) $\langle u_i (t + \Delta t) u_t(t)\rangle_{t,i}$ (averaging over all particles and time, but wait until the system has relaxed until a steady state)
    - Fit with $e^{-t/\tau}$. $\tau(\rho, \Gamma, L)$
- Mean Squared Displacement



- Plot the mean velocity over time and see if it stays at 0 or if we see the flips






## Potential thing to look at
- Analytical solution of 3-particle case with alignment interactions: Probability of entire system to flip

- Usually can look at 1, 2, 3, infinite particle systems


## Make a file writing down the info from the papers and what I want to test from them so that I have it one place.



## Conceptual Questions:
Active matter community usually using assumption that particles are in an overdamped system (like a thick fluid) where the passive particles will want slow down to $v=0$. But active particles are accelerating instantaneously to $v_0$




# 23 June

- Would be interesting to plot one particle in a simulation and see if it is staying in a band or flipping

- Look at average spin (polar order with sign)

- Try making boxwidth bigger to see how that effects simulations when spacial correlation is less

- Also look at mean squared displacement

- Longer simulations

- also account for box edge for nearest neighbor search



## Paper

- Write a general literature overview and send to Till next week. Talk about what people are doing in the field and what people are looking for. Talk about what system I am trying to look at

- Till can meet between 10:20 and 13:00 next week

# 27 June

![Notes from meeting](./Images/june27_meetingnotes.jpg)




# 7 July
Make histogram with bigger bin size
and then do it with weighting the densities by number of particles


Correlation with one particle with itself at different times

Correlation with one particle with its neighbors at the same time

Do these on longer data with yes/no interaction


Make a presentation with plots and summaries of what I did each week



$\log(n)$




# 11 July


- Sweep interactoin strength
    - e.g. keep `boxwidth=100`, `Numparticles=1000` and sweep `interactionfliprate`


- Orientation self-correlation plots:
    - interesting that they decay and don't stay in their flock. Look at animation to see if 1 particle goes back and forth or stays

    - only need to plot until timescale of ~1-10, after that is just noise

    - order lines in correct order

- Scale parameters semi-logarithmically to get a good look at different magnitudes
    - e.g. [1, 2, 5, 10, 20, 50, 100, 200, 500]



- Animation with 1 particle highlighted its own color to see if it's bounding back and forth or staying in the flock


## July 23

- Look at 1D run-and-tumble particles


- Density plot -- histogram of densities by position with colors based on polarization


- Eventually look at turning away interactions





- Introduce active matter
- talk about brownian motion
- talk about 2D cases (flocking, Vicsek)
- then move to 1D cases


- 1-2 pages on what the correlation functions are
- ~1 page on how I implemented my code


- Start writing introduction
- Make an outline of methods


- Write email to Patrick and Till to discuss meeting for next Friday



# 30 July

- Run turnaway interactions with lower interaction and long time scales. Try plotting without y axis.
    - Try to find transition between steady state and moving state,
    - find functions that might characterize this


- run alignsimple sims at longer timescales (~1000s with N=1000 and N=10000 depending on how long it takes to run)
    - Redo MSD and OC on these timescales, can say settled at t=50s
    - see if number of bands stays steady throughout sim
    - also investigate mean polarization over all particles over time.


- look at MSD of alignsimple and make a log-log plot


- write down analytic definition of interactions



# Week of 6 Aug

### Questions:

- Citing software packages?? How do I do that
- What style guide should I use? AIP?
- Should I use `I` or `we` in my paper? I had help from you guys but should it me `me` because it is a Masters' thesis?

- What is the measure of uncertainty? Vicsek mentions that the `"uncertainty is within the width of the data points"` but what does that mean with particle simulation and average velocity? Is he saying the uncertainty over many runs?


# Week of 17 Aug

### Questions:

- Error bars: should I use standard error or 95% confidence (1.96 * \sigma)

- Capitalisation of words like 'figure', 'equation', 'section', or 'chapter' when referring to other parts of my dissertation


- Add diagrams for
    - showing Vicsek model
    - Explaining expected motion of brownian vs run and tumble particles
    - Explaining what systems with high vs low polar order look like
    - explaining ballistic vs diffusive behavior

- Other notes
    - dont need v(t) equation it is trivial
    - make descriptions less redundant and make them in order of complexity
    - add labels to plot to show $\rho$ or $\nu$ so that you don't have to read the caption
    - Talk less about exponents
        - $\nu$ sweep can't even argue that it's polynomial
        - focus on qualitative results

    - Only do MSD of turn-away
    - talk about qualitative resutls
    - show that becomes subdiffusive at high $\nu$ (Slope < 1)


    - MSD tells you the time where ballistic becomes diffusive. Also tells you about long range motion. i.e. the coefficient of the long range motion after becomes diffusive.



# 19 August Meeting
### Questions:
    - ask whether plot subcaption looks better on bottom or in corner
    


# 19 August Meeting
    - talk about caging
    -add intro sentences to each section and subsection to introduce what I am going to go into


    - talk about Equation 32 in Zöttl and Stark paper with MSD diffusion times

    - look at turnaround time equation and make sure it is right
        - might need to take absolute value of F(t)