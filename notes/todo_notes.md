# TODO Notes


## June 18

Changed Python code to load data files that have `snapshot_dt != dt`. But still have to check that it works




- Check that can animate data files that have `snapshot_dt != dt` and that fps of animations work as expected.



## June 23

- Fix extending of simulations to combine SimulationData objects to efficiently read and do correlation functions on extended sims.

- Run longer sims on larger boxes and do MSD and polar order on them



## July 14

Convert old serialized files to new serialization method.-- DONE


## July 23
Fix plotting Orientation Correlation data from .txt files to work with files that aren't sweeps (i.e. for data files that are just 1 line).



## August 2
Am seeing that average reversal time varies logarithmically with `interactoinfliprate` for `alignsimple`. The fliprate collapses to around ~2.0 at low and 0 interaction fliprate. I am running simulations with `I=[0.5, 1.5, 2.5, 3, 3.5, 4, 4.5]` to see if there is a sudden phase transition from collapsed to varying. Maybe it itsn't log it just is an illusion because of the values close to 2.