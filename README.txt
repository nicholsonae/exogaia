To run the code it requires 4 integer arguments:
1. determines the insulating / reflecting properties of the chemical set
2. wires the geochemistry
3. initialises microbe metabolisms 
4. data file number

For the data in the paper the numbers corresponding to each chemical set:
- chemical set A - 120
- chemical set B - 320
- chemical set C - 520
- chemical set D - 920
- chemical set E - 720

e.g. example run:

./exogaia 120 1234564 1234560 0 &

The parameter link_probability determines the connectivity of the geochemical network (between 0 and 1)

The most important data file it spits out is the exogaia_macro_data one which contains: 

timestep | total population | population of the largest species | number of species | planet temperature | total insulating effect of atmosphere (%) | total reflecting effect of atmosphere (%)
