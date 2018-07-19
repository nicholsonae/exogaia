# ExoGaia
This is the ExoGaia code, and details could be seen from the following paper:

Nicholson, A. E., Wilkinson, D. M., Williams, H. T., & Lenton, T. M. (2018). Gaian bottlenecks and planetary habitability maintained by evolving model biospheres: The ExoGaia model. *Monthly Notices of the Royal Astronomical Society*, *477*(1), 727-740. [https://doi.org/10.1093/mnras/sty658](https://doi.org/10.1093/mnras/sty658).

### Getting Started

To begin you'll need a copy of the source code. Either clone directly from the GitHub repository or download it as a zip file. To compile the code on linux machine, you need  `-std=c++11` or `-std=gnu++11` flag to support for the ISO C++ 2011 standard. On Mac, you can compile it with `g++` directly.

```{bash}
$ git clone https://github.com/nicholsonae/exogaia.git
$ cd exogaia
$ g++ exogaia.cpp -std=c++11 -o exogaia 
```

### Usage 
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
