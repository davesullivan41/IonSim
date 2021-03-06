IonSim
======

Python based simulation of an ion engine accelerating through the solar system. This program was partially motivated by an attempt to prove the technical viability of the business idea in [this elevator pitch](https://www.youtube.com/watch?v=m1IhtxvTHQ8), delivered by myself at the 2014 Rice Undergraduate Elevator Pitch Competition.

======

The purpose of this project was to demonstrate the ability of tiny ion engines to accelerate small satellites outside the solar system within a reasonable time frame. The results can be seen most clearly seen in the following plots:
![Jupiter fly-by](Results/1.7kg_4Engines_Assist3/figure_1.png)
This shows the path of the satellite as it exits the solar system. In this case, 4 engines are used to accelerate a 1kg payload. The blue dot highlights a gravity assist from Jupiter.
![Energy plot](Results/1.7kg_4Engines_Assist3/figure_5.png)
Here you can see the kinetic and potential energies of the satellite as it accelerates. The kinetic energy surpasses the potential energy around the Jupiter fly-by which demonstrates that the satellite has reached escape velocity. 

The program works by pre-building a model of the solar system using an Euler-Cromer method, which allows for more flexibility and accuracy than a simple Kepler's laws derivation. Then the satellite is inserted into this system near geosync and acceleration begins. The entire trajectory of the satellite is calculated, saved to a data file, and then plots are produced. Data processing for the example case shown above took approximately 20 minutes on a mid-range desktop computer. 

A careful observer will notice that no planets further away than Jupiter are considered in this model. This approximation was made in order to decrease the complexity of the simulation and because the satellite is clearly able to escape the solar system without influence from any of the planets beyond Jupiter. I ran tests that included Saturn and its effect on the trajectory of the satellite was only significant when the satellite came within close proximity to the planet. 

