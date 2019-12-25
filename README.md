# Arduino_MPC_Library
This is a compact Model Predictive Control (MPC) library for Teensy 4.0 system (or Arduino system in general).

For demonstration purpose, I made the MPC control a state-space model (2 input, 2 output, 4 state) of Jet Transport Aircraft (ref: https://www.mathworks.com/help/control/ug/mimo-state-space-models.html#buv3tp8-1). The result for Hp = 7 and Hu = 4 is:

![Alt text](Result.png "Result.png")

(Computation time: 187 us)
