# Simulator of the Attitude of a Satellite with Magnetorquers and Reaction Wheels

Created by **Loup BESSARD**

---

## Objective
This simulator was created for the nanosatellite **NIMPH** (*Nanosatellite to Investigate Microwave Photonics Hardware*).  
NIMPH has a specific configuration with **four actuators**: one reaction wheel and three magnetorquers.  

The main objectives of this project are:
1. To **model** the satellite’s attitude (orientation) using differential equations and simulate its evolution over time.  
2. To **test and compare different control laws** designed to stabilize or guide the satellite.

---

## Assumptions
- The satellite is **cubic**, **rigid**, and has a **uniform mass distribution**. It is possible to define its dimensions, mass, and inertia tensor.  
- **External perturbations** (friction, solar radiation, etc.) are neglected.  
- The **orbital motion** follows **Keplerian assumptions**.  
- **Actuators** are aligned with the satellite’s body axes.  

---

## Model Capabilities

This simulator allows you to:
- Define the **inertia matrix** of the satellite.  
- Configure between **0 and 3 magnetorquers** and **0 and 3 reaction wheels**.  
- Set the **inertia**, **torque saturation**, and **velocity saturation** of each reaction wheel.  
- Define the **magnetic saturation** limits of magnetorquers.  
- Specify **orbital parameters** (assuming Keplerian motion).  
- Choose between several **reference frames** for magnetic field computation: **NED**, **ECEF**, **ECI**, or **LVLH**.  

---

## Methods and Techniques

The model is based on **differential equations** developed in the accompanying documentation.  
To solve these equations, the simulator uses the **Runge–Kutta 4th-order (RK4)** integration method.  

The **Earth’s magnetic field** is modeled using the **IGRF** (International Geomagnetic Reference Field) model, obtained from the following [repository](https://github.com/wb-bgs/m_IGRF).

Once the dynamic model is established, users can implement and test various **attitude control laws** on top of it.

---

## Control and Command Features

You can define your own **control law** to stabilize or orient the satellite.  
Three examples are provided to illustrate different approaches:

- **Simulation 1**: Uses an **anonymous function** based on **quaternions**.  
- **Simulation 2**: Uses a control function defined in a separate file (`command_exemple2.m`) and expressed in **Cardan angles**. It includes **six additional adaptive states**, whose derivatives are defined in `derivativeAdditionState_exemple2.m`.  
- **Simulation 3**: Extends Simulation 2 by adding a **projection function** (`projection_exemple2.m`), also linked to the adaptive control structure.

Your control function can work with either **quaternions** or **Euler/Cardan angles**.  
It is also possible to add **new dynamic states** (e.g., adaptive gains or variable inertia) and include **external parameters** via `varargin`.  
Note that including extra parameters generally **increases computation time**.

---

## File Structure

| File | Test | Description |
|:------|:------|:-------------|
| `satellite.m` | `satelliteTest.m` | Contains the `satellite` class used to define the spacecraft (inertia, actuators, etc.). It is also the main entry point for running simulations and setting initial conditions. |
| `orbit.m` | `orbitTest.m` | Contains the `orbit` class, which computes the satellite’s position in geocentric coordinates. In this project, the orbit is mainly used to determine the Earth’s magnetic field at each position. |
| `reactionWheel.m` | `reactionWheelTest.m` | Defines the `reactionWheel` class representing a single reaction wheel, including its initial speed, torque limits, and automatic saturation handling. |
| `magnetorquer.m` | `magnetorquerTest.m` | Defines the `magnetorquer` class representing one magnetic actuator with its saturation behavior. |
| `planet.m` | `planetTest.m` | Defines the `planet` class, a database used by `orbit.m` to simulate orbits around different celestial bodies. |
| `tools.m` | `toolsTest.m` | Defines the `tools` class, providing utility functions for coordinate conversions, attitude representation transformations, and plotting. |

The directory `earth_magnetic_field/` is obtained from the external IGRF model repository linked above.

To run the tests, make sure you are in the correct directory, then use:
```runtests('file_nameTest')```

To run all tests, use:
```runtests```

![image Info](https://github.com/BLoup81/Simulator-satellite-with-magnetorqer-and-reaction-wheels/blob/main/diagramme_de_classe.png)

---

## Example Files

Two example scripts are provided to help you get started:

- **`exemple.m`**: Demonstrates how to define a custom satellite configuration using the provided classes.  
- **`exemple2.m`**: Demonstrates how to run simulations and implement control laws. It includes three different setups:
  1. 3 reaction wheels + 3 magnetorquers, controlled via a quaternion-based anonymous function.  
  2. Same configuration, but with an adaptive control law defined in an external function.  
  3. Same as (2), with an additional projection function for enhanced control behavior.
  
---

## Documentation
A detailed report is provided in **Rapport_stage.pdf** (in French), which contains:
- The complete satellite model.  
- The actuator models.  
- The orbital equations.  
- Several control laws and their mathematical derivations.  
- A bibliography for further reading.  