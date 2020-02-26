<img src="fivepoint_pynq/notebooks/server/img/5-point_relpose.png" width="288" height="216" />

# 5-point Relative Pose Problem for PYNQ

This project provides overlay and notebooks to run a websocket-based **5-point 
relative pose** application. It requires [PYNQ](http://www.pynq.io/) version 
`2.5.1` and above to work.

The purpose of this application is to identify the possible relative camera 
motions given five matching points from two calibrated views.

This is a porting of the 
[project](https://bitbucket.org/necst/xohw18_5points_public) by the *NECSTLab* 
at *Politecnico di Milano* that won the Xilinx European 
[Open Hardware Design Contest](http://www.openhw.eu) in 2018, for the *AWS EC2 
F1 Category*.

More info on the project can be found at the following links:

 - [Project website](https://bitbucket.org/necst/xohw18_5points_public)
 - [PDF report](https://bitbucket.org/necst/xohw18_5points_public/src/master/report/report.pdf)
 - [YouTube presentation](https://www.youtube.com/watch?v=UDGWGNdglFs)
 - [Open Hardware 2018 winners & finalists](http://www.openhw.eu/2018-finalists.html)

Refer to the 
[README](https://github.com/Xilinx/5point-PYNQ/tree/master/overlay/README.md) 
in the `overlay` folder for more information regarding the used overlay and how 
it is created.

Please notice that the distributed overlay might not be available for all 
target devices. Supported devices are listed in the overlays 
[README](https://github.com/Xilinx/5point-PYNQ/tree/master/overlay/README.md). 
There you may also find instructions on how to synthesize and use overlays for 
a different device.

## Quick Start

Install the `fivepoint-pynq` package using `pip`:
   ```bash
   pip install fivepoint-pynq
   ```

After the package is installed, to get your own copy of all the notebooks 
available run:
   ```bash
   pynq get-notebooks 5point
   ```

You can try things out by running:
   ```bash
   cd pynq-notebooks
   jupyter notebook
   ```

There are a number of additional options for the `pynq get-notebooks` command,
you can list them by typing 
   ```bash
   pynq get-notebooks --help
   ```

You can also refer to the official 
[PYNQ documentation](https://pynq.readthedocs.io/) for more information 
regarding the *PYNQ Command Line Interface* and in particular the 
`get-notebooks` command.

## Licenses

**5point-PYNQ**: [Apache License 2.0](https://github.com/giunatale/Alveo-PYNQ/blob/master/LICENSE)

**Xilinx Open Hardware 2018 - 5 Points to Rule Them All**: [MIT License](https://bitbucket.org/necst/xohw18_5points_public/src/master/LICENSE.md) (3d-party component)
