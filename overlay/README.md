# Overlay for 5point-PYNQ

The included accelerator is taken from the 
[project](https://bitbucket.org/necst/xohw18_5points_public) by the 
*NECSTLab* at *Politecnico di Milano* that won the Xilinx European 
[Open Hardware Design Contest](http://www.openhw.eu) in 2018, for the 
*AWS EC2 F1 Category*.

## Supported Boards/Shells

Currently, we distribute the overlay only for the following Alveo boards and 
shells:

Shell                    | Board             
-------------------------|-----------------
xilinx_u200_xdma_201830_2|Xilinx Alveo U200
xilinx_u250_xdma_201830_2|Xilinx Alveo U250
xilinx_u280_xdma_201920_1|Xilinx Alveo U280
xilinx_u50_xdma_201920_1|Xilinx Alveo U50

Designs are built using Vitis 2019.2.

## Rebuilding Overlay

To build xclbins for a new shell, move to the `overlay` folder and run the 
`make build` command, passing the appropriate `DEVICE`:
   ```bash
   cd overlay
   make build DEVICE=<target-shell>
   ```

You can then do `make install` to copy the overlays in the appropriate 
notebooks folder. 
If you want to do everything in one go (`build` and `install`):
   ```bash
   make DEVICE=<target-shell>
   ```

To install built xclbins in a specific path (different from the default one):
   ```bash
   make install INSTALL_PATH=<target-path>
   ```

It is expected that the target path contains the appropriate notebooks folders.
Otherwise, the install step will not be carried out.

## Use Built Overlay

Once you have synthesized your overlay, you can use them by first getting the
notebooks, disabling device-based overlays resolution (Refer to `pynq
get-notebooks --help` for info on the used options.) 

```bash
pynq get-notebooks 5point --ignore-overlays --path <target-path>
```

And then copy the overlays in the target directory, using the `make install` 
command

```bash
make install INSTALL_PATH=<target-path>/pynq-notebooks/5point
```
