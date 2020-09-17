# Tutorial

## Basics

... To be written

## Drawing Coplanar Waveguides

... To be written

Meandering resonators can be drawn via the CPW_Meander_Resonator function. Note that the class is a valid CPW object and can thus, be strung together with other CPW objects or linked via the function CPW2CPW. The meanders are configured such that the inner and outer path lengths are equal. The number of meanders depends on the specified turn-radius and maximal width parameters.

## Drawing Capacitors

Basic capacitors can be drawn via the classes Capacitor_Interdigitated (interdigitated capacitors) and Capacitor_GapCPW (a basic gap in a CPW line). Both classes have a cpw_params attribute that can be passed onto CPW2CPW. Thus, the capacitors can act like CPW objects and be strung in series with other CPW lines. To see the basic configurations of interdigitated capacitors that can be drawn, consider the following unit-tests (to be run on configuring a layer_photo in the usual Klayout code):

``` python
#Testing the configurations
cap = Capacitor_Interdigitated(DPoint(600e3,1240e3),
                                 3.3e3, #Finger width
                                 3.3e3, #Capacitor Finger Space
                                 100e3, #Finger length
                                 3.3e3, #Capacitor gap
                                 10e3, #Padding width
                                 5e3, #Side gap
                                 5, 'same_N_right',
                                 trans_in = klayout.db.Trans(0, False, 0, 0))
cap.place( cell, layer_photo )
cap = Capacitor_Interdigitated(DPoint(800e3,1240e3),
                                 3.3e3, #Finger width
                                 3.3e3, #Capacitor Finger Space
                                 100e3, #Finger length
                                 10e3, #Capacitor gap
                                 20e3, #Padding width
                                 5e3, #Side gap
                                 5, 'same_N_left',
                                 trans_in = klayout.db.Trans(0, False, 0, 0))
cap.place( cell, layer_photo )
cap = Capacitor_Interdigitated(DPoint(1000e3,1240e3),
                                 3.3e3, #Finger width
                                 10e3, #Capacitor Finger Space
                                 100e3, #Finger length
                                 10e3, #Capacitor gap
                                 10e3, #Padding width
                                 15e3, #Side gap
                                 5, 'diff_N_start',
                                 trans_in = klayout.db.Trans(0, False, 0, 0))
cap.place( cell, layer_photo )
cap = Capacitor_Interdigitated(DPoint(1200e3,1240e3),
                                 10e3, #Finger width
                                 3.3e3, #Capacitor Finger Space
                                 100e3, #Finger length
                                 20e3, #Capacitor gap
                                 10e3, #Padding width
                                 5e3, #Side gap
                                 5, 'diff_N_end',
                                 trans_in = klayout.db.Trans(0, False, 0, 0))
cap.place( cell, layer_photo )
```

## COMSOL

Klayout designs can be ported into fresh COMSOL simulations after installing COMSOL and MPh (https://mph.readthedocs.io/en/latest/installation.html):

```
pip install mph
```

The MPh library takes care of the Java interface to talk to COMSOL. It could perhaps be bypassed as most of the model-building functionality just uses the JPype library to run the Java raw commands (see COMSOL.py). Here is some boilerplate code from the tutorial ExampleCOMSOL.ipynb:

``` python
from ClassLib.COMSOL import COMSOL_Model
import matplotlib.pyplot as plt

chip_len = CHIP.dx
chip_wid = CHIP.dy
chip_thickness = 0.5e-3

cmdl = COMSOL_Model('supercond',24) #24 threads
cmdl.initialize_model(chip_len, chip_wid, chip_thickness)

#Add conductors (layout is a Klayout object)
cmdl.add_metallic_Klayout(layout, layer_photo)
#Add the structures to define the port from CPW to ground
cmdl.create_port_on_CPW(pad_L,True)
cmdl.create_port_on_CPW(pad_R,False)
#Define the CPW to be a fine-structure by passing a point within the polygon
cmdl.register_fine_structure(res_mean.start*1e-9 + DPoint(1e-6,0), 1.2e-6, 10e-6)

cmdl.build_geom_mater_elec_mesh()
cmdl.set_freq_range(7.210e9,7.214e9,30)
cmdl.run_simulation()
cmdl.save("Test_1.mph")

#Results will be stored within the mph file; but one may choose to access them directly here:
s_param_vals = cmdl.get_sparams()
plt.plot(s_param_vals[0]/1e9, s_param_vals[1],'bo-',label='S11')
plt.plot(s_param_vals[0]/1e9, s_param_vals[2],'ro-',label='S21')
plt.grid(True, which="both")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('f(GHz)')
plt.ylabel('S-parameters')
plt.show()
```

The class COMSOL_Model wraps functionality given in MPh to create a silicon chip with metallic structures on its top surface. Currently the code just supports a simple s-parameter frequency sweep (although the mph file will store the E-field data for inspection via the COMSOL GUI).
