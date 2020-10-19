import mph
import jpype.types as jtypes
import klayout.db as db
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib as mpl


class COMSOL_Model:
    def __init__(self, model_name, num_cores = 2):
        self._engine = mph.Client(cores=num_cores)

        self.model_name = model_name
        self._model = self._engine.create(model_name)

    def initialize_model(self, chip_len, chip_wid, chip_thickness, **kwargs):
        '''
        Initialise the default constructs of the COMSOL model.

        Inputs:
            - chip_len, chip_wid, chip_thickness - Length, width and height of the chip in metres
            - pad_x, pad_y, pad_z - (Optional) The padding space (on either side - e.g. total padding along x-axis is 2*pad_x) between the
                                    simulation's exterior boundaries and the chip. Default value for all 3 axes is 0.5e-3.
        '''
        self.chip_len = chip_len
        self.chip_wid = chip_wid
        self.chip_thickness = chip_thickness

        self.pad_x = kwargs.get('pad_x', 0.5e-3)
        self.pad_y = kwargs.get('pad_y', 0.5e-3)
        self.pad_z = kwargs.get('pad_z', 0.5e-3)

        self._ports = []
        self._conds = []
        self._conds_coords = []
        self._fine_mesh = []
        self._num_sel = 0

        self._model.java.component().create("comp1", True)

        self._model.java.component("comp1").geom().create("geom1", 3)
        self._model.java.component("comp1").mesh().create("mesh1")        

        #Create physics and study for RF analysis (e.g. s-parameters)
        self._model.java.component("comp1").physics().create("emw", "ElectromagneticWaves", "geom1")
        self._model.java.study().create("std1")
        self._model.java.study("std1").create("freq", "Frequency")
        self._model.java.study("std1").feature("freq").set("solnum", "auto")
        self._model.java.study("std1").feature("freq").set("notsolnum", "auto")
        self._model.java.study("std1").feature("freq").set("savesolsref", jtypes.JBoolean(False))
        self._model.java.study("std1").feature("freq").set("ngen", "5")
        self._model.java.sol().create('solRFsparams')
        self._model.java.sol('solRFsparams').createAutoSequence('std1') 

        #Create physics and study for electrostatics (e.g. generating capacitance matrices)
        self._model.java.component("comp1").physics().create("es", "Electrostatics", "geom1")
        self._model.java.component("comp1").physics("es").prop("PortSweepSettings").set("useSweep", jtypes.JBoolean(True))
        self._model.java.param().set("PortName", jtypes.JInt(1))    #Just use default port name...
        self._model.java.study().create("stdCapMat")
        self._model.java.study("stdCapMat").create("capMat","Stationary")
        self._model.java.sol().create('solCapMat')
        self._model.java.sol('solCapMat').createAutoSequence('stdCapMat')

        #N.B. The ordering of the datasets here indicates that s-parameters live in dset1 while capMat lives in dset2.
        self._dset_sParams = "dset1"
        self._dset_capMat = "dset2"

        #Activate the appropriate physics (to be solved) for the given studies
        self._model.java.study("std1").feature("freq").activate("emw", jtypes.JBoolean(True))
        self._model.java.study("stdCapMat").feature("capMat").activate("emw", jtypes.JBoolean(False))
        self._model.java.study("std1").feature("freq").activate("es", jtypes.JBoolean(False))
        self._model.java.study("stdCapMat").feature("capMat").activate("es", jtypes.JBoolean(True))

        #Create the main bounding area
        self._model.java.component("comp1").geom("geom1").lengthUnit("m")
        self._create_block_corner('blk_chip', chip_len,chip_wid,chip_thickness, 0,0,-chip_thickness)
        self._create_block_corner('blk_boundary', chip_len+2*self.pad_x,chip_wid+2*self.pad_y,chip_thickness+2*self.pad_z, -self.pad_x,-self.pad_y,-self.chip_thickness-self.pad_z)

        #Create workplane and subsequent metallic geometry...
        self._model.java.component("comp1").geom("geom1").feature().create("wp1", "WorkPlane")
        self._model.java.component("comp1").geom("geom1").feature("wp1").set("unite", jtypes.JBoolean(True)) #Unite all objects...

        # self._model.java.component("comp1").geom("geom1").run()

    def add_metallic_Klayout(self, kLayoutObj, layer_id, **kwargs):
        '''
        Adds metallic conductors from the Klayout design onto the surface layer of the chip simulation. The idea is to supply the final Klayout layout object
        and the metallic polygons are added to the COMSOL simulation.

        Inputs:
            - kLayoutObj - A Klayout object (i.e. the object used when calling the function to save to a GDS file)
            - layer_id - The index of the layer from which to take the metallic polygons (e.g. in the old tutorials, that would be layer_photo)
            - cell_num - (Optional) The cell index in the Klayout object in which the layer resides. Default value is taken to be 0.
        '''
        cell_num = kwargs.get('cell_num', 0)

        self.kLy2metre = kLayoutObj.dbu*1e-6

        polys = [x for x in kLayoutObj.cell(cell_num).shapes(layer_id).each(db.Shapes.SPolygons)]
        for m in range(len(polys)):
            pol_name = "pol"+str(m)
            #Convert coordinates into metres...
            cur_poly = [[p.x*self.kLy2metre,p.y*self.kLy2metre] for p in polys[m].each_point_hull()]
            sel_x, sel_y, sel_r = self._create_poly(pol_name, cur_poly)
            #Get the selection point - any point inside the polygon...
            self._conds += [self._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]
            self._conds_coords += [np.array(cur_poly)]

    def create_port_on_CPW(self, CPW_obj, is_start=True, len_launch = 20e-6):
        '''
        Creates an RF port on a CPW inlet. The elements form fins to ground from the central CPW stripline. Note that the first port will automatically be
        the 50Ohm excitation port in the E-field plots while the second port is a 50Ohm ground. The s-parameters will calculate S11 and S21 if 2 such ports
        are defined.

        Inputs:
            - CPW_obj - A CPW object that has the attributes: start, end, width, gap
            - is_start - If True, then the port is attached to the start of the CPW, while False attaches the port to the end of the CPW
            - len_launch - (Default: 20e-6) Length of the inlet port fins along the the CPW. It is a good idea to keep it thin w.r.t. CPW gap 
        '''
        vec_launch = CPW_obj.end - CPW_obj.start
        vec_launch /= vec_launch.length()
        if not is_start:
            vec_launch = -vec_launch
            vec_ori = CPW_obj.end
        else:
            vec_ori = CPW_obj.start
        vec_ori *= self.kLy2metre
        vec_perp = db.DVector(vec_launch.y,-vec_launch.x)
        vec_launch *= len_launch

        cpw_wid = CPW_obj.width * self.kLy2metre
        cpw_gap = CPW_obj.gap * self.kLy2metre

        pol_name = "port_launch" + str(len(self._ports))
        
        launches = [vec_ori + vec_perp * cpw_wid*0.5, vec_ori + vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch + vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch + vec_perp * cpw_wid*0.5]
        launches = [[p.x,p.y] for p in launches]
        sel_x, sel_y, sel_r = self._create_poly(pol_name + "a", launches)
        cur_launch = [self._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]

        launches = [vec_ori - vec_perp * cpw_wid*0.5, vec_ori - vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch - vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch - vec_perp * cpw_wid*0.5]
        launches = [[p.x,p.y] for p in launches]
        sel_x, sel_y, sel_r = self._create_poly(pol_name + "b", launches)
        cur_launch += [self._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]

        #Each port is defined as: [portA-selection-name, portB-selection-name, vec_CPW2GND_1] where vec_CPW2GND_1 is a db.DVector pointing in the direction
        #of ground from the CPW for portA.
        cur_launch += [vec_perp]
        self._ports += [cur_launch]

    def register_fine_structure(self, dPoint, min_boundary_dist, max_boundary_dist):
        '''
        Registers a metallic surface polygon as fine-structure when meshing. The registration of the polygon is based off a point on the boundary or within
        the polygon noting that the final structure in Klayout will be the amalgamation of all the pieces...

        Inputs:
            - dPoint - any vector-like object with x and y attributes denoting any position residing within or on the boundary of the polygon.
            - min_boundary_dist - Minimum size of meshed-triangles on polygon's boundary
            - max_boundary_dist - Maximum size of meshed-triangles on polygon's boundary
        '''
        #Move points on boundary to within the chip as the boundary selection shouldn't select the side of the chip...
        if dPoint.x <= 10e-9:
            dPoint.x = 12e-9
        if dPoint.x >= self.chip_len-10e-9:
            dPoint.x = self.chip_len-12e-9
        if dPoint.y <= 10e-9:
            dPoint.y = 12e-9
        if dPoint.y >= self.chip_wid-10e-9:
            dPoint.y = self.chip_wid-12e-9
        #Store the fine meshes as a tuple: (selection-poly-name, min_boundary_dist, max_boundary_dist)
        self._fine_mesh += [(self._create_boundary_selection_sphere(10e-9, dPoint.x,dPoint.y), min_boundary_dist, max_boundary_dist)]

    def build_geom_mater_elec_mesh(self):
        '''
        Builds geometry, sets up materials, sets up electromagnetic parameters/ports and builds the mesh.
        '''
        #Create materials
        self._model.java.component("comp1").geom("geom1").run()
        self._create_material('Vacuum', 1.0)
        self._create_material('Si', 11.7, self._model.java.selection('geom1_blk_chip_dom').entities(3))

        #Note that PEC1 is the default exterior boundary condition
        cond_bounds = [self._get_selection_boundaries(x)[0] for x in self._conds]
        self._model.java.component("comp1").physics("emw").create("pec2", "PerfectElectricConductor", 2)
        self._model.java.component("comp1").physics("emw").feature("pec2").selection().set(jtypes.JArray(jtypes.JInt)(cond_bounds))
        #Create terminals for capacitance matrix simulations
        for cur_term in range(len(cond_bounds)):
            term_name = "term"+str(cur_term)
            self._model.java.component("comp1").physics("es").create(term_name, "Terminal", 2)
            self._model.java.component("comp1").physics("es").feature(term_name).selection().set(jtypes.JArray(jtypes.JInt)([cond_bounds[cur_term]]))
            self._model.java.component("comp1").physics("es").feature(term_name).set("TerminalType", "Voltage")
            self._model.java.component("comp1").physics("es").feature(term_name).set("TerminalName", jtypes.JInt(cur_term+1))
        #Create the excitation ports for RF simulations
        for cur_port_id,cur_port in enumerate(self._ports):
            port_name = "lport" + str(cur_port_id)
            self._model.java.component("comp1").physics("emw").create(port_name, "LumpedPort", 2)
            self._model.java.component("comp1").physics("emw").feature(port_name).set('PortType', 'MultiElementUniform')
            port_bndsA = self._get_selection_boundaries(cur_port[0])
            port_bndsB = self._get_selection_boundaries(cur_port[1])
            self._model.java.component("comp1").physics("emw").feature(port_name).selection().set(jtypes.JArray(jtypes.JInt)(port_bndsA+port_bndsB))
            self._model.java.component("comp1").physics("emw").feature(port_name).feature('ue1').selection().set(jtypes.JArray(jtypes.JInt)(port_bndsA))
            self._model.java.component("comp1").physics("emw").feature(port_name).feature('ue2').selection().set(jtypes.JArray(jtypes.JInt)(port_bndsB))
            self._model.java.component("comp1").physics("emw").feature(port_name).feature('ue1').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([cur_port[2].x,cur_port[2].y,0.0]))
            self._model.java.component("comp1").physics("emw").feature(port_name).feature('ue2').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([-cur_port[2].x,-cur_port[2].y,0.0]))

        #Create mesh
        for mesh_ind, cur_fine_struct in enumerate(self._fine_mesh):
            cur_polys = self._get_selection_boundaries(cur_fine_struct[0])
            if len(cur_polys) == 0:
                continue
            mesh_name = "ftri" + str(mesh_ind)
            self._model.java.component("comp1").mesh("mesh1").create(mesh_name, "FreeTri")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).create("size1", "Size")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).selection().set(jtypes.JArray(jtypes.JInt)(cur_polys))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('custom', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmaxactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hminactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmin', jtypes.JDouble(cur_fine_struct[1]))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmax', jtypes.JDouble(cur_fine_struct[2]))
        self._model.java.component("comp1").mesh("mesh1").create("ftet1", "FreeTet")
        self._model.java.component("comp1").mesh("mesh1").feature("ftet1").create("size1", "Size")
        self._model.java.component("comp1").mesh("mesh1").run()

    def set_freq_range(self, freq_start, freq_end, num_points, use_previous_solns = False):
        '''
        Set the frequency range to sweep in the s-parameter and E-field simulation.

        Inputs:
            - freq_start, freq_end - Start and end (inclusive) frequencies in units of Hertz
            - num_points - Number of points to use in the sweep (must be above 1)
            - use_previous_solns - (Default: False) If True, it will use the previous frequency value as a starting point for the current point; it can have the tendency
                                   to falsely find local minima and cause the solution to have 'discrete steps' - so it's better to keep it False if accuracy is important.
        '''
        assert num_points > 1, "Number of points must be greater than 1 when setting up a frequency sweep"
        self._model.java.study("std1").feature("freq").set("errestandadap", "none")
        self._model.java.study("std1").feature("freq").set("plist", "range({0}[GHz],({1}[GHz]-({0}[GHz]))/{2},{1}[GHz])".format(freq_start*1e-9, freq_end*1e-9, num_points-1))
        if use_previous_solns:
            self._model.java.study("std1").feature("freq").set("preusesol", "yes")
        else:
            self._model.java.study("std1").feature("freq").set("preusesol", "no")
        
    def run_simulation_sparams(self, recompute=True):
        '''
        Run simulation to get s-parameters after running the simulation. Returns a 3-row array in which the rows are: frequency values, S11s, S21s.
        
        Inputs:
            - recompute - (Default True) If true, the solution result is recomputed
        '''
        if (recompute):
            self._model.java.sol('solRFsparams').runAll()
        
        self._model.java.result().numerical().create("ev1", "Eval")
        self._model.java.result().numerical("ev1").set("data", self._dset_sParams)
        self._model.java.result().numerical("ev1").set("expr", "freq")
        freqs = self._model.java.result().numerical("ev1").getData()
        self._model.java.result().numerical("ev1").set("expr", "emw.S11dB")
        s11s = self._model.java.result().numerical("ev1").getData()
        self._model.java.result().numerical("ev1").set("expr", "emw.S21dB")
        s21s = self._model.java.result().numerical("ev1").getData()
        self._model.java.result().numerical().remove("ev1")

        #For some reason the returned arrays are rows of the same data apparently repeated across the columns...
        freqs = np.array(freqs[0])[:,1]
        s11s = np.array(s11s[0])[:,1]
        s21s = np.array(s21s[0])[:,1]

        return np.vstack([freqs,s11s,s21s])
        
    def run_simulation_capMat(self):
        '''
        Runs the simulation and returns a capacitance matrix.
        '''
        num_ports = len(self._conds)
        capMatFull = np.zeros([num_ports,num_ports])
        
        #Setup temporary results dataset 
        self._model.java.result().numerical().create("gmev1", "EvalGlobalMatrix")
        self._model.java.result().numerical("gmev1").set("data", self._dset_capMat)
        self._model.java.result().numerical("gmev1").set("expr", "es.C")
        
        for cur_port in range(num_ports):
            self._model.java.param().set("PortName", jtypes.JInt(cur_port+1))
            #Evaluate column of capacitance matrix
            self._model.java.sol('solCapMat').runAll()
            #Extract column from result
            capCol = self._model.java.result().numerical("gmev1").computeResult()
            capCol = np.array(capCol[0])
            capMatFull[:,cur_port] = capCol[:,cur_port]
        
        self._model.java.result().numerical().remove("gmev1")
        return capMatFull

    def display_conductor_indices(self):
        '''
        Plots a coloured visualisation of the metallic conductors and their corresponding row/column indices of the capacitance matrix.
        '''
        fig, ax = plt.subplots()
        fig.set_size_inches(10,10)
        colMap = mpl.cm.jet
        numConds = len(self._conds_coords)
        colMap = plt.get_cmap('jet', numConds)
        #Add metallic conductors as discrete colours
        coll = PolyCollection(self._conds_coords, array=np.arange(numConds)+1,cmap=colMap, edgecolors='none')
        ax.add_collection(coll)
        #Ensure uniform scale
        ax.autoscale_view()
        ax.set_aspect('equal', 'box')
        #Add a custom colorbar to show the colour key
        cax = fig.colorbar(coll, ticks=1+(np.arange(0,numConds)+0.5)*(numConds-1)/numConds, ax=ax)
        cax.ax.set_yticklabels(np.arange(0,numConds)+1)  # vertically oriented colorbar


    def save(self, file_name):
        self._model.save(file_name)

    def _create_block_corner(self, name, size_x,size_y,size_z, pos_x,pos_y,pos_z):
        '''
        Creates a block in the geometry.

        Inputs:
            - name - Unique name of the geometry object
            - size_x,size_y,size_z - Dimensions of the block in the prescribed units
            - pos_x,pos_y,pos_z    - Corner position of the block in the prescribed units
        '''
        self._model.java.component("comp1").geom("geom1").create(name, "Block")
        self._model.java.component("comp1").geom("geom1").feature(name).set("base", "corner")
        self._model.java.component("comp1").geom("geom1").feature(name).set("size", jtypes.JArray(jtypes.JDouble)([size_x,size_y,size_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set("pos", jtypes.JArray(jtypes.JDouble)([pos_x,pos_y,pos_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set('selresult', 'on') #To generate automatic selections...

    def _create_poly(self, name, poly_coords):
        '''
        Creates a polygon on workplane wp1. The return value is (x,y,r) where (x,y) is a point inside the polygon and r is the radius such that a sphere
        fits inside the polygon...

        Inputs:
            - name - Unique name of the geometry object
            - poly_coords - Coordinates (doesn't have to close) given as a list of lists: [[x1,y1], [x2,y2], ...]
        '''
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(name, "Polygon")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('selresult', 'on') #To generate automatic selections...
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('selresultshow', 'bnd')
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set("source", "table")
        min_ind = 0
        for row_no,cur_pt in enumerate(poly_coords):
            #Table is parsed as: value,row_id,col_id
            # self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).setIndex("table", jtypes.JDouble(cur_pt[0]), jtypes.JInt(row_no), 0)
            # self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).setIndex("table", jtypes.JDouble(cur_pt[1]), jtypes.JInt(row_no), 1)
            if cur_pt[1] <= poly_coords[min_ind][1]:
                min_ind = row_no
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('table', jtypes.JArray(jtypes.JDouble,2)(poly_coords))    #This is much faster...
        #Find point inside polygon...
        min_ind2 = (min_ind + 1) % len(poly_coords)
        vec1 = np.array([ poly_coords[min_ind2][0]-poly_coords[min_ind][0], poly_coords[min_ind2][1]-poly_coords[min_ind][1] ])
        vec2 = np.array([ poly_coords[min_ind-1][0]-poly_coords[min_ind][0], poly_coords[min_ind-1][1]-poly_coords[min_ind][1] ])
        vec1 = vec1/np.linalg.norm(vec1)
        vec2 = vec2/np.linalg.norm(vec2)
        if np.dot(vec1,vec2) < -0.999847695: # 179 degrees
            vec3 = np.array([0,1])  #Assumption of being the lower-most point...
        else:
            vec3 = vec1 + vec2
            vec3 = vec3/np.linalg.norm(vec3)
        epsilon_mov = 10e-9
        rad = np.sqrt(min(1-np.dot(vec1,vec3)**2, 1-np.dot(vec2,vec3)**2)) * 0.9 * epsilon_mov
        vec3 *= epsilon_mov
        return (vec3[0]+poly_coords[min_ind][0], vec3[1]+poly_coords[min_ind][1], rad)

    def _create_material(self, name, rel_permit, selected_domain=''):
        '''
        Creates a polygon on workplane wp1.

        Inputs:
            - name - Unique name of the material
            - rel_permit - Relative permittivity of material
            - selected_domain - Selected domain specified via commands like: self._model.java.selection('geom1_blk_chip_dom').entities(3)
        '''
        self._model.java.component("comp1").material().create(name, "Common")
        if (isinstance(selected_domain, str) and  selected_domain == ''):
            self._model.java.component("comp1").material(name).selection().all()
        else:
            #self._model.java.component("comp1").geom("geom1").feature('blk_chip').outputSelection()[4]
            self._model.java.component("comp1").material(name).selection().set(selected_domain)
        self._model.java.component("comp1").material(name).propertyGroup("def").set("relpermittivity", jtypes.JDouble(rel_permit))
        self._model.java.component("comp1").material(name).propertyGroup("def").set("relpermeability", jtypes.JDouble(1))
        self._model.java.component("comp1").material(name).propertyGroup("def").set("electricconductivity", jtypes.JDouble(0))

    def _create_boundary_selection_sphere(self, radius, pos_x,pos_y,pos_z=0.0):
        '''
        Creates a selection in which all boundaries within the sphere are selected (after the geometry has been fully built). Return value is the selection
        name/ID that can be used to query for domains later via the function _get_selection_boundaries.

        Inputs:
            - radius - Radius of sphere
            - pos_x,pos_y,pos_z - Position of sphere's centre
        '''
        sel_name = 'sel' + str(self._num_sel)
        self._num_sel += 1
        self._model.java.selection().create(sel_name, 'Ball')
        self._model.java.selection(sel_name).set('entitydim', '2')
        self._model.java.selection(sel_name).set('condition', 'intersects')
        self._model.java.selection(sel_name).set('posx', jtypes.JDouble(pos_x))
        self._model.java.selection(sel_name).set('posy', jtypes.JDouble(pos_y))
        self._model.java.selection(sel_name).set('posz', jtypes.JDouble(pos_z))
        self._model.java.selection(sel_name).set('r', jtypes.JDouble(radius))
        return sel_name

    def _get_selection_boundaries(self, sel_name):
        '''
        Returns a list of integers pertaining to the boundaries within a spherical-selection of ID given by sel_name.
        '''
        return [x for x in self._model.java.selection(sel_name).entities(2)]

    