import klayout.db
from math import sqrt, cos, sin, tan, atan2, pi, copysign
from klayout.db import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from klayout.db import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from .Coplanars import CPW_RL_Path, CPW
from .Coplanars import Coil_type_1, Coil_air_bridges, CPW_arc, CPW2CPW

from ClassLib.Coplanars import *
from ClassLib.Shapes import *

class CWave2CPW( Element_Base ):
    '''
    Draws a semi-circle coupler from coplanar waveguide to jelly capacitance plates.
    '''
    def __init__( self, c_wave_cap, params, n_pts=50, trans_in=None ):
        self.c_wave_ref = c_wave_cap
        if isinstance(params, dict):
            self.Z1 = params['Z1']
            self.d_alpha1 = params['d_alpha1']
            self.width1 = params['width1']
            self.gap1 = params['gap1']
            self.Z2 = params['Z2']
            self.d_alpha2 = params['d_alpha2']
            self.width2 = params['width2']
            self.gap2 = params['gap2']
        else:
            # not recommended
            self.Z1 = params[0]
            self.d_alpha1 = params[1]
            self.width1 = params[2]
            self.gap1 = params[3]
            self.Z2 = params[4]
            self.d_alpha2 = params[5]
            self.width2 = params[6]
            self.gap2 = params[7]
        self.n_pts = n_pts
        super().__init__(self.c_wave_ref.origin, trans_in)

    def _get_solid_arc( self, center, R, width, alpha_start, alpha_end, n_inner, n_outer ):
        pts = []

        d_alpha_inner = (alpha_end - alpha_start)/(n_inner - 1)
        d_alpha_outer = -(alpha_end - alpha_start)/(n_outer-1)

        for i in range( 0,n_inner ):
            alpha = alpha_start + d_alpha_inner*i
            pts.append( center + DPoint( cos( alpha ), sin( alpha ) )*(R - width/2) )
        for i in range( 0,n_outer ):
            alpha = alpha_end + d_alpha_outer*i
            pts.append( center + DPoint( cos( alpha ), sin( alpha ) )*(R + width/2) )

        return DSimplePolygon( pts )

    def init_regions( self ):
        origin = DPoint(0,0)
        arc_1_solid = self._get_solid_arc( origin, self.c_wave_ref.in_circle.r + self.gap1 + self.width1/2, self.width1,
                                                    pi/2 - self.d_alpha1/2, pi/2 + self.d_alpha1/2,
                                                    self.n_pts,self.n_pts )
        arc_1_empty = self._get_solid_arc( origin, self.c_wave_ref.in_circle.r + self.gap1/2, self.gap1,
                                                    pi/2 - self.d_alpha1/2, pi/2 + self.d_alpha1/2,
                                                    self.n_pts,self.n_pts )

        arc_2_solid = self._get_solid_arc( origin, self.c_wave_ref.in_circle.r + self.gap2 + self.width2/2, self.width2,
                                                    3/2*pi - self.d_alpha2/2, 3/2*pi + self.d_alpha2/2,
                                                    self.n_pts,self.n_pts )
        arc_2_empty = self._get_solid_arc( origin, self.c_wave_ref.in_circle.r + self.gap2/2, self.gap2,
                                                    3/2*pi - self.d_alpha2/2, 3/2*pi + self.d_alpha2/2,
                                                    self.n_pts,self.n_pts )
        self.metal_region.insert( arc_1_solid )
        self.metal_region.insert( arc_2_solid )
        self.empty_region.insert( arc_1_empty )
        self.empty_region.insert( arc_2_empty )

class CWave( Complex_Base ):
    '''
    Draws a condensator from a circle cutted into 2 pieces.
    '''

    def __init__(self, center, r_out, dr, n_segments, s, alpha, r_curve, delta=40e3, n_pts=50, solid=True, trans_in=None ):
        '''
        Parameters:
        center: DPoint
            A center of circle.
        r_out: float
            The outer radius of a circle (along to the edge of the ground).
        dr: float
            The gap in between the common ground and the circle perimeter.
        n_segments: int
            The number of segments (without turns) composed into a condensator gap.
        s: float
            The half-width of the gap.
        alpha: rad
            The angle of single arc of slice.
        r_curve: float
            The radius of single arc.
        delta: float
            length of the horizontal lines on the ends of the cut
        n_pts: int
            The number of points on the perimeter of the circle.
        solid: ???
        trans_in: Bool
            Initial transformation
        '''
        self.r_out = r_out
        self.dr = dr
        self.r_in = self.r_out - self.dr
        self.n_segments = n_segments
        self.s = s
        self.alpha = alpha
        self.r_curve = r_curve
        self.n_pts = n_pts
        self.delta = delta
        self.L_full = 2*(self.r_out - self.dr) - 2*self.delta
        # calculating parameters of the CPW_RL_Path #
        L_full = self.L_full
        alpha = self.alpha
        if abs(self.alpha) != pi:
            self.L1 = L_full/((self.n_segments+1)*cos(alpha))
            self.L0 = self.L1/2
        else:
            raise ValueError("180 degrees turns in CWave are not supported.")

        if( self.L0 < 0 ):
            print("CPW_RL_Path: impossible parameters combination")

        super(). __init__(center,trans_in)

    def init_primitives(self):
        origin = DPoint(0,0)
        #erased line params
        Z = CPWParameters( 0, self.s/2 )
        shapes = ''
        angles = []
        lengths = []
        # placing circle r_out with dr clearance from ground polygon
        self.empt_circle = Circle( origin,self.r_out, n_pts=self.n_pts, solid=False )
        self.in_circle = Circle( origin, self.r_out - self.dr, n_pts=self.n_pts, solid=True )
        self.empt_circle.empty_region -= self.in_circle.metal_region
        self.primitives["empt_circle"] = self.empt_circle
        self.primitives["in_circle"] = self.in_circle

        self.RL_start = origin + DPoint( -self.in_circle.r,0 )
        shapes += 'LRL'
        angles.extend([self.alpha])
        lengths.extend([self.delta, self.L0])
        #rl_path_start = CPW_RL_Path(self.RL_start, "LRLR", Z, self.r_curve, [self.delta, self.L0], [self.alpha,-self.alpha] )
        #self.primitives["rl_start"] = rl_path_start

        # intermidiate RLRs
        for i in range(self.n_segments):
            if( i%2 == 1 ):
                m_x = -1
            else:
                m_x = 1
            shapes += 'RL'
            angles.extend([-2*m_x*self.alpha])
            lengths.extend([self.L1])
            # prev_path = list(self.primitives.values())[-1]
            # rl_path_p = CPW_RL_Path( prev_path.end, "RLR", Z, self.r_curve, [self.L1], [-m_x*self.alpha,m_x*self.alpha])
            # self.primitives["rl_path_" + str(i)] = rl_path_p
            # ending RLR
        if( self.n_segments%2 == 1 ):
            m_x = 1
        else:
            m_x = -1
        shapes += 'RLRL'
        angles.extend([2*m_x*self.alpha,-m_x*self.alpha])
        lengths.extend([self.L0, self.delta])
        cut = CPW_RL_Path(self.RL_start,shapes,Z,self.r_curve,lengths,angles)
        # prev_path = list(self.primitives.values())[-1]
        # rl_path_end = CPW_RL_Path( prev_path.end, "RLRL", Z, self.r_curve, [self.L0, self.delta], [m_x*self.alpha,-m_x*self.alpha])
        self.primitives["cut"] = cut

        
# Large capacitor mod class
# I (Marcus) think I found a small bug in the large cap class,
# I corrected it below (where noted).
class Large_capacitor_mod( Element_Base ):
    def __init__( self, origin, d_1, d_2, N, gap_2, spacing, gap_3, trans_in = None, fingers_same = False ):
        print("Large_capacitor_mod class has been deprecated - consider using Capacitor_Interdigitated")

        # origin corresponds to the middle of the bottom pad
        self.origin = origin
        # Width of the capacitor
        self.d_1 = d_1
        # Space between 2 pads of capacitor 
        self.d_2 = d_2
        # Number of connections between the 2 pads 
        self.N = N 
        self.gap_2 = gap_2
        self.spacing = spacing
        
        # the gap between pad and each finger tip is hard coded to be 1 um, Eric added a new argument for this gap
        self.gap_3 = gap_3 
        
        super().__init__(self.origin, trans_in)
    
    def init_regions( self ):
        
         
        self.d = (self.N -1)*(2*self.spacing + 2*self.d_1) + self.d_1
        
        self.p0 = DPoint(-self.d/2,0)
        self.p1 = DPoint(self.d/2, 5e3)
        self.metal0 = klayout.db.DBox( self.p0, self.p1)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal0) )
        
        for i in range(self.N-1):
            self.p2 = DPoint(i*(2*self.d_1+2*self.spacing) -self.d/2, 5e3)
            self.p3 = DPoint(i*(2*self.d_1+2*self.spacing) + self.d_1 -self.d/2, 5e3+self.d_2)
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal_region.insert( klayout.db.Box().from_dbox(self.metal1) ) 
            
            self.p4 = DPoint(i*(2*self.d_1+2*self.spacing) + self.d_1 + self.spacing -self.d/2, 6e3)
            self.p5 = DPoint(i*(2*self.d_1+2*self.spacing) + 2*self.d_1 +self.spacing -self.d/2, 5e3+self.d_2 + self.gap_3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal_region.insert( klayout.db.Box().from_dbox(self.metal2) )
            
            self.empty0 = klayout.db.DBox( self.p3, self.p3 + DPoint(-self.d_1, self.gap_3))
            self.empty1 = klayout.db.DBox( self.p3 + DPoint(0,self.gap_3), self.p4 + DPoint(0, -self.gap_3))
            self.empty2 = klayout.db.DBox( self.p4 + DPoint(0,-self.gap_3), self.p4 + DPoint(self.d_1, 0))
            self.empty3 = klayout.db.DBox( self.p4 + DPoint(self.d_1,-self.gap_3), self.p5 + DPoint(self.spacing, 0))
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty0) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty1) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty2) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty3) )
            
        # Last Finger 
        #---- Below Modified -----#
        self.p2 = DPoint((self.N-1)*(2*self.d_1+2*self.spacing) -self.d/2, 5e3)
        self.p3 = DPoint((self.N-1)*(2*self.d_1+2*self.spacing) -self.d/2 + self.d_1, 5e3+self.d_2)
        #-----Above Modified------#
        self.metal3 = klayout.db.DBox( self.p3, self.p2)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal3) )
        
        #Top Pad:
        self.p6 = DPoint(-self.d/2, self.d_2 + 1e3 + 5e3)
        self.p7 =DPoint(self.d/2, self.d_2 + 1e3 + 10e3)
        self.metal3 = klayout.db.DBox( self.p6, self.p7)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal3) )
      
      # No need for empty regions, when inserting capacitors in empty regions:
        self.p8 = DPoint(-self.d/2-self.gap_2,0)
        self.p9 = DPoint(-self.d/2, self.d_2 + 1e3 + 10e3)
        self.empty4 = klayout.db.DBox( self.p8, self.p9)
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty4) )
        
        self.p10 = DPoint(self.d/2+self.gap_2,0)
        self.p11= DPoint(self.d/2, self.d_2 + 1e3 + 10e3)
        self.empty5 = klayout.db.DBox( self.p10, self.p11)
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty5) )
       
        # gap at top of last finger:
        self.p12 = DPoint(self.d/2-self.d_1, self.d_2 + 5e3)
        self.p13 = DPoint(self.d/2, self.d_2 + 5e3 +self.gap_3)
        self.empty6 = klayout.db.DBox(self.p12, self.p13 )
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty6) ) 
        

class Large_capacitor( Element_Base ):
    def __init__( self, origin, d_1, d_2, N, gap_2, spacing, trans_in = None ):
        print("Large_capacitor class has been deprecated - consider using Capacitor_Interdigitated")

        # origin corresponds to the middle of the bottom pad
        self.origin = origin
        # Width of the capacitor
        self.d_1 = d_1
        # Space between 2 pads of capacitor 
        self.d_2 = d_2
        # Number of connections between the 2 pads 
        self.N = N
        self.gap_2 = gap_2
        self.spacing = spacing
        
        super().__init__(self.origin, trans_in)
    
    def init_regions( self ):
        
         
        self.d = (self.N -1)*(2*self.spacing + 2*self.d_1) + self.d_1
        
        self.p0 = DPoint(-self.d/2,0)
        self.p1 = DPoint(self.d/2, 5e3)
        self.metal0 = klayout.db.DBox( self.p0, self.p1)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal0) )
        
        for i in range(self.N-1):
            self.p2 = DPoint(i*(2*self.d_1+2*self.spacing) -self.d/2, 5e3)
            self.p3 = DPoint(i*(2*self.d_1+2*self.spacing) + self.d_1 -self.d/2, 5e3+self.d_2)
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal_region.insert( klayout.db.Box().from_dbox(self.metal1) ) 
            
            self.p4 = DPoint(i*(2*self.d_1+2*self.spacing) + self.d_1 + self.spacing -self.d/2, 6e3)
            self.p5 = DPoint(i*(2*self.d_1+2*self.spacing) + 2*self.d_1 +self.spacing -self.d/2, 5e3+self.d_2+1e3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal_region.insert( klayout.db.Box().from_dbox(self.metal2) )
            
            self.empty0 = klayout.db.DBox( self.p3, self.p3 + DPoint(-self.d_1, 1e3))
            self.empty1 = klayout.db.DBox( self.p3 + DPoint(0,1e3), self.p4 + DPoint(0, -1e3))
            self.empty2 = klayout.db.DBox( self.p4 + DPoint(0,-1e3), self.p4 + DPoint(self.d_1, 0))
            self.empty3 = klayout.db.DBox( self.p4 + DPoint(self.d_1,-1e3), self.p5 + DPoint(self.spacing, 0))
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty0) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty1) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty2) )
            self.empty_region.insert( klayout.db.Box().from_dbox(self.empty3) )
            
         #-----------------------------------#
        # Comment from Marcus: I think there's a little bug in the following two lines
        # Copy corresponding code from Large_capacitor_mod to fix
        self.p2 = DPoint((self.N-1)*(2*self.d_1+self.spacing)-self.d/2, 5e3)
        self.p3 = DPoint((self.N-1)*(2*self.d_1+self.spacing) + self.d_1-self.d/2, 5e3+self.d_2)
        #-----------------------------------#
        self.metal1 = klayout.db.DBox( self.p2, self.p3)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal1) )
      
        self.p6 = DPoint(-self.d/2, self.d_2 + 1e3 + 5e3)
        self.p7 =DPoint(self.d/2, self.d_2 + 1e3 + 10e3)
        self.metal3 = klayout.db.DBox( self.p6, self.p7)
        self.metal_region.insert( klayout.db.Box().from_dbox(self.metal3) )
        
        self.p8 = DPoint(-self.d/2-self.gap_2,0)
        self.p9 = DPoint(-self.d/2, self.d_2 + 1e3 + 10e3)
        self.empty4 = klayout.db.DBox( self.p8, self.p9)
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty4) )
        
        self.p10 = DPoint(self.d/2+self.gap_2,0)
        self.p11= DPoint(self.d/2, self.d_2 + 1e3 + 10e3)
        self.empty5 = klayout.db.DBox( self.p10, self.p11)
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty5) )
        
        self.p12 = DPoint(self.d/2-self.d_1, self.d_2 + 5e3)
        self.p13 = DPoint(self.d/2, self.d_2 + 6e3)
        self.empty6 = klayout.db.DBox(self.p12, self.p13 )
        self.empty_region.insert( klayout.db.Box().from_dbox(self.empty6) )
        
        
class Connector_large_capacitor( Complex_Base ):
    def __init__( self, origin, width, L_sep, gap_1, gap_2, spacing, d_1, d_2, N, gap_3, trans_in = None ):
        print("Connector_large_capacitor class has been deprecated - consider using the cpw_params attribute in Capacitor_Interdigitated")
        
        # origin corresponds to the middle of the bottom pad
        self.origin = origin
        # Width of the capacitor
        self.d_1 = d_1
        # Space between 2 pads of capacitor 
        self.d_2 = d_2
        # Number of connections between the 2 pads 
        self.N = N
        # Width of the line connecting the capacitor to circuit
        self.width = width
        # Gap near connecting wire
        self.gap_1 = gap_1
        # Gap around the capacitor 
        self.gap_2 = gap_2 
        # Spacint between lines of capacitor 
        self.spacing = spacing
        # Distance between transmission line and capacitor
        self.L_sep = L_sep
        self.d = (self.N -1)*(2*self.spacing + 2*self.d_1) + self.d_1
        
        #the gap between pad and each finger tip is hard coded to be 1 um, Eric add a new argument for this gap
        self.gap_3 = gap_3

        super().__init__( origin, trans_in )

        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):      
        
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 + DPoint(0, 100e3)
        self.p2 = self.p1 + DPoint(0, 20e3)
        self.p3 = self.p2 + DPoint(0, 2e3)
        self.p4 = self.p3 + DPoint(0, self.d_2 + 1e3 + 10e3)
        self.p5 = self.p4 + DPoint(0, 2e3)
        self.p6 = self.p5 + DPoint(0, 20e3)
        self.p7 = self.p6 + DPoint(0, self.L_sep)
        
        
        
        self.transmission_connection_1 = CPW( self.width, self.gap_1, self.p0, self.p1)
        self.primitives["transmission_connection_1"] = self.transmission_connection_1
        
        self.capa_connex_1 = CPW( self.d, self.gap_2, self.p2, self.p3)
        self.primitives["capa_connex_1"] = self.capa_connex_1
        
        self.joint_1 = CPW2CPW( self.transmission_connection_1, self.capa_connex_1,
                                 self.transmission_connection_1.end, self.capa_connex_1.start)
        self.primitives["joint_1"] = self.joint_1
        
        self.capa = Large_capacitor_mod(self.p3, self.d_1, self.d_2, self.N, self.gap_2, self.spacing, self.gap_3)
        self.primitives["capa"] = self.capa
        
        self.capa_connex_2 = CPW( self.d, self.gap_2, self.p4, self.p5)
        self.primitives["capa_connex_2"] = self.capa_connex_2
        
        self.transmission_connection_2 = CPW( self.width, self.gap_1, self.p6, self.p7)
        self.primitives["transmission_connection_2"] = self.transmission_connection_2
        
        self.joint_2 = CPW2CPW(self.capa_connex_2, self.transmission_connection_2, 
                                self.capa_connex_2.end, self.transmission_connection_2.start)
        self.primitives["joint_2"] = self.joint_2
        
        self.connections = [self.p7, self.p0]
        self.angle_connections = [0,self.transmission_connection_2.alpha_end]


class Capacitor_Interdigitated( Element_Base ):
    def __init__( self, origin, fing_wid, fing_gap, fing_len, inter_dig_gap, pad_width, side_gap, N, fing_config = 'diff_N_start', trans_in = None ):
        '''
        Draws an interdigitated capacitor in which the first pad (thus, the origin) is in the bottom and the second pad is on the top. The fingers thus extend vertically.
        Thus, if using trans_in with an integer rotation, 3 orients the capacitor left-to-right horizontally.

        The cpw_params attribute can be used to get the CPW parameters (e.g. width and gap) required to join the capacitor into another CPW object.

        Inputs: 
            - origin - Origin corresponds to the middle of the bottom pad
            - fing_wid - Width of the a finger
            - fing_gap - Gap between adjacent fingers
            - fing_len - Length of the fingers (from their associated pads)
            - inter_dig_gap - Gap between a given finger and the opposite pad
            - pad_width - Capacitor pad widths (that is, the padding on the top and bottom of the capacitors that sandwich the fingers)
            - side_gap - Empty space to pad the left and right hand sides of the capacitor
            - N - Number of fingers per side of capacitor
            - fing_config - The arrangement of the fingers:
                                - diff_N_start - There are N-1 fingers on the end-pad slot into the N fingers on the starting pad
                                - diff_N_end   - There are N-1 fingers on the start-pad slot into the N fingers on the ending pad
                                - same_N_left  - There are N fingers on both pads with the left-most finger being from the bottom pad
                                - same_N_right - There are N fingers on both pads with the left-most finger being from the top pad
            - trans_in - Klayout transformation - e.g. klayout.db.Trans(0, False, 0, 0)
        '''    
        # origin corresponds to the middle of the bottom pad
        self.origin = origin
        self.fing_wid = fing_wid
        self.fing_gap = fing_gap
        self.fing_len = fing_len
        self.inter_dig_gap = inter_dig_gap
        self.side_gap = side_gap
        self.pad_width = pad_width
        self.N = N
        self.fing_config = fing_config
        self.start = self.origin
        super().__init__(self.origin, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
    
    def init_regions( self ):
        odd_config = self.fing_config == 'diff_N_start' or self.fing_config == 'diff_N_end'
        if odd_config:
            cap_wid = 2*(self.N-1)*(self.fing_gap + self.fing_wid) + self.fing_wid
        else:
            cap_wid = 2*self.N*(self.fing_gap + self.fing_wid) - self.fing_gap
        
        yUpper = self.pad_width + self.fing_len + self.inter_dig_gap

        #Bottom Pad
        self.metal_region.insert(klayout.db.DBox( DPoint(-cap_wid/2,0), DPoint(cap_wid/2,self.pad_width) ))
        #Top Pad
        self.metal_region.insert(klayout.db.DBox( DPoint(-cap_wid/2,yUpper), DPoint(cap_wid/2,yUpper + self.pad_width) ))
        
        self.connections = [DPoint(0,0),DPoint(0, yUpper + self.pad_width)]
        self.cpw_params = CPWParameters(cap_wid, self.side_gap)
        
        #Write the unit cells:
        #                   __________y2
        #   p3  p6OO     p8
        #   .     OO      .
        #p2 .     OOp5    .
        #  OO     ..      OO
        #  OO     ..      OO
        #  OOp1  p4 p7    OO__________y1
        #
        unit_cell_wid = 2*(self.fing_wid+self.fing_gap)
        if self.fing_config == 'diff_N_start' or self.fing_config == 'same_N_left':
            y1 = self.pad_width
            y2 = yUpper
            cy1 = y1 + self.fing_len        #y-value of first finger's edge
            cy2 = y1 + self.inter_dig_gap   #y-value of second finger's edge
        else:
            y1 = yUpper
            y2 = self.pad_width
            cy1 = y1 - self.fing_len
            cy2 = y1 - self.inter_dig_gap
        for i in range(self.N):
            p1 = DPoint(i*unit_cell_wid - cap_wid/2 + self.fing_wid, y1)
            p2 = DPoint(i*unit_cell_wid - cap_wid/2, cy1)
            p3 = DPoint(p1.x, y2)
            p4 = p1 + DPoint(self.fing_gap,0)
            p5 = DPoint(p4.x + self.fing_wid,cy2)
            p6 = DPoint(p4.x, p3.y)
            p7 = DPoint(p5.x, p4.y)
            p8 = DPoint(p7.x + self.fing_gap, p6.y)

            #First finger of unit-cell
            self.metal_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p1, p2 )))
            self.empty_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p2, p3 )))
            self.empty_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p1, p6 )))
            #Second finger of unit-cell
            if not odd_config or i < self.N-1:
                self.metal_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p5, p6 )))
                self.empty_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p4, p5 )))
                if i < self.N-1:    #Don't need empty region on final finger...
                    self.empty_region.insert(klayout.db.Box().from_dbox(klayout.db.DBox( p7, p8 )))

        #Padding on left and right hand sides
        self.empty_region.insert(klayout.db.DBox( DPoint(-cap_wid/2-self.side_gap,0), DPoint(-cap_wid/2,yUpper+self.pad_width) ))
        self.empty_region.insert(klayout.db.DBox( DPoint(cap_wid/2+self.side_gap,0), DPoint(cap_wid/2,yUpper+self.pad_width) ))

class Capacitor_GapCPW( Element_Base ):
    def __init__( self, origin, cpw_params, gap_size, trans_in = None ):
        '''
        Draws a gap-capacitor in which the first pad (thus, the origin) is in the bottom and the second pad is on the top. The gap is a CPW with no central conductor.
        Thus, if using trans_in with an integer rotation, 3 orients the capacitor left-to-right horizontally.

        The cpw_params attribute can be used to get the CPW parameters (e.g. width and gap) required to join the capacitor into another CPW object.

        Inputs: 
            - origin - Origin corresponds to the middle of the bottom pad
            - cpw_params - CPW or CPWParameters holding the parameters to inherit width and side-gap
            - gap_size - actual size of gap across the transmission line (forming the capacitor)
            - trans_in - Klayout transformation - e.g. klayout.db.Trans(0, False, 0, 0)
        '''    
        # origin corresponds to the middle of the bottom pad
        self.origin = origin
        self.cpw_params = CPWParameters(cpw_params.width, cpw_params.gap)
        self.gap_size = gap_size
        self.start = self.origin
        super().__init__(self.origin, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
    
    def init_regions( self ):
        half_width = self.cpw_params.width*0.5 + self.cpw_params.gap
        self.empty_region.insert(klayout.db.DBox( DPoint(-half_width,0), DPoint(half_width,self.gap_size) ))
        self.connections = [DPoint(0,0),DPoint(0, self.gap_size)]

        