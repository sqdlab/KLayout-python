import klayout.db
from math import sqrt, cos, sin, atan2, pi, copysign
from klayout.db import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from klayout.db import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from ClassLib.BaseClasses import *

class QBit_Flux_1( Element_Base ):
    def __init__( self, origin, params, trans_in=None ):
      
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]
        self.B1_width = params[11]
        self.B1_height = params[12]
        self.B2_height = params[13]
        self.B5_height = params[14]
        self.B6_width = params[15]
        self.B6_height = params[16]
        self.B7_width = params[17]
        self.B7_height = params[18]
        self.dCap = params[19]
        self.gap = params[20]        
        
        # calculated parameters
        self.B2_width = (self.gap - self.dCap - self.B6_width/2 - 4*self.w - self.jos1_a - self.jos2_a)/2
        self.B5_width = self.B2_width
        self.B3_width = self.a + self.jos1_a/2 - self.jos2_a/2 - self.w
        self.B4_width = self.a + self.jos2_a/2 - self.jos1_a/2 - self.w
        self._alpha_1 = 0.5
        self._length1 = (self.b + 2*self.w + self.jos1_b)
        self._alpha_2 = 0.5
        self._length2 = (self.b + 2*self.w - 2*self.f2 - self.jos2_b) #length
        
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 + DPoint( self.dCap,0 )
        self.p2 = self.p1 + DPoint( self.B1_width, (self.B1_height - self.B2_height)/2 )
        self.p3 = self.p2 + DPoint( self.B2_width, self.B2_height )
        self.p4 = self.p3 + DPoint( self.w, -self.w )
        self.p5 = self.p4 + DPoint( self.B3_width, self.w )
        self.p6 = self.p5 - DPoint( 0,self.f2 + self._length2*self._alpha_2 ) 
        self.p7 = self.p6 + DPoint( 2*self.w + self.jos2_a, -self.jos2_b )
        self.p8 = self.p7 - DPoint( 0,self.f2 + (1-self._alpha_2)*self._length2 )
        self.p9 = self.p8 - DPoint( self.B4_width + self.w,0 )
        self.p10 = self.p8 + DPoint( self.B5_width, (self.B5_height - self.B6_height)/2 )
        self.p11 = self.p10 + DPoint( self.B6_width, (self.B6_height - self.B7_height)/2 )

        self.B1 = klayout.db.DBox( self.p1, self.p1 + DPoint( self.B1_width,self.B1_height ) )
        self.B2 = klayout.db.DBox( self.p2, self.p3 )
        self.B3 = klayout.db.DBox( self.p4, self.p5 )
        self.B4 = klayout.db.DBox( self.p9, self.p9 + DPoint( self.B4_width, self.w ) )
        self.B5 = klayout.db.DBox( self.p8, self.p8 + DPoint( self.B5_width, self.B5_height ) )
        self.B6 = klayout.db.DBox( self.p10, self.p10 + DPoint( self.B6_width, self.B6_height ) )
        self.B7 = klayout.db.DBox( self.p11, self.p11 + DPoint( self.B7_width, self.B7_height ) )

        self.poly_1 = self._make_polygon( self._length1*self._alpha_1, self.w, self.d1, self.f1, self.jos1_b )
        self.poly_1.transform( DCplxTrans( 1.0, 270, False, self.p3 ) )
        self.poly_2 = self._make_polygon( self._length1*(1-self._alpha_1), self.w, self.d1, self.f1, self.jos1_b )
        self.poly_2.transform( DCplxTrans( 1.0, 90, False, self.p9 ) )
        self.poly_3 = self._make_polygon( self._length2*self._alpha_2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_3.transform( DCplxTrans( 1.0, 270, False, self.p5 ) )
        self.poly_4 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_4.transform( DCplxTrans( 1.0, 90, False, self.p7 ) ) 
        self.poly_5 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_5.transform( DCplxTrans( 1.0, 270, False, self.p6 ) )           
        self.poly_6 = self._make_polygon( (1-self._alpha_2)*self._length2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_6.transform( DCplxTrans( 1.0, 90, False, self.p8 ) )               
        super().__init__( origin, trans_in )
        
    def _make_polygon( self, length, w, d, f, overlapping ):
        polygon = DSimplePolygon
        p1 = DPoint(0,0)
        p2 = p1 + DPoint( length,0 )
        p3 = p2 + DPoint( 0, w )
        p4 = p3 - DPoint( overlapping,0 )
        p5 = p4 - DPoint( 0, d )
        p6 = p5 - DPoint( f,0 )
        p7 = p6 + DPoint( 0, d )
        p8 = p1 + DPoint( 0,w )
        
        polygon = DSimplePolygon( [p1,p2,p3,p4,p5,p6,p7,p8] )
        return polygon
        
    def init_regions( self ):
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B1) ) 
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B2) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B3) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B4) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B5) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B6) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B7) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_1) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_2) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_3) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_4) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_5) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_6) )
            

class QBit_Flux_2( Element_Base ):
    def __init__( self, origin, params, trans_in=None ):
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]
        self.B1_width = params[11]
        self.B1_height = params[12]
        self.B2_width = params[13]
        self.B5_width = params[14]
        self.B6_width = params[15]
        self.B6_height = params[16]
        self.B7_width = params[17]
        self.B7_height = params[18]
        self.dCap = params[19]
        self.gap = params[20]        
        
        # calculated parameters
        self.B2_height = (self.gap - self.dCap - self.B6_width/2 - 4*self.w - self.jos1_a - self.jos2_a)/2
        self.B5_height = self.B2_height
        self.B3_width = self.a + self.jos1_a/2 - self.jos2_a/2 - self.w
        self.B4_width = self.a + self.jos2_a/2 - self.jos1_a/2 - self.w
        self._alpha_1 = 0.5
        self._length1 = (self.b + 2*self.w + self.jos1_b)
        self._alpha_2 = 0.5
        self._length2 = (self.b + 2*self.w - 2*self.f2 - self.jos2_b) #length
        
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 - DPoint( 0,self.dCap + self.B1_height )
        self.p2 = self.p1 + DPoint( (self.B1_width - self.B2_width)/2, -self.B2_height )
        self.p3 = self.p2 + DPoint( (self.B2_width - self.B3_width)/2 - self.w, 0 )
        self.p4 = self.p3 + DPoint( self.w, -self.w )
        self.p5 = self.p2 + DPoint( (self.B2_width + self.B3_width)/2 , 0 )
        self.p6 = self.p5 - DPoint( 0,self.f2 + self._length2*self._alpha_2 ) 
        self.p7 = self.p6 + DPoint( 2*self.w + self.jos2_a, -self.jos2_b )
        self.p8 = self.p7 - DPoint( 0,self.f2 + (1-self._alpha_2)*self._length2 )
        self.p9 = self.p8 - DPoint( self.w + (self.B5_width + self.B4_width)/2,self.B5_height )
        self.p10 = self.p8 - DPoint( self.w + self.B4_width,0 ) 
        self.p11 = self.p9 + DPoint( (self.B5_width - self.B6_width)/2, -self.B6_height )
        self.p12 = self.p11 + DPoint( (self.B6_width - self.B7_width)/2, -self.B7_height )
        
        self.B1 = klayout.db.DBox( self.p1, self.p1 + DPoint( self.B1_width,self.B1_height ) )
        self.B2 = klayout.db.DBox( self.p2, self.p2 + DPoint( self.B2_width,self.B2_height ) )
        self.B3 = klayout.db.DBox( self.p4, self.p5 )
        self.B4 = klayout.db.DBox( self.p10, self.p10 + DPoint( self.B4_width, self.w ) )
        self.B5 = klayout.db.DBox( self.p9, self.p9 + DPoint( self.B5_width, self.B5_height ) )
        self.B6 = klayout.db.DBox( self.p11, self.p11 + DPoint( self.B6_width, self.B6_height ) )
        self.B7 = klayout.db.DBox( self.p12, self.p12 + DPoint( self.B7_width, self.B7_height ) )

        self.poly_1 = self._make_polygon( self._length1*self._alpha_1, self.w, self.d1, self.f1, self.jos1_b )
        self.poly_1.transform( DCplxTrans( 1.0, 270, False, self.p3 ) )
        self.poly_2 = self._make_polygon( self._length1*(1-self._alpha_1), self.w, self.d1, self.f1, self.jos1_b )
        self.poly_2.transform( DCplxTrans( 1.0, 90, False, self.p10 ) )
        self.poly_3 = self._make_polygon( self._length2*self._alpha_2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_3.transform( DCplxTrans( 1.0, 270, False, self.p5 ) )
        self.poly_4 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_4.transform( DCplxTrans( 1.0, 90, False, self.p7 ) ) 
        self.poly_5 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_5.transform( DCplxTrans( 1.0, 270, False, self.p6 ) )           
        self.poly_6 = self._make_polygon( (1-self._alpha_2)*self._length2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_6.transform( DCplxTrans( 1.0, 90, False, self.p8 ) )               
        super().__init__( origin, trans_in )
        
    def _make_polygon( self, length, w, d, f, overlapping ):
        polygon = DSimplePolygon
        p1 = DPoint(0,0)
        p2 = p1 + DPoint( length,0 )
        p3 = p2 + DPoint( 0, w )
        p4 = p3 - DPoint( overlapping,0 )
        p5 = p4 - DPoint( 0, d )
        p6 = p5 - DPoint( f,0 )
        p7 = p6 + DPoint( 0, d )
        p8 = p1 + DPoint( 0,w )
        
        polygon = DSimplePolygon( [p1,p2,p3,p4,p5,p6,p7,p8] )
        return polygon
        
    def init_regions( self ):
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B1) ) 
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B2) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B3) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B4) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B5) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B6) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B7) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_1) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_2) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_3) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_4) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_5) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_6) )


class QBit_Flux_3( Element_Base ):
    def __init__( self, origin, params, trans_in=None ):
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]
        self.B1_width = params[11]
        self.B1_height = params[12]
        self.B4_width = params[14]
        self.B5_width = params[15]
        self.B5_height = params[16]
        self.B6_width = params[17]
        self.B6_height = params[18]
        self.dCap = params[19]
        self.gap = params[20]        
        
        # calculated parameters
        self.B4_height = (self.gap - self.dCap - self.B6_width/2 - 4*self.w - self.jos1_a - self.jos2_a)/2
        self.B2_width = self.a + self.jos1_a/2 - self.jos2_a/2 - self.w
        self.B1_width = 2*self.w + self.B2_width
        self.B3_width = self.a + self.jos2_a/2 - self.jos1_a/2 - self.w
        self._alpha_1 = 0.5
        self._length1 = (self.b + 2*self.w + self.jos1_b)
        self._alpha_2 = 0.5
        self._length2 = (self.b + 2*self.w - 2*self.f2 - self.jos2_b) #length
        
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 - DPoint( 0,self.dCap + self.B1_height )
        self.p2 = self.p1 + DPoint( self.B1_width - self.w , 0 )
        self.p3 = self.p1 + DPoint( self.w, -self.w )
        self.p4 = self.p2 - DPoint( 0,self.f2 + self._length2*self._alpha_2 ) 
        self.p5 = self.p4 + DPoint( 2*self.w + self.jos2_a, -self.jos2_b )
        self.p6 = self.p5 - DPoint( 0,self.f2 + (1-self._alpha_2)*self._length2 )
        self.p7 = self.p6 - DPoint( self.w + self.B3_width,0 ) 
        self.p8 = self.p6 - DPoint( self.w + (self.B4_width + self.B3_width)/2,self.B4_height )
        self.p9 = self.p8 + DPoint( (self.B4_width - self.B5_width)/2, -self.B5_height )
        self.p10 = self.p9 + DPoint( (self.B5_width - self.B6_width)/2, -self.B6_height )
        
        self.B1 = klayout.db.DBox( self.p1, self.p1 + DPoint( self.B1_width,self.B1_height ) )
        self.B2 = klayout.db.DBox( self.p3, self.p2 )
        self.B3 = klayout.db.DBox( self.p7, self.p7 + DPoint( self.B3_width, self.w ) )
        self.B4 = klayout.db.DBox( self.p8, self.p8 + DPoint( self.B4_width, self.B4_height ) )
        self.B5 = klayout.db.DBox( self.p9, self.p9 + DPoint( self.B5_width, self.B5_height ) )
        self.B6 = klayout.db.DBox( self.p10, self.p10 + DPoint( self.B6_width, self.B6_height ) )

        self.poly_1 = self._make_polygon( self._length1*self._alpha_1, self.w, self.d1, self.f1, self.jos1_b )
        self.poly_1.transform( DCplxTrans( 1.0, 270, False, self.p1 ) )
        self.poly_2 = self._make_polygon( self._length1*(1-self._alpha_1), self.w, self.d1, self.f1, self.jos1_b )
        self.poly_2.transform( DCplxTrans( 1.0, 90, False, self.p7 ) )
        self.poly_3 = self._make_polygon( self._length2*self._alpha_2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_3.transform( DCplxTrans( 1.0, 270, False, self.p2 ) )
        self.poly_4 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_4.transform( DCplxTrans( 1.0, 90, False, self.p5 ) ) 
        self.poly_5 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_5.transform( DCplxTrans( 1.0, 270, False, self.p4 ) )           
        self.poly_6 = self._make_polygon( (1-self._alpha_2)*self._length2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_6.transform( DCplxTrans( 1.0, 90, False, self.p6 ) )               
        super().__init__( origin, trans_in )
        
    def _make_polygon( self, length, w, d, f, overlapping ):
        polygon = DSimplePolygon
        p1 = DPoint(0,0)
        p2 = p1 + DPoint( length,0 )
        p3 = p2 + DPoint( 0, w )
        p4 = p3 - DPoint( overlapping,0 )
        p5 = p4 - DPoint( 0, d )
        p6 = p5 - DPoint( f,0 )
        p7 = p6 + DPoint( 0, d )
        p8 = p1 + DPoint( 0,w )
        
        polygon = DSimplePolygon( [p1,p2,p3,p4,p5,p6,p7,p8] )
        return polygon
        
    def init_regions( self ):
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B1) ) 
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B2) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B3) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B4) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B5) )
        self.metal_region.insert( klayout.db.Box().from_dbox(self.B6) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_1) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_2) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_3) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_4) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_5) )
        self.metal_region.insert( SimplePolygon().from_dpoly(self.poly_6) )


class QBit_Flux_Сshunted( Element_Base ):
    def __init__( self, origin, params, trans_in=None ):
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]
        self.dCap = params[11]
        self.gap = params[12]        
        self.square_a = params[13]
        self.dSquares = params[14]
        self.alum_over = params[15]
        self.B1_width = params[16]
        
        # calculated parameters
        self.B2_width = self.a + self.jos1_a/2 - self.jos2_a/2 - self.w
        self.qbit_width = 2*self.w + self.B2_width
        self.B3_width = self.a + self.jos2_a/2 - self.jos1_a/2 - self.w
        self.B1_height = (self.dSquares - 2*self.w - self.b)/2
        self._alpha_1 = 0.5
        self._length1 = (self.b + 2*self.w + self.jos1_b)
        self._alpha_2 = 0.5
        self._length2 = (self.b + 2*self.w - 2*self.f2 - self.jos2_b) #length
        
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 - DPoint( 0, self.dCap + self.square_a )
        self.p2 = self.p1 + DPoint( self.square_a/2 - self.qbit_width/2, -(self.dCap + self.B1_height) )
        self.p3 = self.p2 + DPoint( self.qbit_width - self.w , 0 )
        self.p4 = self.p2 + DPoint( self.w, -self.w )
        self.p5 = self.p3 - DPoint( 0,self.f2 + self._length2*self._alpha_2 ) 
        self.p6 = self.p5 + DPoint( 2*self.w + self.jos2_a, -self.jos2_b )
        self.p7 = self.p6 - DPoint( 0,self.f2 + (1-self._alpha_2)*self._length2 )
        self.p8 = self.p7 - DPoint( self.w + self.B3_width,0 ) 
        self.p9 = self.p7 - DPoint( self.w + (self.qbit_width + self.B3_width)/2,self.B1_height )
        self.p10 = self.p1 - DPoint( 0, self.square_a + self.b + 2*self.B1_height + 2*self.w )
        
        self.SQ1 = klayout.db.DBox( self.p1, self.p1 + DPoint( self.square_a, self.square_a ) )
        self._B1p1 = self.p2 + DPoint( self.qbit_width/2 - self.B1_width/2,0 )
        self.B1 = klayout.db.DBox( self._B1p1, self._B1p1+ DPoint( self.B1_width,self.B1_height + self.alum_over ) )
        self.B2 = klayout.db.DBox( self.p4, self.p3 )
        self.B3 = klayout.db.DBox( self.p8, self.p8 + DPoint( self.B3_width, self.w ) )
        self._B4p1  = self.p9 + DPoint( self.qbit_width/2 - self.B1_width/2,-self.alum_over )
        self.B4 = klayout.db.DBox( self._B4p1, self._B4p1 + DPoint( self.B1_width, self.B1_height + self.alum_over ) )
        self.SQ2 = klayout.db.DBox( self.p10, self.p10 + DPoint( self.square_a, self.square_a ) )
    
        self.poly_1 = self._make_polygon( self._length1*self._alpha_1, self.w, self.d1, self.f1, self.jos1_b )
        self.poly_1.transform( DCplxTrans( 1.0, 270, False, self.p2 ) )
        self.poly_2 = self._make_polygon( self._length1*(1-self._alpha_1), self.w, self.d1, self.f1, self.jos1_b )
        self.poly_2.transform( DCplxTrans( 1.0, 90, False, self.p8 ) )
        self.poly_3 = self._make_polygon( self._length2*self._alpha_2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_3.transform( DCplxTrans( 1.0, 270, False, self.p3 ) )
        self.poly_4 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_4.transform( DCplxTrans( 1.0, 90, False, self.p6 ) ) 
        self.poly_5 = self._make_polygon( 2*self.jos2_b + self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_5.transform( DCplxTrans( 1.0, 270, False, self.p5 ) )           
        self.poly_6 = self._make_polygon( (1-self._alpha_2)*self._length2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_6.transform( DCplxTrans( 1.0, 90, False, self.p7 ) )               
        super().__init__( origin, trans_in )
        
    def _make_polygon( self, length, w, d, f, overlapping ):
        polygon = DSimplePolygon
        p1 = DPoint(0,0)
        p2 = p1 + DPoint( length,0 )
        p3 = p2 + DPoint( 0, w )
        p4 = p3 - DPoint( overlapping,0 )
        p5 = p4 - DPoint( 0, d )
        p6 = p5 - DPoint( f,0 )
        p7 = p6 + DPoint( 0, d )
        p8 = p1 + DPoint( 0,w )
        
        polygon = DSimplePolygon( [p1,p2,p3,p4,p5,p6,p7,p8] )
        return polygon
        
    def init_regions( self ):
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        # electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()
    
        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.SQ1) ) 
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B1) ) 
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B2) )
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B3) )
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B4) )
        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.SQ2) )
        
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_1) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_2) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_3) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_4) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_5) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_6) )

        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
    
    # overwritig parent method "place" to be able to draw on 2 different layers simultaneously
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, "photo" )
        
        # placing electron-beam litography
        super().place( cell, layer_el, "el" )
        

class QBit_Flux_Сshunted_3JJ( Element_Base ):
    def __init__( self, origin, params, trans_in=None ):
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]
        self.dCap = params[11]
        self.gap = params[12]        
        self.square_a = params[13]
        self.dSquares = params[14]
        self.alum_over = params[15]
        self.B1_width = params[16]
        
        # calculated parameters
        self.B2_width = self.a + self.jos1_a/2 - self.jos2_a/2 - self.w
        self.qbit_width = 2*self.w + self.B2_width
        self.B3_width = self.b - self.jos2_a/2 - self.jos1_a/2 - 2*self.w
        self.B1_height = (self.dSquares - 2*self.w - self.b)/2
        self._alpha_1 = 0.5
        self._length1 = (self.b + 2*self.w + self.jos1_b)
        self._alpha_2 = 0.5
        self._length2 = (self.b + 2*self.w - 2*self.f2 - self.jos2_b) #length
        
        self._length_right = (self.b + 2*self.w - 2*self.jos2_b)/3
        self.p0 = DPoint(0,0)
        self.p1 = self.p0 - DPoint( 0, self.dCap + self.square_a )
        self.p2 = self.p1 + DPoint( self.square_a/2 - self.qbit_width/2, -(self.dCap + self.B1_height) )
        self.p3 = self.p2 + DPoint( self.qbit_width - self.w , 0 )
        self.p4 = self.p2 + DPoint( self.w, -self.w )
        self.p5 = self.p3 + DPoint( 2*self.w + self.jos2_a, -(2*self._length_right + 2*self.jos2_b)) 
        self.p6 = self.p3 - DPoint( 0, self.b + 2*self.w)
        self.p7 = self.p6 - DPoint( self.B3_width,0 )
        self.p8 = self.p7 - DPoint( self.w, self.B1_height ) 
        self.p9 = self.p1 - DPoint( 0, self.square_a + self.b + 2*self.B1_height + 2*self.w )
        
        self.SQ1 = klayout.db.DBox( self.p1, self.p1 + DPoint( self.square_a, self.square_a ) )
        self._B1p1 = self.p2 + DPoint( (self.B2_width + 2*self.w)/2 - self.B1_width/2,0 )
        self.B1 = klayout.db.DBox( self._B1p1, self._B1p1+ DPoint( self.B1_width,self.B1_height + self.alum_over ) )
        self.B2 = klayout.db.DBox( self.p4, self.p3 )
        self.B3 = klayout.db.DBox( self.p7, self.p7 + DPoint( self.B3_width, self.w ) )
        self._B4p1  = self.p8 + DPoint( (self.B3_width + 2*self.w)/2 - self.B1_width/2,-self.alum_over )
        self.B4 = klayout.db.DBox( self._B4p1, self._B4p1 + DPoint( self.B1_width, self.B1_height + self.alum_over ) )
        self.SQ2 = klayout.db.DBox( self.p9, self.p9 + DPoint( self.square_a, self.square_a ) )
    
        self.poly_1 = self._make_polygon( self._length1*self._alpha_1, self.w, self.d1, self.f1, self.jos1_b )
        self.poly_1.transform( DCplxTrans( 1.0, 270, False, self.p2 ) )
        self.poly_2 = self._make_polygon( self._length1*(1-self._alpha_1), self.w, self.d1, self.f1, self.jos1_b )
        self.poly_2.transform( DCplxTrans( 1.0, 90, False, self.p7 ) )
        self.poly_3 = self._make_polygon( self._length_right + self.jos2_b, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_3.transform( DCplxTrans( 1.0, 270, False, self.p3 ) )
        
        self.poly_4 = self._make_polygon( self._length_right + self.jos2_b - self.f2, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_4.transform( DCplxTrans( 1.0, 270, True, self.p5 + DPoint(0,self._length_right + self.jos2_b - self.f2) ) ) 
        _poly_tmp = self._make_polygon( self._length_right + self.jos2_b - self.f2, self.w, self.d2, self.f2, self.jos2_b )
        _poly_tmp.transform( DCplxTrans( 1.0, 90, False, self.p5 + DPoint(0,self.jos2_b + self.f2) ) )
        _reg_tmp4 = Region()
        _reg_tmp4 .insert( SimplePolygon().from_dpoly(self.poly_4) ) 
        _reg_tmp = Region()
        _reg_tmp.insert( SimplePolygon().from_dpoly(_poly_tmp) )
        self._reg_tmp_to_metal = (_reg_tmp + _reg_tmp4).merged()
        
        
        self.poly_5 = self._make_polygon( self._length_right + self.jos2_b, self.w, self.d2, self.f2, self.jos2_b )
        self.poly_5.transform( DCplxTrans( 1.0, 90, True, self.p6 ) )                       
        super().__init__( origin, trans_in )
        
    def _make_polygon( self, length, w, d, f, overlapping ):
        polygon = DSimplePolygon
        p1 = DPoint(0,0)
        p2 = p1 + DPoint( length,0 )
        p3 = p2 + DPoint( 0, w )
        p4 = p3 - DPoint( overlapping,0 )
        p5 = p4 - DPoint( 0, d )
        p6 = p5 - DPoint( f,0 )
        p7 = p6 + DPoint( 0, d )
        p8 = p1 + DPoint( 0,w )
        
        polygon = DSimplePolygon( [p1,p2,p3,p4,p5,p6,p7,p8] )
        return polygon
    
    def init_regions( self ):
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        # electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()
    
        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.SQ1) ) 
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B1) ) 
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B2) )
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B3) )
        self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.B4) )
        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.SQ2) )
        
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_1) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_2) )
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_3) )
        self.metal_regions["el"] = self.metal_regions["el"] +  self._reg_tmp_to_metal
        self.metal_regions["el"].insert( SimplePolygon().from_dpoly(self.poly_5) )

        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
    
    # overwritig parent method "place" to be able to draw on 2 different layers simultaneously
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, region_name="photo" )
        
        # placing electron-beam litography
        super().place( cell, layer_el, region_name="el" )

class Test_Qubit_Box( Element_Base ):
    def __init__( self, origin, length, trans_in=None ):
        self.origin = origin
        self.length = length                      
        super().__init__( origin, trans_in )
    
    def init_regions( self ):
        ##TODO    
        # If using metal_regions, when we place klayout.db.Box().from_dbox() objects, the relative origin is not modified
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        #electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()

        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(klayout.db.DBox( klayout.db.DPoint(0, 0), klayout.db.DPoint(self.length, self.length) ) ) )
        
        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
        #self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(klayout.db.DBox( klayout.db.DPoint(0, 0), klayout.db.DPoint(self.length, self.length) ).move(klayout.db.Vector(self.origin.x, self.origin.y) ) ) )
    # Overwrite place from Element_Base
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, "photo" )
        
        # placing electron-beam litography
        super().place( cell, layer_el, "el" )
        
        
class QBit_UQ( Element_Base ):
    def __init__( self, origin, params, trans_in=None , Loop=True):
        
        
        self.params = params
        # Dimensions between the 2 pads
        self.pad_width = params[0]
        self.pad_length = params[1]
        # Distance between the 2 pads 
        self.d = params[2]
        # This version has a loop whose dimensions are precised here, There should also 
        # be a version without any loop, next thing to do
        self.loop_size_1 = params[3]
        self.loop_size_2 = params[4]
        # width of the metal line connected to the pads
        self.w_1 = params[5]
        # width of the U shape connected to the previous metal line 
        self.w_2 = params[6]
        # Dimensions of the Josephson junction 
        self.JJ_w_1 = params[7]
        self.JJ_w_2 = params[8]
        # Distance on the gds file between the 2 parts of the Josephson junction 
        self.JJ_d = params[9]  
        self.Loop = Loop
        
        # calculated parameters
        self.length_1=(self.d-self.JJ_w_2-self.loop_size_1-self.w_2)/2 +20e3
        self.length_2=(self.loop_size_1-self.JJ_d)/2 
        
        if self.Loop:
        
            self.p0 = DPoint(0,0)
            self.p1 = self.p0 + DPoint( self.pad_length, self.pad_width )
            self.p2 = self.p1 + DPoint( (-self.pad_length/2 - self.w_1/2), -20e3 )
            self.p3 = self.p2 + DPoint( self.w_1, self.length_1 )
            self.p4 = self.p3 + DPoint( -(self.w_1+self.loop_size_2)/2, 0 )
            self.p5 = self.p4 + DPoint( self.loop_size_2, self.w_2) 
            self.p6 = self.p4 + DPoint( 0, self.w_2)
            self.p7 = self.p6 + DPoint( self.w_2, self.length_2 )
            self.p8 = self.p5 - DPoint( self.w_2, 0) 
            self.p9 = self.p8 + DPoint( self.w_2, self.length_2 )
            self.p10 = self.p7 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p11 = self.p10 + DPoint( self.JJ_w_1, self.length_2)
            self.p12 = self.p9 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p13 = self.p12 + DPoint( self.JJ_w_1, self.length_2)
            self.p14 = self.p6 + DPoint( 0, self.loop_size_1)
            self.p15 = self.p14 + DPoint( self.loop_size_2, self.JJ_w_2)
            self.p16 = self.p2 + DPoint( 0, self.length_1 + self.w_2 + self.loop_size_1 + self.JJ_w_2)
            self.p17 = self.p16 + DPoint( self.w_1, self.length_1 )
            self.p18 = self.p0 + DPoint( 0, self.d + self.pad_width )
            self.p19 = self.p18 + DPoint( self.pad_length, self.pad_width )
            
            self.Photo1 = klayout.db.DBox( self.p0, self.p1)
            self.Photo2 = klayout.db.DBox( self.p18, self.p19)
            
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal3 = klayout.db.DBox( self.p6, self.p7)
            self.metal4 = klayout.db.DBox( self.p8, self.p9)
            self.metal5 = klayout.db.DBox( self.p10, self.p11)
            self.metal6 = klayout.db.DBox( self.p12, self.p13)
            self.metal7 = klayout.db.DBox( self.p14, self.p15)
            self.metal8 = klayout.db.DBox( self.p16, self.p17)
        
        else:
        
            self.p0 = DPoint(0,0)
            self.p1 = self.p0 + DPoint( self.pad_length, self.pad_width )
            self.p2 = self.p1 + DPoint( (-self.pad_length/2 - self.w_1/2), -20e3 )
            self.p3 = self.p2 + DPoint( self.w_1, self.length_1 )
            self.p4 = self.p3 + DPoint( -self.w_1, 0 )
            self.p5 = self.p4 + DPoint( self.loop_size_2/2 + self.w_1/2, self.w_2) 
            self.p6 = self.p5 - DPoint( self.w_2, 0) 
            self.p7 = self.p6 + DPoint( self.w_2, self.length_2 )
            self.p8 = self.p7 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p9 = self.p8 + DPoint( self.JJ_w_1, self.length_2)
            self.p10 = self.p2 + DPoint(0, self.length_1 + self.w_2 + self.loop_size_1)
            self.p11 = self.p10 + DPoint( self.loop_size_2/2 + self.w_1/2, self.JJ_w_2)
            self.p12 = self.p10 + DPoint( 0 , self.JJ_w_2)
            self.p13 = self.p12 + DPoint( self.w_1, self.length_1 )
            self.p14 = self.p0 + DPoint( 0, self.d + self.pad_width )
            self.p15 = self.p14 + DPoint( self.pad_length, self.pad_width )
            
            self.Photo1 = klayout.db.DBox( self.p0, self.p1)
            self.Photo2 = klayout.db.DBox( self.p14, self.p15)
            
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal3 = klayout.db.DBox( self.p6, self.p7)
            self.metal4 = klayout.db.DBox( self.p8, self.p9)
            self.metal5 = klayout.db.DBox( self.p10, self.p11)
            self.metal6 = klayout.db.DBox( self.p12, self.p13)
            
        super().__init__( origin, trans_in )
        
        self.start = self.connections[0]
        self.end = self.connections[1]

    
    def init_regions( self ):
        
        self.connections = [DPoint(0,0),DPoint(0, self.d + 2*self.pad_width)]
        
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        # electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()
    
        if self.Loop:
        
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal2) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal3) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal4) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal5) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal6) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal7) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal8) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo2) )
        
        else:
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal2) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal3) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal4) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal5) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal6) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo2) )
        

        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
    
    # overwritig parent method "place" to be able to draw on 2 different layers simultaneously
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, region_name="photo" )  
        # placing electron-beam litography
        super().place( cell, layer_el, region_name="el" )
        
        
class IBM_Transmon( Element_Base ):
    def __init__( self, origin_ini, params, trans_in=None , Loop=True):
        
        
        self.params = params
        # Dimensions between the 2 pads
        self.pad_width = params[0]
        self.pad_length = params[1]
        # Distance between the 2 pads 
        self.d = params[2]
        # This version has a loop whose dimensions are precised here, There should also 
        # be a version without any loop, next thing to do
        self.loop_size_1 = params[3]
        self.loop_size_2 = params[4]
        # width of the metal line connected to the pads
        self.w_1 = params[5]
        # width of the U shape connected to the previous metal line 
        self.w_2 = params[6]
        # Dimensions of the Josephson junction 
        self.JJ_w_1 = params[7]
        self.JJ_w_2 = params[8]
        # Distance on the gds file between the 2 parts of the Josephson junction 
        self.JJ_d = params[9]  
        self.Loop = Loop
        # Caracteristic of the connection between the qubit and the CPW
        self.connections_spacing = params[10]
        self.connectors_Length = params[11]
        self.connectors_width = params[12]
        self.CPW_width = params[13]
        
        # calculated parameters
        self.length_1=(self.d-self.JJ_w_2-self.loop_size_1-self.w_2)/2 +20e3
        self.length_2=(self.loop_size_1-self.JJ_d)/2 
        self.origin = origin_ini - DPoint(self.pad_length/2, self.d + 2*self.pad_width + self.connections_spacing)
        
        if self.Loop:
        
            self.p0 = DPoint(0,0)
            self.p1 = self.p0 + DPoint( self.pad_length, self.pad_width )
            self.p2 = self.p1 + DPoint( (-self.pad_length/2 - self.w_1/2), -20e3 )
            self.p3 = self.p2 + DPoint( self.w_1, self.length_1 )
            self.p4 = self.p3 + DPoint( -(self.w_1+self.loop_size_2)/2, 0 )
            self.p5 = self.p4 + DPoint( self.loop_size_2, self.w_2) 
            self.p6 = self.p4 + DPoint( 0, self.w_2)
            self.p7 = self.p6 + DPoint( self.w_2, self.length_2 )
            self.p8 = self.p5 - DPoint( self.w_2, 0) 
            self.p9 = self.p8 + DPoint( self.w_2, self.length_2 )
            self.p10 = self.p7 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p11 = self.p10 + DPoint( self.JJ_w_1, self.length_2)
            self.p12 = self.p9 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p13 = self.p12 + DPoint( self.JJ_w_1, self.length_2)
            self.p14 = self.p6 + DPoint( 0, self.loop_size_1)
            self.p15 = self.p14 + DPoint( self.loop_size_2, self.JJ_w_2)
            self.p16 = self.p2 + DPoint( 0, self.length_1 + self.w_2 + self.loop_size_1 + self.JJ_w_2)
            self.p17 = self.p16 + DPoint( self.w_1, self.length_1 )
            self.p18 = self.p0 + DPoint( 0, self.d + self.pad_width )
            self.p19 = self.p18 + DPoint( self.pad_length, self.pad_width )
            self.p20 = self.p18 + DPoint(self.pad_length/2 - self.CPW_width/2 -self.connections_spacing , self.pad_width)
            self.p21 = self.p20 + DPoint(-self.connectors_width, self.connectors_Length)
            self.p22 = self.p18 + DPoint(self.pad_length/2 + self.CPW_width/2 +self.connections_spacing , self.pad_width)
            self.p23 = self.p22 + DPoint(self.connectors_width, self.connectors_Length)
            
            self.Photo1 = klayout.db.DBox( self.p0, self.p1)
            self.Photo2 = klayout.db.DBox( self.p18, self.p19)
            self.Photo3 = klayout.db.DBox( self.p20, self.p21)
            self.Photo4 = klayout.db.DBox( self.p22, self.p23)
            
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal3 = klayout.db.DBox( self.p6, self.p7)
            self.metal4 = klayout.db.DBox( self.p8, self.p9)
            self.metal5 = klayout.db.DBox( self.p10, self.p11)
            self.metal6 = klayout.db.DBox( self.p12, self.p13)
            self.metal7 = klayout.db.DBox( self.p14, self.p15)
            self.metal8 = klayout.db.DBox( self.p16, self.p17)
        
        else:
        
            self.p0 = DPoint(0,0)
            self.p1 = self.p0 + DPoint( self.pad_length, self.pad_width )
            self.p2 = self.p1 + DPoint( (-self.pad_length/2 - self.w_1/2), -20e3 )
            self.p3 = self.p2 + DPoint( self.w_1, self.length_1 )
            self.p4 = self.p3 + DPoint( -self.w_1, 0 )
            self.p5 = self.p4 + DPoint( self.loop_size_2/2 + self.w_1/2, self.w_2) 
            self.p6 = self.p5 - DPoint( self.w_2, 0) 
            self.p7 = self.p6 + DPoint( self.w_2, self.length_2 )
            self.p8 = self.p7 - DPoint( (self.w_2 + self.JJ_w_1)/2, 0)
            self.p9 = self.p8 + DPoint( self.JJ_w_1, self.length_2)
            self.p10 = self.p2 + DPoint(0, self.length_1 + self.w_2 + self.loop_size_1)
            self.p11 = self.p10 + DPoint( self.loop_size_2/2 + self.w_1/2, self.JJ_w_2)
            self.p12 = self.p10 + DPoint( 0 , self.JJ_w_2)
            self.p13 = self.p12 + DPoint( self.w_1, self.length_1 )
            self.p14 = self.p0 + DPoint( 0, self.d + self.pad_width )
            self.p15 = self.p14 + DPoint( self.pad_length, self.pad_width )
            self.p16 = self.p14 + DPoint(self.pad_length/2 - self.CPW_width/2 -self.connections_spacing , self.pad_width)
            self.p17 = self.p16 + DPoint(-self.connectors_width, self.connectors_Length)
            self.p18 = self.p14 + DPoint(self.pad_length/2 + self.CPW_width/2 +self.connections_spacing , self.pad_width)
            self.p19 = self.p18 + DPoint(self.connectors_width, self.connectors_Length)
            
            self.Photo1 = klayout.db.DBox( self.p0, self.p1)
            self.Photo2 = klayout.db.DBox( self.p14, self.p15)
            self.Photo3 = klayout.db.DBox( self.p16, self.p17)
            self.Photo4 = klayout.db.DBox( self.p18, self.p19)
            
            self.metal1 = klayout.db.DBox( self.p2, self.p3)
            self.metal2 = klayout.db.DBox( self.p4, self.p5)
            self.metal3 = klayout.db.DBox( self.p6, self.p7)
            self.metal4 = klayout.db.DBox( self.p8, self.p9)
            self.metal5 = klayout.db.DBox( self.p10, self.p11)
            self.metal6 = klayout.db.DBox( self.p12, self.p13)
            
        super().__init__( self.origin, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
    
    def init_regions( self ):
        
        self.connections = [DPoint(0,0),DPoint(0, self.d + 2*self.pad_width)]
        
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        # electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()
    
        if self.Loop:
        
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal2) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal3) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal4) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal5) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal6) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal7) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal8) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo2) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo3) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo4) )
        
        else:
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal1) ) 
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal2) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal3) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal4) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal5) )
            self.metal_regions["el"].insert( klayout.db.Box().from_dbox(self.metal6) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo2) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo3) )
            self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(self.Photo4) )
        

        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
    
    # overwritig parent method "place" to be able to draw on 2 different layers simultaneously
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, region_name="photo" )  
        # placing electron-beam litography
        super().place( cell, layer_el, region_name="el" )        


class Test_Qubit_Box( Element_Base ):
    def __init__( self, origin, length, trans_in=None ):
        self.origin = origin
        self.length = length                      
        super().__init__( origin, trans_in )
    
    def init_regions( self ):
        ##TODO    
        # If using metal_regions, when we place klayout.db.Box().from_dbox() objects, the relative origin is not modified
        # photolitography regions 
        self.metal_regions["photo"] = Region()
        self.empty_regions["el"] = Region()
        #electron-beam litography regions
        self.metal_regions["el"] = Region() 
        self.empty_regions["photo"] = Region()

        self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(klayout.db.DBox( klayout.db.DPoint(0, 0), klayout.db.DPoint(self.length, self.length) ) ) )
        
        for tags in ['photo', 'el']:
            self.metal_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
            self.empty_regions[tags].move(klayout.db.Vector(self.origin.x, self.origin.y))
        #self.metal_regions["photo"].insert( klayout.db.Box().from_dbox(klayout.db.DBox( klayout.db.DPoint(0, 0), klayout.db.DPoint(self.length, self.length) ).move(klayout.db.Vector(self.origin.x, self.origin.y) ) ) )
    # Overwrite place from Element_Base
    def place( self, cell, layer_photo, layer_el ):
        # placing photolitography
        super().place( cell, layer_photo, "photo" )
        
        # placing electron-beam litography
        super().place( cell, layer_el, "el" )