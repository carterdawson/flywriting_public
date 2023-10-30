"""
A flight drawing class
"""
import numpy as np
import math

def sort_second(val):
    return val[1]


class Drawing:

    max_points = 0

    def __init__(self, str_p, wpt_o, scale_x=1.0, scale_y=1.0, deg_rotate=0, deg_cull_angle=0.0):
        self.str_points = str_p
        self.wpt_origin = wpt_o
        self.scale_x = scale_x
        self.scale_y = scale_y
        self.num_culled = 0
        self.lat_scale = 0.000284757
        self.lon_scale = 0.000388889
        self.deg_rotate = deg_rotate

        #parse my points
        # npt_source : original points, as parsed
        # npt_full   : npt_source, scaled and rotated
        # npt_points : npt_full, culled (maybe)
        # npt_culled : any culled points
        self.npt_source = self.npt_full = self.npt_points = np.array(self._parse_points())

        #scale the X points
        self.npt_full[:,0] = self.scale_x * self.npt_full[:,0]
        #and the Y points
        self.npt_full[:,1] = -self.scale_y * self.npt_full[:,1]

        #if we need to rotate...
        if self.deg_rotate != 0:
            #create a rotation matrix
            theta = np.radians(self.deg_rotate)
            c, s = np.cos(theta), np.sin(theta)
            R = np.array(((c,-s), (s,c)))
            #... and apply it to npt_full
            self.npt_full = (self.npt_full.dot(R))

        #if there is a cull points angle...
        self.deg_cull_angle = deg_cull_angle
        if deg_cull_angle > 0:
            #cull based on that angle... this sets npt_culled and overwrites npt_points
            self._cull_points_angle(deg_cull_angle)
        else:
            #cull based on number of points... this sets npt_culled and overwrites npt_points
            self._cull_points_count()



        #calculate the center and the origin
        x_min_max = (self.npt_full[:,0].min(),self.npt_full[:,0].max())
        y_min_max = (self.npt_full[:,1].min(),self.npt_full[:,1].max())

        x_center = (x_min_max[1]-x_min_max[0])/2+x_min_max[0]
        y_center = (y_min_max[1]-y_min_max[0])/2+y_min_max[0]
        self.pt_center = (x_center, y_center)
        self.pt_origin = self.npt_full[0]

    def _cull_points_angle(self, deg_cull_angle):
        #create a culling mask
        num_points = np.size(self.npt_full, 0)
        self.culling_mask = np.ones(num_points, dtype=bool)
        
        #calculate the angle that we need to cull by
        rad_cull_angle = self.deg_cull_angle * 2 * math.pi / 360.0

        #initialize idx_a, b, and c
        idx_a = 0
        idx_b = 1
        idx_c = 2

        #while I have not gotten to the end...
        while idx_c < num_points:
            #calculate the angles for ab and ac, and the difference between them
            rad_ab = self._rad_between(self.npt_source[idx_a], self.npt_source[idx_b])
            rad_ac = self._rad_between(self.npt_source[idx_a], self.npt_source[idx_c])
            rad_diff = abs(rad_ac - rad_ab)
            #if that difference is sufficiently small
            if (rad_diff < rad_cull_angle):
                #cull the one before C (which, in the trivial case, is B)
                self.culling_mask[idx_c-1] = False
                #increment idx_c
                idx_c += 1
                #we're now proceeding with A and B the same, but C advanced by 1
            #otherwise...
            else:
                #not culling here, so reset a to the current c-1, and b and c accordingly
                idx_a = idx_c-1
                idx_b = idx_a + 1
                idx_c = idx_b + 1
                #we're now proceeding with A where B was, B one above that, and C one above B
        # we've now gone through all the points, and culled any whose angles are smaller than
        # our reference

        #apply our culling mask
        self.npt_points = self.npt_full[self.culling_mask]
        self.npt_culled = self.npt_full[np.logical_not(self.culling_mask)]
        self.num_culled = np.size(self.npt_culled, 0)

    def _cull_points_count(self):
        #if I have a max_points...
        if self.max_points > 0:
            #and I have more points than that...
            num_points = np.size(self.npt_points, 0)
            num_to_cull = num_points - self.max_points
            if num_to_cull > 0:
                #get an iterator
                #it = np.nditer(self.npt_points, flags=['f_index'])
                #for each point
                angles = []
                for i in range(0,num_points):
                    #if I'm almost at the end...
                    if i > (num_points-3):
                        #we're done
                        break
                    #get this point, and the next two
                    a = self.npt_points[i]
                    b = self.npt_points[i+1]
                    c = self.npt_points[i+2]

                    #compute the angle for A and C
                    ac = self._rad_between(a,c)
                    #compute the angle for A and B
                    ab = self._rad_between(a,b)
                    #save the index of this point along with the angle difference
                    angles.append((i+1, math.fabs(ac-ab)))
                #end for each point (mostly)
                #sort the resulting list in ascending order by the angle difference
                angles.sort(key=sort_second)
                #delete everything past the number of points I need to cull
                angles = np.array(angles[0:num_to_cull], dtype=int)
                #create a mask from that
                mask = np.ones(num_points, dtype=bool)
                for a in angles:
                    mask[a[0]] = False
                #save my new npt_points as masked by the above list
                #print(f"before:{self.npt_points}")
                self.npt_points = self.npt_full[mask]
                self.npt_culled = self.npt_full[np.logical_not(mask)]
                self.num_culled = num_to_cull
                #print(f"after:{self.npt_points}")

    def _rad_between(self, pt_a, pt_b):
        line = pt_b - pt_a
        return math.atan(line[0]/line[1])


    def get_x(self, culled=True):
        if culled:
            return self.npt_points[:,0]
        else:
            return self.npt_full[:,0]

    def get_y(self, culled=True):
        if culled:
            return self.npt_points[:,1]
        else:
            return self.npt_full[:,1]

    def get_lat(self, culled=True):
        return self.lat_from_y(self.get_y(culled))
    
    def get_lon(self, culled=True):
        return self.lon_from_x(self.get_x(culled))



    def set_scale(self, lat_scale, lon_scale):
        self.lat_scale = lat_scale
        self.lon_scale = lon_scale

    def _parse_points(self):
        apt_points = []
        #get all of the points
        astr_points = self.str_points.split("|")
        #for each point in our string...
        for str_point in astr_points:
            #convert it into an array of strings
            astr_point = str_point.split(",")
            #and then into integers, adding it into our list of points
            apt_points.append((int(astr_point[0]), int(astr_point[1])))

        return apt_points

    def lat_from_y(self, y):
        lat_scale = self.lat_scale
        return self.wpt_origin[0] + lat_scale*(y-self.pt_origin[1])

    def lon_from_x(self, x):
        lon_scale = self.lon_scale
        return self.wpt_origin[1] + lon_scale*(x-self.pt_origin[0])

    def wpt_from_pt(self, pt):
        return (self.lat_from_y(pt[1]), self.lon_from_x(pt[0]))



    
if __name__ == "__main__":


    test_drawing = Drawing("1,0|2,0|3,0|4,1|5,0|6,0|7,0",
        (36.81, -121.96),
        deg_cull_angle=10.0
    )

    Drawing.max_points = 5


    test_drawing = Drawing(
        "1018,254|981,267|947,285|917,301|893,318|868,337|844,359",
        (36.81, -121.96),
        0.5, 0.5,
        -30
    )

    print(f"Full: {test_drawing.npt_full}")
    print(f"Culled: {test_drawing.npt_culled}")
    print(f"current: {test_drawing.npt_points}")


    Drawing.max_points = 0


    penguin = Drawing(
        "150,1007|202,1017|254,1018|299,1010|343,996|392,985|449,979|485,982|510,989|519,1000|510,1007|490,1018|450,1024|411,1027|368,1024|325,1017|278,1000|240,978|204,942|181,908|163,868|152,820|146,765|147,713|155,649|165,601|179,546|195,491|210,446|227,411|244,400|262,418|273,453|281,488|288,526|294,572|301,617|306,657|310,704|310,750|305,795|299,833|288,866|274,903|255,932|246,942|234,946|223,937|217,921|212,878|211,822|213,780|215,742|214,702|210,669|201,636|196,586|194,552|196,521|197,470|206,432|221,381|237,333|264,270|295,222|324,187|363,169|397,164|433,172|475,181|518,177|561,184|582,196|560,203|535,207|503,216|471,229|444,245|428,261|421,292|439,320|465,358|485,409|492,454|497,532|498,596|497,648|494,712|490,750|488,792|485,831",
        (37.70338, -121.668091)
    )
    print(penguin.pt_center)
    y = penguin.get_y()

    #lat = penguin.lat_from_y(y)
    #print(lat)
    #lat2 = penguin.get_lat()
    #print(lat2)
    



