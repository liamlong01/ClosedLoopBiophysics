from glumpy import app, gloo, gl, glm

import numpy as np
from numpy.matlib import repmat

class NeuronRender(object):

    def __init__(self, cell, conn, elec = None):
        self.cell = cell
        self.elec = elec

        self.conn = conn

        self.maxcolor = 0.0
        self.mincolor = 0.0

    def show(self):


        cubesize = 1e-2
        elecsize = 2e-2

        window = app.Window()


        vertex = """
        uniform mat4   u_model;         // Model matrix
        uniform mat4   u_view;          // View matrix
        uniform mat4   u_projection;    // Projection matrix
        attribute vec3 a_position;      // Vertex position
        attribute vec4 a_color;

        varying vec4 v_color;           //varying imemcolor
        

        void main()
        {
            v_color = a_color;
            gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
        } """



        fragment = """
        varying vec4 v_color;
        void main()
        {
            gl_FragColor = v_color;
        } """

        xnorm = np.max(self.cell.xstart)
        ynorm = np.max(self.cell.ystart)
        znorm = np.max(self.cell.zstart)
        norm = np.max([xnorm,ynorm,znorm])/3 

        cubes = []
        elecs = []
        print(len(self.cell.xstart))
        for i in range(len(self.cell.xstart)):

            centre1 = np.array([self.cell.xstart[i]/norm , self.cell.ystart[i]/norm, self.cell.zstart[i]/norm])
            centre2 = np.array([self.cell.xend[i]/norm , self.cell.yend[i]/norm, self.cell.zend[i]/norm])
            mid = np.mean([centre1,centre2], axis = 0)
            # print(mid)
         
            print("seg:", mid)
            V = np.zeros(8, [("a_position", np.float32, 3),("a_color", np.float32, 4)])
            V["a_position"] = [mid + [cubesize,cubesize,cubesize], mid+[0	, +cubesize	, +cubesize	], mid+[0	,0	, +cubesize	], mid+[ +cubesize	,0	, +cubesize	], mid+[ +cubesize	,0,0	], mid+[ +cubesize	, +cubesize	,0	], mid+[0	, +cubesize	,0	], mid+[0	,0	,0	]] 
            V["a_color"]    = [[0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 0, 1], [1, 1, 0, 1], [1, 1, 1, 1], [1, 0, 1, 1], [1, 0, 0, 1]]

            I = np.array([0,1,2, 0,2,3,  0,3,4, 0,4,5,  0,5,6, 0,6,1,
              1,6,7, 1,7,2,  7,4,3, 7,3,2,  4,7,6, 4,6,5], dtype=np.uint32)

            V = V.view(gloo.VertexBuffer)
            I = I.view(gloo.IndexBuffer)

            cube = gloo.Program(vertex, fragment)
            cube["a_position"] = V["a_position"]
            cube["a_color"] = V["a_color"]

            view = np.eye(4,dtype=np.float32)
            model = np.eye(4,dtype=np.float32)
            projection = np.eye(4,dtype=np.float32)
            glm.translate(view, 0,0,-5)
            cube['u_model'] = model
            cube['u_view'] = view
            cube['u_projection'] = projection

            cubes.append(cube)
        
        if self.elec:

            
            cubesize = elecsize

            for i in range(round(len(self.elec.X.flatten())/123)):
                i = i*123
         
                mid = np.array( [self.elec.X.flatten()[i], self.elec.Y.flatten()[i], self.elec.Z.flatten()[i]])/norm
                
                V = np.zeros(8, [("a_position", np.float32, 3),("a_color", np.float32, 4)])
                V["a_position"] = [mid + [cubesize,cubesize,cubesize], mid+[0   , +cubesize , +cubesize ], mid+[0   ,0  , +cubesize ], mid+[ +cubesize  ,0  , +cubesize ], mid+[ +cubesize  ,0,0    ], mid+[ +cubesize  , +cubesize ,0  ], mid+[0   , +cubesize ,0  ], mid+[0   ,0  ,0  ]] 
                V["a_color"]    = repmat([0,1,0,1], 8,1).tolist()

                I = np.array([0,1,2, 0,2,3,  0,3,4, 0,4,5,  0,5,6, 0,6,1,
                  1,6,7, 1,7,2,  7,4,3, 7,3,2,  4,7,6, 4,6,5], dtype=np.uint32)

                V = V.view(gloo.VertexBuffer)
                I = I.view(gloo.IndexBuffer)

                cube = gloo.Program(vertex, fragment)
                cube["a_position"] = V["a_position"]
                cube["a_color"] = V["a_color"]

                view = np.eye(4,dtype=np.float32)
                model = np.eye(4,dtype=np.float32)
                projection = np.eye(4,dtype=np.float32)
                glm.translate(view, 0,0,-5)
                cube['u_model'] = model
                cube['u_view'] = view
                cube['u_projection'] = projection

                elecs.append(cube)
        print("cubes loaded")








        self.phi, self.theta = 0,0
        self.cubes = cubes
        self.elecs =elecs


        @window.event
        def on_resize(width, height):
            ratio = width / float(height)
            for cube in self.cubes:
               cube['u_projection'] = glm.perspective(100, ratio, 2.0, 100.0)
            for cube in self.elecs:
               cube['u_projection'] = glm.perspective(100, ratio, 2.0, 100.0)
         


        @window.event
        def on_draw(dt):
           
            window.clear()
            for cube in self.cubes:
                cube.draw(gl.GL_TRIANGLES, I)
            for cube in self.elecs:
                cube.draw(gl.GL_TRIANGLES, I)
             

            # Make cube rotate
            self.theta += 10 # degrees
            self.phi += 0.5 # degrees
         
            model = np.eye(4, dtype=np.float32)
            glm.rotate(model, self.theta, 0, 1, 0)
            
            
            colordata = self.conn.recv()
            self.maxcolor = np.max([self.maxcolor, np.max(colordata)])
            self.mincolor = np.min([self.mincolor, np.min(colordata)])
            for cube, col in zip(self.cubes, colordata):
           
                cube['u_model'] = model
                cube['a_color'] = repmat([(col-self.mincolor)/(self.maxcolor-self.mincolor), 0, 1 - ((col-self.mincolor)/(self.maxcolor-self.mincolor)), .3], 8, 1).tolist()

     
            for cube in self.elecs:
               
                cube['u_model'] = model
                cube['a_color'] = repmat([0,1,0,1], 8, 1).tolist()


        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)

        app.run()

def start(celldata, conn, elec):
    print(elec)
    renderer = NeuronRender(celldata, conn, elec=elec)

    renderer.show()


if __name__ == '__main__':

    from LFPylite import CellData
    from aberraAxon import MyelinatedCell 

    mcell = MyelinatedCell()

    mcell.loadcell('L1_DAC_bNAC219_1')  

    renderer = NeuronRender(CellData(mcell))

    renderer.show()
