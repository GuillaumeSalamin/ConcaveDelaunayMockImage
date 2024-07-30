# Last file who compile together result of all PartialImage.py result to build the final Image
import numpy as np
import matplotlib.pyplot as plt
import yaml
import ctypes
import glfw
import pickle
from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader
import os
os.environ["SDL_VIDEO_X11_FORCE_EGL"] = "1"

class App:
    def __init__(self, data, width=640, height=400,alpha=1):
        self.width = width
        self.height = height

        # Initialize GLFW
        if not glfw.init():
            raise RuntimeError("Failed to initialize GLFW")

        # Create an off-screen context
        glfw.window_hint(glfw.VISIBLE, glfw.FALSE)
        glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
        glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
        glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
        self.offscreen_window = glfw.create_window(self.width, self.height, "", None, None)
        if not self.offscreen_window:
            glfw.terminate()
            raise RuntimeError("Failed to create GLFW off-screen window")

        # Make the off-screen OpenGL context current
        glfw.make_context_current(self.offscreen_window)

        # Initialize OpenGL
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        # Initialize the triangle
        self.triangle = MultripleTriangle(data,alpha)
        self.shader = self.triangle.createShader()
        glUseProgram(self.shader)

    def render_to_image(self):
        # Clear the screen
        glClear(GL_COLOR_BUFFER_BIT)

        # Render the triangle
        glUseProgram(self.shader)
        glBindVertexArray(self.triangle.vao)
        glDrawArrays(GL_TRIANGLES,0,self.triangle.num_triangles * 3)
        #glDrawElements(GL_TRIANGLES, self.triangle.num_triangles * 3, GL_UNSIGNED_INT, None)

        # Read the pixels from the framebuffer
        #glReadBuffer(GL_FRONT)
        #pixels = glReadPixels(0, 0, self.width, self.height, GL_RGBA, GL_UNSIGNED_BYTE)
        #image = np.frombuffer(pixels, dtype=np.uint8).reshape(self.height, self.width, 4)
        #image = np.flipud(image)  # Flip the image vertically
        pixels = np.empty((self.height, self.width, 4), dtype=np.uint8)
        #glReadPixels(0, 0, self.width, self.height, GL_RGBA, GL_FLOAT, pixels)
        glReadPixels(0, 0, self.width, self.height, GL_RGBA, GL_UNSIGNED_BYTE, pixels)

        return pixels[:, :, 0]#image

    def destroy(self):
        self.triangle.destroy()
        glDeleteProgram(self.shader)
        glfw.destroy_window(self.offscreen_window)
        glfw.terminate()

class MultripleTriangle:
    def __init__(self, data,alpha):

        

        self.num_triangles = len(data)//3

        # Create vertex data for all triangles
        self.vertices = []
        self.indices = []

        for ii in range(self.num_triangles):
            base_index = ii * 3

            self.vertices.extend([
                data[3*ii][0], data[3*ii][1], 0.0, data[3*ii][2]*alpha,
                data[3*ii+1][0], data[3*ii+1][1], 0.0, data[3*ii+1][2]*alpha,
                data[3*ii+2][0], data[3*ii+2][1], 0.0, data[3*ii+2][2]*alpha 
            ])

            self.indices.extend([
                base_index, base_index + 1, base_index + 2
            ])

        self.vertices = np.array(self.vertices, dtype=np.float32)
        self.indices = np.array(self.indices, dtype=np.uint32)

        # Create VAO, VBO, and EBO
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)

        self.vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, self.vertices.nbytes, self.vertices, GL_STATIC_DRAW)

        self.ebo = glGenBuffers(1)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ebo)
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, self.indices.nbytes, self.indices, GL_STATIC_DRAW)

        glEnableVertexAttribArray(0)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 16, ctypes.c_void_p(0))

        glEnableVertexAttribArray(1)
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 16, ctypes.c_void_p(12))

        self.shader = self.createShader()

    def createShader(self):
        vertex_shader = """
        #version 330 core

        layout(location = 0) in vec3 aPos; // The position variable has attribute position 0
        layout(location = 1) in float aAlpha; // The alpha variable has attribute position 1

        out float vertexAlpha; // Specify an alpha output to the fragment shader

        void main()
        {
            gl_Position = vec4(aPos, 1.0); // Assign the position of the vertex
            vertexAlpha = aAlpha; // Pass the alpha value to the fragment shader
        }
        """

        fragment_shader = """
        #version 330 core

        in float vertexAlpha; // Input alpha from the vertex shader

        out vec4 FragColor; // Output fragment color

        void main()
        {
            FragColor = vec4(1.0, 1.0, 1.0, vertexAlpha); // Set the fragment color to white with varying alpha
        }
        """

        shader = compileProgram(
            compileShader(vertex_shader, GL_VERTEX_SHADER),
            compileShader(fragment_shader, GL_FRAGMENT_SHADER)
        )
        return shader

    def destroy(self):
        glDeleteVertexArrays(1, (self.vao,))
        glDeleteBuffers(1, (self.vbo,))
        glDeleteBuffers(1, (self.ebo,))

def adapt_satured_image(matrix,alpha):
    ans = np.zeros((1000,1000))
    for ii in range(1000):
        for jj in range(1000):
            tmp = matrix[ii,jj]
            if tmp>=256*alpha:
                ans[ii,jj] = 0
            else:
                ans[ii,jj] = tmp/alpha
    return ans


if __name__ == "__main__":
    with open("paraImage.yaml") as stream:
        try:
            para = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    triangleFolder = para["triangleFolder"]
    ImagePath = para["ImagePath"]


    matrix = np.zeros((1000,1000))
    for ii in range(100):#[10,11,12,13,14,15,16]:
        print(f'{ii+1}/100')
        file_path = f'{ImagePath}/Partial_Image_{ii}_data.pickle'
        with open(file_path, 'rb') as file:
            # Deserialize and retrieve the variable from the file
           data = pickle.load(file)

        #print(f'The variable data_{ii} has been loaded successfully.')

        app = App(data,width=1000, height=1000,alpha=1)
        image = app.render_to_image()
        matrix+=image
        
        del image
        app.destroy()

        app = App(data,width=1000, height=1000,alpha=256)
        image = adapt_satured_image(app.render_to_image(),256)
        matrix+=image
        del data
        del image


    plt.matshow(matrix, cmap='gray')  # Display only alpha channel
    plt.axis('off')  # Turn off axes
    plt.title('Alpha Channel')
    plt.gca().invert_yaxis()
    plt.savefig(f'{ImagePath}/img_linear.png')
    plt.show()

    plt.matshow(matrix, cmap='gray',norm='log')  # Display only alpha channel
    plt.axis('off')  # Turn off axes
    plt.title('Alpha Channel')
    plt.gca().invert_yaxis()
    #plt.gca().set_facecolor('black')
    plt.savefig(f'{ImagePath}/img_log.png')
    plt.show()

    file_path_data = f'{ImagePath}/matrix.pickle' 
    with open(file_path_data, 'wb') as file:
        # Serialize and write the variable to the file
        pickle.dump(matrix, file)

