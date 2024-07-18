import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLUT import *
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib
import os
from pyvirtualdisplay import Display

# Start virtual display
display = Display(visible=0, size=(800, 600))
display.start()


def draw_triangles_with_opacity(vertices, save_image=False, image_filename="output.jpg"):
    # Initialize Pygame and set up the OpenGL context
    pygame.init()
    screen = pygame.display.set_mode((2000, 2000), DOUBLEBUF | OPENGL)
    pygame.display.set_caption('OpenGL Triangle Drawing')

    # Basic OpenGL setup
    glClearColor(0.0, 0.0, 0.0, 1.0)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE)

    # Convert vertices list to numpy array
    vertices = np.array(vertices, dtype=np.float32).flatten()

    # Vertex shader source code
    vertex_shader_source = """
    #version 330 core
    layout(location = 0) in vec2 position;
    layout(location = 1) in float vertex_opacity;
    out float frag_opacity;

    void main()
    {
        gl_Position = vec4(position, 0.0, 1.0);
        frag_opacity = vertex_opacity;
    }
    """
    
    # Fragment shader source code
    fragment_shader_source = """
    #version 330 core
    in float frag_opacity;
    out vec4 fragColor;

    void main()
    {
        fragColor = vec4(1.0, 1.0, 1.0, frag_opacity);
    }
    """
    
    # Compile vertex shader
    vertex_shader = glCreateShader(GL_VERTEX_SHADER)
    glShaderSource(vertex_shader, vertex_shader_source)
    glCompileShader(vertex_shader)
    if not glGetShaderiv(vertex_shader, GL_COMPILE_STATUS):
        print(glGetShaderInfoLog(vertex_shader).decode())
        raise RuntimeError("Vertex shader compilation failed")

    # Compile fragment shader
    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER)
    glShaderSource(fragment_shader, fragment_shader_source)
    glCompileShader(fragment_shader)
    if not glGetShaderiv(fragment_shader, GL_COMPILE_STATUS):
        print(glGetShaderInfoLog(fragment_shader).decode())
        raise RuntimeError("Fragment shader compilation failed")

    # Link shaders into a program
    shader_program = glCreateProgram()
    glAttachShader(shader_program, vertex_shader)
    glAttachShader(shader_program, fragment_shader)
    glLinkProgram(shader_program)
    if not glGetProgramiv(shader_program, GL_LINK_STATUS):
        print(glGetProgramInfoLog(shader_program).decode())
        raise RuntimeError("Shader program linking failed")

    glUseProgram(shader_program)

    # Set up vertex array object and vertex buffer object
    VAO = glGenVertexArrays(1)
    VBO = glGenBuffers(1)
    glBindVertexArray(VAO)

    glBindBuffer(GL_ARRAY_BUFFER, VBO)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_STATIC_DRAW)

    position = glGetAttribLocation(shader_program, 'position')
    glEnableVertexAttribArray(position)
    glVertexAttribPointer(position, 2, GL_FLOAT, GL_FALSE, 3 * ctypes.sizeof(ctypes.c_float), ctypes.c_void_p(0))

    vertex_opacity = glGetAttribLocation(shader_program, 'vertex_opacity')
    glEnableVertexAttribArray(vertex_opacity)
    glVertexAttribPointer(vertex_opacity, 1, GL_FLOAT, GL_FALSE, 3 * ctypes.sizeof(ctypes.c_float), ctypes.c_void_p(2 * ctypes.sizeof(ctypes.c_float)))

    # Event loop
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == QUIT:
                running = False

        glClear(GL_COLOR_BUFFER_BIT)

        # Draw the triangle
        glDrawArrays(GL_TRIANGLES, 0, len(vertices) // 3)

        size = screen.get_size()
        buffer = glReadPixels(0, 0, *size, GL_RGBA,  GL_UNSIGNED_BYTE)

        pygame.display.flip()


        # ajout tardiff pour recupere la matrix
        
        # Convert the buffer to a NumPy array
        pixel_matrix = np.frombuffer(buffer, dtype=np.uint8).reshape(size[1], size[0], 4)

        if save_image:
            save_image = False
            screen_surf = pygame.image.fromstring(buffer, size, "RGBA")
            pygame.image.save(screen_surf, image_filename)
            print(f"Image saved as {image_filename}")

    glDeleteVertexArrays(1, [VAO])
    glDeleteBuffers(1, [VBO])
    glDeleteProgram(shader_program)
    pygame.quit()

    return pixel_matrix



def Download_simplices(para_a):

    file_path = f'/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/totalHalo/triangle/Full_Image_{para_a}_data.pickle'
    with open(file_path, 'rb') as file:
        # Deserialize and retrieve the variable from the file
        loaded_data = pickle.load(file)

        print(f'The variable data_{para_a} has been loaded successfully.')

    return loaded_data

def adaptOpacity(vertices,opacity_factor=1.0):
    for ii in range(len(vertices)):
        vertices[ii][2] = vertices[ii][2]*opacity_factor/7
    return vertices


para_a_list = [80,81,82,83]
greyMat_tot = np.zeros((2000,2000))
for para_a in para_a_list:
# Generate a 2D triangle with opacity
#vertices = generate_2d_triangles_with_opacity()
    vertices = Download_simplices(para_a)
    vertices  = adaptOpacity(vertices,opacity_factor=5.0)
    # Draw the generated triangle
    pixelMatrix = draw_triangles_with_opacity(vertices, save_image=True, image_filename='testMatlpolib.jpg')

    #import matplotlib
    #matplotlib.use('Agg')
    greyMat = pixelMatrix[:,:,0]
    greyMat_tot += greyMat

plt.matshow(greyMat_tot,cmap='grey',norm='linear')
plt.savefig('/home/astro/ggsalami/ggsalami/TP4b/pythonAnalysis/pythonScript/totalHalo/testMatlpolib_plt.jpg')