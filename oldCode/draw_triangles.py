import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLUT import *
import numpy as np

def draw_triangles(vertices, opacity_factor = 1.0):
    # Initialize Pygame and set up the OpenGL context
    pygame.init()
    screen = pygame.display.set_mode((800, 600), DOUBLEBUF | OPENGL)
    pygame.display.set_caption('OpenGL Triangle Drawing')

    # Basic OpenGL setup
    glClearColor(0.0, 0.0, 0.0, 1.0)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE)

    # Vertices must be flattened into a single list
    vertices = np.array(vertices, dtype=np.float32).flatten()

    # Vertex shader source code
    vertex_shader_source = """
    #version 330 core
    layout(location = 0) in vec3 position;
    void main()
    {
        gl_Position = vec4(position, 1.0);
    }
    """
    
    # Fragment shader source code
    fragment_shader_source = """
    #version 330 core
    out vec4 fragColor;
    void main()
    {
        fragColor = vec4(1.0, 1.0, 1.0, """ + str(opacity_factor / (len(vertices) // 9)) + """);
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
    glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(0))
    glEnableVertexAttribArray(position)

    # Event loop
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == QUIT:
                running = False

        # Clear the screen
        glClear(GL_COLOR_BUFFER_BIT)

        # Draw the triangles
        glDrawArrays(GL_TRIANGLES, 0, len(vertices) // 3)

        # Swap buffers
        pygame.display.flip()

    # Clean up
    glDeleteVertexArrays(1, [VAO])
    glDeleteBuffers(1, [VBO])
    glDeleteProgram(shader_program)
    pygame.quit()

# Example usage:
vertices = [
    # Triangle 1
    [-0.5, -0.5, 0.0],
    [ 0.5, -0.5, 0.0],
    [ 0.0,  0.5, 0.0],
    # Triangle 2
    [ 0.5,  0.5, 0.0],
    [-0.5,  0.5, 0.0],
    [ 0.0, -0.5, 0.0]
]

draw_triangles(vertices)

def generate_random_triangles(num_triangles):
    """
    Generate an arbitrary number of random triangles.
    
    Parameters:
        num_triangles (int): The number of triangles to generate.

    Returns:
        list: A list of vertices for the generated triangles.
    """
    vertices = []
    for _ in range(num_triangles):
        # Each triangle has 3 vertices, and each vertex has 3 coordinates (x, y, z)
        triangle = np.random.uniform(-1.0, 1.0, (3, 3)).tolist()
        vertices.extend(triangle)
    return vertices

num_triangles = 10  # Specify the number of random triangles
vertices = generate_random_triangles(num_triangles)

# Draw the generated triangles
draw_triangles(vertices)
