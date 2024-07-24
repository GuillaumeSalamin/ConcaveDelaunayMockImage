import glfw
from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader
import numpy as np



def drawOneTriangle(points):
    # Initialize GLFW
    if not glfw.init():
        raise RuntimeError("Failed to initialize GLFW")

    # Disable window creation hints
    glfw.window_hint(glfw.VISIBLE, glfw.FALSE)
    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)

    # Create an off-screen context
    offscreen_width, offscreen_height = 800, 600
    offscreen_window = glfw.create_window(offscreen_width, offscreen_height, "", None, None)
    if not offscreen_window:
        glfw.terminate()
        raise RuntimeError("Failed to create GLFW off-screen window")

    # Make the off-screen OpenGL context current
    glfw.make_context_current(offscreen_window)

    # Vertex shader code
    vertex_shader = """
    #version 330 core
    layout (location = 0) in vec2 position;
    void main()
    {
        gl_Position = vec4(position, 0.0, 1.0);
    }
    """

    # Fragment shader code
    fragment_shader = """
    #version 330 core
    out vec4 outColor;
    void main()
    {
        outColor = vec4(1.0, 1.0, 1.0, 1.0);  // white color
    }
    """

    # Compile shaders and create shader program
    shader = compileProgram(
        compileShader(vertex_shader, GL_VERTEX_SHADER),
        compileShader(fragment_shader, GL_FRAGMENT_SHADER)
    )

    # Define the triangle vertices
    #vertices = np.array([
    #    -0.5, -0.5,
    #    0.5, -0.5,
    #    0.0, 0.5
    #], dtype=np.float32)
    vertices = points

    # Create a Vertex Buffer Object (VBO)
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_STATIC_DRAW)

    # Create a Vertex Array Object (VAO)
    vao = glGenVertexArrays(1)
    glBindVertexArray(vao)

    # Specify the layout of the vertex data
    position = glGetAttribLocation(shader, 'position')
    glEnableVertexAttribArray(position)
    glVertexAttribPointer(position, 2, GL_FLOAT, GL_FALSE, 0, None)

    # Unbind VAO and VBO
    glBindVertexArray(0)
    glBindBuffer(GL_ARRAY_BUFFER, 0)

    # Create a framebuffer object (FBO)
    fbo = glGenFramebuffers(1)
    glBindFramebuffer(GL_FRAMEBUFFER, fbo)

    # Create a texture to render to
    texture = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, texture)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, offscreen_width, offscreen_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, None)
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0)

    # Check framebuffer status
    if glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE:
        print("Framebuffer is not complete!")
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        glDeleteFramebuffers(1, [fbo])
        glfw.terminate()
        exit()

    # Render to FBO
    glBindFramebuffer(GL_FRAMEBUFFER, fbo)
    glViewport(0, 0, offscreen_width, offscreen_height)
    glClearColor(0.0, 0.0, 0.0, 1.0)
    glClear(GL_COLOR_BUFFER_BIT)

    glUseProgram(shader)
    glBindVertexArray(vao)
    glDrawArrays(GL_TRIANGLES, 0, 3)
    glBindVertexArray(0)
    glUseProgram(0)

    # Read pixels from the framebuffer
    pixels = glReadPixels(0, 0, offscreen_width, offscreen_height, GL_RGBA, GL_UNSIGNED_BYTE)
    pixels = np.frombuffer(pixels, dtype=np.uint8).reshape(offscreen_height, offscreen_width, 4)

    # Clean up resources
    glDeleteProgram(shader)
    glDeleteBuffers(1, [vbo])
    glDeleteVertexArrays(1, [vao])
    glDeleteTextures(1, [texture])
    glDeleteFramebuffers(1, [fbo])
    glfw.terminate()
    return pixels

def triangleDraw(triangle,alpha):
    tri = np.array(drawOneTriangle(np.array([triangle[0][0],triangle[0][1],triangle[1][0],triangle[1][1],triangle[2][0],triangle[2][1]],dtype=np.float32)))
    mat = tri[:,:,0]*alpha
    return mat