R"(
#version 150

uniform mat4 Projection;
uniform vec4 pos_color;
uniform vec4 neg_color;
in vec4 quad_attr;
in float px;
in float py;
in float r;
in float sx;
out vec4 base_color;
out vec2 txcoord;
out float strength;

void main() {
  // flip between base colors based on whether strength points in or out of the screen
  base_color = pos_color*step(0.0f, sx) + neg_color*step(0.0f, -sx);

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // magnitude of strength density scales fragment color
  strength = abs(sx)/(r*r);

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px + 3.5f*r*quad_attr.x, py + 3.5f*r*quad_attr.y, 0.f, 1.f);
}
)"
