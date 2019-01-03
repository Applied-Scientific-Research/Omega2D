R"(
#version 150

uniform mat4 Projection;
uniform vec4 this_color;
uniform float rad;
in vec4 quad_attr;
in float px;
in float py;
out vec4 base_color;
out vec2 txcoord;

void main() {
  // color pass-through
  base_color = this_color;

  // pass through texture coordinates
  txcoord = quad_attr.xy;

  // make 4 verts as a single primitive and set texture coords - see other shaders
  gl_Position = Projection * vec4(px + rad*quad_attr.x, py + rad*quad_attr.y, 0.f, 1.f);
}
)"
